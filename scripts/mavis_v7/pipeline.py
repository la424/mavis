"""
Top-level pipeline orchestration.

Usage:
    from mavis_v7.config import build_benchmark_config
    from mavis_v7.pipeline import run_pipeline

    results = run_pipeline(
        configs=build_benchmark_config(),
        variants_csv="benchmark_variants_v5.csv",
        structure_dir="structures/",
        preprocessed_dir="processed/",
        foldx_binary="/path/to/foldx",
        output_dir="results/",
        dry_run=True,     # skip real FoldX, use mock values
    )
"""

from pathlib import Path
from typing import Dict, Optional
import pandas as pd

from .config import SystemConfig
from .preprocessing import preprocess_all
from .variant_loading import load_variants
from .metrics import (
    collect_partner_labels, compute_monomer_metrics, compute_multimer_metrics,
    populate_ddg_confidence,
)
from .foldx_runner import compute_three_axis_ddg, _REPAIRED_CACHE
from .mechanism import (
    compute_disruption_score, assign_tier, compute_evaluability, classify_mechanism,
)


def run_pipeline(
    configs: Dict[str, SystemConfig],
    variants_csv: Path,
    structure_dir: Path,
    preprocessed_dir: Path,
    foldx_binary: Optional[Path] = None,
    rotabase: Optional[Path] = None,
    output_dir: Path = Path("results"),
    n_runs: int = 5,
    dry_run: bool = False,
    verbose: bool = True,
) -> Dict:
    """
    Run the full pipeline. Returns a dict with:
      'df'           — the variant dataframe with all computed columns
      'preprocessing_provenance'
      'partner_labels'
      'foldx_summary' — counts of successes/failures per axis
    """
    structure_dir = Path(structure_dir)
    preprocessed_dir = Path(preprocessed_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    foldx_work_dir = output_dir / "foldx_runs"
    foldx_work_dir.mkdir(parents=True, exist_ok=True)

    if verbose:
        print("=" * 70)
        print("MAVIS v7 — Pipeline Run")
        print("=" * 70)

    # ---------------- 1. Preprocessing ----------------
    if verbose:
        print("\n[1/6] Preprocessing crystal/NMR structures")
    provenance = preprocess_all(configs, structure_dir, preprocessed_dir, verbose=verbose)

    # ---------------- 2. Variant loading ----------------
    if verbose:
        print("\n[2/6] Loading variants")
    df = load_variants(variants_csv, configs, verbose=verbose)

    # ---------------- 3. Monomer structural metrics ----------------
    if verbose:
        print("\n[3/6] Computing monomer structural metrics")
    df = compute_monomer_metrics(df, configs, structure_dir, verbose=verbose)

    # ---------------- 4. Multimer structural metrics ----------------
    if verbose:
        print("\n[4/6] Computing multimer structural metrics")
    df = compute_multimer_metrics(
        df, configs, preprocessed_dir, structure_dir, verbose=verbose
    )
    partner_labels = collect_partner_labels(configs)

    # ---------------- 5. FoldX three-axis DDG ----------------
    if verbose:
        mode = "DRY_RUN (mock values)" if dry_run else "LIVE"
        print(f"\n[5/6] Computing FoldX three-axis DDG ({mode}, n_runs={n_runs})")

    # Initialize DDG columns
    df["ddg_monomer"] = pd.NA
    df["ddg_monomer_sd"] = pd.NA
    df["ddg_monomer_runs"] = [[] for _ in range(len(df))]
    for pl in partner_labels:
        df[f"ddg_fold_{pl}"] = pd.NA
        df[f"ddg_fold_{pl}_sd"] = pd.NA
        df[f"ddg_binding_{pl}"] = pd.NA
        df[f"ddg_binding_{pl}_sd"] = pd.NA
    df["foldx_errors"] = ""

    n_ok_mono = 0
    n_ok_multi = 0
    n_fail = 0

    for idx, row in df.iterrows():
        sys_name = row["system"]
        cfg = configs[sys_name]
        multi = cfg.multimer

        variant = row.to_dict()
        try:
            ddg = compute_three_axis_ddg(
                variant, cfg,
                structure_dir=structure_dir,
                preprocessed_dir=preprocessed_dir,
                foldx_work_dir=foldx_work_dir,
                foldx_binary=foldx_binary,
                rotabase=rotabase,
                n_runs=n_runs,
                dry_run=dry_run,
                verbose=False,
            )
        except Exception as e:
            ddg = {"ddg_monomer": None, "ddg_monomer_sd": None, "ddg_monomer_runs": [],
                   "ddg_fold_by_partner": {}, "ddg_binding_by_partner": {},
                   "ddg_binding_sd_by_partner": {}, "ddg_binding_runs_by_partner": {},
                   "foldx_errors": [f"exception: {e}"]}

        df.at[idx, "ddg_monomer"] = ddg["ddg_monomer"]
        df.at[idx, "ddg_monomer_sd"] = ddg["ddg_monomer_sd"]
        df.at[idx, "ddg_monomer_runs"] = ddg["ddg_monomer_runs"]
        if ddg["ddg_monomer"] is not None:
            n_ok_mono += 1

        # Map partner chain results → partner labels
        for p_chain, (fold_m, fold_sd, _) in ddg["ddg_fold_by_partner"].items():
            # Find the partner label for this chain under this gene
            pl = None
            for gene_name, ch in multi.chain_map.items():
                if ch == p_chain:
                    pl = gene_name
                    break
            if pl is None:
                pl = f"chain_{p_chain}"
            df.at[idx, f"ddg_fold_{pl}"] = fold_m
            df.at[idx, f"ddg_fold_{pl}_sd"] = fold_sd

        for p_chain, bind_val in ddg["ddg_binding_by_partner"].items():
            pl = None
            for gene_name, ch in multi.chain_map.items():
                if ch == p_chain:
                    pl = gene_name
                    break
            if pl is None:
                pl = f"chain_{p_chain}"
            df.at[idx, f"ddg_binding_{pl}"] = bind_val
            # Also store the SD across the 5 mutant replicates
            sd_val = ddg.get("ddg_binding_sd_by_partner", {}).get(p_chain)
            if sd_val is not None:
                df.at[idx, f"ddg_binding_{pl}_sd"] = sd_val

        if ddg["ddg_binding_by_partner"]:
            n_ok_multi += 1
        if ddg["foldx_errors"]:
            df.at[idx, "foldx_errors"] = ";".join(ddg["foldx_errors"])
            if ddg["ddg_monomer"] is None and not ddg["ddg_binding_by_partner"]:
                n_fail += 1

    if verbose:
        print(f"  Monomer DDG:  {n_ok_mono}/{len(df)}")
        print(f"  Multimer DDG: {n_ok_multi}/{len(df)}")
        print(f"  Total FoldX failures: {n_fail}")

    # Populate per-partner confidence flags
    df = populate_ddg_confidence(df, partner_labels)

    # ---------------- 6. Scoring, tier, mechanism ----------------
    if verbose:
        print("\n[6/6] Applying scoring, tier, and mechanism classification")

    score_df = df.apply(lambda r: compute_disruption_score(r, partner_labels), axis=1)
    for c in score_df.columns:
        df[c] = score_df[c]

    df["mavis_tier"] = df["mavis_score"].apply(assign_tier)
    df = compute_evaluability(df, partner_labels)

    # v5: multi-threshold mechanism classification
    # Compute mech calls at five thresholds:
    #   - 4 uniform (1.0, 1.5, 2.0, 2.5)
    #   - 1 Sapozhnikov per-axis (mono=2.9, fold=2.9, bind=3.5)
    # Default `mavis_mechanism` (no suffix) aliases the t10 column for
    # backward compatibility with downstream scripts.
    threshold_specs = [
        ('t10',  1.0),
        ('t15',  1.5),
        ('t20',  2.0),
        ('t25',  2.5),
        ('tSAP', {'monomer': 2.9, 'fold': 2.9, 'binding': 3.5}),
    ]

    for tag, thr in threshold_specs:
        mech_results = df.apply(
            lambda r, _t=thr: classify_mechanism(r, partner_labels, ddg_destab=_t), axis=1
        )
        df[f"mavis_mechanism_{tag}"] = [r[0] for r in mech_results]
        df[f"mavis_mechanism_partner_{tag}"] = [r[1] for r in mech_results]
        df[f"mavis_external_evidence_flag_{tag}"] = [r[2] for r in mech_results]

    # Backward-compatible aliases at t=1.0 (the historical default)
    df["mavis_mechanism"] = df["mavis_mechanism_t10"]
    df["mavis_mechanism_partner"] = df["mavis_mechanism_partner_t10"]
    df["mavis_external_evidence_flag"] = df["mavis_external_evidence_flag_t10"]

    if verbose:
        print(f"\n  Mechanism distribution at t=1.0 (mavis_mechanism):")
        print(df["mavis_mechanism"].value_counts().to_string())
        print(f"\n  Mechanism distribution at Sapozhnikov per-axis:")
        print(df["mavis_mechanism_tSAP"].value_counts().to_string())
        print(f"\n  Tier distribution:")
        print(df["mavis_tier"].value_counts().to_string())

    return {
        "df": df,
        "preprocessing_provenance": provenance,
        "partner_labels": partner_labels,
        "foldx_summary": {
            "n_ok_monomer": n_ok_mono,
            "n_ok_multimer": n_ok_multi,
            "n_fail": n_fail,
        },
    }
