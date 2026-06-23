#!/usr/bin/env python3
"""Generic MAVIS runner - structural scoring for ANY genes via a YAML systems config.

  export FOLDX_BINARY=/path/to/foldx
  python run.py --config configs/my_systems.yaml --variants my_variants.csv \
                --structures ./structures --out results/my_run

Variant CSV needs columns: gene,ref_aa,position,alt_aa
A 'system' column is filled in automatically from the config (one row per system the
gene participates in). Pass --no-fanout if your CSV already has 'system'. Extra columns
(e.g. AlphaMissense, franklin) are carried through for the concordance step.
"""
import argparse
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "scripts"))

import pandas as pd
from mavis_v7.config import validate_config
from mavis_v7.config_io import load_config, autofanout


def main():
    ap = argparse.ArgumentParser(description="Run MAVIS structural scoring from a YAML systems config.")
    ap.add_argument("--config", required=True, help="YAML systems config (see configs/chd_systems.yaml)")
    ap.add_argument("--variants", required=True, help="Variant CSV: gene,ref_aa,position,alt_aa[,system]")
    ap.add_argument("--structures", default="structures", help="Directory of AlphaFold/PDB structures")
    ap.add_argument("--out", default="results/run", help="Output directory")
    ap.add_argument("--foldx", default=os.environ.get("FOLDX_BINARY", "foldx"),
                    help="FoldX binary (default: $FOLDX_BINARY)")
    ap.add_argument("--n-runs", type=int, default=5, help="FoldX replicate runs (default 5)")
    ap.add_argument("--no-fanout", action="store_true",
                    help="CSV already has a 'system' column; do not auto-expand")
    ap.add_argument("--dry-run", action="store_true",
                    help="Load config, validate structures, fan out variants - but do not run FoldX")
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()

    structures = Path(args.structures)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    configs = load_config(args.config)
    print(f"Loaded {len(configs)} system(s) from {args.config}")

    ok, issues = validate_config(configs, structures, out / "_preprocessed")
    if issues:
        print(f"\n{len(issues)} structure/config issue(s):")
        for i in issues:
            print("  -", i)
        if not ok and not args.dry_run:
            print("\nResolve the missing/invalid structures above (or re-run with --dry-run to just "
                  "plan). Aborting before FoldX.")
            sys.exit(2)
    else:
        print(f"All referenced structures found under {structures}/")

    df = pd.read_csv(args.variants, encoding="utf-8-sig")
    df.columns = [c.strip() for c in df.columns]
    for col in ("gene", "ref_aa", "position", "alt_aa"):
        if col not in df.columns:
            print(f"Variant CSV missing required column: {col}")
            sys.exit(2)

    if args.no_fanout or "system" in df.columns:
        if "system" not in df.columns:
            print("--no-fanout given but the CSV has no 'system' column.")
            sys.exit(2)
        expanded = df
        print(f"{len(expanded)} variant rows (using existing 'system' column)")
    else:
        expanded, unmatched = autofanout(df, configs)
        if unmatched:
            print(f"\nWARNING: {len(unmatched)} gene(s) match no configured system and were skipped: "
                  f"{', '.join(unmatched)}")
        print(f"Auto-fanned {len(df)} variant(s) -> {len(expanded)} (gene x system) rows")

    if len(expanded) == 0:
        print("No variant rows to run (nothing matched a system). Aborting.")
        sys.exit(2)

    expanded_path = out / "variants_expanded.csv"
    expanded.to_csv(expanded_path, index=False)
    print(f"Wrote expanded variants -> {expanded_path}")

    if args.dry_run:
        print("\n[dry-run] Config valid and variants expanded. Stopping before FoldX.")
        return

    from mavis_v7.pipeline import run_pipeline
    res = run_pipeline(
        configs=configs,
        variants_csv=expanded_path,
        structure_dir=structures,
        preprocessed_dir=out / "_preprocessed",
        output_dir=out,
        foldx_binary=Path(args.foldx),
        rotabase=None,
        n_runs=args.n_runs,
        dry_run=False,
        verbose=not args.quiet,
    )
    final = out / "structural_results.csv"
    res["df"].to_csv(final, index=False)
    print(f"\nDone: {len(res['df'])} rows -> {final}")
    print("Four-way concordance (AlphaMissense + Franklin/ClinVar) is the next, separate step: run "
          "scripts/apply_concordance_v5.py on this output (see docs/adding_your_own_genes.md).")


if __name__ == "__main__":
    main()
