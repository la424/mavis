#!/usr/bin/env python3
"""
MAVIS v7 — Neighborhood (Pipeline 2) extraction.

Ports the ±3 neighborhood contact extraction from CHD_pipeline_v6_0_A3-6
Cells 4 and 5 to MAVIS v7. Produces a second set of structural footprint
columns ('nbhd_*') that run alongside the existing single-residue Pipeline 1.

CONCEPT
-------
Pipeline 1 uses the variant position only:
  - monomer_n_contacts                  (intra-chain contacts at variant_pos)
  - multi_{partner}_inter_contacts       (inter-chain contacts at variant_pos)

Pipeline 2 extends to a ±3 residue window around the variant site:
  - nbhd_mono_contacts_weighted          (distance-weighted intra-chain contacts)
  - multi_{partner}_nbhd_inter_weighted   (distance-weighted inter-chain contacts)
  - multi_{partner}_nbhd_has_interface    (interface anywhere in ±3 window?)

Distance weights (Spec 2.2.1, CHD v6.0):
  offset=0  weight=1.00
  offset=±1 weight=0.75
  offset=±2 weight=0.50
  offset=±3 weight=0.25

Flanking residues with pLDDT < 50 are zeroed out (do not contribute), but do
NOT invalidate the window. The variant position ITSELF must have pLDDT >= 50
for the window to be evaluable at all.

DDG values are NOT recomputed. Neighborhood scoring is purely a structural
footprint metric; the same FoldX DDG columns feed both pipelines.

OUTPUTS
-------
New columns added to results CSV:
  nbhd_mono_contacts_weighted
  nbhd_mono_contacts_raw
  nbhd_mono_evaluable
  nbhd_mono_n_eval_positions
  multi_{partner}_nbhd_contacts_weighted
  multi_{partner}_nbhd_inter_weighted
  multi_{partner}_nbhd_evaluable
  multi_{partner}_nbhd_has_interface

Concordance columns and pipeline_agreement are added by
apply_concordance_v3.py (separate step).

USAGE
-----
    cd ~/mavis_v7
    python3 mavis_v7_neighborhood.py

Reads:   results/mavis_v7_results_corrected.csv + structure files
Writes:  results/mavis_v7_results_with_nbhd.csv

Requires the mavis_v7 package importable from the working directory (i.e.
run from ~/mavis_v7/ or with PYTHONPATH set).
"""
from __future__ import annotations
import argparse
import sys
from pathlib import Path
from typing import Dict

import pandas as pd

# mavis_v7 package imports — expected to be importable from CWD
try:
    from mavis_v7.config import build_benchmark_config
    from mavis_v7.structure_loading import (
        load_monomer, load_multimer, get_plddt,
        count_intra_contacts, count_interface,
    )
    from mavis_v7.metrics import collect_partner_labels, _partner_label_for_chain
except ImportError as e:
    print(f"ERROR: could not import mavis_v7 package: {e}", file=sys.stderr)
    print("       Run this script from ~/mavis_v7/ (with mavis_v7/ as a subdirectory)", file=sys.stderr)
    sys.exit(1)


# ============================================================================
# Neighborhood weights (CHD Spec 2.2.1)
# ============================================================================
NBHD_WEIGHTS = {0: 1.00, 1: 0.75, 2: 0.50, 3: 0.25}
PLDDT_GATE = 50.0


# ============================================================================
# Extraction functions (ported verbatim from CHD Cell 4 / Cell 5)
# ============================================================================
def extract_neighborhood_monomer(contact_map: Dict[int, int],
                                 plddt_map: Dict[int, float],
                                 variant_pos: int) -> Dict:
    """
    Extract ±3 neighborhood contacts with distance-weighted sum.

    Flanking residues with pLDDT < 50 are zeroed (do not contribute), but do
    not invalidate the window. Variant position must have pLDDT >= 50 for
    evaluability.

    Returns dict with keys nbhd_mono_contacts_weighted,
    nbhd_mono_contacts_raw, nbhd_mono_evaluable, nbhd_mono_n_eval_positions.
    """
    result = {
        "nbhd_mono_contacts_weighted": None,
        "nbhd_mono_contacts_raw": None,
        "nbhd_mono_evaluable": False,
        "nbhd_mono_n_eval_positions": 0,
    }
    if not contact_map or not plddt_map:
        return result

    variant_plddt = plddt_map.get(variant_pos)
    if variant_plddt is None or variant_plddt < PLDDT_GATE:
        return result

    result["nbhd_mono_evaluable"] = True
    weighted_sum = 0.0
    raw_sum = 0.0
    n_eval = 0
    for offset in range(-3, 4):
        pos = variant_pos + offset
        weight = NBHD_WEIGHTS.get(abs(offset), 0)
        pos_plddt = plddt_map.get(pos)
        if pos_plddt is not None and pos_plddt >= PLDDT_GATE:
            contacts = contact_map.get(pos, 0)
            weighted_sum += weight * contacts
            raw_sum += contacts
            n_eval += 1

    result["nbhd_mono_contacts_weighted"] = round(weighted_sum, 3)
    result["nbhd_mono_contacts_raw"] = raw_sum
    result["nbhd_mono_n_eval_positions"] = n_eval
    return result


def extract_neighborhood_multimer(contact_map: Dict[int, int],
                                  inter_map: Dict[int, int],
                                  plddt_map: Dict[int, float],
                                  variant_pos: int) -> Dict:
    """
    Extract ±3 neighborhood for a multimer chain. Returns weighted intra +
    inter contacts, and whether ANY residue in the window sits at the
    interface.
    """
    result = {
        "nbhd_multi_contacts_weighted": None,
        "nbhd_multi_inter_weighted": None,
        "nbhd_multi_evaluable": False,
        "nbhd_multi_has_interface": False,
    }
    if not contact_map or not plddt_map:
        return result

    variant_plddt = plddt_map.get(variant_pos)
    if variant_plddt is None or variant_plddt < PLDDT_GATE:
        return result

    result["nbhd_multi_evaluable"] = True
    weighted_contacts = 0.0
    weighted_inter = 0.0
    for offset in range(-3, 4):
        pos = variant_pos + offset
        weight = NBHD_WEIGHTS.get(abs(offset), 0)
        pos_plddt = plddt_map.get(pos)
        if pos_plddt is not None and pos_plddt >= PLDDT_GATE:
            weighted_contacts += weight * contact_map.get(pos, 0)
            inter = inter_map.get(pos, 0) if inter_map else 0
            weighted_inter += weight * inter
            if inter > 0:
                result["nbhd_multi_has_interface"] = True

    result["nbhd_multi_contacts_weighted"] = round(weighted_contacts, 3)
    result["nbhd_multi_inter_weighted"] = round(weighted_inter, 3)
    return result


# ============================================================================
# Monomer neighborhood pass
# ============================================================================
def compute_monomer_nbhd(df: pd.DataFrame, configs, structure_dir: Path,
                         verbose: bool = True) -> pd.DataFrame:
    """Extract nbhd_mono_* columns for every variant. Cache per (system, gene)."""
    cache = {}
    out_cols = {
        "nbhd_mono_contacts_weighted": [],
        "nbhd_mono_contacts_raw": [],
        "nbhd_mono_evaluable": [],
        "nbhd_mono_n_eval_positions": [],
    }

    for _, row in df.iterrows():
        sys_name = row["system"]
        gene = row["gene"]
        pos = int(row["position_mono"])
        key = (sys_name, gene)

        if key not in cache:
            cfg = configs[sys_name]
            struct, plddt, mspec = load_monomer(gene, cfg, structure_dir)
            if struct is not None and mspec is not None:
                # Compute contacts over the full chain so ±3 neighbors are covered.
                # count_intra_contacts with positions=None does this.
                contacts = count_intra_contacts(struct, mspec.chain_id, positions=None)
            else:
                contacts = {}
                plddt = {}
            cache[key] = {"contacts": contacts, "plddt": plddt}

        c = cache[key]
        nbhd = extract_neighborhood_monomer(c["contacts"], c["plddt"], pos)
        for k in out_cols:
            out_cols[k].append(nbhd[k])

    for k, v in out_cols.items():
        df[k] = v

    if verbose:
        n_eval = df["nbhd_mono_evaluable"].sum()
        print(f"  Monomer nbhd: {n_eval}/{len(df)} evaluable")
    return df


# ============================================================================
# Multimer neighborhood pass
# ============================================================================
def compute_multimer_nbhd(df: pd.DataFrame, configs, preprocessed_dir: Path,
                          source_dir: Path, verbose: bool = True) -> pd.DataFrame:
    """
    Extract multi_{partner}_nbhd_* columns for every variant × partner.
    Cache per system, then per (chain_a, chain_b) pair within a system.
    """
    all_partners = collect_partner_labels(configs)

    # Initialize nbhd columns for every partner label
    for pl in all_partners:
        df[f"multi_{pl}_nbhd_contacts_weighted"] = pd.NA
        df[f"multi_{pl}_nbhd_inter_weighted"] = pd.NA
        df[f"multi_{pl}_nbhd_evaluable"] = False
        df[f"multi_{pl}_nbhd_has_interface"] = False

    # Cache structure per system
    struct_cache = {}
    load_failures = []  # track systems where load_multimer returned None

    for idx, row in df.iterrows():
        sys_name = row["system"]
        gene = row["gene"]
        pos = int(row["position_multi"])
        cfg = configs[sys_name]
        multi = cfg.multimer
        my_chain = multi.chain_map.get(gene)
        if my_chain is None:
            continue

        if sys_name not in struct_cache:
            struct = load_multimer(cfg, preprocessed_dir, source_dir)
            if struct is None:
                # Build the path load_multimer would have tried, for diagnostics
                expected_dir = source_dir if multi.structure_type == "AF" else preprocessed_dir
                expected_path = Path(expected_dir) / multi.pdb_file
                msg = (f"  WARNING: load_multimer returned None for system {sys_name!r} "
                       f"(expected PDB at {expected_path}, exists={expected_path.exists()}); "
                       f"all variants in this system will have nbhd_multi_evaluable=False")
                print(msg, file=sys.stderr)
                load_failures.append((sys_name, str(expected_path)))
            struct_cache[sys_name] = {"struct": struct, "chains": {}}
        struct = struct_cache[sys_name]["struct"]
        if struct is None:
            continue

        # pLDDT of variant chain (used for gating the window)
        chain_cache = struct_cache[sys_name]["chains"]
        if my_chain not in chain_cache:
            chain_cache[my_chain] = {
                "plddt": get_plddt(struct, my_chain, plddt_gate=multi.plddt_gate),
                "contacts": count_intra_contacts(struct, my_chain, positions=None),
                "inter_by_partner": {},
            }
        my_plddt = chain_cache[my_chain]["plddt"]
        my_contacts = chain_cache[my_chain]["contacts"]

        for p_chain in multi.pairwise_partners.get(gene, []):
            pl = _partner_label_for_chain(multi, p_chain)
            if p_chain not in chain_cache[my_chain]["inter_by_partner"]:
                inter, iface = count_interface(struct, my_chain, p_chain)
                chain_cache[my_chain]["inter_by_partner"][p_chain] = (inter, iface)
            inter, _ = chain_cache[my_chain]["inter_by_partner"][p_chain]

            nbhd = extract_neighborhood_multimer(my_contacts, inter, my_plddt, pos)
            df.at[idx, f"multi_{pl}_nbhd_contacts_weighted"] = nbhd["nbhd_multi_contacts_weighted"]
            df.at[idx, f"multi_{pl}_nbhd_inter_weighted"]    = nbhd["nbhd_multi_inter_weighted"]
            df.at[idx, f"multi_{pl}_nbhd_evaluable"]          = nbhd["nbhd_multi_evaluable"]
            df.at[idx, f"multi_{pl}_nbhd_has_interface"]      = nbhd["nbhd_multi_has_interface"]

    if verbose:
        # Count variants with any partner evaluable (true diagnostic of pipeline health)
        eval_cols = [c for c in df.columns if c.startswith("multi_") and c.endswith("_nbhd_evaluable")]
        n_any_eval = df[eval_cols].apply(lambda r: any(bool(v) for v in r), axis=1).sum() if eval_cols else 0
        # Count variants with any partner at nbhd interface
        iface_cols = [c for c in df.columns if c.startswith("multi_") and c.endswith("_nbhd_has_interface")]
        n_any_iface = df[iface_cols].apply(lambda r: any(bool(v) for v in r), axis=1).sum() if iface_cols else 0
        print(f"  Multimer nbhd: {n_any_eval}/{len(df)} evaluable, "
              f"{n_any_iface}/{len(df)} have ≥1 partner at nbhd interface")
        if load_failures:
            print(f"  NOTE: {len(load_failures)} system(s) failed to load; see warnings above")

    return df


# ============================================================================
# Main driver
# ============================================================================
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input",        default="results/mavis_v7_results_corrected.csv")
    ap.add_argument("--output",       default="results/mavis_v7_results_with_nbhd.csv")
    ap.add_argument("--structures",   default="structures",
                    help="Directory with monomer AF PDBs and AF multimer PDBs")
    ap.add_argument("--preprocessed", default="processed",
                    help="Directory with preprocessed xray/NMR multimer PDBs (*_processed.pdb)")
    ap.add_argument("--multimer-src", default="structures",
                    help="Directory with AF multimer PDBs (fold_*_model_0.pdb)")
    ap.add_argument("--config",       default="mavis_v7/system_configs.yaml",
                    help="System config file (YAML). Default matches mavis_v7 layout.")
    args = ap.parse_args()

    in_path = Path(args.input); out_path = Path(args.output)
    if not in_path.exists():
        print(f"ERROR: input CSV not found: {in_path}", file=sys.stderr)
        sys.exit(1)

    # Fail fast on obvious path-mismatch, preserving the silent-0/44 lesson
    structures_dir = Path(args.structures)
    preprocessed_dir = Path(args.preprocessed)
    multimer_src = Path(args.multimer_src)
    missing = [str(p) for p in [structures_dir, preprocessed_dir, multimer_src] if not p.exists()]
    if missing:
        print(f"WARNING: the following path(s) do not exist: {missing}", file=sys.stderr)
        print(f"         If this is wrong, rerun with --structures/--preprocessed/--multimer-src", file=sys.stderr)

    df = pd.read_csv(in_path)
    print(f"Loaded {len(df)} variants, {len(df.columns)} columns")

    # position_mono / position_multi may be the same as position, but
    # safe-handle: if those columns aren't present, fall back to 'position'.
    if "position_mono" not in df.columns:
        df["position_mono"] = df["position"]
    if "position_multi" not in df.columns:
        df["position_multi"] = df["position"]

    # Load system configs
    configs = build_benchmark_config()
    print(f"Loaded {len(configs)} system configs")

    # Monomer nbhd
    print("\n[1/2] Extracting monomer neighborhoods...")
    df = compute_monomer_nbhd(df, configs, structures_dir)

    # Multimer nbhd
    print("\n[2/2] Extracting multimer neighborhoods...")
    df = compute_multimer_nbhd(df, configs, preprocessed_dir, multimer_src)

    # Save
    df.to_csv(out_path, index=False)
    print(f"\n✓ Wrote {out_path}  ({len(df)} rows, {len(df.columns)} cols)")

    # Summary
    print("\nNeighborhood summary:")
    print(f"  nbhd_mono_evaluable:          {df['nbhd_mono_evaluable'].sum()}/{len(df)}")
    nbhd_eval_cols = [c for c in df.columns if c.startswith("multi_") and c.endswith("_nbhd_evaluable")]
    n_any_eval = df[nbhd_eval_cols].apply(lambda r: any(bool(v) for v in r), axis=1).sum() if nbhd_eval_cols else 0
    nbhd_iface_cols = [c for c in df.columns if c.startswith("multi_") and c.endswith("_nbhd_has_interface")]
    n_any_iface = df[nbhd_iface_cols].apply(lambda r: any(bool(v) for v in r), axis=1).sum() if nbhd_iface_cols else 0
    print(f"  ≥1 nbhd partner evaluable:    {n_any_eval}/{len(df)}")
    print(f"  ≥1 nbhd interface detected:   {n_any_iface}/{len(df)}")


if __name__ == "__main__":
    main()
