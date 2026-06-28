"""
prepare_chd_input.py — turn the raw CHD variant CSV into the pipeline's per-system input.

The raw file (variants_with_alphamissense_and_franklin_expanded.csv) has columns:
    gene, ref_aa, position, alt_aa, AlphaMissense, AlphaMissense_pathogenicity, franklin
It has NO `system` column. The package pipeline routes each variant to a system via that
column, and CHD scoring is BIDIRECTIONAL — a variant in gene G is scored against EVERY kept
system in which G appears (as hub OR as partner). So one input variant fans out to multiple
rows, one per (variant, system) pair.

This script:
  1. Reads the raw CSV (utf-8-sig — the file has a BOM).
  2. Excludes gene in {zic5, zic2} (no complexes).  163 -> 145.
  3. For each variant, finds every kept system whose chain_map contains that gene, and emits
     one row per system with the `system` column set.
  4. Writes the fanned-out CSV the pipeline consumes.

Dedup note: the supplied 163-file already reflects the DVL2 R367G/R367Q dedup (165->163), so no
further dedup is applied here. (If you ever start from the 165 pre-dedup file, dedup on
[gene,position,ref_aa,alt_aa] first.)

Usage:
    python3 prepare_chd_input.py \
        --in  variants_with_alphamissense_and_franklin_expanded.csv \
        --out chd_input_per_system.csv
"""

import argparse
import csv
import sys

EXCLUDE_GENES = {"zic5", "zic2"}


def gene_to_systems(configs):
    """Map each gene -> list of system names whose chain_map contains it."""
    g2s = {}
    for sysname, sc in configs.items():
        for gene in sc.multimer.chain_map.keys():
            g2s.setdefault(gene.lower(), []).append(sysname)
    return g2s


def prepare(in_path, out_path, configs):
    with open(in_path, encoding="utf-8-sig") as f:
        rows = list(csv.DictReader(f))

    # Guard: AlphaMissense must be a class string and AlphaMissense_pathogenicity a float.
    # Catches the column-transposition seen in the kpna1/kpna6/tcf7l1 annotation batch,
    # which otherwise silently zeroes am_hit in downstream concordance.
    def _isnum(x):
        try:
            float(x); return True
        except (TypeError, ValueError):
            return False
    _VALID_AM = {"likely_pathogenic", "likely_benign", "ambiguous", "pathogenic", "benign"}
    _bad = [(r.get("gene"), f"{r.get('ref_aa','')}{r.get('position','')}{r.get('alt_aa','')}",
             r.get("AlphaMissense"), r.get("AlphaMissense_pathogenicity"))
            for r in rows
            if str(r.get("AlphaMissense", "")).strip().lower() not in _VALID_AM
            or not _isnum(r.get("AlphaMissense_pathogenicity"))]
    if _bad:
        raise ValueError(
            f"AlphaMissense columns look transposed/garbled in {len(_bad)} row(s) "
            f"(want AlphaMissense=class string, AlphaMissense_pathogenicity=float). "
            f"First few: {_bad[:5]}")

    g2s = gene_to_systems(configs)

    # Report genes in the input that map to NO kept system (will be dropped)
    input_genes = {r["gene"].lower() for r in rows}
    mapped = set(g2s.keys())
    unmapped = sorted(input_genes - mapped - EXCLUDE_GENES)
    if unmapped:
        print(f"  NOTE: input genes with no kept system (dropped): {unmapped}", file=sys.stderr)

    out_rows = []
    n_excluded = 0
    n_unmapped = 0
    for r in rows:
        g = r["gene"].lower()
        if g in EXCLUDE_GENES:
            n_excluded += 1
            continue
        systems = g2s.get(g, [])
        if not systems:
            n_unmapped += 1
            continue
        for sysname in systems:
            nr = dict(r)
            nr["gene"] = g                     # normalize lowercase
            nr["system"] = sysname
            out_rows.append(nr)

    # Write — original columns + system
    fieldnames = list(rows[0].keys())
    if "system" not in fieldnames:
        fieldnames = fieldnames + ["system"]
    with open(out_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(out_rows)

    # Summary
    uniq_variants = {(r["gene"].lower(), r["position"], r["ref_aa"], r["alt_aa"])
                     for r in rows if r["gene"].lower() not in EXCLUDE_GENES
                     and g2s.get(r["gene"].lower())}
    print(f"  input rows:            {len(rows)}")
    print(f"  excluded (zic5/zic2):  {n_excluded}")
    print(f"  unmapped (no system):  {n_unmapped}")
    print(f"  unique analyzed variants: {len(uniq_variants)}")
    print(f"  fanned-out (variant,system) rows written: {len(out_rows)} -> {out_path}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="in_path",
                    default="variants_with_alphamissense_and_franklin_expanded.csv")
    ap.add_argument("--out", dest="out_path", default="chd_input_per_system.csv")
    args = ap.parse_args()

    # import the package config (adjust path/import as needed in your env)
    from mavis_v7.config import build_benchmark_config  # noqa
    try:
        from mavis_v7.build_chd_config import build_chd_config
    except ImportError:
        # if build_chd_config was pasted into config.py instead:
        from mavis_v7.config import build_chd_config

    prepare(args.in_path, args.out_path, build_chd_config())
