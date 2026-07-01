#!/usr/bin/env python3
"""
apply_chd_concordance.py - Versioned CHD four-way concordance (re)assembly.

Re-derives the CHD four-way concordance (structure tier . ddG . AlphaMissense .
Franklin) on top of the LOCKED structural layer, using the shared engine
mavis_v7.concordance.add_concordance - the exact engine that produced the
reference (verified to round-trip on the reference byte-for-byte).

INPUT MODEL (why it consumes the concordance reference, not chd_structural_results.csv):
  The canonical CHD structural layer embeds a pLDDT reconciliation applied AFTER
  FoldX and BEFORE the final concordance; that reconciled intermediate was not
  persisted, so reference_outputs/chd_concordance_results_FIXED.csv is the
  authoritative carrier of the locked structural columns. This script treats those
  structural columns as fixed (NO FoldX, NO structural re-derivation), repairs the
  known AlphaMissense / AlphaMissense_pathogenicity column transposition
  (kpna1/kpna6/tcf7l1) - cross-checked against the committed annotation table -
  recomputes the four-way concordance via the engine, and appends a pathogenic-only
  Franklin threshold (no VUS). On already-correct inputs it is an exact identity.

Usage:
  python scripts/apply_chd_concordance.py \
      --in  reference_outputs/chd_concordance_results_FIXED.csv \
      --am  variants_with_alphamissense_and_franklin_expanded.csv \
      --out reference_outputs/chd_concordance_results_FIXED.csv
"""
import argparse, sys
from pathlib import Path
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from mavis_v7.concordance import add_concordance, std_franklin

FRANKLIN_PATHONLY = {"pathogenic", "likely pathogenic"}   # NO VUS counted
KEYS = ["gene", "ref_aa", "position", "alt_aa"]


def _is_float(x):
    try:
        float(x); return True
    except (ValueError, TypeError):
        return False


def repair_am_transposition(df, am_path):
    """Swap AlphaMissense <-> AlphaMissense_pathogenicity on rows where they are
    transposed (numeric score in the class column AND non-numeric class in the
    score column). Cross-check repaired class values against the committed table."""
    cls, sco = "AlphaMissense", "AlphaMissense_pathogenicity"
    mask = df[cls].map(_is_float) & ~df[sco].map(_is_float)
    n = int(mask.sum())
    if n:
        tmp = df.loc[mask, cls].copy()
        df.loc[mask, cls] = df.loc[mask, sco].values
        df.loc[mask, sco] = tmp.values
    am = pd.read_csv(am_path)[KEYS + [cls]].copy()
    for k in KEYS:
        am[k] = am[k].astype(str)
    chk = df[KEYS + [cls]].copy()
    for k in KEYS:
        chk[k] = chk[k].astype(str)
    chk = chk.merge(am.rename(columns={cls: cls + "_src"}), on=KEYS, how="left")
    bad = int((chk[cls].astype(str).str.lower()
               != chk[cls + "_src"].astype(str).str.lower()).sum())
    return df, n, bad


def add_pathonly(df):
    """Pathogenic-only Franklin threshold: mirror the engine's strict/relaxed/t3
    assembly, substituting a path/likely-path-only Franklin vote (no VUS)."""
    def row(r):
        fr = std_franklin(r.get("franklin"))
        fp = 1 if (isinstance(fr, str) and fr.lower() in FRANKLIN_PATHONLY) else 0
        se = bool(r.get("structure_evaluable", False)); de = bool(r.get("ddg_evaluable", False))
        ae = bool(r.get("am_evaluable", False));        fe = bool(r.get("franklin_evaluable", False))

        def asm(sv, dv, av):
            s = d = 0
            if se: s += sv; d += 1
            if de: s += dv; d += 1
            if ae: s += av; d += 1
            if fe: s += fp; d += 1
            return s, max(d, 1)

        sv_s = int(r.get("structure_vote_strict", 0)); sv_t3 = int(r.get("structure_vote_t3", 0))
        ss, sd = asm(sv_s,  int(r.get("ddg_vote_strict", 0)),  int(r.get("am_vote_strict", 0)))
        rs, rd = asm(sv_s,  int(r.get("ddg_vote_relaxed", 0)), int(r.get("am_vote_relaxed", 0)))
        ts, td = asm(sv_t3, int(r.get("ddg_vote_strict", 0)),  int(r.get("am_vote_strict", 0)))
        return pd.Series({
            "franklin_vote_pathonly": fp,
            "concordance_pathonly_strict": f"{ss}/{sd}", "concordance_pathonly_strict_n": ss, "concordance_pathonly_strict_denom": sd,
            "concordance_pathonly_relaxed": f"{rs}/{rd}", "concordance_pathonly_relaxed_n": rs, "concordance_pathonly_relaxed_denom": rd,
            "concordance_pathonly_t3": f"{ts}/{td}", "concordance_pathonly_t3_n": ts, "concordance_pathonly_t3_denom": td,
        })
    add = df.apply(row, axis=1)
    for c in add.columns:
        df[c] = add[c]
    return df


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--in", dest="inp", required=True)
    ap.add_argument("--am", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    df = pd.read_csv(a.inp)
    partners = sorted({c[len("ddg_fold_"):] for c in df.columns
                       if c.startswith("ddg_fold_") and not c.endswith("_sd")})
    df, n_rep, n_bad = repair_am_transposition(df, a.am)
    print(f"AM transposition: repaired {n_rep} rows | cross-check mismatches vs source: {n_bad}")
    if n_bad:
        sys.exit(f"ERROR: {n_bad} AlphaMissense values disagree with committed table after repair")
    df = add_concordance(df, partners)   # overwrites engine concordance cols in place
    df = add_pathonly(df)                # appends pathonly cols
    df.to_csv(a.out, index=False)
    print(f"WROTE {len(df)} rows x {df.shape[1]} cols -> {a.out}")


if __name__ == "__main__":
    main()
