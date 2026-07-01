#!/usr/bin/env python3
"""
collapse_chd.py - Surgical per-variant collapse update for the CHD concordance.

The collapsed deliverable (reference_outputs/chd_concordance_collapsed.{csv,xlsx})
carries CURATED columns - cohort, scope_note, control_citation (literature PMIDs/
DOIs) - that exist in NO other file. This script PRESERVES the existing collapsed
file and updates it surgically from the regenerated full results:

  * strict per-variant aggregation recomputed from the full file (max evidence
    across a variant's systems) - integer/string columns, so it changes only
    variants whose full-file four-way actually changed (the AlphaMissense-repair
    variants);
  * AlphaMissense class/score corrected ONLY for the repaired transposed variants
    (detected as rows whose collapsed class value is numeric);
  * NEW relaxed and pathogenic-only concordance columns appended (same
    max-across-systems rule on the relaxed / pathonly four-way);
  * every other column - structural hits, best_plddt, n_systems, min_tier, and all
    curated columns - preserved byte-for-byte.

The strict aggregation rule was verified to reproduce the existing collapse exactly
(structure/ddg/am/franklin hits; conc_n = max four-way; conc_denom = max denom).

Usage:
  python scripts/collapse_chd.py \
      --full      /tmp/chd_FIXED_regen.csv \
      --collapsed reference_outputs/chd_concordance_collapsed.csv \
      --out       reference_outputs/chd_concordance_collapsed.csv \
      --xlsx      reference_outputs/MAVIS_CHD_concordance_collapsed.xlsx
"""
import argparse
import pandas as pd

KEYS = ["gene", "ref_aa", "position", "alt_aa"]


def _is_float(x):
    try:
        float(x); return True
    except (ValueError, TypeError):
        return False


def aggregate(full):
    """Per-variant max-across-systems aggregation for strict/relaxed/pathonly."""
    f = full.copy()
    f["g"] = f["gene"].str.lower()
    recs = []
    for (g, r, p, a), grp in f.groupby(["g", "ref_aa", "position", "alt_aa"]):
        amc = grp["AlphaMissense"].dropna()
        ams = grp["AlphaMissense_pathogenicity"].dropna()
        recs.append({
            "k": (g, r, p, a),
            "am_hit_s": int(grp["am_vote_strict"].max()),
            "cn_s": int(grp["four_way_strict"].max()),
            "cd_s": int(grp["four_way_strict_denom"].max()),
            "cn_r": int(grp["four_way_relaxed"].max()),
            "cd_r": int(grp["four_way_relaxed_denom"].max()),
            "cn_p": int(grp["concordance_pathonly_strict_n"].max()),
            "cd_p": int(grp["concordance_pathonly_strict_denom"].max()),
            "am_class": (amc.iloc[0] if len(amc) else None),
            "am_score": (float(ams.iloc[0]) if len(ams) else None),
        })
    return pd.DataFrame(recs)


def _canon(s):
    def f(x):
        if pd.isna(x): return "<NA>"
        try: return f"{float(x):.6f}"
        except (ValueError, TypeError): return str(x)
    return s.map(f)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--full", required=True)
    ap.add_argument("--collapsed", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--xlsx", default=None, help="also write an .xlsx copy")
    a = ap.parse_args()

    C = pd.read_csv(a.collapsed).copy()
    F = pd.read_csv(a.full)
    A = aggregate(F)
    maps = {c: dict(zip(A["k"], A[c]))
            for c in ["am_hit_s", "cn_s", "cd_s", "cn_r", "cd_r", "cn_p", "cd_p", "am_class", "am_score"]}

    C["k"] = list(zip(C["gene"].str.lower(), C["ref_aa"], C["position"], C["alt_aa"]))
    pull = lambda col: [maps[col].get(k) for k in C["k"]]
    orig = C.copy()

    # --- strict aggregation recompute (integer/string -> churn-safe) ---
    C["am_hit"] = pull("am_hit_s")
    C["conc_n"] = pull("cn_s")
    C["conc_denom"] = pull("cd_s")
    C["concordance"] = [f"{n}/{d}" for n, d in zip(C["conc_n"], C["conc_denom"])]

    # --- AlphaMissense class/score: repair ONLY transposed rows (class is numeric) ---
    m_am = C["AlphaMissense_class"].map(_is_float)
    C.loc[m_am, "AlphaMissense_class"] = [maps["am_class"].get(k) for k in C.loc[m_am, "k"]]
    C.loc[m_am, "AlphaMissense_score"] = [maps["am_score"].get(k) for k in C.loc[m_am, "k"]]

    # --- NEW: relaxed + pathogenic-only concordance columns ---
    C["conc_relaxed_n"] = pull("cn_r"); C["conc_relaxed_denom"] = pull("cd_r")
    C["concordance_relaxed"] = [f"{n}/{d}" for n, d in zip(C["conc_relaxed_n"], C["conc_relaxed_denom"])]
    C["conc_pathonly_n"] = pull("cn_p"); C["conc_pathonly_denom"] = pull("cd_p")
    C["concordance_pathonly"] = [f"{n}/{d}" for n, d in zip(C["conc_pathonly_n"], C["conc_pathonly_denom"])]

    # --- diff report vs input collapsed ---
    changed = {c: int((_canon(orig[c]) != _canon(C[c])).sum())
               for c in orig.columns if c in C.columns and (_canon(orig[c]) != _canon(C[c])).any()}
    print("CHANGED existing cols (col: nrows):", changed)
    print("NEW cols:", [c for c in C.columns if c not in orig.columns and c != "k"])
    for c in ["am_hit", "concordance"]:
        mm = (_canon(orig[c]) != _canon(C[c])).values
        if mm.any():
            print(f"  {c}: {orig.loc[mm,'variant'].tolist()} -> {C.loc[mm,c].tolist()}")

    C = C.drop(columns=["k"])
    C.to_csv(a.out, index=False)
    print(f"WROTE {len(C)} rows x {C.shape[1]} cols -> {a.out}")
    if a.xlsx:
        C.to_excel(a.xlsx, index=False)
        print(f"WROTE xlsx -> {a.xlsx}")


if __name__ == "__main__":
    main()
