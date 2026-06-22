#!/usr/bin/env python3
"""
MAVIS v7 — mech_consistency pLDDT-reconciliation patch + runner.

PURPOSE
-------
The Track-B grader `grade_mechanism_consistency` operates on the raw pipeline
mechanism call and does NOT honor the per-variant pLDDT interface-exclusions that
`compute_structural_agreement` already applies. This causes variants whose ONLY
inconsistency is a fold/interface false-fire AT A LOW-pLDDT INTERFACE (which the
batch decisions explicitly excluded from scoring) to be graded `inconsistent`,
double-standarding the two paired headline metrics.

This module applies a minimal, auditable reconciliation: a variant is dropped
from the mech_consistency denominator (graded "N/A (pLDDT-excluded)") iff
  (a) it carries a logged pLDDT interface-exclusion flag, AND
  (b) the inconsistency is DRIVEN BY the excluded axis — i.e. the pipeline's
      mechanism call fires the very axis that was flagged low-confidence.
Variants whose inconsistency stands on a confidently-modeled axis are NOT excused.

This mirrors exactly the structural_agreement convention (omit the excluded axis),
so both paired metrics now use one gating rule.

HOW TO USE LOCALLY
------------------
Two integration options for `apply_concordance_v5.py`:

  OPTION 1 (wrapper, no edit to the rubric body — recommended):
    After the existing per-threshold grading loop in main() Stage 6, post-process
    the mech_consistency_{suf} columns with `reconcile_grade(...)` below, using the
    PLDDT_INTERFACE_EXCLUSIONS table. Re-derive mech_consistency_summary from t25.

  OPTION 2 (inline guard): paste `PLDDT_INTERFACE_EXCLUSIONS` and the guard block
    (see GUARD_SNIPPET string) at the TOP of grade_mechanism_consistency, right
    after the `empty = []` line. The guard returns ("N/A (pLDDT-excluded)", [], [])
    before the rubric runs, when the call fires the excluded axis.

The runner at the bottom reproduces the reconciled headline against the corrected CSV
using OPTION 1 (so the uploaded script body is untouched and the number is auditable).
"""
import sys, types, argparse
from pathlib import Path

# ---------------------------------------------------------------------------
# Logged pLDDT interface-exclusions (variant -> (partner, plddt, excluded_axes))
# excluded_axes: which axis/axes are flagged low-confidence and thus must not be
# counted as a scored false-positive/disruptor for this variant.
#   - 'fold'    -> the ddg_fold_{partner} / multimer-fold interface call is excluded
#   - 'binding' -> the ddg_binding_{partner} interface call is excluded
# Source: corrections-log Batches 9, 10, 11.
# ---------------------------------------------------------------------------
PLDDT_INTERFACE_EXCLUSIONS = {
    ('tnni3', 'R145G'): {'partner': 'tnnc1',   'plddt': 54.93, 'axes': {'fold', 'binding'}},  # B9
    ('tnni3', 'R145Q'): {'partner': 'tnnc1',   'plddt': 54.93, 'axes': {'fold', 'binding'}},  # B9
    ('calm1', 'D96V'):  {'partner': 'cacna1c', 'plddt': 69.45, 'axes': {'fold', 'binding'}},  # B10 (borderline)
    ('smad4', 'D351H'): {'partner': 'smad3',   'plddt': 57.37, 'axes': {'binding'}},          # B11 (bind only; fold->unknown already)
}

# Which mechanism labels assert which axes (mirror of MECH_LABEL_ASSERTIONS, fold collapsed).
# Used to decide whether the inconsistency is DRIVEN BY an excluded axis.
def _call_fires_excluded_axis(mech_call, excluded_axes, mech_assertions):
    info = mech_assertions.get(mech_call)
    if info is None:
        return False
    fold_m_a, fold_c_a, bind_a, iface_a, direction = info
    fold_a = fold_m_a or fold_c_a
    if 'fold' in excluded_axes and fold_a:
        return True
    if 'binding' in excluded_axes and bind_a:
        return True
    # interface-only DNG fire at an excluded interface also counts as excluded
    if iface_a and ('fold' in excluded_axes or 'binding' in excluded_axes) and not (fold_a or bind_a):
        return True
    return False


def reconcile_grade(gene, variant, raw_grade, mech_call, mech_assertions):
    """
    OPTION-1 post-processor. Returns the reconciled grade.
    Only DOWNGRADES an 'inconsistent'/'partial' to 'N/A (pLDDT-excluded)' when the
    variant is flagged AND the call fires the excluded axis. 'consistent' stays.
    """
    key = (str(gene).lower(), str(variant))
    excl = PLDDT_INTERFACE_EXCLUSIONS.get(key)
    if excl is None:
        return raw_grade
    if raw_grade in ('consistent', 'N/A'):
        return raw_grade  # nothing to excuse
    if _call_fires_excluded_axis(mech_call, excl['axes'], mech_assertions):
        return 'N/A (pLDDT-excluded)'
    return raw_grade


# Snippet for OPTION-2 inline integration (documentation only).
GUARD_SNIPPET = '''
    # --- pLDDT interface-exclusion guard (mech_consistency <-> structural_agreement parity) ---
    _excl = PLDDT_INTERFACE_EXCLUSIONS.get((str(row.get("gene","")).lower(),
              str(row.get("ref_aa",""))+str(row.get("position",""))+str(row.get("alt_aa",""))))
    if _excl is not None:
        _info = MECH_LABEL_ASSERTIONS.get(mech_call)
        if _info is not None:
            _fold_a = _info[0] or _info[1]; _bind_a = _info[2]; _iface_a = _info[3]
            if (("fold" in _excl["axes"] and _fold_a) or
                ("binding" in _excl["axes"] and _bind_a) or
                (_iface_a and not (_fold_a or _bind_a))):
                return "N/A (pLDDT-excluded)", [], []
'''


def _weighted(grades):
    s = [g for g in grades if g in ('consistent', 'partial', 'inconsistent')]
    pts = sum(1.0 if g == 'consistent' else 0.5 if g == 'partial' else 0.0 for g in s)
    return pts, len(s), (pts / len(s) if s else 0.0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--results', default='results/mavis_v7_concordance_v5_out.csv',
                    help='CSV that already has raw mech_consistency_{suf} columns')
    ap.add_argument('--script-dir', default='/mnt/user-data/uploads',
                    help='dir containing apply_concordance_v5.py (for MECH_LABEL_ASSERTIONS)')
    args = ap.parse_args()

    # stub the missing module so we can import the assertions table + rubric
    m = types.ModuleType('mavis_v7'); ms = types.ModuleType('mavis_v7.mechanism')
    ms.classify_mechanism = lambda *a, **k: ('No structural effect detected', '', '')
    m.mechanism = ms; sys.modules['mavis_v7'] = m; sys.modules['mavis_v7.mechanism'] = ms
    sys.path.insert(0, args.script_dir)
    import importlib, pandas as pd
    A = importlib.import_module('apply_concordance_v5')

    df = pd.read_csv(args.results)
    TAGS = ['t10', 't15', 't20', 't25', 'tSAP']

    print('=== mech_consistency: raw vs pLDDT-reconciled ===')
    print(f'{"thr":5s}  {"raw":>22s}   {"reconciled":>22s}')
    for suf in TAGS:
        raw_col = f'mech_consistency_{suf}'
        if raw_col not in df.columns:
            print(f'  {suf}: column {raw_col} absent — run the Track-B grader first'); continue
        recon = []
        for _, r in df.iterrows():
            v = f"{r['ref_aa']}{r['position']}{r['alt_aa']}"
            recon.append(reconcile_grade(r['gene'], v, r[raw_col],
                                         r.get(f'mavis_mechanism_corrected_{suf}'),
                                         A.MECH_LABEL_ASSERTIONS))
        df[f'mech_consistency_{suf}_plddt_reconciled'] = recon
        rp, rn, rf = _weighted(df[raw_col])
        cp, cn, cf = _weighted(pd.Series(recon))
        print(f'  {suf:5s}  {rp:.1f}/{rn}={rf:.3f}        {cp:.1f}/{cn}={cf:.3f}')

    df['mech_consistency_summary_plddt_reconciled'] = df['mech_consistency_t25_plddt_reconciled']
    out = Path(args.results).with_name('mavis_v7_concordance_v5_reconciled.csv')
    df.to_csv(out, index=False)
    print(f'\nwrote {out}')
    print('\nExcused (N/A pLDDT-excluded) at t25:')
    for _, r in df.iterrows():
        if r['mech_consistency_t25_plddt_reconciled'] == 'N/A (pLDDT-excluded)':
            print(f"  {r['gene']} {r['ref_aa']}{r['position']}{r['alt_aa']}  "
                  f"(raw was {r['mech_consistency_t25']}; call={r['mavis_mechanism_corrected_t25']})")


if __name__ == '__main__':
    main()
