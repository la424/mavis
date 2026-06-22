#!/usr/bin/env python3
"""
MAVIS v7 — tier/ddG corroboration column patch.

Adds per-variant structural-evidence columns derived ENTIRELY from columns already present in the
pipeline output. No new computation, no ground-truth dependency. Run as a post-processing stage on
the concordance output (or the corrected results CSV).

WHAT IT ADDS
------------
1. `structural_evidence_strength` (per variant) — a legible relabeling of mavis_tier encoding the
   calibrated gradient (T1 strongest -> T4 minimal). Makes the tier self-documenting for downstream
   (CHD) use.

2. `axis_corroboration_{mono,fold,bind}` (per axis) — the orthogonal-channel agreement flag, applied
   ONLY where the literature ground-truth token is `unknown` (elsewhere it's 'has_ground_truth').
   Compares two INDEPENDENT channels:
     - contact-severity channel: mavis_tier (Grantham severity x contact count; ddG-independent)
     - thermodynamic channel:    FoldX ddG on that axis, vs a 1.0 kcal/mol detectability floor
   Values:
     corroborated_disruption          - both channels fire (elevated confidence of real disruption)
     corroborated_disruption_redundant- both fire BUT evidence is buried_core/partially_buried on the
                                         monomer axis (channels may key off the same embedding fact;
                                         corroboration is weaker / partly redundant)
     corroborated_silent              - both quiet (the 'unknown' is likely a true neutral)
     channels_conflict                - disagree (lowered confidence; often non-thermodynamic mechanism)
     ddg_unevaluable                  - no partner / no ddG available for that axis
     has_ground_truth                 - axis was literature-grounded; corroboration not applicable

DESIGN NOTES
------------
- The redundancy sub-flag (DECIDED - Luke) fires when BOTH channels fire AND the monomer evidence
  token is buried_core or partially_buried: a large Grantham change at a densely-packed buried
  position lights up both channels without being two independent observations of fold destabilization.
- Strength and mechanism are kept SEPARATE (the framework's axis separation): this patch adds an
  evidence-STRENGTH column; the mechanism call lives in mavis_mechanism_*. A downstream user reads
  them together.
- Independence basis: mavis_score = Grantham severity x (mono + max inter contacts); contains NO ddG
  term (confirmed: apply_concordance_v5 / mechanism.py line 71; v13 removed the nbhd ddG term).
"""
import argparse
import pandas as pd

FLOOR = 1.0  # kcal/mol detectability floor for the FoldX channel

STRENGTH_MAP = {
    '1': 'T1: strong structural evidence',
    '2': 'T2: moderate structural evidence',
    '3': 'T3: uncertain structural evidence',
    '4': 'T4: minimal structural evidence',
}

def _tn(t):
    return str(t).split()[-1] if pd.notna(t) else '?'

def _fv(x):
    try: return float(x)
    except: return None

def _tier_fires(row):
    return _tn(row.get('mavis_tier')) in ('1', '2')

def _monomer_is_buried(row):
    ev = str(row.get('mavis_score_evidence', ''))
    return ('buried_core(monomer)' in ev) or ('partially_buried(monomer)' in ev)

def _foldx_axis(row, axis):
    bp = row.get('best_inter_partner')
    bp = bp if (isinstance(bp, str) and bp) else None
    if axis == 'mono':
        return _fv(row.get('ddg_monomer'))
    if axis == 'fold':
        return _fv(row.get(f'ddg_fold_{bp}')) if bp else None
    if axis == 'bind':
        return _fv(row.get(f'ddg_binding_{bp}')) if bp else None
    return None

GT_COL = {'mono': 'expected_ddg_monomer',
          'fold': 'expected_ddg_fold_complex',
          'bind': 'expected_ddg_binding'}

def corroboration_for_axis(row, axis):
    gt = row.get(GT_COL[axis])
    if gt != 'unknown':
        return 'has_ground_truth'
    val = _foldx_axis(row, axis)
    if val is None:
        return 'ddg_unevaluable'
    tier_fires = _tier_fires(row)
    ddg_fires = abs(val) >= FLOOR
    if tier_fires and ddg_fires:
        if axis == 'mono' and _monomer_is_buried(row):
            return 'corroborated_disruption_redundant'
        return 'corroborated_disruption'
    if (not tier_fires) and (not ddg_fires):
        return 'corroborated_silent'
    return 'channels_conflict'

def add_columns(df):
    df['structural_evidence_strength'] = df['mavis_tier'].map(lambda t: STRENGTH_MAP.get(_tn(t), 'unknown_tier'))
    for axis in ('mono', 'fold', 'bind'):
        df[f'axis_corroboration_{axis}'] = df.apply(lambda r: corroboration_for_axis(r, axis), axis=1)
    return df

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--in', dest='infile', default='results/p2_corrected/mavis_v7_concordance.csv')
    ap.add_argument('--out', dest='outfile', default='results/p2_corrected/mavis_v7_concordance_annotated.csv')
    args = ap.parse_args()

    df = pd.read_csv(args.infile)
    df = add_columns(df)
    df.to_csv(args.outfile, index=False)

    print(f'Wrote {args.outfile}')
    print('\n=== structural_evidence_strength distribution ===')
    print(df['structural_evidence_strength'].value_counts().to_string())
    print('\n=== axis_corroboration (unknown axes only) ===')
    for axis in ('mono', 'fold', 'bind'):
        col = f'axis_corroboration_{axis}'
        vc = df[df[col] != 'has_ground_truth'][col].value_counts()
        if len(vc):
            print(f'  [{axis}]')
            for k, v in vc.items():
                print(f'     {k}: {v}')
    print('\n=== corroborated_disruption cases (the strongest "unmeasured but doubly-suggested") ===')
    for axis in ('mono', 'fold', 'bind'):
        col = f'axis_corroboration_{axis}'
        sub = df[df[col].isin(['corroborated_disruption', 'corroborated_disruption_redundant'])]
        for _, r in sub.iterrows():
            tag = 'REDUNDANT' if r[col].endswith('redundant') else 'independent'
            print(f"   {r['gene']} {r['ref_aa']}{r['position']}{r['alt_aa']} [{axis}] role={r['role']} ({tag})")

if __name__ == '__main__':
    main()
