"""
concordance.py — Four-way concordance layer (OPTIONAL module).

Ported faithfully from CHD_pipeline_v6_0_A3-6.ipynb Cell 10 (`compute_concordance_v6`)
into the mavis_v7 package conventions.

DESIGN (Path A, scope-fenced):
  - This is a COLUMN-CONSUMING voting layer. It computes nothing structural.
  - It reads engine outputs (mavis_tier, per-axis ddg + confidence) plus EXTERNAL
    annotations the user supplies as columns (AlphaMissense, franklin).
  - It is OPTIONAL: gated on the external columns being present. If AlphaMissense
    and franklin are both absent, the four-way concordance degrades gracefully to
    the internal structural evidence only (tier + ddG), which is exactly the
    public-tool "structural-only" mode.

This module is what the BENCHMARK does NOT use (benchmark = structural eval only)
and what CHD + the public-tool's optional-external mode DO use.

Faithfulness notes (differences from notebook, made explicit):
  1. The notebook used a single `ddg_confidence` string column ('high'/'low'/...).
     The package uses per-axis `ddg_{partner}_confident` booleans + `ddg_monomer_confident`.
     `derive_ddg_confidence()` below reconstructs a single ddg_confidence string from the
     per-axis flags so the ported voting logic is unchanged. If the input already has a
     `ddg_confidence` column, it is used as-is.
  2. Thresholds (DDG_DESTAB=1.0, DDG_HIGHLY=2.0) imported from config — identical to notebook.
  3. AM/Franklin thresholds and the strict/relaxed/t3 vote construction are copied verbatim.

Usage:
    from .concordance import add_concordance
    df = add_concordance(df, partner_labels)   # adds concordance_* columns
"""

import pandas as pd
from typing import List

from .constants import ss, sf
from .config import DDG_DESTAB, DDG_HIGHLY


# ============================================================================
# Franklin label normalizer (ported verbatim from notebook Cell 9)
# ============================================================================
def std_franklin(v) -> str:
    """Normalize Franklin classification strings to canonical labels."""
    if pd.isna(v):
        return 'No data'
    v = str(v).strip()
    vl = v.lower()
    if 'pathogenic' in vl and 'likely' in vl:
        return 'Likely Pathogenic'
    elif 'pathogenic' in vl:
        return 'Pathogenic'
    elif 'benign' in vl and 'likely' in vl:
        return 'Likely Benign'
    elif 'benign' in vl:
        return 'Benign'
    elif 'vus' in vl and 'high' in vl:
        return 'VUS (high)'
    elif 'vus' in vl and 'mid' in vl:
        return 'VUS (mid)'
    elif 'vus' in vl and 'low' in vl:
        return 'VUS (low)'
    elif 'vus' in vl:
        return 'VUS'
    return v


# ============================================================================
# Reconcile per-axis confidence flags → single ddg_confidence string
# ============================================================================
def derive_ddg_confidence(row, partner_labels: List[str]) -> str:
    """
    The notebook concordance used a single 'ddg_confidence' string ('high'/'low'/'').
    The package carries per-axis booleans. Reconstruct the string so the ported
    voting logic is byte-for-byte unchanged.

    Rule (matches notebook semantics):
      - 'high'  if ANY evaluated axis is confident (monomer or any partner)
      - 'low'   if axes exist but NONE is confident
      - ''      if no ddG axis was evaluated at all
    If the row already has a non-empty 'ddg_confidence', that value is returned.
    """
    existing = ss(row.get('ddg_confidence'))
    if existing:
        return existing.lower()

    any_axis_present = False
    any_axis_confident = False

    if pd.notna(row.get('ddg_monomer')):
        any_axis_present = True
        if row.get('ddg_monomer_confident', False):
            any_axis_confident = True

    for pl in partner_labels:
        for col in (f"ddg_fold_{pl}", f"ddg_binding_{pl}"):
            if col in row.index and pd.notna(row.get(col)):
                any_axis_present = True
                if row.get(f"ddg_{pl}_confident", False):
                    any_axis_confident = True

    if not any_axis_present:
        return ''
    return 'high' if any_axis_confident else 'low'


# ============================================================================
# Confidence-gated max |ddG| across the three axes (ported from notebook)
# ============================================================================
def _gated_max_abs_ddg(row, partner_labels: List[str]) -> float:
    """
    Three-axis confidence-gated max_abs_ddg:
      |ddg_monomer| (if ddg_monomer_confident),
      |ddg_fold_{pl}| and |ddg_binding_{pl}| (if ddg_{pl}_confident).
    Falls back to aggregate ddg_multimer_{max,min} only when NO per-partner
    binding/fold columns exist. Identical logic to notebook Cell 10.
    """
    ddg_vals = []
    if pd.notna(row.get('ddg_monomer')) and row.get('ddg_monomer_confident', False):
        ddg_vals.append(abs(float(row['ddg_monomer'])))

    has_pp = False
    for pl in partner_labels:
        bind_col = f"ddg_binding_{pl}"
        fold_col = f"ddg_fold_{pl}"
        conf = row.get(f"ddg_{pl}_confident", False)
        if bind_col in row.index and pd.notna(row.get(bind_col)):
            has_pp = True
            if conf:
                ddg_vals.append(abs(float(row[bind_col])))
        if fold_col in row.index and pd.notna(row.get(fold_col)):
            has_pp = True
            if conf:
                ddg_vals.append(abs(float(row[fold_col])))

    if not has_pp:
        for col in ('ddg_multimer_max', 'ddg_multimer_min'):
            v = row.get(col)
            if pd.notna(v):
                ddg_vals.append(abs(float(v)))

    return max(ddg_vals) if ddg_vals else 0.0


# ============================================================================
# Four-way concordance (ported verbatim from notebook compute_concordance_v6)
# ============================================================================
def compute_concordance(row, partner_labels: List[str]) -> pd.Series:
    """
    Four-way concordance: strict, relaxed, T3-inclusive, with adjusted denominators.
    Votes over: structure tier, ddG, AlphaMissense, Franklin — each gated on its
    own evaluability. Denominator = number of evaluable axes (min 1).
    """
    tier = ss(row.get('mavis_tier'))                  # package col (was v6_tier)
    am = ss(row.get('AlphaMissense')).lower()
    fr = std_franklin(row.get('franklin'))
    fr_lower = fr.lower() if isinstance(fr, str) else ''
    ddg_conf = derive_ddg_confidence(row, partner_labels)

    max_abs_ddg = _gated_max_abs_ddg(row, partner_labels)

    struct_eval = row.get('structure_evaluable', False)
    ddg_eval = row.get('ddg_evaluable', False)
    am_eval = row.get('am_evaluable', False)
    fr_eval = row.get('franklin_evaluable', False)

    tier_t12 = tier in ('Tier 1', 'Tier 2')
    tier_t123 = tier_t12 or tier == 'Tier 3'

    ddg_strict = 1 if (ddg_conf == 'high' and max_abs_ddg >= DDG_HIGHLY) else 0
    ddg_relaxed = 1 if (ddg_conf not in ('low', '') and max_abs_ddg >= DDG_DESTAB) else 0
    am_strict = 1 if am == 'likely_pathogenic' else 0
    am_relaxed = 1 if am in ('likely_pathogenic', 'ambiguous') else 0
    fr_strict = 1 if fr_lower in ('pathogenic', 'likely pathogenic', 'vus (high)') else 0
    fr_relaxed = 1 if fr_lower in ('pathogenic', 'likely pathogenic', 'vus (high)', 'vus (mid)') else 0

    def build(tier_v, ddg_v, am_v, fr_v):
        s, d = 0, 0
        if struct_eval: s += tier_v; d += 1
        if ddg_eval:    s += ddg_v;  d += 1
        if am_eval:     s += am_v;   d += 1
        if fr_eval:     s += fr_v;   d += 1
        return s, max(d, 1)

    s_s, s_d = build(1 if tier_t12 else 0, ddg_strict, am_strict, fr_strict)
    r_s, r_d = build(1 if tier_t12 else 0, ddg_relaxed, am_relaxed, fr_relaxed)
    t3_s, t3_d = build(1 if tier_t123 else 0, ddg_strict, am_strict, fr_strict)

    return pd.Series({
        'four_way_strict': s_s, 'four_way_strict_denom': s_d,
        'concordance_strict': f"{s_s}/{s_d}",
        'four_way_relaxed': r_s, 'four_way_relaxed_denom': r_d,
        'concordance_relaxed': f"{r_s}/{r_d}",
        'four_way_t3': t3_s, 'four_way_t3_denom': t3_d,
        'concordance_t3': f"{t3_s}/{t3_d}",
        'structure_vote_strict': 1 if tier_t12 else 0,
        'structure_vote_t3': 1 if tier_t123 else 0,
        'ddg_vote_strict': ddg_strict,
        'ddg_vote_relaxed': ddg_relaxed,
        'am_vote_strict': am_strict,
        'am_vote_relaxed': am_relaxed,
        'franklin_vote_strict': fr_strict,
        'franklin_vote_relaxed': fr_relaxed,
        'max_abs_ddg_gated': round(max_abs_ddg, 4),
        'franklin_std': fr,
    })


# ============================================================================
# Public entry point
# ============================================================================
def add_concordance(df: pd.DataFrame, partner_labels: List[str]) -> pd.DataFrame:
    """
    Add four-way concordance columns to df. OPTIONAL layer.

    Evaluability columns are set here if absent:
      - am_evaluable       = AlphaMissense column present and non-null per row
      - franklin_evaluable = franklin column present and non-null per row
    (structure_evaluable / ddg_evaluable come from mechanism.compute_evaluability.)

    If neither AlphaMissense nor franklin columns exist, the external votes are
    simply never evaluable and concordance reduces to the internal structural
    evidence (tier + ddG) — the public-tool structural-only behavior. No error.
    """
    df = df.copy()

    if 'AlphaMissense' in df.columns:
        df['am_evaluable'] = df['AlphaMissense'].notna()
    else:
        df['am_evaluable'] = False

    if 'franklin' in df.columns:
        df['franklin_evaluable'] = df['franklin'].notna()
    else:
        df['franklin_evaluable'] = False

    conc = df.apply(lambda r: compute_concordance(r, partner_labels), axis=1)
    for c in conc.columns:
        df[c] = conc[c]
    return df
