#!/usr/bin/env python3
"""
MAVIS v7 — Concordance post-processor v4.

Layered on top of v3. New in v4:

  (1) Derived `expected_mech_class` column using the v2 grading rubric —
      structural-expectation classification derived from the four expected_*
      axis columns plus evidence_axes. Role does NOT determine mechanism
      class. Values:
        structurally_uncommitted       — all axes not-tested
        interface_uncommitted_magnitude — topology positive, DDG axes silent
        structurally_silent             — all committed axes negative
        fold_mechanism                  — fold positive, binding neg/not_tested
        ppi_destab_mechanism / ppi_stab_mechanism
        mixed_structural                — both fold and binding positive

  (2) Mechanism-consistency grading per rubric v2 decision tree. Columns:
        mech_consistency_t10 / _t15 / _t20           (Pipeline 1)
        nbhd_mech_consistency_t10 / _t15 / _t20      (Pipeline 2)
        mech_consistency_summary                     (conservative across thresholds)
        nbhd_mech_consistency_summary
        mech_false_positive_axes_t10 / _t15 / _t20   (which tested-negative axes fired)
        mech_missed_positive_axes_t10 / _t15 / _t20  (which tested-positive axes missed)
        (and nbhd_ versions)
      Grades: consistent / partial / inconsistent / N/A.

  (3) Structural-signal / external-consensus split of concordance:
        structural_signal_strict / relaxed    (tier + DDG only; per pipeline)
        external_consensus_strict / relaxed   (AM + Franklin; pipeline-independent)
        nbhd_structural_signal_strict / relaxed
      Preserves existing concordance_* columns; additions let downstream
      readers distinguish STRUCTURAL-DISRUPTION detection from
      PATHOGENICITY-PREDICTION evidence.

  (4) Evidence-aware mild_destab handling (Approach 3):
      - mild_destab + structural evidence (binding/monomer/both) → positive axis
      - mild_destab + soft evidence (functional/population)     → not-tested
      This reflects annotator confidence: structural evidence means the
      mild_destab annotation came from measured ΔΔG; functional evidence
      means it was an inferred guess.

  (5) Per-variant 95% confidence intervals on each DDG axis. For each of
      ddg_monomer, ddg_fold_{partner}, ddg_binding_{partner}, two CI flavors:
        ci95_internal_low/high — using FoldX's 5-run SD * 1.96. Tight,
                                  captures run-to-run reproducibility only.
        ci95_sapozhnikov_low/high — using global ±2.9 (fold) / ±3.5 (binding)
                                  prediction intervals from Sapozhnikov et al.
                                  2023 (BMC Bioinformatics). Realistic,
                                  captures structure quality and conformational
                                  uncertainty.
      For each CI, distinguishable_*_from_{0,t10,t15,t20,t25} flags indicate
      whether the CI excludes the corresponding threshold magnitude.

  (6) Four-threshold mechanism sweep: t=1.0, 1.5, 2.0, 2.5 kcal/mol. Reflects
      the full range of thresholds used in the literature: 1.0 (CAGI5
      frataxin convention, prior CHD continuity), 1.5 (Caldararu/Guerois
      classification optimum), 2.0 (intermediate), 2.5 (Sapozhnikov
      confident-detection threshold). No single threshold designated as
      "primary" — the rubric reports all four side-by-side.

  (7) Role-independent grading. The `role` field does NOT enter
      mechanism-consistency grading or mech_class derivation. Every variant
      is graded against its structural annotations uniformly. Pathogenicity
      is evaluated separately via the external_consensus track.

  (8) Mode detection. If the upstream benchmark CSV lacks role or the
      expected_* columns, the mechanism-consistency and mech_class
      computations skip gracefully ("Mode B" — structural-signal only).
      CI columns are emitted in both modes since they are pipeline-internal.

Vote rules, evaluability, Pipeline 2 scoring, per-threshold mechanism
classification, and pipeline_agreement are unchanged from v3.

USAGE
-----
    cd ~/mavis_v7
    python3 apply_concordance_v4.py

Reads:  results/mavis_v7_results_with_nbhd.csv
        AM_variants_mavis_mechanism_test.xlsx
Writes: results/mavis_v7_concordance.csv
"""
from __future__ import annotations
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd


# ============================================================================
# Configuration (unchanged from v3)
# ============================================================================
DDG_DESTAB = 1.0
DDG_HIGHLY = 2.0
AM_STRICT  = {"likely_pathogenic"}
AM_RELAXED = {"likely_pathogenic", "ambiguous"}
FRANKLIN_STRICT  = {"pathogenic", "likely_pathogenic", "vus_high"}
FRANKLIN_RELAXED = {"pathogenic", "likely_pathogenic", "vus_high", "vus_mid"}
FOOTPRINT_TIERS = {"Tier 1", "Tier 2"}
LOW_PLDDT_CUTOFF = 70.0
HIGH_SD_CUTOFF = 1.0
DDG_SWEEP_THRESHOLDS = (1.0, 1.5, 2.0, 2.5)
# v5: Sapozhnikov-confident per-axis threshold (mono=2.9, fold=2.9, bind=3.5)
# Used as the 5th threshold in mechanism and structural_agreement evaluation.
SAPOZHNIKOV_PER_AXIS = {'monomer': 2.9, 'fold': 2.9, 'binding': 3.5}
# Tag used in column suffixes for the Sapozhnikov per-axis threshold.
SAPOZHNIKOV_TAG = "tSAP"
# Full list of (tag, threshold_or_dict) pairs for the v5 sweep
THRESHOLD_SPECS = [
    ('t10', 1.0),
    ('t15', 1.5),
    ('t20', 2.0),
    ('t25', 2.5),
    (SAPOZHNIKOV_TAG, SAPOZHNIKOV_PER_AXIS),
]

# ----------------------------------------------------------------------------
# CI uncertainty model parameters (Sapozhnikov et al. 2023 BMC Bioinformatics)
# ----------------------------------------------------------------------------
# Sapozhnikov 2023 reports global 95% prediction intervals of ±2.9 kcal/mol for
# folding stability and ±3.5 kcal/mol for binding stability when comparing
# FoldX predictions to experimental ΔΔG. These are realistic upper bounds on
# the disagreement between FoldX and experiment, accounting for structure
# quality, conformational variability, and biochemical factors.
SAPOZHNIKOV_FOLD_CI95 = 2.9
SAPOZHNIKOV_BIND_CI95 = 3.5

# Internal FoldX run-to-run uncertainty uses Z=1.96 (95% CI two-sided) on the
# 5-run SD reported by BuildModel/AnalyseComplex. This is OPTIMISTIC — it only
# captures run-to-run reproducibility, not the structure-quality and
# conformational uncertainty that Sapozhnikov bounds capture.
INTERNAL_CI_Z = 1.96

NBHD_TIER_THRESHOLDS = [(5.0, "Tier 1"), (3.0, "Tier 2"), (1.5, "Tier 3")]
DEFAULT_TIER = "Tier 4"

DISRUPTION_POINTS = [(20, 4.0), (10, 3.0), (4, 2.0), (1, 1.0)]
CONTACT_DRIVEN_THRESHOLD = 6
BURIAL_RANK = {"unknown": 0, "surface_exposed": 1, "partially_buried": 2, "buried_core": 3}
RANK_TO_BURIAL = {v: k for k, v in BURIAL_RANK.items()}

# ----------------------------------------------------------------------------
# Mechanism-consistency rubric (v4)
# ----------------------------------------------------------------------------
# Columns required to run mech-consistency and derive expected_mech_class.
# If any are missing, the script falls back to "Mode B" (structural-only).
GROUND_TRUTH_COLS = [
    "role",
    "expected_ddg_monomer",
    "expected_ddg_fold_complex",
    "expected_ddg_binding",
    "expected_topology",
]

# Axis-status classification per expected column value
_POSITIVE_AXIS_TOKENS = {"destab", "mild_destab", "stab"}
_NEGATIVE_AXIS_TOKENS = {"neutral"}
_NOTTEST_AXIS_TOKENS  = {"unknown", ""}

_AT_INTERFACE_TOKENS = {"at_interface"}
_AWAY_INTERFACE_TOKENS = {"away_from_interface"}

# MAVIS mechanism labels → structural axes that label asserts POSITIVE.
# Each label is a dict of booleans for (fold_monomer_affected,
# fold_complex_affected, binding_affected, interface_asserted) and a direction
# tag. Derived from classify_mechanism_at's label space.
MECH_LABEL_ASSERTIONS = {
    # Label                                                                  fold_m, fold_c, bind, iface, dir
    # ---- Pure fold destab (mono/multi/both, with/without interface) ----
    "Monomer fold destabilization":                                        (True,  False, False, False, "destab"),
    "Monomer fold destabilization at interface":                           (True,  False, False, True,  "destab"),
    "Multimer fold destabilization":                                       (False, True,  False, False, "destab"),
    "Multimer fold destabilization at interface":                          (False, True,  False, True,  "destab"),
    "Both fold destabilization":                                           (True,  True,  False, False, "destab"),
    "Both fold destabilization at interface":                              (True,  True,  False, True,  "destab"),
    # ---- Pure fold stab ----
    "Monomer fold stabilization":                                          (True,  False, False, False, "stab"),
    "Monomer fold stabilization at interface":                             (True,  False, False, True,  "stab"),
    "Multimer fold stabilization":                                         (False, True,  False, False, "stab"),
    "Multimer fold stabilization at interface":                            (False, True,  False, True,  "stab"),
    "Both fold stabilization":                                             (True,  True,  False, False, "stab"),
    "Both fold stabilization at interface":                                (True,  True,  False, True,  "stab"),
    # ---- Fold + PPI same-direction destab ----
    "Monomer fold + PPI destabilization":                                  (True,  False, True,  False, "destab"),
    "Multimer fold + PPI destabilization":                                 (False, True,  True,  False, "destab"),
    "Both fold + PPI destabilization":                                     (True,  True,  True,  False, "destab"),
    # ---- Fold + PPI same-direction stab ----
    "Monomer fold + PPI stabilization":                                    (True,  False, True,  False, "stab"),
    "Multimer fold + PPI stabilization":                                   (False, True,  True,  False, "stab"),
    "Both fold + PPI stabilization":                                       (True,  True,  True,  False, "stab"),
    # ---- Fold destab + PPI stab (conflicting) ----
    "Monomer fold destab. + PPI stabilization (conflicting)":              (True,  False, True,  False, "conflicting"),
    "Multimer fold destab. + PPI stabilization (conflicting)":             (False, True,  True,  False, "conflicting"),
    "Both fold destab. + PPI stabilization (conflicting)":                 (True,  True,  True,  False, "conflicting"),
    # ---- Fold stab + PPI destab (conflicting) ----
    "Monomer fold stab. + PPI destabilization (conflicting)":              (True,  False, True,  False, "conflicting"),
    "Multimer fold stab. + PPI destabilization (conflicting)":             (False, True,  True,  False, "conflicting"),
    "Both fold stab. + PPI destabilization (conflicting)":                 (True,  True,  True,  False, "conflicting"),
    # ---- Pure PPI ----
    "PPI destabilization":                                                 (False, False, True,  False, "destab"),
    "PPI stabilization":                                                   (False, False, True,  False, "stab"),
    "PPI conflicting (mixed partner signals)":                             (False, False, True,  False, "conflicting"),
    # ---- Structural without DDG ----
    "Interface variant (DDG neutral)":                                     (False, False, False, True,  "neutral"),
    "Structural variant — contact-driven (DDG neutral)":                   (False, False, False, False, "neutral"),
    "Structural variant — burial-driven (DDG neutral)":                    (False, False, False, False, "neutral"),
    # ---- Catch-all ----
    "No structural effect detected":                                       (False, False, False, False, "none"),
    "Structure unevaluable":                                               (False, False, False, False, "unevaluable"),
    "Structure low-confidence at variant site":                            (False, False, False, False, "unevaluable"),
}


# ============================================================================
# Helpers (unchanged)
# ============================================================================
def ss(v):
    if v is None or (isinstance(v, float) and np.isnan(v)): return ""
    if pd.isna(v): return ""
    return str(v)


def sf(v, default=0.0):
    if v is None: return default
    if pd.isna(v): return default
    try: return float(v)
    except (ValueError, TypeError): return default


def discover_partners(df):
    return sorted({c.replace("ddg_binding_", "")
                   for c in df.columns
                   if c.startswith("ddg_binding_")
                   and not c.endswith("_sd")
                   and not c.endswith("_indistinguishable")
                   and not c.endswith("_confident")})


def normalize_franklin(v):
    if pd.isna(v): return ""
    s = str(v).strip().lower()
    canon = {
        "benign": "benign", "likely benign": "likely_benign", "likely_benign": "likely_benign",
        "vus (low)": "vus_low", "vus_low": "vus_low", "vus low": "vus_low",
        "vus (mid)": "vus_mid", "vus_mid": "vus_mid", "vus mid": "vus_mid",
        "vus (high)": "vus_high", "vus_high": "vus_high", "vus high": "vus_high",
        "pathogenic": "pathogenic",
        "likely pathogenic": "likely_pathogenic", "likely_pathogenic": "likely_pathogenic",
    }
    return canon.get(s, s)


def normalize_am_class(v):
    if pd.isna(v): return ""
    return str(v).strip().lower()


def derive_ddg_confidence(row, partners):
    site_status = str(row.get("site_plddt_status", "")).lower()
    if site_status == "crystal":
        return "high"
    if not bool(row.get("structure_evaluable", True)):
        return "low"
    mono_plddt = row.get("monomer_plddt")
    mono_high = pd.notna(mono_plddt) and float(mono_plddt) >= 70
    mono_med  = pd.notna(mono_plddt) and float(mono_plddt) >= 50
    partner_plddts = [row.get(f"multi_{p}_plddt") for p in partners]
    partner_present = any(pd.notna(x) for x in partner_plddts)
    partner_high = any(pd.notna(x) and float(x) >= 70 for x in partner_plddts)
    partner_med  = any(pd.notna(x) and float(x) >= 50 for x in partner_plddts)
    if mono_high and partner_high: return "high"
    if mono_high and not partner_present: return "high"
    if mono_med or partner_med: return "medium"
    return "low"


def compute_max_abs_ddg(row, partners):
    vals = []
    if pd.notna(row.get("ddg_monomer")) and bool(row.get("ddg_monomer_confident", False)):
        vals.append(abs(float(row["ddg_monomer"])))
    for p in partners:
        if not bool(row.get(f"ddg_{p}_confident", False)):
            continue
        fv = row.get(f"ddg_fold_{p}")
        if pd.notna(fv):
            vals.append(abs(float(fv)))
        bv = row.get(f"ddg_binding_{p}")
        if pd.notna(bv) and not bool(row.get(f"ddg_binding_{p}_indistinguishable", False)):
            vals.append(abs(float(bv)))
    return max(vals) if vals else 0.0


# ============================================================================
# NEW in v4: per-variant confidence intervals on each DDG axis
# ============================================================================
def compute_ddg_cis(row, partners):
    """
    Compute two flavors of 95% CI on each DDG axis:

    Internal CI: ddg ± 1.96 * SD, where SD is the FoldX run-to-run standard
      deviation across the n=5 BuildModel replicates. Captures run-to-run
      reproducibility only — does NOT capture structure-quality uncertainty.
      Tight, optimistic.

    Sapozhnikov CI: ddg ± 2.9 (fold) or ± 3.5 (binding) kcal/mol, the global
      95% prediction intervals from Sapozhnikov et al. 2023 (BMC Bioinformatics)
      for FoldX-vs-experimental ΔΔG agreement. Captures structure quality,
      conformational variability, and biochemical factors. Wide, realistic.
      Independent of per-variant SD.

    For each CI, also computes distinguishable_from_threshold flags at each
    threshold in DDG_SWEEP_THRESHOLDS. A variant is "distinguishable from
    threshold T" if the absolute value of the DDG, minus the CI half-width,
    exceeds T — meaning the CI excludes the [-T, +T] range.

    Returns a dict of new column values to merge into the row.
    """
    out = {}

    # Helper to compute CI bounds and flags for one DDG axis
    def _ci_for_axis(ddg_val, ddg_sd, axis_kind, prefix):
        """
        ddg_val   — point estimate (mean across runs)
        ddg_sd    — run-to-run SD (None if missing)
        axis_kind — 'fold' (use ±2.9) or 'binding' (use ±3.5)
        prefix    — column prefix, e.g. 'ddg_monomer' or 'ddg_fold_brca1'
        """
        result = {}
        sapoz_w = SAPOZHNIKOV_FOLD_CI95 if axis_kind == "fold" else SAPOZHNIKOV_BIND_CI95

        if pd.isna(ddg_val) or ddg_val is None:
            # No DDG value → all CI columns empty
            result[f"{prefix}_ci95_internal_low"]    = None
            result[f"{prefix}_ci95_internal_high"]   = None
            result[f"{prefix}_ci95_sapozhnikov_low"] = None
            result[f"{prefix}_ci95_sapozhnikov_high"]= None
            for thr in (0.0,) + DDG_SWEEP_THRESHOLDS:
                tag = "0" if thr == 0.0 else f"t{int(thr*10)}"
                result[f"{prefix}_distinguishable_internal_from_{tag}"]    = None
                result[f"{prefix}_distinguishable_sapozhnikov_from_{tag}"] = None
            return result

        ddg_val = float(ddg_val)
        abs_ddg = abs(ddg_val)

        # Internal CI: needs SD. If SD missing or zero (n=1 case), CI = point estimate
        if pd.notna(ddg_sd) and ddg_sd is not None:
            half = INTERNAL_CI_Z * float(ddg_sd)
            internal_low  = round(ddg_val - half, 4)
            internal_high = round(ddg_val + half, 4)
        else:
            internal_low  = internal_high = round(ddg_val, 4)
            half = 0.0
        result[f"{prefix}_ci95_internal_low"]  = internal_low
        result[f"{prefix}_ci95_internal_high"] = internal_high

        # Sapozhnikov CI: fixed half-width per axis kind
        result[f"{prefix}_ci95_sapozhnikov_low"]  = round(ddg_val - sapoz_w, 4)
        result[f"{prefix}_ci95_sapozhnikov_high"] = round(ddg_val + sapoz_w, 4)

        # Distinguishable flags — at threshold T, the CI excludes [-T, +T] iff
        # |ddg| - half_width > T (the CI's near edge is beyond the threshold).
        # We tag the threshold with its scaled int form: 0, t10, t15, t20, t25.
        for thr in (0.0,) + DDG_SWEEP_THRESHOLDS:
            tag = "0" if thr == 0.0 else f"t{int(thr*10)}"
            result[f"{prefix}_distinguishable_internal_from_{tag}"]    = bool((abs_ddg - half)    > thr)
            result[f"{prefix}_distinguishable_sapozhnikov_from_{tag}"] = bool((abs_ddg - sapoz_w) > thr)
        return result

    # Monomer axis (fold kind)
    out.update(_ci_for_axis(
        row.get("ddg_monomer"), row.get("ddg_monomer_sd"),
        axis_kind="fold", prefix="ddg_monomer",
    ))

    # Per-partner: fold-in-complex and binding axes
    for p in partners:
        out.update(_ci_for_axis(
            row.get(f"ddg_fold_{p}"), row.get(f"ddg_fold_{p}_sd"),
            axis_kind="fold", prefix=f"ddg_fold_{p}",
        ))
        out.update(_ci_for_axis(
            row.get(f"ddg_binding_{p}"), row.get(f"ddg_binding_{p}_sd"),
            axis_kind="binding", prefix=f"ddg_binding_{p}",
        ))

    return out


def low_structural_confidence(row, sd_cols):
    if str(row.get("site_plddt_status", "")).lower() == "crystal":
        return False
    plddt = row.get("best_plddt")
    if pd.notna(plddt) and float(plddt) < LOW_PLDDT_CUTOFF:
        return True
    for c in sd_cols:
        v = row.get(c)
        if pd.notna(v) and float(v) > HIGH_SD_CUTOFF:
            return True
    return False


# ============================================================================
# Pipeline 2 scoring (unchanged from v3)
# ============================================================================
def assign_tier(score):
    if pd.isna(score): return "Tier 4 (unevaluable)"
    for threshold, tier in NBHD_TIER_THRESHOLDS:
        if score >= threshold: return tier
    return DEFAULT_TIER


def compute_nbhd_score(row, partners):
    # (unchanged from v3; omitted here for brevity — identical to v3 line-for-line)
    # ---- BEGIN v3-IDENTICAL BLOCK ----
    mono_nbhd_w = row.get("nbhd_mono_contacts_weighted")
    has_mono_nbhd = pd.notna(mono_nbhd_w)
    mono_eval = bool(row.get("nbhd_mono_evaluable", False))

    best_inter_w = 0.0
    has_any_partner_nbhd = False
    best_partner = None
    for p in partners:
        inter_w = row.get(f"multi_{p}_nbhd_inter_weighted")
        if pd.notna(inter_w):
            has_any_partner_nbhd = True
            if float(inter_w) > best_inter_w:
                best_inter_w = float(inter_w)
                best_partner = p

    if not has_mono_nbhd and not has_any_partner_nbhd:
        return {
            "nbhd_mono_score": None,
            "nbhd_inter_score": None,
            "nbhd_ddg_score": None,
            "nbhd_context_score": None,
            "nbhd_final_score": None,
            "nbhd_best_partner": None,
            "nbhd_evaluable": False,
        }

    mono_score = 0.0
    if has_mono_nbhd and mono_eval:
        mw = float(mono_nbhd_w)
        for thr, pts in DISRUPTION_POINTS:
            if mw >= thr:
                mono_score = pts
                break

    inter_score = 0.0
    if best_partner is not None:
        for thr, pts in DISRUPTION_POINTS:
            if best_inter_w >= thr:
                inter_score = pts
                break

    # v13 P2 BUGFIX: removed spurious ddg_score term from neighborhood score.
    # Pipeline 2 (neighborhood) scoring should reflect neighborhood-level
    # structural disruption (mono + inter + context), NOT the variant's own
    # DDG which is already represented in Pipeline 1. The previous formula
    # double-counted the DDG signal and inflated nbhd_final_score for any
    # variant with high |DDG|. Now the neighborhood score is a pure
    # structural-context signal.

    context_score = 0.0
    burial = ss(row.get("monomer_burial")).lower()
    if burial == "buried_core":
        context_score += 0.5
    any_iface_nbhd = False
    for p in partners:
        if bool(row.get(f"multi_{p}_nbhd_has_interface", False)):
            any_iface_nbhd = True
            break
    if any_iface_nbhd:
        context_score += 0.5

    final = mono_score + inter_score + context_score
    return {
        "nbhd_mono_score": round(mono_score, 2),
        "nbhd_inter_score": round(inter_score, 2),
        "nbhd_ddg_score": 0.0,  # v13: no longer contributes; kept for backward column compat
        "nbhd_context_score": round(context_score, 2),
        "nbhd_final_score": round(final, 2),
        "nbhd_best_partner": best_partner,
        "nbhd_evaluable": True,
    }


# ============================================================================
# Mechanism classification (v5: thin wrapper around modular classify_mechanism)
# ============================================================================
# Add the modular package to the path so we can import the patched classifier
# with per-axis threshold dispatch and 33-category fold split.
_THIS_DIR = Path(__file__).resolve().parent
if str(_THIS_DIR) not in sys.path:
    sys.path.insert(0, str(_THIS_DIR))
from mavis_v7.mechanism import classify_mechanism as _modular_classify_mechanism


def classify_mechanism_at(row, partners, threshold, tier_col="mavis_tier"):
    """
    Classify mechanism at a given threshold (scalar or per-axis dict).

    Wraps the modular classify_mechanism with two adaptations:
      1. Masks binding ΔΔG to 0 for partners flagged as indistinguishable
         from the per-system baseline (preserving v4 behavior in this script).
      2. Substitutes the requested tier column (mavis_tier vs nbhd_tier) so
         that contact-driven / burial-driven decisions reflect the correct
         pipeline's tier assignment.

    Returns just the mechanism label (the modular classifier returns a tuple;
    we discard partner and external_evidence_flag here for column-shape
    compatibility with the rest of this script).
    """
    # Build a row copy with binding masked for indistinguishable partners
    row_view = row.copy()
    for pl in partners:
        flag_col = f"ddg_binding_{pl}_indistinguishable"
        bd_col = f"ddg_binding_{pl}"
        if flag_col in row_view.index and bool(row_view[flag_col]):
            if bd_col in row_view.index and pd.notna(row_view[bd_col]):
                row_view[bd_col] = 0.0

    # Substitute the requested tier (modular classifier reads 'mavis_tier' by default)
    if tier_col != "mavis_tier":
        row_view = row_view.copy()
        row_view["mavis_tier"] = row_view.get(tier_col)

    mech, _partner, _ext = _modular_classify_mechanism(row_view, partners, ddg_destab=threshold)
    return mech


# ============================================================================
# Concordance voting (unchanged from v3)
# ============================================================================
def footprint_vote(tier):
    return 1 if str(tier) in FOOTPRINT_TIERS else 0


def ddg_votes(row):
    conf = str(row.get("ddg_confidence_derived", ""))
    mxd = float(row.get("max_abs_ddg", 0.0))
    strict = 1 if (conf == "high" and mxd >= DDG_HIGHLY) else 0
    relaxed = 1 if (conf in ("high", "medium") and mxd >= DDG_DESTAB) else 0
    return strict, relaxed


def am_votes(row):
    cls = normalize_am_class(row.get("AM class"))
    return (1 if cls in AM_STRICT else 0, 1 if cls in AM_RELAXED else 0)


def franklin_votes(row):
    fr = normalize_franklin(row.get("franklin"))
    return (1 if fr in FRANKLIN_STRICT else 0, 1 if fr in FRANKLIN_RELAXED else 0)


def axis_evaluable(row):
    return {
        "struct":   bool(row.get("structure_evaluable", True)),
        "ddg":      str(row.get("ddg_confidence_derived", "")) in ("high", "medium"),
        "am":       pd.notna(row.get("AM pathogenicity")),
        "franklin": pd.notna(row.get("franklin")) and str(row.get("franklin", "")).strip() != "",
    }


def compute_concordance(row, tier_val, prefix="", include_external=True):
    fp_v = footprint_vote(tier_val)
    ddg_s, ddg_r = ddg_votes(row)
    am_s, am_r = am_votes(row)
    fr_s, fr_r = franklin_votes(row)
    ev = axis_evaluable(row)
    if not include_external:
        ev = {**ev, "am": False, "franklin": False}

    def assemble(tv, dv, av, fv):
        s, d = 0, 0
        if ev["struct"]:   s += int(bool(tv)); d += 1
        if ev["ddg"]:      s += dv; d += 1
        if ev["am"]:       s += av; d += 1
        if ev["franklin"]: s += fv; d += 1
        return s, d

    strict_n,  strict_d  = assemble(fp_v, ddg_s, am_s, fr_s)
    relaxed_n, relaxed_d = assemble(fp_v, ddg_r, am_r, fr_r)
    suffix = "full" if include_external else "struct"
    return {
        f"{prefix}concordance_strict_{suffix}":  f"{strict_n}/{strict_d}" if strict_d else "NA",
        f"{prefix}concordance_relaxed_{suffix}": f"{relaxed_n}/{relaxed_d}" if relaxed_d else "NA",
        f"{prefix}concordance_strict_{suffix}_n":     strict_n,
        f"{prefix}concordance_strict_{suffix}_denom": strict_d,
        f"{prefix}concordance_relaxed_{suffix}_n":    relaxed_n,
        f"{prefix}concordance_relaxed_{suffix}_denom":relaxed_d,
    }


# ============================================================================
# NEW in v4: structural-signal / external-consensus split
# ============================================================================
def compute_signal_consensus_split(row, tier_val, prefix=""):
    """
    Split the existing 4-axis concordance into two orthogonal tracks:
      - structural_signal: tier-footprint + DDG magnitude. Measures PIPELINE
        detection of structural disruption. Independent of pathogenicity.
      - external_consensus: AlphaMissense + Franklin. Measures whether EXTERNAL
        predictors agree on pathogenicity.

    Denominators are per-row based on axis evaluability, so unevaluable axes
    don't count. By construction: concordance_full_n = signal_n + consensus_n
    (within matching strict/relaxed flavor). Preserved for reconstructibility.
    """
    fp_v = footprint_vote(tier_val)
    ddg_s, ddg_r = ddg_votes(row)
    am_s, am_r = am_votes(row)
    fr_s, fr_r = franklin_votes(row)
    ev = axis_evaluable(row)

    # Structural signal = struct tier + DDG axes only
    sig_strict_n = (int(bool(fp_v)) if ev["struct"] else 0) + (ddg_s if ev["ddg"] else 0)
    sig_strict_d = (1 if ev["struct"] else 0) + (1 if ev["ddg"] else 0)
    sig_relax_n  = (int(bool(fp_v)) if ev["struct"] else 0) + (ddg_r if ev["ddg"] else 0)
    sig_relax_d  = sig_strict_d

    # External consensus = AM + Franklin only
    con_strict_n = (am_s if ev["am"] else 0) + (fr_s if ev["franklin"] else 0)
    con_strict_d = (1 if ev["am"] else 0) + (1 if ev["franklin"] else 0)
    con_relax_n  = (am_r if ev["am"] else 0) + (fr_r if ev["franklin"] else 0)
    con_relax_d  = con_strict_d

    return {
        f"{prefix}structural_signal_strict":  f"{sig_strict_n}/{sig_strict_d}" if sig_strict_d else "NA",
        f"{prefix}structural_signal_relaxed": f"{sig_relax_n}/{sig_relax_d}"  if sig_relax_d  else "NA",
        f"{prefix}structural_signal_strict_n":     sig_strict_n,
        f"{prefix}structural_signal_strict_denom": sig_strict_d,
        f"{prefix}structural_signal_relaxed_n":    sig_relax_n,
        f"{prefix}structural_signal_relaxed_denom":sig_relax_d,
        # External consensus is not prefixed — it's pipeline-independent and
        # identical for P1 and P2. We only emit it once (prefix="").
        f"external_consensus_strict":  f"{con_strict_n}/{con_strict_d}" if con_strict_d else "NA",
        f"external_consensus_relaxed": f"{con_relax_n}/{con_relax_d}"  if con_relax_d  else "NA",
        f"external_consensus_strict_n":     con_strict_n,
        f"external_consensus_strict_denom": con_strict_d,
        f"external_consensus_relaxed_n":    con_relax_n,
        f"external_consensus_relaxed_denom":con_relax_d,
    }


# ============================================================================
# NEW in v4: axis-status classification for mechanism-consistency grading
# ============================================================================
# Evidence-aware rule: mild_destab is graded positive only when evidence_axes
# is structural (binding/monomer/both); otherwise treated as not-tested.
STRUCTURAL_EVIDENCE = {"binding", "monomer", "both"}
SOFT_EVIDENCE       = {"functional", "population"}


def _tokenize_axis(val, evidence=""):
    """Classify an expected_ddg_* value into positive / negative / not_tested.
    Evidence-aware: mild_destab → positive only under structural evidence."""
    s = ss(val).strip().lower()
    e = ss(evidence).strip().lower()
    if s == "mild_destab":
        return "positive" if e in STRUCTURAL_EVIDENCE else "not_tested"
    if s in _POSITIVE_AXIS_TOKENS: return "positive"
    if s in _NEGATIVE_AXIS_TOKENS: return "negative"
    return "not_tested"  # unknown, empty, or any unrecognized value


def _tokenize_topology(val):
    """Classify expected_topology as positive (at_interface), negative (away_from_interface), or not_tested."""
    s = ss(val).strip().lower()
    if s in _AT_INTERFACE_TOKENS: return "positive"
    if s in _AWAY_INTERFACE_TOKENS: return "negative"
    return "not_tested"


def classify_axis_status(row):
    """
    Map the four ground-truth axis columns to {positive, negative, not_tested},
    using evidence_axes to disambiguate mild_destab.
    Returns None if any required column is missing from the row (Mode B).
    """
    if any(c not in row.index for c in ["expected_ddg_monomer",
                                        "expected_ddg_fold_complex",
                                        "expected_ddg_binding",
                                        "expected_topology"]):
        return None
    evidence = row.get("evidence_axes", "")
    return {
        "fold_monomer": _tokenize_axis(row["expected_ddg_monomer"], evidence),
        "fold_complex": _tokenize_axis(row["expected_ddg_fold_complex"], evidence),
        "binding":      _tokenize_axis(row["expected_ddg_binding"], evidence),
        "topology":     _tokenize_topology(row["expected_topology"]),
    }


def _axis_direction(row):
    """
    For each positive axis, return the direction tag (destab / stab) so we can
    detect direction flips. A positive axis with value 'destab' or 'mild_destab'
    implies direction=destab; 'stab' implies direction=stab.
    """
    out = {}
    for axis_name, col in [
        ("fold_monomer", "expected_ddg_monomer"),
        ("fold_complex", "expected_ddg_fold_complex"),
        ("binding",      "expected_ddg_binding"),
    ]:
        s = ss(row.get(col)).strip().lower()
        if s in ("destab", "mild_destab"): out[axis_name] = "destab"
        elif s == "stab":                  out[axis_name] = "stab"
        else:                               out[axis_name] = None
    return out


# ============================================================================
# NEW in v4: derived expected_mech_class (role-independent, rubric v2)
# ============================================================================
def derive_expected_mech_class(row):
    """
    Derive a mechanism class label from structural annotations alone
    (role does NOT enter). Returns one of:

      structurally_uncommitted        — all four axes not-tested (N/A case)
      interface_uncommitted_magnitude — topology positive, all DDG axes not-tested
      structurally_silent              — all committed axes negative (or
                                         pathogenic-with-functional-evidence
                                         and no firm positive axis)
      fold_mechanism                   — fold axis positive, binding neg/not-tested
      ppi_destab_mechanism             — binding positive (destab direction),
                                         fold neg/not-tested
      ppi_stab_mechanism               — binding positive (stab direction),
                                         fold neg/not-tested
      mixed_structural                 — both fold and binding positive
      NA                                — required columns missing (Mode B)
    """
    # Mode B guard
    if any(c not in row.index for c in GROUND_TRUTH_COLS):
        return "NA"
    axes = classify_axis_status(row)
    if axes is None:
        return "NA"

    evidence = ss(row.get("evidence_axes")).strip().lower()
    role = ss(row.get("role")).strip().lower()

    # Commitment flags
    fold_m_c = axes["fold_monomer"] != "not_tested"
    fold_c_c = axes["fold_complex"] != "not_tested"
    bind_c   = axes["binding"] != "not_tested"
    iface_c  = axes["topology"] != "not_tested"
    any_committed = fold_m_c or fold_c_c or bind_c or iface_c
    fold_committed = fold_m_c or fold_c_c

    # N/A path: nothing committed
    if not any_committed:
        return "structurally_uncommitted"

    # Positive axis flags
    fold_pos = axes["fold_monomer"] == "positive" or axes["fold_complex"] == "positive"
    bind_pos = axes["binding"] == "positive"
    iface_pos = axes["topology"] == "positive"

    # Soft-evidence override for pathogenic variants with no firm structural axis.
    # If evidence is functional/population and no committed axis is destab/stab
    # at the raw annotation level, classify as structurally_silent even if
    # mild_destab would otherwise have been positive (the evidence-aware
    # tokenizer has already demoted it to not_tested in that case, but we
    # re-check the raw annotations here to be explicit).
    if role.startswith("pathogenic") and evidence in SOFT_EVIDENCE:
        raw_strong = any(
            ss(row.get(c)).strip().lower() in ("destab", "stab")
            for c in ("expected_ddg_monomer", "expected_ddg_fold_complex", "expected_ddg_binding")
        )
        if not raw_strong:
            return "structurally_silent"

    # Positive-axis branches
    if fold_pos and bind_pos:
        return "mixed_structural"
    if fold_pos:
        return "fold_mechanism"
    if bind_pos:
        bind_dir = _axis_direction(row).get("binding")
        return "ppi_stab_mechanism" if bind_dir == "stab" else "ppi_destab_mechanism"

    # No positive ΔΔG axes
    if iface_pos and not fold_committed and not bind_c:
        # Interface-only case: topology positive but ΔΔG axes uncommitted
        return "interface_uncommitted_magnitude"

    # All committed axes negative (or only topology committed negative)
    return "structurally_silent"


# ============================================================================
# NEW in v4: mechanism-consistency grading per rubric v2
# ============================================================================
def grade_mechanism_consistency(row, mech_call, expected_class, axes):
    """
    Apply the rubric v2 decision tree.

    Returns a tuple (grade, false_positive_axes, missed_positive_axes) where:
      - grade in {'consistent', 'partial', 'inconsistent', 'N/A'}
      - false_positive_axes: list of tested-negative axes MAVIS fired on
      - missed_positive_axes: list of tested-positive axes MAVIS did not fire

    Parameters
    ----------
    row : pd.Series
    mech_call : str — the MAVIS mechanism label
    expected_class : str — output of derive_expected_mech_class
    axes : dict|None — output of classify_axis_status
    """
    empty = []

    # Step 1: unevaluable structure
    if mech_call in ("Structure unevaluable", "Structure low-confidence at variant site"):
        return "N/A", empty, empty

    # Mode B or uncommitted / interface-only-uncommitted
    if expected_class in ("NA", "structurally_uncommitted",
                          "interface_uncommitted_magnitude") or axes is None:
        return "N/A", empty, empty

    # Resolve MAVIS label assertions
    label_info = MECH_LABEL_ASSERTIONS.get(mech_call)
    if label_info is None:
        return "N/A", empty, empty
    fold_m_a, fold_c_a, bind_a, iface_a, direction = label_info
    # v5 Rubric B: track fold subtypes separately. fold_asserted (legacy
    # "any fold subtype") is retained for the structurally_silent branch
    # below, where the question is just "did MAVIS fire any fold axis".
    fold_asserted = fold_m_a or fold_c_a

    mavis_no_effect = (mech_call == "No structural effect detected")
    mavis_dng = mech_call in (
        "Interface variant (DDG neutral)",
        "Structural variant — contact-driven (DDG neutral)",
        "Structural variant — burial-driven (DDG neutral)",
    )

    # Axis flags from ground truth — v5 Rubric B keeps fold subtypes separate
    fold_m_pos = axes["fold_monomer"] == "positive"
    fold_m_neg = axes["fold_monomer"] == "negative"
    fold_m_unknown = axes["fold_monomer"] == "untested"
    fold_c_pos = axes["fold_complex"] == "positive"
    fold_c_neg = axes["fold_complex"] == "negative"
    fold_c_unknown = axes["fold_complex"] == "untested"
    # Aggregate flags retained for branches that don't need subtype distinction
    fold_pos = fold_m_pos or fold_c_pos
    fold_neg = fold_m_neg and not fold_c_pos  # both negative (or one neg + other unknown), neither positive
    bind_pos = axes["binding"] == "positive"
    bind_neg = axes["binding"] == "negative"
    iface_pos = axes["topology"] == "positive"
    iface_neg = axes["topology"] == "negative"
    any_pos = fold_pos or bind_pos

    # Compute false positives and missed positives
    dirs = _axis_direction(row)
    exp_fold_dir = dirs.get("fold_monomer") or dirs.get("fold_complex")
    exp_bind_dir = dirs.get("binding")

    fp, mp, dflip = [], [], []

    # ---- v5 Rubric B fold-subtype checks (Q1=A, Q2=partial, Q3=strict) ----
    # Q1: axes annotated 'unknown' (untested) are excluded from grading.
    # Q2: missing one fold subtype while matching the other → counted in mp,
    #     producing a partial grade in the aggregation step below.
    # Q3: applies uniformly, including to conflicting calls.

    # Missed monomer-fold: ground truth expects monomer destab/stab, call
    # doesn't assert monomer fold axis.
    if fold_m_pos and not fold_m_a:
        mp.append("fold_monomer")
    # False-positive monomer-fold: ground truth says monomer neutral, call
    # asserts monomer fold axis.
    if fold_m_neg and fold_m_a:
        fp.append("fold_monomer")
    # Missed multimer-fold: ground truth expects fold_complex destab/stab,
    # call doesn't assert multimer fold axis.
    if fold_c_pos and not fold_c_a:
        mp.append("fold_complex")
    # False-positive multimer-fold: ground truth says fold_complex neutral,
    # call asserts multimer fold axis.
    if fold_c_neg and fold_c_a:
        fp.append("fold_complex")

    # Direction flip on fold (only check if at least one fold subtype is
    # asserted AND at least one is annotated positive AND direction is concrete).
    if (fold_pos and fold_asserted and exp_fold_dir
            and direction not in ("conflicting", "none", "neutral")):
        if ((exp_fold_dir == "destab" and direction == "stab") or
            (exp_fold_dir == "stab" and direction == "destab")):
            dflip.append("fold")

    # ---- Binding axis (unchanged from v4 — no subtype distinction needed) ----
    if bind_pos and not bind_a:
        mp.append("binding")
    if bind_neg and bind_a:
        fp.append("binding")
    if (bind_pos and bind_a and exp_bind_dir
            and direction not in ("conflicting", "none", "neutral")):
        if ((exp_bind_dir == "destab" and direction == "stab") or
            (exp_bind_dir == "stab" and direction == "destab")):
            dflip.append("binding")

    # Topology false-fire only if MAVIS asserts interface AND ground truth is
    # explicitly away, AND MAVIS isn't already fully asserting via fold/bind.
    if iface_neg and iface_a and not (fold_asserted or bind_a):
        fp.append("topology")

    # Step 3: structurally_silent expects silence
    if expected_class == "structurally_silent":
        if mavis_no_effect or mavis_dng:
            return "consistent", empty, empty
        # Any specific-axis firing on a structurally-silent variant is a
        # false positive. Per rubric v2, ANY false-fire → inconsistent
        # (no partial credit on structurally_silent).
        # v5: report fold subtype granularly when applicable
        axes_fired = []
        if fold_m_a: axes_fired.append("fold_monomer")
        if fold_c_a: axes_fired.append("fold_complex")
        if bind_a:   axes_fired.append("binding")
        if iface_a and not fold_asserted and not bind_a and iface_neg:
            axes_fired.append("topology")
        return "inconsistent", axes_fired, empty

    # Step 4: direction flip = inconsistent
    if dflip:
        return "inconsistent", fp, mp

    # Step 5: "No structural effect" when positive axes exist = inconsistent
    if mp and mavis_no_effect:
        return "inconsistent", fp, mp

    # Step 6: all positive axes correctly called?
    if not mp:
        if not fp:
            return "consistent", empty, empty
        return "partial", fp, empty

    # Step 7: primary positive axis missed
    if fp:
        # Missed-and-wrong → inconsistent
        return "inconsistent", fp, mp
    # Magnitude miss case: topology positive, MAVIS gives DNG-only
    if mavis_dng and iface_pos:
        return "partial", empty, mp
    # Axis miss without false-fire → partial
    return "partial", empty, mp


def summarize_mech_consistency(grades):
    """Most-conservative grade across thresholds. Ranks:
    inconsistent > partial > consistent > N/A."""
    rank = {"inconsistent": 3, "partial": 2, "consistent": 1, "N/A": 0}
    if not grades:
        return "N/A"
    return max(((rank.get(g, 0), g) for g in grades), key=lambda x: x[0])[1]


# ============================================================================
# Pipeline agreement (unchanged from v3)
# ============================================================================
def classify_agreement(p1_tier, p2_tier):
    t1 = ss(p1_tier); t2 = ss(p2_tier)
    if not t2 or "Awaiting" in t2 or "unevaluable" in t2.lower():
        return "Partially unevaluable"
    t1_high = t1 in FOOTPRINT_TIERS
    t2_high = t2 in FOOTPRINT_TIERS
    if t1_high and t2_high: return "Concordant high"
    if not t1_high and not t2_high: return "Concordant low"
    if t2_high and not t1_high: return "Neighborhood-elevated"
    if t1_high and not t2_high: return "Neighborhood-depressed"
    return "Partially unevaluable"


# ============================================================================
# v5 NEW: structural_agreement, directional_agreement, per-axis votes,
#         tier_structural_signal_type diagnostic, evaluation_note divergence
# ============================================================================

# Constants for sub-threshold direction credit
DDG_MIN_DETECTABLE = 0.5  # below this, |ΔΔG| is treated as "no signal"


def _operative_axis(row, partners, axis_kind):
    """Return (value, partner) for the operative partner on this axis,
    selected by max |ΔΔG| among confident partners with internal-CI
    distinguishability from 0. Returns (None, None) if no qualifying partner."""
    best_val, best_p, best_abs = None, None, -1.0
    for p in partners:
        if not bool(row.get(f"ddg_{p}_confident", False)): continue
        v = row.get(f"ddg_{axis_kind}_{p}")
        if not pd.notna(v): continue
        # Use internal CI95-from-0 distinguishability as the gate
        ci_dist_col = f"ddg_{axis_kind}_{p}_distinguishable_internal_from_0"
        ci_dist = bool(row.get(ci_dist_col, False))
        if not ci_dist:
            continue
        if abs(float(v)) > best_abs:
            best_val = float(v); best_p = p; best_abs = abs(float(v))
    return best_val, best_p


def _axis_check(predicted_value, fires, ground_truth):
    """
    Per-axis grading rule (called only after CI95 distinguishability gate
    has been applied upstream).
    Returns (in_denom, agrees).
      GT unknown → not in denom
      GT destab/mild_destab → fires AND positive direction
      GT stab → fires AND negative direction
      GT neutral → not fires
    """
    gt = str(ground_truth).strip().lower() if pd.notna(ground_truth) else 'unknown'
    if gt in ('unknown', '', 'nan'):
        return False, False
    if gt in ('destab', 'mild_destab'):
        return True, fires and predicted_value > 0
    if gt == 'stab':
        return True, fires and predicted_value < 0
    return True, not fires  # neutral


def compute_structural_agreement(row, partners, mono_thr, fold_thr, bind_thr):
    """
    Per-variant structural_agreement: count gradeable axes (tier + 3 DDG)
    against per-axis ground truth.

    Tier is conditioned on expected_mech_class (Approach 3): the tier is
    expected to fire iff expected_mech_class != 'structurally_silent'.

    DDG axes use symmetric internal-CI95 gating: an axis is excluded if 0
    falls within the internal CI of the prediction (axis is internally
    indistinguishable from 0).

    Returns (n_correct, n_gradeable).
    """
    n, d = 0, 0

    # Tier (Approach 3: conditioned on expected_mech_class)
    if bool(row.get("structure_evaluable", True)):
        emc = str(row.get('expected_mech_class', '')).strip().lower()
        if emc in ('structurally_silent', 'mixed_structural', 'fold_mechanism',
                   'ppi_destab_mechanism', 'ppi_stab_mechanism'):
            d += 1
            tier_fires = str(row.get('mavis_tier')) in FOOTPRINT_TIERS
            tier_expected = emc != 'structurally_silent'
            if tier_fires == tier_expected:
                n += 1

    # Monomer ΔΔG
    if bool(row.get('ddg_monomer_confident', False)) and pd.notna(row.get('ddg_monomer')):
        ci_dist = bool(row.get('ddg_monomer_distinguishable_internal_from_0', False))
        if ci_dist:
            v = float(row['ddg_monomer'])
            fires = abs(v) >= mono_thr
            in_d, ok = _axis_check(v, fires, row.get('expected_ddg_monomer'))
            if in_d:
                d += 1
                if ok: n += 1

    # Fold ΔΔG (operative partner — selects partner with the largest
    # CI-distinguishable signal, which is what should be graded against the
    # complex-context expectation)
    fv, fp = _operative_axis(row, partners, 'fold')
    if fv is not None:
        fires = abs(fv) >= fold_thr
        in_d, ok = _axis_check(fv, fires, row.get('expected_ddg_fold_complex'))
        if in_d:
            d += 1
            if ok: n += 1

    # Binding ΔΔG (operative partner)
    bv, bp = _operative_axis(row, partners, 'binding')
    if bv is not None:
        fires = abs(bv) >= bind_thr
        in_d, ok = _axis_check(bv, fires, row.get('expected_ddg_binding'))
        if in_d:
            d += 1
            if ok: n += 1

    return n, d


def compute_directional_agreement(row, partners, mono_thr, fold_thr, bind_thr):
    """
    Per-variant directional_agreement: like structural_agreement but awards
    half-credit for sub-threshold-but-correct-direction axis predictions.

    Specifically, a DDG axis where:
      - The ground truth expects destab/stab (in either direction),
      - The predicted value is in the "detectable" range
        (|value| >= DDG_MIN_DETECTABLE),
      - The direction matches ground truth (positive for destab, negative for stab),
      - But |value| < threshold (sub-threshold magnitude)
    receives a partial credit of 0.5 instead of 0.

    Returns (n_full_credit, n_half_credit, n_gradeable).
    """
    n_full, n_half, d = 0, 0, 0

    # Tier (binary, no partial credit on tier)
    if bool(row.get("structure_evaluable", True)):
        emc = str(row.get('expected_mech_class', '')).strip().lower()
        if emc in ('structurally_silent', 'mixed_structural', 'fold_mechanism',
                   'ppi_destab_mechanism', 'ppi_stab_mechanism'):
            d += 1
            tier_fires = str(row.get('mavis_tier')) in FOOTPRINT_TIERS
            tier_expected = emc != 'structurally_silent'
            if tier_fires == tier_expected:
                n_full += 1

    def _grade_axis(value, threshold, ground_truth):
        """Returns (in_denom, full_credit_bool, half_credit_bool)."""
        gt = str(ground_truth).strip().lower() if pd.notna(ground_truth) else 'unknown'
        if gt in ('unknown', '', 'nan'):
            return False, False, False
        v = float(value)
        fires = abs(v) >= threshold
        if gt in ('destab', 'mild_destab'):
            if fires and v > 0: return True, True, False
            if abs(v) >= DDG_MIN_DETECTABLE and v > 0: return True, False, True  # sub-threshold but right direction
            return True, False, False
        if gt == 'stab':
            if fires and v < 0: return True, True, False
            if abs(v) >= DDG_MIN_DETECTABLE and v < 0: return True, False, True
            return True, False, False
        # neutral
        return True, (not fires), False

    # Monomer
    if bool(row.get('ddg_monomer_confident', False)) and pd.notna(row.get('ddg_monomer')):
        ci_dist = bool(row.get('ddg_monomer_distinguishable_internal_from_0', False))
        if ci_dist:
            in_d, full, half = _grade_axis(row['ddg_monomer'], mono_thr, row.get('expected_ddg_monomer'))
            if in_d:
                d += 1
                if full: n_full += 1
                elif half: n_half += 1

    # Fold (operative partner)
    fv, fp = _operative_axis(row, partners, 'fold')
    if fv is not None:
        in_d, full, half = _grade_axis(fv, fold_thr, row.get('expected_ddg_fold_complex'))
        if in_d:
            d += 1
            if full: n_full += 1
            elif half: n_half += 1

    # Binding (operative partner)
    bv, bp = _operative_axis(row, partners, 'binding')
    if bv is not None:
        in_d, full, half = _grade_axis(bv, bind_thr, row.get('expected_ddg_binding'))
        if in_d:
            d += 1
            if full: n_full += 1
            elif half: n_half += 1

    return n_full, n_half, d


def compute_axis_votes(row, partners):
    """
    Per-axis vote columns: did each axis fire at strict (≥2.0) and relaxed (≥1.0) thresholds?
    Returns dict with 6 keys: {monomer,fold,binding} × {strict,relaxed}.
    """
    out = {}
    # Monomer
    mv = row.get('ddg_monomer')
    if pd.notna(mv) and bool(row.get('ddg_monomer_confident', False)):
        v = float(mv)
        out['ddg_monomer_vote_strict'] = abs(v) >= 2.0
        out['ddg_monomer_vote_relaxed'] = abs(v) >= 1.0
    else:
        out['ddg_monomer_vote_strict'] = None
        out['ddg_monomer_vote_relaxed'] = None
    # Fold (operative)
    fv, _ = _operative_axis(row, partners, 'fold')
    if fv is not None:
        out['ddg_fold_vote_strict'] = abs(fv) >= 2.0
        out['ddg_fold_vote_relaxed'] = abs(fv) >= 1.0
    else:
        out['ddg_fold_vote_strict'] = None
        out['ddg_fold_vote_relaxed'] = None
    # Binding (operative)
    bv, _ = _operative_axis(row, partners, 'binding')
    if bv is not None:
        out['ddg_binding_vote_strict'] = abs(bv) >= 2.0
        out['ddg_binding_vote_relaxed'] = abs(bv) >= 1.0
    else:
        out['ddg_binding_vote_strict'] = None
        out['ddg_binding_vote_relaxed'] = None
    return out


def compute_tier_structural_signal_type(row, partners):
    """
    Diagnostic column: categorize the tier's underlying structural signal
    as fold_related / binding_related / both / ambiguous / none.

    fold_related: high monomer_n_contacts AND non-surface burial AND
                  low/zero inter-chain contacts
    binding_related: any partner with high inter-chain contacts AND
                     surface-exposed
    both: both patterns present (e.g., HBB W37 — buried at tetramer interface)
    ambiguous: tier fires but neither pattern dominates
    none: tier doesn't fire (Tier 3-4)
    """
    if str(row.get('mavis_tier')) not in FOOTPRINT_TIERS:
        return 'none'

    mono_c = sf(row.get('monomer_n_contacts'), 0)
    burial = ss(row.get('monomer_burial')).lower()
    HIGH_MONO = 5  # threshold for "high monomer contacts" — calibrated against benchmark
    HIGH_INTER = 3  # threshold for "high inter-chain contacts"

    high_mono_contacts = mono_c >= HIGH_MONO
    non_surface_burial = burial in ('buried_core', 'partially_buried')

    max_inter = 0.0
    for p in partners:
        ic = row.get(f'multi_{p}_inter_contacts')
        if pd.notna(ic) and float(ic) > max_inter:
            max_inter = float(ic)
    high_inter = max_inter >= HIGH_INTER

    fold_pattern = high_mono_contacts and non_surface_burial and not high_inter
    binding_pattern = high_inter and burial == 'surface_exposed'
    both_pattern = high_mono_contacts and non_surface_burial and high_inter

    if both_pattern: return 'both'
    if fold_pattern: return 'fold_related'
    if binding_pattern: return 'binding_related'
    return 'ambiguous'


def compute_evaluation_note(row, sa_score_d_25, mc_grade_25):
    """
    Flag variants where structural_agreement and mech_consistency disagree
    directionally at t=2.5, with a short reason. Helps identify case studies
    for the Discussion section without manual hunting.

    Returns a string note, or empty string if no notable divergence.
    """
    sa_n, sa_d = sa_score_d_25
    if sa_d == 0:
        return ""
    sa = sa_n / sa_d

    mc = str(mc_grade_25) if pd.notna(mc_grade_25) else ""

    # Notable divergence: high SA but inconsistent MC
    if sa >= 0.75 and mc == 'inconsistent':
        return "high_axis_agreement_but_inconsistent_synthesis"
    # Or: low SA but consistent MC (right answer for wrong reasons)
    if sa <= 0.25 and mc == 'consistent':
        return "consistent_synthesis_but_low_axis_agreement"
    # Or: borderline SA with partial MC
    if 0.4 <= sa <= 0.6 and mc == 'partial':
        return "borderline_axis_agreement_partial_synthesis"
    return ""


# ============================================================================
# Main
# ============================================================================
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", default="results/mavis_v7_results_with_nbhd.csv")
    ap.add_argument("--external", default="AM_variants_mavis_mechanism_test.xlsx")
    ap.add_argument("--outdir", default="results")
    args = ap.parse_args()

    rp = Path(args.results); ep = Path(args.external); od = Path(args.outdir)
    od.mkdir(parents=True, exist_ok=True)
    if not rp.exists():
        print(f"ERROR: results CSV not found: {rp}", file=sys.stderr); sys.exit(1)
    if not ep.exists():
        print(f"ERROR: external xlsx not found: {ep}", file=sys.stderr); sys.exit(1)

    df = pd.read_csv(rp)
    ext = pd.read_excel(ep)
    required = {"gene", "variant", "AM pathogenicity", "AM class", "franklin"}
    missing = required - set(ext.columns)
    if missing:
        print(f"ERROR: external xlsx missing columns: {sorted(missing)}", file=sys.stderr); sys.exit(1)

    df["gene"] = df["gene"].astype(str).str.lower()
    df["variant"] = df["variant"].astype(str)
    ext["gene"] = ext["gene"].astype(str).str.lower()
    ext["variant"] = ext["variant"].astype(str)
    merged = df.merge(ext, on=["gene", "variant"], how="left", validate="one_to_one")

    # Detect operating mode
    gt_missing = [c for c in GROUND_TRUTH_COLS if c not in merged.columns]
    mode_A = len(gt_missing) == 0
    print(f"\nOperating mode: {'A (benchmark eval, ground truth present)' if mode_A else 'B (structural-signal only)'}")
    if not mode_A:
        print(f"  Missing ground-truth columns: {gt_missing}")
        print(f"  Will skip mech_consistency and expected_mech_class derivation.")

    partners = discover_partners(merged)
    sd_cols = [c for c in merged.columns
               if c.endswith("_sd") and c.startswith(("ddg_monomer", "ddg_fold", "ddg_binding"))]

    # Stage 1: ddg confidence + magnitude
    merged["ddg_confidence_derived"] = merged.apply(
        lambda r: derive_ddg_confidence(r, partners), axis=1)
    merged["max_abs_ddg"] = merged.apply(
        lambda r: compute_max_abs_ddg(r, partners), axis=1)
    merged["low_structural_confidence"] = merged.apply(
        lambda r: low_structural_confidence(r, sd_cols), axis=1)

    # Stage 1b NEW: per-variant 95% CIs on each DDG axis
    # (internal SD-based + Sapozhnikov-bound, with distinguishable_from_threshold flags)
    ci_records = [compute_ddg_cis(r, partners) for _, r in merged.iterrows()]
    ci_df = pd.DataFrame(ci_records)
    for c in ci_df.columns:
        merged[c] = ci_df[c].values

    # Stage 2: Pipeline 1 mechanism at each threshold (v5: 5 thresholds incl. Sapozhnikov per-axis)
    # Column naming: mech_t10, mech_t15, mech_t20, mech_t25, mech_tSAP
    for tag, t in THRESHOLD_SPECS:
        merged[f"mech_{tag}"] = merged.apply(
            lambda r, _t=t: classify_mechanism_at(r, partners, _t, tier_col="mavis_tier"), axis=1)

    # Stage 3: Pipeline 2 — nbhd_score + nbhd_tier
    nbhd_records = [compute_nbhd_score(r, partners) for _, r in merged.iterrows()]
    nbhd_df = pd.DataFrame(nbhd_records)
    for c in nbhd_df.columns:
        merged[c] = nbhd_df[c].values
    merged["nbhd_tier"] = merged["nbhd_final_score"].apply(assign_tier)
    merged["nbhd_structure_evaluable"] = merged["nbhd_final_score"].notna()

    # Stage 4: Pipeline 2 mechanism at each threshold (v5: 5 thresholds)
    # Columns: nbhd_mech_t10, nbhd_mech_t15, nbhd_mech_t20, nbhd_mech_t25, nbhd_mech_tSAP
    for tag, t in THRESHOLD_SPECS:
        merged[f"nbhd_mech_{tag}"] = merged.apply(
            lambda r, _t=t: classify_mechanism_at(r, partners, _t, tier_col="nbhd_tier"), axis=1)

    # Stage 5: concordance for both pipelines × both modes (unchanged)
    p1_full = merged.apply(
        lambda r: compute_concordance(r, r["mavis_tier"], prefix="", include_external=True),
        axis=1, result_type="expand")
    p1_struct = merged.apply(
        lambda r: compute_concordance(r, r["mavis_tier"], prefix="", include_external=False),
        axis=1, result_type="expand")
    p2_full = merged.apply(
        lambda r: compute_concordance(r, r["nbhd_tier"], prefix="nbhd_", include_external=True),
        axis=1, result_type="expand")
    p2_struct = merged.apply(
        lambda r: compute_concordance(r, r["nbhd_tier"], prefix="nbhd_", include_external=False),
        axis=1, result_type="expand")
    for src in [p1_full, p1_struct, p2_full, p2_struct]:
        for c in src.columns:
            merged[c] = src[c].values

    # Stage 5b NEW: structural-signal / external-consensus split (both pipelines)
    p1_split = merged.apply(
        lambda r: compute_signal_consensus_split(r, r["mavis_tier"], prefix=""),
        axis=1, result_type="expand")
    p2_split = merged.apply(
        lambda r: compute_signal_consensus_split(r, r["nbhd_tier"], prefix="nbhd_"),
        axis=1, result_type="expand")
    # external_consensus is pipeline-independent — take it from p1_split only
    # to avoid duplicate columns
    external_cols = [c for c in p1_split.columns if c.startswith("external_consensus")]
    for c in p1_split.columns:
        if c.startswith("external_consensus"):
            merged[c] = p1_split[c].values
        elif c.startswith("structural_signal"):
            merged[c] = p1_split[c].values
    for c in p2_split.columns:
        if c.startswith("nbhd_structural_signal"):
            merged[c] = p2_split[c].values

    # Stage 6 NEW: derived expected_mech_class + mechanism consistency grading
    if mode_A:
        merged["expected_mech_class"] = merged.apply(derive_expected_mech_class, axis=1)

        # Precompute axes once per row, keyed by DataFrame index for safe lookup
        axes_by_idx = {idx: classify_axis_status(r) for idx, r in merged.iterrows()}

        def _grade_one(r, mech_col):
            g, fp, mp = grade_mechanism_consistency(
                r, r.get(mech_col), r.get("expected_mech_class"),
                axes_by_idx.get(r.name))
            return g, ",".join(fp) if fp else "", ",".join(mp) if mp else ""

        # v5: iterate over all 5 thresholds (4 uniform + Sapozhnikov per-axis)
        all_tags = [tag for tag, _ in THRESHOLD_SPECS]  # t10, t15, t20, t25, tSAP
        for suf in all_tags:
            # Pipeline 1
            p1 = [_grade_one(r, f"mech_{suf}") for _, r in merged.iterrows()]
            merged[f"mech_consistency_{suf}"] = [x[0] for x in p1]
            merged[f"mech_false_positive_axes_{suf}"] = [x[1] for x in p1]
            merged[f"mech_missed_positive_axes_{suf}"] = [x[2] for x in p1]
            # Pipeline 2
            p2 = [_grade_one(r, f"nbhd_mech_{suf}") for _, r in merged.iterrows()]
            merged[f"nbhd_mech_consistency_{suf}"] = [x[0] for x in p2]
            merged[f"nbhd_mech_false_positive_axes_{suf}"] = [x[1] for x in p2]
            merged[f"nbhd_mech_missed_positive_axes_{suf}"] = [x[2] for x in p2]

        # v5 Headline: use t=2.5 grade as the canonical 'summary'. Other
        # thresholds are also reported. Threshold-stable flag now considers
        # all 5 thresholds.
        merged["mech_consistency_summary"] = merged["mech_consistency_t25"]
        merged["nbhd_mech_consistency_summary"] = merged["nbhd_mech_consistency_t25"]

        def _stable(r, prefix=""):
            # v5: stable means same grade across all 5 thresholds
            grades = [r.get(f"{prefix}mech_consistency_{s}") for s in all_tags]
            return len(set(grades)) == 1

        merged["mech_consistency_threshold_stable"] = merged.apply(
            lambda r: _stable(r, prefix=""), axis=1)
        merged["nbhd_mech_consistency_threshold_stable"] = merged.apply(
            lambda r: _stable(r, prefix="nbhd_"), axis=1)

        # ==========================================================
        # v5 Stage 6b: structural_agreement, directional_agreement,
        #              per-axis votes, tier diagnostic, evaluation_note
        # ==========================================================
        # structural_agreement at all 5 thresholds
        for tag, t in THRESHOLD_SPECS:
            if isinstance(t, dict):
                mt, ft, bt = t['monomer'], t['fold'], t['binding']
            else:
                mt = ft = bt = float(t)
            sa_results = merged.apply(
                lambda r, _m=mt, _f=ft, _b=bt:
                    compute_structural_agreement(r, partners, _m, _f, _b), axis=1)
            merged[f"structural_agreement_n_{tag}"] = [r[0] for r in sa_results]
            merged[f"structural_agreement_d_{tag}"] = [r[1] for r in sa_results]
            merged[f"structural_agreement_{tag}"] = [
                r[0]/r[1] if r[1] > 0 else None for r in sa_results
            ]

        # directional_agreement at all 5 thresholds (sub-threshold-direction credit)
        for tag, t in THRESHOLD_SPECS:
            if isinstance(t, dict):
                mt, ft, bt = t['monomer'], t['fold'], t['binding']
            else:
                mt = ft = bt = float(t)
            da_results = merged.apply(
                lambda r, _m=mt, _f=ft, _b=bt:
                    compute_directional_agreement(r, partners, _m, _f, _b), axis=1)
            merged[f"directional_agreement_full_{tag}"] = [r[0] for r in da_results]
            merged[f"directional_agreement_half_{tag}"] = [r[1] for r in da_results]
            merged[f"directional_agreement_d_{tag}"] = [r[2] for r in da_results]
            merged[f"directional_agreement_{tag}"] = [
                (r[0] + 0.5*r[1])/r[2] if r[2] > 0 else None for r in da_results
            ]

        # Per-axis vote columns (threshold-independent strict/relaxed flags)
        vote_records = merged.apply(lambda r: compute_axis_votes(r, partners), axis=1)
        for col in ['ddg_monomer_vote_strict', 'ddg_monomer_vote_relaxed',
                    'ddg_fold_vote_strict', 'ddg_fold_vote_relaxed',
                    'ddg_binding_vote_strict', 'ddg_binding_vote_relaxed']:
            merged[col] = vote_records.apply(lambda d, _c=col: d.get(_c))

        # Tier diagnostic column
        merged["tier_structural_signal_type"] = merged.apply(
            lambda r: compute_tier_structural_signal_type(r, partners), axis=1)

        # Evaluation note (divergence flag at t=2.5)
        merged["evaluation_note"] = merged.apply(
            lambda r: compute_evaluation_note(
                r,
                (r.get("structural_agreement_n_t25", 0) or 0,
                 r.get("structural_agreement_d_t25", 0) or 0),
                r.get("mech_consistency_t25")
            ),
            axis=1
        )

        # v14 confidence columns: p1/p2 ddg_concordance per threshold
        # Labels: concordant_disruption / structural_only / ddg_only / concordant_silent
        def _ddg_concordance_label(tier_high, ddg_fires):
            if tier_high and ddg_fires: return "concordant_disruption"
            if tier_high and not ddg_fires: return "structural_only"
            if (not tier_high) and ddg_fires: return "ddg_only"
            return "concordant_silent"

        for tag, t in THRESHOLD_SPECS:
            if isinstance(t, dict):
                # Use the most stringent (binding) for the label cutoff at Sapozhnikov
                cutoff = max(t['monomer'], t['fold'], t['binding'])
            else:
                cutoff = float(t)
            # P1 (mavis_tier)
            def _p1_concordance(r, _c=cutoff):
                tier_high = str(r.get("mavis_tier")) in FOOTPRINT_TIERS
                ddg_fires = sf(r.get("max_abs_ddg"), 0.0) >= _c
                return _ddg_concordance_label(tier_high, ddg_fires)
            merged[f"p1_ddg_concordance_{tag}"] = merged.apply(_p1_concordance, axis=1)
            # P2 (nbhd_tier)
            def _p2_concordance(r, _c=cutoff):
                tier_high = str(r.get("nbhd_tier")) in FOOTPRINT_TIERS
                ddg_fires = sf(r.get("max_abs_ddg"), 0.0) >= _c
                return _ddg_concordance_label(tier_high, ddg_fires)
            merged[f"p2_ddg_concordance_{tag}"] = merged.apply(_p2_concordance, axis=1)
    else:
        print("  (Mode B: skipping expected_mech_class and mech_consistency)")

    # Stage 7: pipeline agreement
    merged["pipeline_agreement"] = merged.apply(
        lambda r: classify_agreement(r.get("mavis_tier"), r.get("nbhd_tier")), axis=1)

    # Save
    master_path = od / "mavis_v7_concordance.csv"
    merged.to_csv(master_path, index=False)
    print(f"\nWrote {master_path}  ({len(merged)} rows, {len(merged.columns)} cols)")

    # Summary
    print("\n=== Pipeline agreement ===")
    print(merged["pipeline_agreement"].value_counts().to_string())
    print("\n=== Pipeline 1 vs Pipeline 2 tier ===")
    print(pd.crosstab(merged["mavis_tier"], merged["nbhd_tier"]).to_string())

    if mode_A:
        print("\n=== expected_mech_class distribution ===")
        print(merged["expected_mech_class"].value_counts(dropna=False).to_string())

        print("\n=== Mechanism consistency by threshold (per rubric v2) ===")
        print("  Per-threshold (P1 / P2):")
        for suf, label in (("t10", "t=1.0"), ("t15", "t=1.5"),
                           ("t20", "t=2.0"), ("t25", "t=2.5")):
            for pipeline, prefix in [("P1", ""), ("P2", "nbhd_")]:
                col = f"{prefix}mech_consistency_{suf}"
                counts = merged[col].value_counts()
                total_graded = sum(counts.get(g, 0) for g in ["consistent", "partial", "inconsistent"])
                if total_graded:
                    score = (counts.get("consistent", 0) + 0.5 * counts.get("partial", 0)) / total_graded
                    print(f"    {pipeline} {label}: consistent={counts.get('consistent', 0):2d}, "
                          f"partial={counts.get('partial', 0):2d}, "
                          f"inconsistent={counts.get('inconsistent', 0):2d}, "
                          f"N/A={counts.get('N/A', 0):2d}, "
                          f"score={score:.3f}")
        print("  Conservative summary (t=2.5 grade — Sapozhnikov-confident threshold):")
        for pipeline, col in [("P1 (single-residue)", "mech_consistency_summary"),
                              ("P2 (neighborhood)",    "nbhd_mech_consistency_summary")]:
            counts = merged[col].value_counts()
            total_graded = sum(counts.get(g, 0) for g in ["consistent", "partial", "inconsistent"])
            if total_graded:
                score = (counts.get("consistent", 0) + 0.5 * counts.get("partial", 0)) / total_graded
                print(f"    {pipeline}: consistent={counts.get('consistent', 0)}, "
                      f"partial={counts.get('partial', 0)}, "
                      f"inconsistent={counts.get('inconsistent', 0)}, "
                      f"N/A={counts.get('N/A', 0)}, "
                      f"score={score:.3f}")
        n_stable_p1 = int(merged["mech_consistency_threshold_stable"].sum())
        n_stable_p2 = int(merged["nbhd_mech_consistency_threshold_stable"].sum())
        print(f"  Threshold-stable variants (all four thresholds agree): "
              f"P1={n_stable_p1}/{len(merged)}, P2={n_stable_p2}/{len(merged)}")

        # Distinguishability summary using monomer axis as a quick gauge
        print("\n=== DDG distinguishability (monomer axis, for orientation) ===")
        for tag, label in (("0", "from 0"), ("t10", "from 1.0"), ("t15", "from 1.5"),
                           ("t20", "from 2.0"), ("t25", "from 2.5")):
            int_col = f"ddg_monomer_distinguishable_internal_from_{tag}"
            sap_col = f"ddg_monomer_distinguishable_sapozhnikov_from_{tag}"
            if int_col in merged.columns and sap_col in merged.columns:
                n_int = merged[int_col].sum()
                n_sap = merged[sap_col].sum()
                print(f"  monomer distinguishable {label} kcal/mol: "
                      f"internal CI={int(n_int)}/{len(merged)}, "
                      f"Sapozhnikov CI={int(n_sap)}/{len(merged)}")

        print("\n=== Concordance means (primary n=33, mech controls excluded) ===")
        primary = merged[merged["role"].isin(
            ["pathogenic", "benign", "pathogenic_lof", "pathogenic_gof"])].copy()
        primary["truth_binary"] = primary["role"].apply(
            lambda r: "pathogenic" if r != "benign" else "benign")
        for pipeline, prefix in [("P1 (single-residue)", ""), ("P2 (neighborhood)", "nbhd_")]:
            for mode in ("full", "struct"):
                for flavor in ("strict", "relaxed"):
                    coln = f"{prefix}concordance_{flavor}_{mode}_n"
                    cold = f"{prefix}concordance_{flavor}_{mode}_denom"
                    if coln not in primary.columns: continue
                    msg = f"  {pipeline} {flavor}/{mode}:"
                    for tc in ("pathogenic", "benign"):
                        sub = primary[primary["truth_binary"] == tc]
                        if len(sub):
                            msg += f"  {tc}={sub[coln].mean():.2f}/{sub[cold].mean():.2f}"
                    print(msg)

        print("\n=== Signal / consensus split means (primary) ===")
        for pipeline, prefix in [("P1", ""), ("P2", "nbhd_")]:
            for flavor in ("strict", "relaxed"):
                sig_n = f"{prefix}structural_signal_{flavor}_n"
                sig_d = f"{prefix}structural_signal_{flavor}_denom"
                con_n = f"external_consensus_{flavor}_n"
                con_d = f"external_consensus_{flavor}_denom"
                if sig_n not in primary.columns: continue
                msg = f"  {pipeline} {flavor}:"
                for tc in ("pathogenic", "benign"):
                    sub = primary[primary["truth_binary"] == tc]
                    if len(sub):
                        msg += (f"  {tc}: signal={sub[sig_n].mean():.2f}/{sub[sig_d].mean():.2f}, "
                                f"consensus={sub[con_n].mean():.2f}/{sub[con_d].mean():.2f}")
                print(msg)


if __name__ == "__main__":
    main()
