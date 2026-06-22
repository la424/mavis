"""
Mechanism classification — 17 categories (16 from v6 + 1 new low-confidence).

v7 changes from v6-A3:
  1. Renamed labels (epistemic, not phenotypic):
       "Benign (structurally evaluated)"            → "No structural effect detected"
       "Stabilizing (potential GoF)"                → "Fold stabilization"
       "Stabilizing at interface (potential GoF)"   → "Fold stabilization at interface"
       "Stabilizing + PPI (potential GoF)"          → "Fold + PPI stabilization"
       "PPI stabilization (potential GoF)"          → "PPI stabilization"
  2. NEW category: "Structure low-confidence at variant site"
       — fires when structure is evaluable but NO axis passes the confidence gate
  3. pLDDT gate added to contact-driven and burial-driven branches
       (previously ungated — v7 fix)
  4. Crystal/NMR structures bypass pLDDT gating entirely via the structure
     loading sentinel (pLDDT = 100 when plddt_gate=False).

Category list (renamed categories marked with *):

  Fold disrupted:
    1  Fold + PPI destabilization
    2  Fold destabilization at interface
    3  Fold destabilization
    4  Fold destab. + PPI stabilization (conflicting)

  Fold stabilized *:
    5  Fold + PPI stabilization *
    6  Fold stabilization at interface *
    7  Fold stabilization *
    8  Fold stab. + PPI destabilization (conflicting)

  Fold neutral:
    9  PPI destabilization (mono neutral)
   10  PPI stabilization *
   11  PPI conflicting (mixed partner signals)
   12  Interface variant (DDG neutral)
   13  Structural variant — contact-driven (DDG neutral)
   14  Structural variant — burial-driven (DDG neutral)
   15  No structural effect detected *

  Meta:
   16  Structure unevaluable
   17  Structure low-confidence at variant site  *** NEW ***
"""

import pandas as pd
from typing import Tuple

from .constants import sf, si, ss, sb, BURIAL_RANK, RANK_TO_BURIAL
from .config import DDG_DESTAB, DDG_HIGHLY, DISRUPTION_POINTS, CONTACT_DRIVEN_THRESHOLD, TIER_THRESHOLDS, DEFAULT_TIER


# ============================================================================
# Tier assignment
# ============================================================================
def assign_tier(score):
    if pd.isna(score):
        return DEFAULT_TIER
    s = float(score)
    for thresh, label in TIER_THRESHOLDS:
        if s >= thresh:
            return label
    return DEFAULT_TIER


# ============================================================================
# Pipeline 1 structural disruption score
# ============================================================================
def compute_disruption_score(row, partner_labels):
    """
    Pipeline 1 score — severity × (mono + max inter-chain contacts).
    pLDDT gate applied to partner contacts (multi_{pl}_plddt >= 50).
    """
    score = 0.0
    ev = []

    sev = sf(row.get("substitution_severity"), 0)
    mono_c = sf(row.get("monomer_n_contacts"), 0)

    max_inter = 0.0
    best_inter_partner = None
    gated_iface_partners = []

    for pl in partner_labels:
        plddt_col = f"multi_{pl}_plddt"
        ic_col = f"multi_{pl}_inter_contacts"
        iface_col = f"multi_{pl}_is_interface"

        pl_plddt = sf(row.get(plddt_col), 0) if plddt_col in row.index else 0
        if pl_plddt < 50:
            continue

        if iface_col in row.index and sb(row.get(iface_col)):
            gated_iface_partners.append(pl)

        if ic_col in row.index and pd.notna(row[ic_col]):
            ic_val = float(row[ic_col])
            if ic_val > max_inter:
                max_inter = ic_val
                best_inter_partner = pl

    total_contacts = mono_c + max_inter
    disruption = round(sev * total_contacts, 2)

    for thresh, pts in DISRUPTION_POINTS:
        if disruption >= thresh:
            score += pts
            ev.append(f"disruption({disruption:.1f})")
            break
    else:
        ev.append(f"no_disruption({disruption:.1f})")

    # Interface bonus
    n_gated = len(gated_iface_partners)
    if n_gated >= 2:
        score += 2.0
        ev.append(f"multi_interface({n_gated})")
    elif n_gated == 1:
        score += 1.5
        ev.append(f"interface({gated_iface_partners[0]})")

    # Burial (pLDDT-gated)
    mono_burial = ss(row.get("monomer_burial"))
    best_rank = BURIAL_RANK.get(mono_burial, 0)
    best_source = "monomer"
    for pl in partner_labels:
        burial_col = f"multi_{pl}_burial"
        plddt_col = f"multi_{pl}_plddt"
        if burial_col in row.index and pd.notna(row.get(burial_col)):
            if sf(row.get(plddt_col), 0) >= 50:
                pl_rank = BURIAL_RANK.get(ss(row[burial_col]), 0)
                if pl_rank > best_rank:
                    best_rank = pl_rank
                    best_source = pl
    best_burial = RANK_TO_BURIAL.get(best_rank, "unknown")
    if best_burial == "buried_core":
        score += 2.0
        ev.append(f"buried_core({best_source})")
    elif best_burial == "partially_buried":
        score += 1.0
        ev.append(f"partially_buried({best_source})")

    # pLDDT multiplier (uses best_plddt from data layer — already sentinel-filled for crystal)
    bp = row.get("best_plddt")
    if bp is not None and not pd.isna(bp):
        bp_val = float(bp)
        if bp_val < 50:
            score *= 0.4
            ev.append(f"very_low_plddt({int(bp_val)})")
        elif bp_val < 70:
            score *= 0.7
            ev.append(f"low_plddt({int(bp_val)})")

    return pd.Series({
        "mavis_score": round(score, 2),
        "mavis_score_evidence": ";".join(ev),
        "contact_disruption": disruption,
        "total_contacts": total_contacts,
        "max_inter_contacts": max_inter,
        "best_inter_partner": best_inter_partner or "",
        "best_burial": best_burial,
        "best_burial_source": best_source,
        "interface_partners_gated": ";".join(gated_iface_partners) if gated_iface_partners else "",
        "n_interface_partners_gated": n_gated,
    })


# ============================================================================
# Evaluability
# ============================================================================
def compute_evaluability(df, partner_labels):
    """
    structure_evaluable — best_plddt >= 50 (or crystal/NMR sentinel)
    multimer_evaluable  — at least one partner pLDDT >= 50
    ddg_evaluable       — ddg_confidence is not low/empty (if column present)
    """
    df["structure_evaluable"] = df["best_plddt"].apply(
        lambda x: pd.notna(x) and float(x) >= 50
    )

    def any_partner_confident(row):
        for pl in partner_labels:
            col = f"multi_{pl}_plddt"
            if col in row.index and pd.notna(row[col]) and float(row[col]) >= 50:
                return True
        return False
    df["multimer_evaluable"] = df.apply(any_partner_confident, axis=1)

    # ddg_evaluable: True when any per-axis confidence flag is True.
    # Per-axis flags are ddg_monomer_confident + ddg_{partner}_confident for
    # each partner label; these are populated by populate_ddg_confidence().
    conf_cols = ["ddg_monomer_confident"] + [f"ddg_{pl}_confident" for pl in partner_labels]
    conf_cols = [c for c in conf_cols if c in df.columns]
    if conf_cols:
        df["ddg_evaluable"] = df[conf_cols].any(axis=1)
    else:
        df["ddg_evaluable"] = False
    return df


# ============================================================================
# Mechanism classifier — 17 categories
# ============================================================================
def _resolve_thresholds(ddg_destab):
    """
    Normalize the threshold parameter into a per-axis dict.
    Accepts:
      - None → use config.DDG_DESTAB for all axes (backward compat)
      - scalar (float) → uniform threshold for all axes
      - dict with keys 'monomer', 'fold', 'binding' → per-axis thresholds
    Returns dict {'monomer': float, 'fold': float, 'binding': float}.
    """
    if ddg_destab is None:
        return {'monomer': DDG_DESTAB, 'fold': DDG_DESTAB, 'binding': DDG_DESTAB}
    if isinstance(ddg_destab, (int, float)):
        v = float(ddg_destab)
        return {'monomer': v, 'fold': v, 'binding': v}
    if isinstance(ddg_destab, dict):
        return {
            'monomer': float(ddg_destab.get('monomer', DDG_DESTAB)),
            'fold': float(ddg_destab.get('fold', DDG_DESTAB)),
            'binding': float(ddg_destab.get('binding', DDG_DESTAB)),
        }
    raise TypeError(f"ddg_destab must be None, scalar, or dict; got {type(ddg_destab)}")


def classify_mechanism(row, partner_labels, ddg_destab=None) -> Tuple[str, str, bool]:
    """
    Return (mechanism_label, partner_name, external_evidence_flag).

    The fold axis is split: monomer-driven and multimer-driven fold
    disruption are reported as distinct categories. This produces 32
    possible mechanism strings (vs 17 in the v6/v7-pre-v5 scheme).

    Args:
      row: pd.Series — the variant row from the pipeline DataFrame.
      partner_labels: list[str] — partner gene labels (e.g. ['bard1']).
      ddg_destab: None, scalar, or dict {'monomer','fold','binding'}.
        Defaults to config.DDG_DESTAB for all axes (backward compat).

    Per-axis thresholds enable Sapozhnikov-confident grading
    (mono=2.9, fold=2.9, bind=3.5 kcal/mol).

    Contact/burial branches are pLDDT-gated: fire only if best_plddt >= 50.
    """
    thr = _resolve_thresholds(ddg_destab)
    mono_thr = thr['monomer']
    fold_thr = thr['fold']
    bind_thr = thr['binding']

    tier = ss(row.get("mavis_tier"))
    is_high_tier = tier in ("Tier 1", "Tier 2")
    structure_eval = row.get("structure_evaluable", False)

    ddg_m = row.get("ddg_monomer")
    mono_conf = row.get("ddg_monomer_confident", False)

    mono_c = sf(row.get("monomer_n_contacts"), 0)
    gated_iface = ss(row.get("interface_partners_gated"))
    is_iface = len(gated_iface) > 0

    # ----- Structure unevaluable -----
    if not structure_eval:
        return "Structure unevaluable", "", False

    # ----- Fold axis: split into monomer-driven vs multimer-driven -----
    has_dm = pd.notna(ddg_m)
    dm_val = float(ddg_m) if has_dm else 0.0
    mono_fold_destab = has_dm and mono_conf and dm_val > mono_thr
    mono_fold_stab = has_dm and mono_conf and dm_val < -mono_thr

    partner_fold_destab = False
    partner_fold_stab = False
    any_axis_confident = mono_conf
    for pl in partner_labels:
        fv = row.get(f"ddg_fold_{pl}")
        bv = row.get(f"ddg_binding_{pl}")
        conf = row.get(f"ddg_{pl}_confident", False)
        if conf and (pd.notna(fv) or pd.notna(bv)):
            any_axis_confident = True
        if pd.notna(fv) and conf:
            f_val = float(fv)
            if f_val > fold_thr:
                partner_fold_destab = True
            if f_val < -fold_thr:
                partner_fold_stab = True

    # Build fold-split category indicator. None means no fold signal.
    # Categories: 'mono_destab', 'multi_destab', 'both_destab',
    #             'mono_stab',   'multi_stab',   'both_stab'
    fold_cat = None
    if mono_fold_destab and partner_fold_destab:
        fold_cat = 'both_destab'
    elif mono_fold_destab:
        fold_cat = 'mono_destab'
    elif partner_fold_destab:
        fold_cat = 'multi_destab'
    elif mono_fold_stab and partner_fold_stab:
        fold_cat = 'both_stab'
    elif mono_fold_stab:
        fold_cat = 'mono_stab'
    elif partner_fold_stab:
        fold_cat = 'multi_stab'

    # ----- PPI axis -----
    ppi_destab = False
    ppi_stab = False
    ppi_partner_destab = None
    ppi_partner_stab = None
    has_per_partner = False

    for pl in partner_labels:
        bv = row.get(f"ddg_binding_{pl}")
        if pd.notna(bv):
            has_per_partner = True
            conf = row.get(f"ddg_{pl}_confident", False)
            if not conf:
                continue
            bv_f = float(bv)
            if bv_f > bind_thr:
                ppi_destab = True
                if ppi_partner_destab is None:
                    ppi_partner_destab = pl
            if bv_f < -bind_thr:
                ppi_stab = True
                if ppi_partner_stab is None:
                    ppi_partner_stab = pl

    ppi_partner = ppi_partner_destab or ppi_partner_stab or ""

    # ----- Structure low-confidence at variant site -----
    # Fires when structure exists but NO axis (monomer or any partner) is confident.
    if not any_axis_confident:
        return "Structure low-confidence at variant site", "", False

    # ----- Fold + PPI destabilization (same direction) -----
    if fold_cat in ('mono_destab', 'multi_destab', 'both_destab') and ppi_destab:
        prefix = {'mono_destab': 'Monomer fold', 'multi_destab': 'Multimer fold',
                  'both_destab': 'Both fold'}[fold_cat]
        return f"{prefix} + PPI destabilization", ppi_partner_destab or "", False

    # ----- Fold destab + PPI stab (conflicting) -----
    if fold_cat in ('mono_destab', 'multi_destab', 'both_destab') and ppi_stab:
        prefix = {'mono_destab': 'Monomer fold destab.', 'multi_destab': 'Multimer fold destab.',
                  'both_destab': 'Both fold destab.'}[fold_cat]
        return f"{prefix} + PPI stabilization (conflicting)", ppi_partner_stab or "", False

    # ----- Fold + PPI stabilization (same direction) -----
    if fold_cat in ('mono_stab', 'multi_stab', 'both_stab') and ppi_stab:
        prefix = {'mono_stab': 'Monomer fold', 'multi_stab': 'Multimer fold',
                  'both_stab': 'Both fold'}[fold_cat]
        return f"{prefix} + PPI stabilization", ppi_partner_stab or "", False

    # ----- Fold stab + PPI destab (conflicting) -----
    if fold_cat in ('mono_stab', 'multi_stab', 'both_stab') and ppi_destab:
        prefix = {'mono_stab': 'Monomer fold stab.', 'multi_stab': 'Multimer fold stab.',
                  'both_stab': 'Both fold stab.'}[fold_cat]
        return f"{prefix} + PPI destabilization (conflicting)", ppi_partner_destab or "", False

    # ----- Pure fold (no PPI signal) -----
    if fold_cat in ('mono_destab', 'multi_destab', 'both_destab'):
        prefix = {'mono_destab': 'Monomer fold destabilization',
                  'multi_destab': 'Multimer fold destabilization',
                  'both_destab': 'Both fold destabilization'}[fold_cat]
        suffix = " at interface" if is_iface else ""
        return f"{prefix}{suffix}", "", False

    if fold_cat in ('mono_stab', 'multi_stab', 'both_stab'):
        prefix = {'mono_stab': 'Monomer fold stabilization',
                  'multi_stab': 'Multimer fold stabilization',
                  'both_stab': 'Both fold stabilization'}[fold_cat]
        suffix = " at interface" if is_iface else ""
        return f"{prefix}{suffix}", "", False

    # ----- Fold neutral -----
    if ppi_destab and not ppi_stab:
        return "PPI destabilization", ppi_partner_destab or "", False
    if ppi_stab and not ppi_destab:
        return "PPI stabilization", ppi_partner_stab or "", False
    if ppi_destab and ppi_stab:
        return "PPI conflicting (mixed partner signals)", ppi_partner, False
    if is_iface:
        return "Interface variant (DDG neutral)", "", False

    # Contact-driven / burial-driven — pLDDT-gated
    best_plddt = sf(row.get("best_plddt"), 0)
    if is_high_tier and best_plddt >= 50:
        if mono_c >= CONTACT_DRIVEN_THRESHOLD:
            return "Structural variant — contact-driven (DDG neutral)", "", False
        else:
            return "Structural variant — burial-driven (DDG neutral)", "", False

    return "No structural effect detected", "", False


# ============================================================================
# Canonical label list (v5 fold-split scheme — 32 categories)
# ============================================================================
ALL_MECHANISMS = [
    # Pure fold destab (6)
    "Monomer fold destabilization",
    "Monomer fold destabilization at interface",
    "Multimer fold destabilization",
    "Multimer fold destabilization at interface",
    "Both fold destabilization",
    "Both fold destabilization at interface",
    # Pure fold stab (6)
    "Monomer fold stabilization",
    "Monomer fold stabilization at interface",
    "Multimer fold stabilization",
    "Multimer fold stabilization at interface",
    "Both fold stabilization",
    "Both fold stabilization at interface",
    # Fold + PPI same direction destab (3)
    "Monomer fold + PPI destabilization",
    "Multimer fold + PPI destabilization",
    "Both fold + PPI destabilization",
    # Fold + PPI same direction stab (3)
    "Monomer fold + PPI stabilization",
    "Multimer fold + PPI stabilization",
    "Both fold + PPI stabilization",
    # Fold destab + PPI stab conflicting (3)
    "Monomer fold destab. + PPI stabilization (conflicting)",
    "Multimer fold destab. + PPI stabilization (conflicting)",
    "Both fold destab. + PPI stabilization (conflicting)",
    # Fold stab + PPI destab conflicting (3)
    "Monomer fold stab. + PPI destabilization (conflicting)",
    "Multimer fold stab. + PPI destabilization (conflicting)",
    "Both fold stab. + PPI destabilization (conflicting)",
    # Pure PPI (3)
    "PPI destabilization",
    "PPI stabilization",
    "PPI conflicting (mixed partner signals)",
    # Structural without DDG (3)
    "Interface variant (DDG neutral)",
    "Structural variant — contact-driven (DDG neutral)",
    "Structural variant — burial-driven (DDG neutral)",
    # Catch-all (3)
    "No structural effect detected",
    "Structure unevaluable",
    "Structure low-confidence at variant site",
]
