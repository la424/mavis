"""
MAVIS v7 — Phase 3 evaluation module (revised).

FOUR levels of benchmark evaluation:
  Level 1: Binary pathogenic vs benign (TPR/TNR/accuracy + baselines)
  Level 2: Three-class pathogenic_lof/pathogenic_gof/benign confusion matrix
  Level 3: Detailed axis-level mechanism matching + HBB quantitative
  Level 4: GoF mechanism detection (per-variant + summary, headline)

Primary threshold: 1.0 kcal/mol (matches the MAVIS pipeline's production
operating threshold used in the CHD paper). Results at 1.5 and 2.0 reported
side-by-side for transparency across reasonable FoldX literature thresholds.

Sensitivity analyses (mechanism controls, pLDDT) reported at threshold 1.0.

No scipy dependency — bootstrap and correlation via pure numpy.
"""
from __future__ import annotations

import random
from dataclasses import dataclass
from typing import Callable, Dict, List, Tuple

import numpy as np
import pandas as pd


# ============================================================================
# Pure-numpy statistics
# ============================================================================
def pearson_r(x: np.ndarray, y: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mx, my = x.mean(), y.mean()
    num = ((x - mx) * (y - my)).sum()
    den = float(np.sqrt(((x - mx) ** 2).sum() * ((y - my) ** 2).sum()))
    return float(num / den) if den > 0 else float("nan")


def spearman_rho(x: np.ndarray, y: np.ndarray) -> float:
    def ranks(a):
        a = np.asarray(a, dtype=float)
        order = np.argsort(a, kind="mergesort")
        ranks_out = np.empty_like(order, dtype=float)
        n = len(a)
        i = 0
        while i < n:
            j = i
            while j + 1 < n and a[order[j + 1]] == a[order[i]]:
                j += 1
            avg = (i + j) / 2.0 + 1
            for k in range(i, j + 1):
                ranks_out[order[k]] = avg
            i = j + 1
        return ranks_out
    return pearson_r(ranks(x), ranks(y))


def bootstrap_metric(
    predictions: np.ndarray,
    labels: np.ndarray,
    metric_fn: Callable[[np.ndarray, np.ndarray], float],
    n_iter: int = 1000,
    seed: int = 42,
) -> Tuple[float, float, float]:
    rng = random.Random(seed)
    n = len(predictions)
    if n == 0:
        return float("nan"), float("nan"), float("nan")
    point = metric_fn(predictions, labels)
    values = []
    for _ in range(n_iter):
        idx = [rng.randrange(n) for _ in range(n)]
        try:
            v = metric_fn(predictions[idx], labels[idx])
            if not (np.isnan(v) or np.isinf(v)):
                values.append(v)
        except Exception:
            continue
    if not values:
        return point, float("nan"), float("nan")
    values = np.array(sorted(values))
    return float(point), float(np.percentile(values, 2.5)), float(np.percentile(values, 97.5))


def tpr(predictions: np.ndarray, labels: np.ndarray) -> float:
    pos = labels == 1
    if pos.sum() == 0:
        return float("nan")
    return float((predictions[pos] == 1).sum() / pos.sum())


def tnr(predictions: np.ndarray, labels: np.ndarray) -> float:
    neg = labels == 0
    if neg.sum() == 0:
        return float("nan")
    return float((predictions[neg] == 0).sum() / neg.sum())


def accuracy(predictions: np.ndarray, labels: np.ndarray) -> float:
    if len(predictions) == 0:
        return float("nan")
    return float((predictions == labels).sum() / len(predictions))


# ============================================================================
# Classifier config and helpers
# ============================================================================
@dataclass
class ClassifierConfig:
    name: str
    plddt_gate: int = 50
    ddg_destab_threshold: float = 1.0


def _monomer_confident(row, plddt_gate: int) -> bool:
    p = row.get("monomer_plddt")
    if pd.isna(p):
        return False
    return float(p) >= plddt_gate


def _partner_confident(row, pl: str, plddt_gate: int) -> bool:
    p = row.get(f"multi_{pl}_plddt")
    if pd.isna(p):
        return False
    return float(p) >= plddt_gate


# ============================================================================
# Level 1 classifiers
# ============================================================================
def mavis_full_positive(row, cfg: ClassifierConfig) -> bool:
    """fold OR binding OR interface (all pLDDT-gated at cfg.plddt_gate)."""
    # Fold
    fold_disrupted = False
    dm = row.get("ddg_monomer")
    if pd.notna(dm) and _monomer_confident(row, cfg.plddt_gate):
        if abs(float(dm)) > cfg.ddg_destab_threshold:
            fold_disrupted = True
    if not fold_disrupted:
        for col in row.index:
            if col.startswith("ddg_fold_") and not col.endswith("_sd"):
                pl = col.replace("ddg_fold_", "")
                fv = row.get(col)
                if pd.notna(fv) and _partner_confident(row, pl, cfg.plddt_gate):
                    if abs(float(fv)) > cfg.ddg_destab_threshold:
                        fold_disrupted = True
                        break
    # Binding
    binding_disrupted = False
    for col in row.index:
        if (col.startswith("ddg_binding_")
                and not col.endswith("_sd")
                and not col.endswith("_indistinguishable")):
            pl = col.replace("ddg_binding_", "")
            bv = row.get(col)
            flag = bool(row.get(f"{col}_indistinguishable", False))
            if pd.notna(bv) and _partner_confident(row, pl, cfg.plddt_gate) and not flag:
                if abs(float(bv)) > cfg.ddg_destab_threshold:
                    binding_disrupted = True
                    break
    # Interface
    iface = str(row.get("interface_partners_gated", "")).strip()
    is_iface = False
    if iface and iface.lower() != "nan":
        for pl in iface.split(";"):
            ic = row.get(f"multi_{pl}_inter_contacts")
            if pd.notna(ic) and float(ic) > 0:
                disrupt = row.get("contact_disruption", 0)
                if pd.notna(disrupt) and float(disrupt) >= 1.5:
                    is_iface = True
                    break
    return fold_disrupted or binding_disrupted or is_iface


def monomer_only_positive(row, cfg: ClassifierConfig) -> bool:
    dm = row.get("ddg_monomer")
    if pd.notna(dm) and _monomer_confident(row, cfg.plddt_gate):
        return abs(float(dm)) > cfg.ddg_destab_threshold
    return False


def structural_score_only_positive(row, cfg: ClassifierConfig) -> bool:
    tier = str(row.get("mavis_tier", "")).strip()
    return tier in ("Tier 1", "Tier 2")


# ============================================================================
# Level 1 — single-threshold engine + multi-threshold merged view
# ============================================================================
def level1_single_threshold(
    df: pd.DataFrame,
    threshold: float,
    plddt_gate: int,
    include_mechanism_controls: bool,
    n_bootstrap: int,
) -> pd.DataFrame:
    cfg = ClassifierConfig(name=f"t{threshold}",
                           plddt_gate=plddt_gate,
                           ddg_destab_threshold=threshold)
    df_eval = df.copy()
    if not include_mechanism_controls:
        df_eval = df_eval[df_eval["role"] != "mechanism_control"]

    def label_fn(r):
        phen = str(r.get("phenotype", "")).lower()
        if phen.startswith("pathogenic"):
            return 1
        if phen == "benign":
            return 0
        return -1

    df_eval = df_eval.copy()
    df_eval["_label"] = df_eval.apply(label_fn, axis=1)
    df_eval = df_eval[df_eval["_label"] != -1]
    labels = df_eval["_label"].values.astype(int)

    classifiers = {
        "MAVIS_full": mavis_full_positive,
        "monomer_only": monomer_only_positive,
        "structural_score_only": structural_score_only_positive,
    }

    rows = []
    for name, fn in classifiers.items():
        preds = df_eval.apply(lambda r: 1 if fn(r, cfg) else 0, axis=1).values.astype(int)
        tpr_pt, tpr_lo, tpr_hi = bootstrap_metric(preds, labels, tpr, n_bootstrap, seed=42)
        tnr_pt, tnr_lo, tnr_hi = bootstrap_metric(preds, labels, tnr, n_bootstrap, seed=43)
        acc_pt, acc_lo, acc_hi = bootstrap_metric(preds, labels, accuracy, n_bootstrap, seed=44)
        n_pos = int((labels == 1).sum())
        n_neg = int((labels == 0).sum())
        tp = int(((preds == 1) & (labels == 1)).sum())
        fn = int(((preds == 0) & (labels == 1)).sum())
        tn = int(((preds == 0) & (labels == 0)).sum())
        fp = int(((preds == 1) & (labels == 0)).sum())
        rows.append({
            "classifier": name, "threshold": threshold, "plddt_gate": plddt_gate,
            "n_pos": n_pos, "n_neg": n_neg,
            "TP": tp, "FN": fn, "TN": tn, "FP": fp,
            "TPR": round(tpr_pt, 3),
            "TPR_CI_lo": round(tpr_lo, 3), "TPR_CI_hi": round(tpr_hi, 3),
            "TNR": round(tnr_pt, 3),
            "TNR_CI_lo": round(tnr_lo, 3), "TNR_CI_hi": round(tnr_hi, 3),
            "accuracy": round(acc_pt, 3),
            "accuracy_CI_lo": round(acc_lo, 3), "accuracy_CI_hi": round(acc_hi, 3),
        })
    return pd.DataFrame(rows)


def level1_merged_table(
    df: pd.DataFrame,
    thresholds: List[float] = (1.0, 1.5, 2.0),
    plddt_gate: int = 50,
    include_mechanism_controls: bool = False,
    n_bootstrap: int = 1000,
) -> Tuple[pd.DataFrame, Dict]:
    per_threshold = {}
    for t in thresholds:
        per_threshold[t] = level1_single_threshold(
            df, threshold=t, plddt_gate=plddt_gate,
            include_mechanism_controls=include_mechanism_controls,
            n_bootstrap=n_bootstrap,
        )
    classifiers = list(per_threshold[thresholds[0]]["classifier"])
    rows = []
    for cls in classifiers:
        row = {"classifier": cls}
        for t in thresholds:
            sub = per_threshold[t][per_threshold[t]["classifier"] == cls].iloc[0]
            if "n_pos" not in row:
                row["n_pos"] = sub["n_pos"]
                row["n_neg"] = sub["n_neg"]
            tkey = f"t{int(t*10):02d}"
            row[f"TP_{tkey}"] = sub["TP"]
            row[f"FN_{tkey}"] = sub["FN"]
            row[f"TN_{tkey}"] = sub["TN"]
            row[f"FP_{tkey}"] = sub["FP"]
            row[f"TPR_{tkey}"] = sub["TPR"]
            row[f"TPR_{tkey}_CI"] = f"[{sub['TPR_CI_lo']}, {sub['TPR_CI_hi']}]"
            row[f"TNR_{tkey}"] = sub["TNR"]
            row[f"TNR_{tkey}_CI"] = f"[{sub['TNR_CI_lo']}, {sub['TNR_CI_hi']}]"
            row[f"acc_{tkey}"] = sub["accuracy"]
            row[f"acc_{tkey}_CI"] = f"[{sub['accuracy_CI_lo']}, {sub['accuracy_CI_hi']}]"
        rows.append(row)
    return pd.DataFrame(rows), per_threshold


# ============================================================================
# Level 2 — three-class (with mechanism column chosen per threshold)
# ============================================================================
# Mechanism → predicted class mapping.
#
# Per discussion in Phase 3, "Fold destab. + PPI stabilization (conflicting)"
# is mapped to pred_gof rather than pred_lof. Rationale: the signature of
# tighter-than-WT binding combined with fold destabilization in the complex
# context is the structural fingerprint of released auto-inhibition —
# a classical GoF mechanism (e.g., PIK3CA E545K: nSH2 contact loosens while
# kinase domain gains activity; Zhao & Vogt 2008).
#
# "Fold stab. + PPI destab. (conflicting)" remains in pred_lof pending data
# to support a symmetric remap. This category was empty in the current
# benchmark at both thresholds.
# ============================================================================
# v5 NOTE: LoF/GoF mech-to-pred-class mapping has been DROPPED.
# ============================================================================
# Previous versions mapped mechanism strings to pred_lof / pred_gof / pred_benign
# using the destabilization → LoF / stabilization → GoF convention. This
# convention is leaky in three known ways (stabilization can drive LoF via
# trapping inactive conformations; destabilization can drive GoF via
# auto-inhibition release as in PIK3CA E545K; many GoF mechanisms are
# structurally invisible — kinetic, post-translational, allosteric).
#
# Under v5 the framework reports per-phenotype detection rates instead:
# "structural disruption detected in X/Y variants of phenotype Z" — a
# convention-independent measure that doesn't claim the pipeline can
# distinguish GoF from LoF based on structural signal alone.
#
# The function `is_structural_disruption_detected` replaces `mech_to_pred_class`
# as the canonical "did the pipeline see structural perturbation" check.

FOOTPRINT_TIERS = frozenset({"Tier 1", "Tier 2"})


def is_structural_disruption_detected(row, ddg_threshold: float = 1.0) -> bool:
    """
    v5: detection = (mavis_tier in {Tier 1, Tier 2}) OR any DDG axis fires.
    Convention-independent: doesn't care about destab vs stab direction,
    just whether MAVIS produced a non-silent prediction.

    Args:
      row: pd.Series with mavis_tier and ddg_* columns
      ddg_threshold: scalar threshold (a single value applied to all axes)
    """
    if str(row.get("mavis_tier")) in FOOTPRINT_TIERS:
        return True
    # Monomer
    if pd.notna(row.get("ddg_monomer")) and bool(row.get("ddg_monomer_confident", False)):
        if abs(float(row["ddg_monomer"])) >= ddg_threshold:
            return True
    # Per-partner fold + binding (any partner firing is detection)
    for col in row.index:
        if not isinstance(col, str):
            continue
        if col.startswith("ddg_fold_") or col.startswith("ddg_binding_"):
            if col.endswith("_sd") or col.endswith("_indistinguishable") \
               or "_distinguishable" in col or col.endswith("_confident"):
                continue
            # Extract partner label
            if col.startswith("ddg_fold_"):
                partner = col[len("ddg_fold_"):]
            else:
                partner = col[len("ddg_binding_"):]
            conf_col = f"ddg_{partner}_confident"
            if conf_col not in row.index or not bool(row.get(conf_col, False)):
                continue
            v = row.get(col)
            if pd.notna(v) and abs(float(v)) >= ddg_threshold:
                return True
    return False


def level2_per_phenotype_detection(
    df: pd.DataFrame,
    ddg_threshold: float = 1.0,
    include_mechanism_controls: bool = True,
) -> pd.DataFrame:
    """
    v5 replacement for level2_three_class.

    Reports detection rate per phenotype:
      - pathogenic / pathogenic_lof / pathogenic_gof: fraction with structural
        disruption detected
      - benign: fraction correctly silent (no structural disruption detected)

    Detection is convention-independent (any structural signal counts).
    No LoF/GoF prediction mapping is used.

    By default, mechanism_control variants are included (v5 dissolves the
    mechanism_control category into the regular pool — see methods §2).
    """
    df_eval = df.copy()
    if not include_mechanism_controls:
        df_eval = df_eval[df_eval["role"] != "mechanism_control"]

    df_eval = df_eval.copy()
    df_eval["_detected"] = df_eval.apply(
        lambda r: is_structural_disruption_detected(r, ddg_threshold), axis=1
    )

    rows = []
    # Group by phenotype
    for phenotype, sub in df_eval.groupby("phenotype", dropna=False):
        n = len(sub)
        if pd.isna(phenotype) or phenotype == "":
            continue
        n_detected = int(sub["_detected"].sum())
        if str(phenotype).lower() == "benign":
            # For benigns, correctly_silent is the meaningful number
            n_correctly_silent = n - n_detected
            rows.append({
                "phenotype": phenotype,
                "n": n,
                "n_detected_structural_disruption": n_detected,
                "n_correctly_silent": n_correctly_silent,
                "rate": round(n_correctly_silent / n, 3) if n else float("nan"),
                "rate_meaning": "fraction_correctly_silent",
            })
        else:
            rows.append({
                "phenotype": phenotype,
                "n": n,
                "n_detected_structural_disruption": n_detected,
                "n_correctly_silent": n - n_detected,
                "rate": round(n_detected / n, 3) if n else float("nan"),
                "rate_meaning": "fraction_with_disruption_detected",
            })

    return pd.DataFrame(rows)


# ============================================================================
# Level 3 — axis matching + HBB
# ============================================================================
def discretize_ddg(val: float, threshold: float = 1.0) -> str:
    if pd.isna(val):
        return "unknown"
    v = float(val)
    if v > threshold:
        return "destab"
    if v > 0.5:
        return "mild_destab"
    if v < -threshold:
        return "stab"
    return "neutral"


def level3_mechanism_axis(
    df: pd.DataFrame,
    threshold: float = 1.0,
    include_mechanism_controls: bool = False,
) -> pd.DataFrame:
    df_eval = df.copy()
    if not include_mechanism_controls:
        df_eval = df_eval[df_eval["role"] != "mechanism_control"]

    def fold_pred(row):
        vals = []
        dm = row.get("ddg_monomer")
        if pd.notna(dm) and row.get("ddg_monomer_confident", False):
            vals.append(float(dm))
        for col in row.index:
            if col.startswith("ddg_fold_") and not col.endswith("_sd"):
                pl = col.replace("ddg_fold_", "")
                if row.get(f"ddg_{pl}_confident", False):
                    fv = row.get(col)
                    if pd.notna(fv):
                        vals.append(float(fv))
        if not vals:
            return "unknown"
        return discretize_ddg(max(vals, key=abs), threshold)

    def binding_pred(row):
        vals = []
        for col in row.index:
            if (col.startswith("ddg_binding_")
                    and not col.endswith("_sd")
                    and not col.endswith("_indistinguishable")):
                pl = col.replace("ddg_binding_", "")
                bv = row.get(col)
                flag = bool(row.get(f"{col}_indistinguishable", False))
                conf = row.get(f"ddg_{pl}_confident", False)
                if pd.notna(bv) and conf and not flag:
                    vals.append(float(bv))
        if not vals:
            return "neutral"
        return discretize_ddg(max(vals, key=abs), threshold)

    def topology_pred(row):
        iface = str(row.get("interface_partners_gated", "")).strip()
        return "at_interface" if iface and iface.lower() != "nan" else "away_from_interface"

    df_eval = df_eval.copy()
    df_eval["_fold_pred"] = df_eval.apply(fold_pred, axis=1)
    df_eval["_binding_pred"] = df_eval.apply(binding_pred, axis=1)
    df_eval["_topology_pred"] = df_eval.apply(topology_pred, axis=1)

    rows = []
    for axis, pred_col, truth_col in [
        ("fold", "_fold_pred", "expected_ddg_fold_complex"),
        ("binding", "_binding_pred", "expected_ddg_binding"),
        ("topology", "_topology_pred", "expected_topology"),
    ]:
        if truth_col not in df_eval.columns:
            continue
        sub = df_eval[df_eval[truth_col].astype(str).str.lower().isin(
            {"destab", "mild_destab", "neutral", "stab",
             "at_interface", "away_from_interface"}
        )]
        total = len(sub)
        if total == 0:
            continue
        def hit(r):
            p, t = str(r[pred_col]), str(r[truth_col]).lower()
            if p == t:
                return 1
            if axis in ("fold", "binding"):
                if p == "mild_destab" and t == "destab":
                    return 1
                if p == "destab" and t == "mild_destab":
                    return 1
            return 0
        hits = int(sub.apply(hit, axis=1).sum())
        rows.append({
            "axis": axis, "threshold": threshold,
            "n_evaluable": total, "n_hit": hits,
            "accuracy": round(hits / total, 3),
        })
    return pd.DataFrame(rows)


def level3_hbb_quantitative(df: pd.DataFrame) -> Dict:
    expected = {"W37Y": 2.0, "W37A": 5.0, "W37G": 7.0, "W37E": 9.0}
    exp_vals, obs_vals, variant_obs = [], [], {}
    for variant, exp_val in expected.items():
        sub = df[(df["system"] == "hemoglobin_tetramer") & (df["variant"] == variant)]
        if len(sub) == 0:
            continue
        obs = sub.iloc[0].get("ddg_binding_hba1_2")
        if pd.isna(obs):
            continue
        exp_vals.append(exp_val)
        obs_vals.append(float(obs))
        variant_obs[variant] = float(obs)

    if len(exp_vals) < 2:
        return {"n": 0, "pearson_r": float("nan"), "spearman_rho": float("nan"),
                "mae": float("nan"), "rmse": float("nan"), "variants": {}}
    exp_arr = np.array(exp_vals)
    obs_arr = np.array(obs_vals)
    return {
        "n": len(exp_vals),
        "pearson_r": round(pearson_r(exp_arr, obs_arr), 4),
        "spearman_rho": round(spearman_rho(exp_arr, obs_arr), 4),
        "mae": round(float(np.mean(np.abs(exp_arr - obs_arr))), 3),
        "rmse": round(float(np.sqrt(np.mean((exp_arr - obs_arr) ** 2))), 3),
        "variants": variant_obs,
    }


# ============================================================================
# Level 4 — GoF mechanism detection (NEW)
# ============================================================================
def level4_gof_detection(
    df: pd.DataFrame,
    thresholds: List[float] = (1.0, 1.5, 2.0),
    plddt_gate: int = 50,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    v5: Per-variant analysis across all 8 pathogenic_gof-phenotype variants.

    Reframed under v5 to report "structural disruption detected in pathogenic_gof
    variants" rather than "GoF mechanism predicted." Same numerator (any structural
    signal counts), different framing — does NOT use the leaky LoF/GoF mapping.

    Stabilization signals are still surfaced per-variant for interpretation
    (auto-inhibition release biology context, e.g. PIK3CA E545K), but no
    automatic mapping to pred_gof is performed.
    """
    gofs = df[df["phenotype"] == "pathogenic_gof"].copy()

    per_variant_rows = []
    for _, r in gofs.iterrows():
        # Surface stabilization signals per variant (interpretive context, not a prediction)
        stab_signals = []
        dm = r.get("ddg_monomer")
        if pd.notna(dm) and float(dm) < -1.0:
            stab_signals.append(f"monomer={dm:+.2f}")
        for col in df.columns:
            if col.startswith("ddg_fold_") and not col.endswith("_sd"):
                pl = col.replace("ddg_fold_", "")
                fv = r.get(col)
                if pd.notna(fv) and float(fv) < -1.0:
                    stab_signals.append(f"fold_{pl}={fv:+.2f}")
            elif (col.startswith("ddg_binding_")
                  and not col.endswith("_sd")
                  and not col.endswith("_indistinguishable")):
                pl = col.replace("ddg_binding_", "")
                bv = r.get(col)
                flag = bool(r.get(f"{col}_indistinguishable", False))
                if pd.notna(bv) and float(bv) < -1.0 and not flag:
                    stab_signals.append(f"bind_{pl}={bv:+.2f}")

        row_info = {
            "system": r["system"],
            "variant": r["variant"],
            "role": r["role"],
            "phenotype": r["phenotype"],
            "stabilization_signals": "; ".join(stab_signals) if stab_signals else "(none)",
        }
        for t in thresholds:
            cfg = ClassifierConfig(name=f"t{t}", plddt_gate=plddt_gate,
                                   ddg_destab_threshold=t)
            tkey = f"t{int(t*10):02d}"
            row_info[f"mavis_positive_{tkey}"] = bool(mavis_full_positive(r, cfg))
            row_info[f"monomer_positive_{tkey}"] = bool(monomer_only_positive(r, cfg))
            row_info[f"struct_positive_{tkey}"] = bool(structural_score_only_positive(r, cfg))
            # v5: per-threshold structural-disruption detection (convention-independent)
            row_info[f"structural_disruption_detected_{tkey}"] = bool(
                is_structural_disruption_detected(r, ddg_threshold=t)
            )
        # v5: report mechanism string at t10 and t15 (the patched columns)
        for mech_col in ("mavis_mechanism_corrected_t10", "mavis_mechanism_corrected_t15"):
            if mech_col in r.index:
                row_info[mech_col] = r.get(mech_col, "")
        per_variant_rows.append(row_info)

    per_variant_df = pd.DataFrame(per_variant_rows)

    n_gof = len(gofs)
    summary_rows = []
    for t in thresholds:
        tkey = f"t{int(t*10):02d}"
        mavis_hits = int(per_variant_df[f"mavis_positive_{tkey}"].sum())
        mono_hits = int(per_variant_df[f"monomer_positive_{tkey}"].sum())
        struct_hits = int(per_variant_df[f"struct_positive_{tkey}"].sum())
        struct_dis_hits = int(per_variant_df[f"structural_disruption_detected_{tkey}"].sum())
        summary_rows.append({
            "classifier_threshold": f"structural_disruption_detected_t{int(t*10):02d}",
            "n_gof_total": n_gof,
            "MAVIS_full_hits": mavis_hits,
            "MAVIS_full_recall": round(mavis_hits / n_gof, 3),
            "monomer_only_hits": mono_hits,
            "monomer_only_recall": round(mono_hits / n_gof, 3),
            "structural_score_hits": struct_hits,
            "structural_score_recall": round(struct_hits / n_gof, 3),
            "any_structural_signal_hits": struct_dis_hits,
            "any_structural_signal_recall": round(struct_dis_hits / n_gof, 3),
        })
    return per_variant_df, pd.DataFrame(summary_rows)


# ============================================================================
# Top-level runner
# ============================================================================
def run_full_evaluation(
    df: pd.DataFrame,
    thresholds: List[float] = (1.0, 1.5, 2.0),
    n_bootstrap: int = 1000,
) -> Dict:
    out = {"thresholds": list(thresholds)}

    # Level 1
    out["level1_merged"], out["level1_per_threshold"] = level1_merged_table(
        df, thresholds=thresholds, plddt_gate=50,
        include_mechanism_controls=False, n_bootstrap=n_bootstrap,
    )
    # Sensitivities at primary (t=1.0)
    out["level1_sensitivity_incl_controls"] = level1_single_threshold(
        df, threshold=1.0, plddt_gate=50,
        include_mechanism_controls=True, n_bootstrap=n_bootstrap,
    )
    out["level1_sensitivity_strict_plddt"] = level1_single_threshold(
        df, threshold=1.0, plddt_gate=70,
        include_mechanism_controls=False, n_bootstrap=n_bootstrap,
    )

    # Level 2 — v5: replaced three-class confusion (LoF/GoF mapping was leaky)
    # with per-phenotype detection rates. Computed at all sweep thresholds.
    # Detection = high tier OR any DDG axis fires (convention-independent).
    out["level2_per_phenotype_detection"] = {}
    for t in thresholds:
        out["level2_per_phenotype_detection"][t] = level2_per_phenotype_detection(
            df, ddg_threshold=t, include_mechanism_controls=True
        )

    # Merged level 2 — per-phenotype detection rate at each threshold
    # (replaces former level2_merged_recall)
    if out["level2_per_phenotype_detection"]:
        sample = next(iter(out["level2_per_phenotype_detection"].values()))
        phenotypes = sorted(sample["phenotype"].unique().tolist())
    else:
        phenotypes = []
    l2_rows = []
    for phen in phenotypes:
        row = {"phenotype": phen}
        for t, det_df in out["level2_per_phenotype_detection"].items():
            tkey = f"t{int(t*10):02d}"
            m = det_df[det_df["phenotype"] == phen]
            if not m.empty:
                row[f"n_{tkey}"] = int(m.iloc[0]["n"])
                row[f"n_detected_{tkey}"] = int(m.iloc[0]["n_detected_structural_disruption"])
                row[f"n_correctly_silent_{tkey}"] = int(m.iloc[0]["n_correctly_silent"])
                row[f"rate_{tkey}"] = m.iloc[0]["rate"]
                row[f"rate_meaning_{tkey}"] = m.iloc[0]["rate_meaning"]
        l2_rows.append(row)
    out["level2_merged_detection"] = pd.DataFrame(l2_rows)

    # Level 3 axis
    out["level3_axis_per_threshold"] = {}
    for t in thresholds:
        out["level3_axis_per_threshold"][t] = level3_mechanism_axis(
            df, threshold=t, include_mechanism_controls=False)

    axes = ["fold", "binding", "topology"]
    l3_rows = []
    for axis in axes:
        row = {"axis": axis}
        for t, ax_df in out["level3_axis_per_threshold"].items():
            tkey = f"t{int(t*10):02d}"
            m = ax_df[ax_df["axis"] == axis]
            if not m.empty:
                row[f"n_{tkey}"] = m.iloc[0]["n_evaluable"]
                row[f"accuracy_{tkey}"] = m.iloc[0]["accuracy"]
        l3_rows.append(row)
    out["level3_axis_merged"] = pd.DataFrame(l3_rows)

    # Level 3 HBB
    out["level3_hbb_quantitative"] = level3_hbb_quantitative(df)

    # Level 4
    out["level4_per_variant"], out["level4_summary"] = level4_gof_detection(
        df, thresholds=thresholds, plddt_gate=50)

    return out
