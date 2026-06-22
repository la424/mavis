"""
MAVIS v7 baseline-correction post-processor.

Adds per-partner `binding_indistinguishable_from_baseline` flags to the results
CSV and re-runs mechanism classification with flagged values treated as neutral.

Baseline logic (per system, per partner):
  - Eligible variants: multi_{partner}_inter_contacts == 0
  - If n_eligible >= 3: threshold = max(0.5, 1.0 * MAD(eligible_binding_values))
                        center = median(eligible_binding_values)
  - If n_eligible < 3:  center = 0.0, threshold = 0.5  (fallback)

A variant's binding_{partner} is flagged as indistinguishable if:
  abs(ddg_binding_{partner} - center) <= threshold

During reclassification, flagged binding values are treated as 0.0 for the
classify_mechanism function, which correctly routes them to non-PPI categories.

Usage:
  cd ~/mavis_v7
  python3 mavis_v7_baseline_correct.py

Writes:
  results/mavis_v7_results_corrected.csv  (with baseline flags + new mechanism)
  results/baseline_summary.csv            (per-system-partner baseline stats)
"""
from pathlib import Path
import json
import sys
import pandas as pd
import numpy as np

from mavis_v7.config import build_benchmark_config
from mavis_v7.mechanism import classify_mechanism


# ============================================================================
# Baseline computation
# ============================================================================
def compute_baseline(df: pd.DataFrame, partner_label: str) -> dict:
    """
    Compute per-system baseline for one partner label.
    Returns dict: {system: (center, threshold, n_eligible, n_used)}
    """
    baselines = {}
    ic_col = f"multi_{partner_label}_inter_contacts"
    bd_col = f"ddg_binding_{partner_label}"

    if ic_col not in df.columns or bd_col not in df.columns:
        return baselines

    for sys_name in df["system"].unique():
        sys_rows = df[df["system"] == sys_name]
        if ic_col not in sys_rows.columns:
            continue

        # Eligible: structural contacts == 0 (and not NaN) and binding DDG not NaN
        eligible_mask = (
            sys_rows[ic_col].fillna(-1).astype(float) == 0.0
        ) & sys_rows[bd_col].notna()
        eligible = sys_rows[eligible_mask][bd_col].astype(float).values

        n = len(eligible)
        if n >= 3:
            center = float(np.median(eligible))
            mad = float(np.median(np.abs(eligible - center)))
            threshold = max(0.5, 1.0 * mad)
            baselines[sys_name] = {
                "center": round(center, 4),
                "threshold": round(threshold, 4),
                "n_eligible": n,
                "method": "per_system_median",
                "source_values": [round(v, 4) for v in eligible.tolist()],
            }
        else:
            # Fallback
            baselines[sys_name] = {
                "center": 0.0,
                "threshold": 0.5,
                "n_eligible": n,
                "method": "universal_fallback",
                "source_values": [round(v, 4) for v in eligible.tolist()],
            }

    return baselines


def apply_baseline_flags(df: pd.DataFrame, partner_labels: list) -> tuple:
    """
    Add ddg_binding_{partner}_indistinguishable columns to df.
    Returns (df_with_flags, baseline_summary_df).
    """
    summary_rows = []

    # Pre-allocate all flag columns in one operation to avoid pandas
    # performance warnings from insertion one-at-a-time.
    new_cols = {f"ddg_binding_{pl}_indistinguishable": False for pl in partner_labels}
    df = pd.concat([df, pd.DataFrame(new_cols, index=df.index)], axis=1)

    for pl in partner_labels:
        baselines = compute_baseline(df, pl)
        flag_col = f"ddg_binding_{pl}_indistinguishable"

        bd_col = f"ddg_binding_{pl}"
        if bd_col not in df.columns:
            continue

        for sys_name, b in baselines.items():
            center = b["center"]
            thresh = b["threshold"]
            ic_col = f"multi_{pl}_inter_contacts"

            # Baseline-indistinguishable flag applies only to variants that
            # do NOT contact the partner (inter_contacts == 0). An interface
            # variant with a small binding DDG is real signal, not baseline
            # noise; zeroing it would suppress a legitimate low-magnitude
            # interaction change. This mirrors the eligibility filter used
            # during baseline computation.
            mask = (df["system"] == sys_name) & df[bd_col].notna()
            if ic_col in df.columns:
                mask &= (df[ic_col].fillna(-1).astype(float) == 0.0)

            for idx in df[mask].index:
                val = float(df.at[idx, bd_col])
                if abs(val - center) <= thresh:
                    df.at[idx, flag_col] = True

            summary_rows.append({
                "system": sys_name,
                "partner": pl,
                "baseline_center": b["center"],
                "baseline_threshold": b["threshold"],
                "n_eligible": b["n_eligible"],
                "method": b["method"],
                "source_values": ";".join(str(v) for v in b["source_values"]),
            })

    summary_df = pd.DataFrame(summary_rows)
    return df, summary_df


# ============================================================================
# Reclassification with flag-aware binding values
# ============================================================================
def reclassify_with_flags(df: pd.DataFrame, partner_labels: list) -> pd.DataFrame:
    """
    Re-run classify_mechanism with flag-masked binding values, at all five
    v5 thresholds. Writes one set of corrected columns per threshold:
      - mavis_mechanism_corrected_t10
      - mavis_mechanism_corrected_t15
      - mavis_mechanism_corrected_t20
      - mavis_mechanism_corrected_t25
      - mavis_mechanism_corrected_tSAP

    The legacy column `mavis_mechanism_corrected` (no suffix) is kept as an
    alias for `_t10` to preserve backward compatibility with downstream
    scripts (run_evaluate.py, build_report.py).
    """
    threshold_specs = [
        ('t10',  1.0),
        ('t15',  1.5),
        ('t20',  2.0),
        ('t25',  2.5),
        ('tSAP', {'monomer': 2.9, 'fold': 2.9, 'binding': 3.5}),
    ]

    # Build masked-row generator once — masking is threshold-independent
    masked_rows = []
    for _, row in df.iterrows():
        row_masked = row.copy()
        for pl in partner_labels:
            flag_col = f"ddg_binding_{pl}_indistinguishable"
            bd_col = f"ddg_binding_{pl}"
            if flag_col in row_masked.index and bool(row_masked[flag_col]):
                if bd_col in row_masked.index and pd.notna(row_masked[bd_col]):
                    row_masked[bd_col] = 0.0
        masked_rows.append(row_masked)

    new_cols_dict = {}
    for tag, thr in threshold_specs:
        mechs, partners, exts = [], [], []
        for r in masked_rows:
            m, p, e = classify_mechanism(r, partner_labels, ddg_destab=thr)
            mechs.append(m); partners.append(p); exts.append(e)
        new_cols_dict[f"mavis_mechanism_corrected_{tag}"] = mechs
        new_cols_dict[f"mavis_mechanism_partner_corrected_{tag}"] = partners
        new_cols_dict[f"mavis_external_evidence_flag_corrected_{tag}"] = exts

    # Backward-compatible aliases
    new_cols_dict["mavis_mechanism_corrected"] = new_cols_dict["mavis_mechanism_corrected_t10"]
    new_cols_dict["mavis_mechanism_partner_corrected"] = new_cols_dict["mavis_mechanism_partner_corrected_t10"]
    new_cols_dict["mavis_external_evidence_flag_corrected"] = new_cols_dict["mavis_external_evidence_flag_corrected_t10"]

    new_cols = pd.DataFrame(new_cols_dict, index=df.index)
    df = pd.concat([df, new_cols], axis=1)

    return df


# ============================================================================
# Main
# ============================================================================
def main():
    results_path = Path("results/mavis_v7_results.csv")
    if not results_path.exists():
        print(f"✗ Results file not found: {results_path}")
        sys.exit(1)

    df = pd.read_csv(results_path)
    print(f"Loaded {len(df)} variants, {len(df.columns)} columns")

    # Discover partner labels from ddg_binding_* columns (excluding _sd and _indistinguishable)
    partner_labels = sorted({
        c.replace("ddg_binding_", "")
        for c in df.columns
        if c.startswith("ddg_binding_")
        and not c.endswith("_sd")
        and not c.endswith("_indistinguishable")
    })
    print(f"Partner labels: {partner_labels}")

    # Apply baseline flags
    print("\n[1/2] Computing per-system baselines and applying flags...")
    df, summary_df = apply_baseline_flags(df, partner_labels)

    summary_path = Path("results/baseline_summary.csv")
    summary_df.to_csv(summary_path, index=False)
    print(f"    Baseline summary: {summary_path}")

    # Count flags
    n_flagged_total = 0
    for pl in partner_labels:
        flag_col = f"ddg_binding_{pl}_indistinguishable"
        if flag_col in df.columns:
            n = df[flag_col].sum()
            n_flagged_total += n
    print(f"    Total flag hits across all variant×partner combos: {n_flagged_total}")

    # Reclassify
    print("\n[2/2] Re-running mechanism classification with flag-masked binding...")
    df = reclassify_with_flags(df, partner_labels)

    # Compare mechanism distributions
    print("\n  Original mechanism distribution:")
    print(df["mavis_mechanism"].value_counts().to_string())
    print("\n  Corrected mechanism distribution:")
    print(df["mavis_mechanism_corrected"].value_counts().to_string())

    # Diff: variants where mechanism CHANGED
    changed = df[df["mavis_mechanism"] != df["mavis_mechanism_corrected"]]
    print(f"\n  {len(changed)} variants had mechanism change after correction:")
    for _, r in changed.iterrows():
        print(f"    {r['system']}/{r['variant']}: "
              f"'{r['mavis_mechanism']}' → '{r['mavis_mechanism_corrected']}'")

    out_path = Path("results/mavis_v7_results_corrected.csv")
    df.to_csv(out_path, index=False)
    print(f"\n✓ Wrote corrected CSV: {out_path}")


if __name__ == "__main__":
    main()
