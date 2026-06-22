#!/usr/bin/env python3
"""
verify_v5.py — Verify the v5 pipeline downstream against cached intermediate.

Runs the downstream chain (baseline_correct → apply_concordance_v5 → run_evaluate)
against the cached `inputs/intermediate/mavis_v7_results.csv` from the prior
FoldX run. This lets us validate the v5 framework code without re-running FoldX.

The headline numbers are checked against expected v5 values from the
populated methods sketch (documentation/methods_metrics_sketch_v3.md).
If implementation matches the sketch: validation passes.

Usage:
    cd PIPELINE_CURRENT/
    PYTHONPATH=scripts python3 verification/verify_v5.py
"""
from __future__ import annotations
import os
import sys
import shutil
import subprocess
from pathlib import Path
import pandas as pd

# Resolve paths relative to PIPELINE_CURRENT/
HERE = Path(__file__).resolve().parent
PIPE = HERE.parent
SCRIPTS = PIPE / "scripts"
INPUTS_INTERMEDIATE = PIPE / "inputs" / "intermediate"
INPUTS_AM = PIPE / "inputs" / "AM_variants_mavis_mechanism_test.xlsx"
WORK = PIPE / "verification" / "_work"

EXPECTED = {
    # Headlines from sketch v3 (Rubric B)
    "structural_agreement_tSAP": 0.709,
    "structural_agreement_t25": 0.718,
    "mech_consistency_tSAP": 0.750,
    "mech_consistency_t25": 0.761,
    "hbb_pearson_r": 0.894,
    "level1_tpr_t10": 0.913,
    "phenotype_detection_pathogenic_t10": 24,  # of 26
    "phenotype_detection_pathogenic_gof_t10": 5,  # of 8
    "phenotype_detection_benign_silent_t10": 4,  # of 10
}


def _run(cmd, cwd):
    print(f"  $ {' '.join(cmd)}")
    r = subprocess.run(cmd, cwd=str(cwd), capture_output=True, text=True)
    if r.returncode != 0:
        print(f"  STDERR: {r.stderr}")
        raise RuntimeError(f"Command failed: {cmd}")
    return r.stdout


def main():
    if WORK.exists():
        shutil.rmtree(WORK)
    WORK.mkdir(parents=True)
    (WORK / "results").mkdir()
    print(f"Working dir: {WORK}")

    # Stage 1: copy cached pre-correction results into the workspace
    src = INPUTS_INTERMEDIATE / "mavis_v7_results.csv"
    if not src.exists():
        print(f"ERROR: cached intermediate {src} not found.")
        return 1
    shutil.copy(src, WORK / "results" / "mavis_v7_results.csv")
    print(f"\nStage 1: staged cached input")

    # Stage 2: baseline correction (writes results/mavis_v7_results_corrected.csv
    # with v5 multi-threshold corrected mechanism columns)
    print(f"\nStage 2: baseline_correct")
    env = os.environ.copy()
    env["PYTHONPATH"] = str(SCRIPTS)
    subprocess.run(
        ["python3", str(SCRIPTS / "mavis_v7_baseline_correct.py")],
        cwd=str(WORK), env=env, check=True
    )

    # Stage 3: neighborhood (writes results/mavis_v7_results_with_nbhd.csv)
    print(f"\nStage 3: neighborhood")
    # If the pre-cached _with_nbhd file exists from prior run, use it (v13 P2
    # bugfix only matters in apply_concordance_v5, not the _with_nbhd file
    # which retains the same nbhd_* columns format).
    nbhd_cached = INPUTS_INTERMEDIATE / "mavis_v7_results_with_nbhd.csv"
    if nbhd_cached.exists():
        print(f"  Using cached nbhd output: {nbhd_cached}")
        shutil.copy(nbhd_cached, WORK / "results" / "mavis_v7_results_with_nbhd.csv")
    else:
        subprocess.run(
            ["python3", str(SCRIPTS / "mavis_v7_neighborhood.py")],
            cwd=str(WORK), env=env, check=True
        )

    # Stage 4: apply_concordance_v5
    print(f"\nStage 4: apply_concordance_v5")
    subprocess.run(
        ["python3", str(SCRIPTS / "apply_concordance_v5.py"),
         "--results", str(WORK / "results" / "mavis_v7_results_with_nbhd.csv"),
         "--external", str(INPUTS_AM),
         "--outdir", str(WORK / "results")],
        env=env, check=True
    )

    # Stage 5: run_evaluate (Track A under v5)
    print(f"\nStage 5: run_evaluate")
    subprocess.run(
        ["python3", str(SCRIPTS / "run_evaluate.py")],
        cwd=str(WORK), env=env, check=True
    )

    # Stage 6: verify outputs against EXPECTED
    print(f"\n{'='*70}\nVERIFICATION\n{'='*70}")
    conc_csv = WORK / "results" / "mavis_v7_concordance.csv"
    if not conc_csv.exists():
        print(f"  FAIL: {conc_csv} not produced")
        return 1
    df = pd.read_csv(conc_csv)
    print(f"  Loaded concordance: {len(df)} variants × {len(df.columns)} cols")

    actual = {}

    # structural_agreement
    for tag in ("t25", "tSAP"):
        n_col = f"structural_agreement_n_{tag}"
        d_col = f"structural_agreement_d_{tag}"
        if n_col in df.columns and d_col in df.columns:
            n_total = df[n_col].fillna(0).sum()
            d_total = df[d_col].fillna(0).sum()
            actual[f"structural_agreement_{tag}"] = round(n_total/d_total, 3) if d_total else 0

    # mech_consistency
    for tag in ("t25", "tSAP"):
        col = f"mech_consistency_{tag}"
        if col in df.columns:
            vc = df[col].value_counts()
            n = len(df) - df[col].isna().sum()
            score = (vc.get("consistent", 0) + 0.5*vc.get("partial", 0)) / n if n else 0
            actual[f"mech_consistency_{tag}"] = round(score, 3)

    # Level 1 TPR @ t=1.0 (from saved CSV)
    l1_csv = WORK / "results" / "evaluation" / "level1_t10.csv"
    if l1_csv.exists():
        l1 = pd.read_csv(l1_csv)
        m = l1[l1["classifier"] == "MAVIS_full"]
        if not m.empty:
            actual["level1_tpr_t10"] = round(float(m.iloc[0]["TPR"]), 3)

    # Per-phenotype detection counts
    det_csv = WORK / "results" / "evaluation" / "level2_detection_t10.csv"
    if det_csv.exists():
        det = pd.read_csv(det_csv)
        for phen, key in [("pathogenic", "phenotype_detection_pathogenic_t10"),
                          ("pathogenic_gof", "phenotype_detection_pathogenic_gof_t10")]:
            sub = det[det["phenotype"] == phen]
            if not sub.empty:
                actual[key] = int(sub.iloc[0]["n_detected_structural_disruption"])
        bsub = det[det["phenotype"] == "benign"]
        if not bsub.empty:
            actual["phenotype_detection_benign_silent_t10"] = int(bsub.iloc[0]["n_correctly_silent"])

    # HBB Pearson r — read from saved level3_hbb_quantitative.csv
    hbb_csv = WORK / "results" / "evaluation" / "level3_hbb_quantitative.csv"
    if hbb_csv.exists():
        # The CSV has columns like {pearson_r, spearman_rho, mae, ...} or rows
        hbb_df = pd.read_csv(hbb_csv)
        # Find pearson_r in either columns or first row
        if "pearson_r" in hbb_df.columns:
            actual["hbb_pearson_r"] = round(float(hbb_df.iloc[0]["pearson_r"]), 3)
        elif "metric" in hbb_df.columns:
            m = hbb_df[hbb_df["metric"] == "pearson_r"]
            if not m.empty:
                actual["hbb_pearson_r"] = round(float(m.iloc[0]["value"]), 3)

    # Compare actual vs expected
    print(f"\n{'Metric':<45} {'Expected':>12} {'Actual':>12} {'Status'}")
    print(f"{'-'*45} {'-'*12} {'-'*12} {'-'*7}")
    all_pass = True
    for k, exp in EXPECTED.items():
        act = actual.get(k)
        if act is None:
            print(f"{k:<45} {exp:>12} {'(missing)':>12} FAIL")
            all_pass = False
            continue
        # Tolerance for floats
        if isinstance(exp, float):
            ok = abs(act - exp) < 0.01
        else:
            ok = act == exp
        status = "OK" if ok else "FAIL"
        print(f"{k:<45} {exp:>12} {act:>12} {status}")
        if not ok:
            all_pass = False

    print(f"\n{'='*70}")
    if all_pass:
        print("ALL HEADLINE NUMBERS MATCH EXPECTED. v5 implementation verified.")
        return 0
    else:
        print("VERIFICATION FAILED — see above.")
        return 1


if __name__ == "__main__":
    sys.exit(main())
