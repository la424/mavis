#!/usr/bin/env python3
"""
verify_stage6.py — Downstream-only verification for the v5 framework.

Re-runs Stage 6 (mechanism reclassification + structural_agreement +
mech_consistency + diagnostic columns) against the cached intermediate
file, without touching FoldX or the structural metrics computation.

Compares output to expected v5 headline numbers.

Usage:
  python verify_stage6.py --intermediate /path/to/mavis_v7_results_with_nbhd.csv
"""
import argparse
import sys
from pathlib import Path
import subprocess
import pandas as pd

# v5 expected headlines (from methods_metrics_sketch_v3 final values)
EXPECTED_V5 = {
    'structural_agreement': {
        't10': 0.658, 't15': 0.718, 't20': 0.709,
        't25': 0.718, 'tSAP': 0.709,
    },
    'mech_consistency': {
        't10': 0.602, 't15': 0.693, 't20': 0.716,
        't25': 0.761, 'tSAP': 0.750,
    },
    'threshold_stable': 28,
    'level1_TPR_t10_MAVIS_full': 0.913,
    'level1_TPR_t10_monomer_only': 0.391,
    'level3_HBB_pearson_r': 0.89,
    'level2_pathogenic_detection_t10': 24,
    'level2_pathogenic_gof_detection_t10': 5,
    'level2_benign_silent_t25': 7,
}

TOLERANCE = 0.005  # within 0.5 percentage points


def check(label, actual, expected, tol=TOLERANCE):
    """Print a checkmark or cross for a numeric assertion."""
    if expected is None:
        print(f"  [SKIP] {label}: no expected value")
        return None
    if actual is None:
        print(f"  [FAIL] {label}: actual is None")
        return False
    delta = abs(actual - expected)
    if delta <= tol:
        print(f"  [ OK ] {label}: {actual:.3f} (expected {expected:.3f})")
        return True
    print(f"  [FAIL] {label}: {actual:.3f} (expected {expected:.3f}, diff {delta:.3f})")
    return False


def run_concordance(input_csv, output_dir, scripts_dir, am_xlsx):
    """Run apply_concordance_v5.py and return path to concordance.csv."""
    cmd = [
        sys.executable,
        str(scripts_dir / 'apply_concordance_v5.py'),
        '--results', str(input_csv),
        '--external', str(am_xlsx),
        '--outdir', str(output_dir),
    ]
    print(f"\n  Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"\n  ERROR: apply_concordance_v5 failed")
        print(result.stdout[-2000:])
        print("STDERR:")
        print(result.stderr[-2000:])
        sys.exit(1)
    return output_dir / 'mavis_v7_concordance.csv'


def run_evaluation(corrected_csv, output_dir, scripts_dir):
    """Run run_evaluate.py to produce Track A outputs."""
    cmd = [
        sys.executable,
        str(scripts_dir / 'run_evaluate.py'),
        '--input', str(corrected_csv),
        '--output-dir', str(output_dir),
    ]
    print(f"\n  Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"\n  ERROR: run_evaluate failed")
        print(result.stdout[-2000:])
        print("STDERR:")
        print(result.stderr[-2000:])
        sys.exit(1)
    return output_dir / 'evaluation'


def verify_concordance(concordance_csv):
    """Verify Track B (concordance) headlines."""
    df = pd.read_csv(concordance_csv)
    results = []
    print("\n" + "=" * 70)
    print("Track B verification (apply_concordance_v5)")
    print("=" * 70)

    # structural_agreement
    print("\nstructural_agreement at all 5 thresholds:")
    for tag in ('t10', 't15', 't20', 't25', 'tSAP'):
        n_col = f'structural_agreement_n_{tag}'
        d_col = f'structural_agreement_d_{tag}'
        if n_col not in df.columns or d_col not in df.columns:
            print(f"  [SKIP] {tag}: columns not present")
            continue
        n = df[n_col].fillna(0).sum()
        d = df[d_col].fillna(0).sum()
        score = float(n / d) if d > 0 else 0.0
        ok = check(f"structural_agreement_{tag}", score,
                   EXPECTED_V5['structural_agreement'].get(tag))
        results.append(ok)

    # mech_consistency
    print("\nmech_consistency at all 5 thresholds:")
    for tag in ('t10', 't15', 't20', 't25', 'tSAP'):
        col = f'mech_consistency_{tag}'
        if col not in df.columns:
            print(f"  [SKIP] {tag}: column not present")
            continue
        vc = df[col].value_counts()
        n = len(df) - df[col].isna().sum()
        score = float((vc.get('consistent', 0) + 0.5 * vc.get('partial', 0)) / n) if n else 0.0
        ok = check(f"mech_consistency_{tag}", score,
                   EXPECTED_V5['mech_consistency'].get(tag))
        results.append(ok)

    # threshold-stable
    if 'mech_consistency_threshold_stable' in df.columns:
        n_stable = int(df['mech_consistency_threshold_stable'].sum())
        ok = check("threshold_stable_count", n_stable,
                   EXPECTED_V5['threshold_stable'], tol=2)  # allow ±2 (small variation OK)
        results.append(ok)

    # Spot-check worked example: PIK3CA E545K
    print("\nWorked example: PIK3CA E545K")
    sub = df[(df['gene'].str.lower() == 'pik3ca') & (df['variant'] == 'E545K')]
    if len(sub) > 0:
        r = sub.iloc[0]
        for tag in ('t10', 'tSAP'):
            mech = r[f'mech_{tag}']
            grade = r[f'mech_consistency_{tag}']
            print(f"  {tag}: mech='{mech}', consistency={grade}")

    return results


def verify_evaluation(eval_dir):
    """Verify Track A headlines."""
    print("\n" + "=" * 70)
    print("Track A verification (run_evaluate)")
    print("=" * 70)
    results = []

    # Level 1
    l1_csv = eval_dir / 'level1_merged.csv'
    if l1_csv.exists():
        l1 = pd.read_csv(l1_csv)
        full = l1[l1['classifier'] == 'MAVIS_full']
        mono = l1[l1['classifier'] == 'monomer_only']
        if len(full) > 0:
            ok = check("Level 1 MAVIS_full TPR @ t=1.0",
                       float(full.iloc[0]['TPR_t10']),
                       EXPECTED_V5['level1_TPR_t10_MAVIS_full'], tol=0.005)
            results.append(ok)
        if len(mono) > 0:
            ok = check("Level 1 monomer_only TPR @ t=1.0",
                       float(mono.iloc[0]['TPR_t10']),
                       EXPECTED_V5['level1_TPR_t10_monomer_only'], tol=0.005)
            results.append(ok)

    # Level 3 HBB
    l3_csv = eval_dir / 'level3_hbb_quantitative.csv'
    if l3_csv.exists():
        l3 = pd.read_csv(l3_csv)
        if len(l3) > 0:
            ok = check("Level 3 HBB Pearson r",
                       float(l3.iloc[0]['pearson_r']),
                       EXPECTED_V5['level3_HBB_pearson_r'], tol=0.01)
            results.append(ok)

    # Level 2 per-phenotype detection
    l2_csv = eval_dir / 'level2_merged_detection.csv'
    if l2_csv.exists():
        l2 = pd.read_csv(l2_csv)
        for phen, key in [
            ('pathogenic', 'level2_pathogenic_detection_t10'),
            ('pathogenic_gof', 'level2_pathogenic_gof_detection_t10'),
        ]:
            sub = l2[l2['phenotype'] == phen]
            if len(sub) > 0:
                actual = int(sub.iloc[0]['n_detected_t10'])
                ok = check(f"Level 2 {phen} detected @ t=1.0",
                           actual, EXPECTED_V5[key], tol=1)
                results.append(ok)

    return results


def main():
    parser = argparse.ArgumentParser(description="Verify v5 framework downstream")
    parser.add_argument('--intermediate', required=True,
                        help='Path to mavis_v7_results_with_nbhd.csv')
    parser.add_argument('--am', required=True,
                        help='Path to AM xlsx')
    parser.add_argument('--scripts-dir', required=True,
                        help='Path to PIPELINE_CURRENT/scripts/')
    parser.add_argument('--corrected',
                        help='Path to mavis_v7_results_corrected.csv (for Track A)')
    parser.add_argument('--output-dir', default='./verification_output',
                        help='Where to write verification artifacts')
    args = parser.parse_args()

    intermediate = Path(args.intermediate).resolve()
    am_xlsx = Path(args.am).resolve()
    scripts_dir = Path(args.scripts_dir).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(exist_ok=True, parents=True)

    print(f"Verifying v5 framework using:")
    print(f"  intermediate: {intermediate}")
    print(f"  AM xlsx:      {am_xlsx}")
    print(f"  scripts:      {scripts_dir}")
    print(f"  output:       {output_dir}")

    # Track B
    concordance_csv = run_concordance(intermediate, output_dir, scripts_dir, am_xlsx)
    track_b_results = verify_concordance(concordance_csv)

    # Track A — needs the corrected CSV (output of baseline_correct, not concordance)
    track_a_results = []
    if args.corrected:
        eval_dir = run_evaluation(Path(args.corrected).resolve(), output_dir, scripts_dir)
        track_a_results = verify_evaluation(eval_dir)
    else:
        print("\n[SKIP Track A] --corrected not provided")

    # Summary
    all_results = [r for r in (track_b_results + track_a_results) if r is not None]
    n_pass = sum(1 for r in all_results if r is True)
    n_fail = sum(1 for r in all_results if r is False)
    print("\n" + "=" * 70)
    print(f"Verification summary: {n_pass}/{n_pass + n_fail} passed")
    print("=" * 70)
    if n_fail > 0:
        print(f"\n[FAIL] {n_fail} checks failed. Review output above.")
        sys.exit(1)
    print("\n[PASS] All v5 framework checks passed.")
    sys.exit(0)


if __name__ == '__main__':
    main()
