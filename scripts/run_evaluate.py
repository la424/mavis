"""
MAVIS v7 — Phase 3 evaluation runner.

Runs the four-level benchmark against the corrected results CSV. Primary
threshold is 1.0 (matches pipeline production); results at 1.0, 1.5, 2.0
reported side-by-side.

Outputs:
  - Merged per-level tables (paper-facing)
  - Per-threshold tables (supplemental / programmatic)
  - Level 4 per-variant GoF detection table

Usage:
  cd ~/mavis_v7
  python3 run_evaluate.py
"""
from pathlib import Path
import sys
import pandas as pd

from mavis_v7.evaluation import run_full_evaluation


def print_hr(title: str, char: str = "=", width: int = 88):
    print()
    print(char * width)
    print(title)
    print(char * width)


def main():
    import argparse
    parser = argparse.ArgumentParser(description="Run Track A (Phase 3) evaluation")
    parser.add_argument('--input', default='results/mavis_v7_results_corrected.csv',
                        help='Path to mavis_v7_results_corrected.csv (output of baseline_correct)')
    parser.add_argument('--output-dir', default=None,
                        help='Where to write per-level CSVs (default: <input dir>/evaluation/)')
    args = parser.parse_args()

    results_path = Path(args.input)
    if not results_path.exists():
        print(f"✗ Corrected results CSV not found: {results_path}")
        print("  Run mavis_v7_baseline_correct.py first.")
        sys.exit(1)

    if args.output_dir is None:
        out_dir = results_path.parent / "evaluation"
    else:
        out_dir = Path(args.output_dir) / "evaluation"

    df = pd.read_csv(results_path)
    print(f"Loaded {len(df)} variants from {results_path}")
    print(f"  Roles: {df['role'].value_counts().to_dict()}")
    print(f"  Phenotypes: {df['phenotype'].value_counts().to_dict()}")

    results = run_full_evaluation(df, thresholds=[1.0, 1.5, 2.0], n_bootstrap=1000)

    # --------------------------------------------------------------------
    # Level 1
    # --------------------------------------------------------------------
    print_hr("LEVEL 1 — Binary pathogenic vs benign (multi-threshold view)")
    print("  Primary threshold: 1.0 (MAVIS production, matches CHD paper)")
    print("  Mechanism controls EXCLUDED; pLDDT gate >= 50; bootstrap 95% CI")
    print()
    with pd.option_context("display.max_columns", None, "display.width", 220):
        print(results["level1_merged"].to_string(index=False))

    print_hr("LEVEL 1 — Per-threshold detail", char="-")
    for t, tbl in results["level1_per_threshold"].items():
        print(f"\n  --- threshold {t} ---")
        print(tbl.to_string(index=False))

    print_hr("LEVEL 1 — Sensitivity: mechanism controls included (threshold 1.0)", char="-")
    print(results["level1_sensitivity_incl_controls"].to_string(index=False))

    print_hr("LEVEL 1 — Sensitivity: strict pLDDT >= 70 (threshold 1.0)", char="-")
    print(results["level1_sensitivity_strict_plddt"].to_string(index=False))

    # --------------------------------------------------------------------
    # Level 2
    # --------------------------------------------------------------------
    print_hr("LEVEL 2 — Per-phenotype detection rate (v5: replaces 3-class confusion)")
    print()
    print("v5 reports detection rate per phenotype rather than LoF/GoF classifier")
    print("recall — see methods §9.1. 'Detection' = high tier OR any DDG axis fires.")
    print()
    print(results["level2_merged_detection"].to_string(index=False))

    # --------------------------------------------------------------------
    # Level 3
    # --------------------------------------------------------------------
    print_hr("LEVEL 3 — Per-axis matching (merged)")
    print()
    print(results["level3_axis_merged"].to_string(index=False))

    print_hr("LEVEL 3 — HBB W37 quantitative recovery (threshold-independent)")
    hbb = results["level3_hbb_quantitative"]
    print(f"  n variants:    {hbb['n']}")
    print(f"  Pearson r:     {hbb['pearson_r']}")
    print(f"  Spearman rho:  {hbb['spearman_rho']}")
    print(f"  MAE:           {hbb['mae']} kcal/mol")
    print(f"  RMSE:          {hbb['rmse']} kcal/mol")
    if hbb.get("variants"):
        print("  Per-variant predictions:")
        for v, obs in hbb["variants"].items():
            print(f"    {v}: {obs:.3f}")

    # --------------------------------------------------------------------
    # Level 4 — Structural disruption detection in pathogenic_gof variants
    # (v5: reframed from "GoF mechanism detection" — same numerator, no LoF/GoF mapping)
    # --------------------------------------------------------------------
    print_hr("LEVEL 4 — Structural disruption detected in pathogenic_gof variants")
    print(f"  All {len(results['level4_per_variant'])} pathogenic_gof-phenotype variants")
    print()
    print("  Summary by classifier × threshold:")
    print(results["level4_summary"].to_string(index=False))

    print_hr("LEVEL 4 — Per-variant detection", char="-")
    pv = results["level4_per_variant"]
    display_cols = ["system", "variant", "role", "phenotype", "stabilization_signals",
                    "mavis_positive_t10", "monomer_positive_t10", "struct_positive_t10",
                    "structural_disruption_detected_t10"]
    if "mavis_mechanism_corrected_t10" in pv.columns:
        display_cols.append("mavis_mechanism_corrected_t10")
    with pd.option_context("display.max_columns", None, "display.width", 200,
                            "display.max_colwidth", 50):
        print(pv[display_cols].to_string(index=False))

    # --------------------------------------------------------------------
    # Write CSV outputs
    # --------------------------------------------------------------------
    out_dir.mkdir(parents=True, exist_ok=True)

    # Level 1
    results["level1_merged"].to_csv(out_dir / "level1_merged.csv", index=False)
    for t, tbl in results["level1_per_threshold"].items():
        tkey = f"t{int(t*10):02d}"
        tbl.to_csv(out_dir / f"level1_{tkey}.csv", index=False)
    results["level1_sensitivity_incl_controls"].to_csv(
        out_dir / "level1_sensitivity_incl_controls.csv", index=False)
    results["level1_sensitivity_strict_plddt"].to_csv(
        out_dir / "level1_sensitivity_strict_plddt.csv", index=False)

    # Level 2 (v5: per-phenotype detection rates, replaces the legacy
    # three-class confusion matrix that depended on the LoF/GoF mapping)
    results["level2_merged_detection"].to_csv(
        out_dir / "level2_merged_detection.csv", index=False)
    for t, det in results["level2_per_phenotype_detection"].items():
        tkey = f"t{int(t*10):02d}"
        det.to_csv(out_dir / f"level2_detection_{tkey}.csv", index=False)

    # Level 3
    results["level3_axis_merged"].to_csv(out_dir / "level3_axis_merged.csv", index=False)
    for t, ax in results["level3_axis_per_threshold"].items():
        tkey = f"t{int(t*10):02d}"
        ax.to_csv(out_dir / f"level3_axis_{tkey}.csv", index=False)
    pd.DataFrame([results["level3_hbb_quantitative"]]).to_csv(
        out_dir / "level3_hbb_quantitative.csv", index=False)

    # Level 4
    results["level4_per_variant"].to_csv(out_dir / "level4_per_variant.csv", index=False)
    results["level4_summary"].to_csv(out_dir / "level4_summary.csv", index=False)

    print_hr("Wrote result tables")
    print(f"  {out_dir}/")
    for f in sorted(out_dir.glob("*.csv")):
        print(f"    {f.name}")


if __name__ == "__main__":
    main()
