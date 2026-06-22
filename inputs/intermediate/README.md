# Intermediate cached files

These files are produced by the pipeline and cached here for fast downstream
verification. A fresh `run_live.py` execution overwrites them.

## Files
- **mavis_v7_results.csv** — output of Stages 1-6 (modular `pipeline.run_pipeline`).
  Contains structural features, FoldX three-axis ΔΔG values, P1 scoring,
  mavis_tier, mavis_mechanism. ~262 columns × 44 rows.
- **mavis_v7_results_corrected.csv** — output of Stage 7 (`baseline_correct.py`).
  Adds per-(system, partner) baseline-indistinguishable flags + reclassified
  mechanism. **Track A consumes this.**
- **mavis_v7_results_with_nbhd.csv** — output of Stage 8 (`neighborhood.py`).
  Adds nbhd_* columns: ±3 residue window aggregation. **Track B consumes this.**

## How they're regenerated
1. Stage 1-6: `python3 scripts/run_live.py` (~24-48h with FoldX)
2. Stage 7:   `python3 scripts/mavis_v7_baseline_correct.py` (~5 sec)
3. Stage 8:   `python3 scripts/mavis_v7_neighborhood.py` (~30 sec, needs structures)

## Cached versions provenance
- These cached files were last regenerated on Apr 22 2026 by Luke
- Pre-Phase-2 patches (the multi-threshold mechanism integration we'll add now)
- After applying the patches, these files will get NEW columns: `mavis_mechanism_t10/t15/t20/t25` and `mavis_mechanism_corrected_t10/t15/t20/t25`
