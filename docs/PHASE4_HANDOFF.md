# MAVIS v7 — Phase 4 Hand-Off

**Date:** 2026-05-05
**State at hand-off:** v5 framework fully implemented and verified end-to-end against cached intermediate. All 16 verification checks pass. Ready for full FoldX rerun on Luke's machine.

---

## What's in PIPELINE_CURRENT/

### Modular package patches (committed)
- `scripts/mavis_v7/mechanism.py` — `classify_mechanism` accepts per-axis threshold dict; 33-category fold split; `ALL_MECHANISMS` updated; `_resolve_thresholds()` helper
- `scripts/mavis_v7/pipeline.py` — writes `mavis_mechanism_t{10,15,20,25,SAP}` plus backward-compat `mavis_mechanism` alias
- `scripts/mavis_v7/evaluation.py` — Track A v5: LoF/GoF mapping dropped; `level2_per_phenotype_detection` replaces `level2_three_class`; Level 4 reframed
- `scripts/mavis_v7_baseline_correct.py` — `reclassify_with_flags` writes `mavis_mechanism_corrected_t{10,15,20,25,SAP}` plus alias

### Pipeline scripts
- `scripts/run_live.py` — full Phase 1–3 runner (8 stages: preprocessing → variant loading → metrics → FoldX → P1 scoring → baseline correction → neighborhood)
- `scripts/mavis_v7_neighborhood.py` — Pipeline 2 (±3 neighborhood)
- `scripts/apply_concordance_v5.py` — Track B v5 concordance assembly with all v5 features:
  - 5-threshold mechanism sweep including Sapozhnikov per-axis
  - Symmetric internal CI95 gating
  - v13 P2 bugfix (DDG term removed from neighborhood score)
  - v14 tier × DDG concordance columns
  - structural_agreement at all 5 thresholds
  - directional_agreement at all 5 thresholds (supplementary)
  - tier_structural_signal_type diagnostic
  - evaluation_note divergence flag
  - threshold-stable flag
  - per-axis vote columns
  - Rubric B fold-subtype-aware mech_consistency grading
- `scripts/run_evaluate.py` — Track A v5 evaluator with --input/--output-dir argparse
- `scripts/build_report.py` — 9-sheet xlsx report builder (unchanged from v4)

### Verification
- `verification/verify_stage6.py` — runs apply_concordance_v5 + run_evaluate against cached intermediate and checks 16 expected v5 headline values

### Documentation
- `documentation/methods_metrics_sketch_v4.md` — populated, verified methods-paper sketch with all v5 numbers
- `documentation/design_decisions.md` — comprehensive v5 design provenance (16 sections)
- `documentation/RUN_PIPELINE.md` — runtime instructions (this file's sibling)

---

## How to run v5 end-to-end on a fresh machine

### Prerequisites
- Python 3.14 virtualenv at `~/mavis_v7/.venv/` with pandas, numpy, openpyxl, biopython, scikit-learn
- FoldX 5.1 binary at `/Users/lukearnce/Downloads/foldx5_Mac_0/foldx_20270131` (path set in `mavis_v7/config.py`; edit if different)
- Structure files in `inputs/raw/structures/` (AlphaFold monomer .pdb files + complex multimer PDBs)

### Phase 4 step-by-step

```bash
cd ~/mavis_v7
source .venv/bin/activate

# Stage 1-7: full pipeline (preprocessing → FoldX → P1 scoring → baseline correction)
# This is the long-running step (5 replicates × 44 variants × N partners × FoldX BuildModel+AnalyseComplex)
python scripts/run_live.py \
  --variants inputs/raw/benchmark_variants_v5.csv \
  --structures-dir inputs/raw/structures/ \
  --output-dir outputs/track_b/

# Stage 7b: neighborhood scoring (Pipeline 2)
python scripts/mavis_v7_neighborhood.py \
  --input outputs/track_b/mavis_v7_results_corrected.csv \
  --output outputs/track_b/mavis_v7_results_with_nbhd.csv

# Track B: apply v5 concordance + diagnostics
python scripts/apply_concordance_v5.py \
  --results outputs/track_b/mavis_v7_results_with_nbhd.csv \
  --external inputs/AM_variants_mavis_mechanism_test.xlsx \
  --outdir outputs/track_b/

# Track A: run binary classifier metrics + HBB + per-phenotype detection
python scripts/run_evaluate.py \
  --input outputs/track_b/mavis_v7_results_corrected.csv \
  --output-dir outputs/track_a/

# Build 9-sheet xlsx report
python scripts/build_report.py \
  --concordance outputs/track_b/mavis_v7_concordance.csv \
  --output outputs/track_b/mavis_v7_concordance_v5.xlsx

# OPTIONAL: verify against expected v5 headlines
python verification/verify_stage6.py \
  --intermediate outputs/track_b/mavis_v7_results_with_nbhd.csv \
  --corrected outputs/track_b/mavis_v7_results_corrected.csv \
  --am inputs/AM_variants_mavis_mechanism_test.xlsx \
  --scripts-dir scripts \
  --output-dir verification_run/
```

### Expected wall-clock time
- run_live.py with 5 FoldX replicates on 44 variants × ~2 partners average: roughly 4–6 hours on a single-core MacBook (FoldX is single-threaded per variant)
- Neighborhood scoring: ~2 minutes
- apply_concordance_v5 + run_evaluate + build_report: <30 seconds combined
- verify_stage6: <30 seconds

---

## Expected output values (sanity check)

If the verify script runs to completion, you should see:

```
Track B verification (apply_concordance_v5)
  [ OK ] structural_agreement_t10: 0.658 (expected 0.658)
  [ OK ] structural_agreement_t15: 0.718 (expected 0.718)
  [ OK ] structural_agreement_t20: 0.709 (expected 0.709)
  [ OK ] structural_agreement_t25: 0.718 (expected 0.718)
  [ OK ] structural_agreement_tSAP: 0.709 (expected 0.709)
  [ OK ] mech_consistency_t10: 0.602 (expected 0.602)
  [ OK ] mech_consistency_t15: 0.693 (expected 0.693)
  [ OK ] mech_consistency_t20: 0.716 (expected 0.716)
  [ OK ] mech_consistency_t25: 0.761 (expected 0.761)
  [ OK ] mech_consistency_tSAP: 0.750 (expected 0.750)
  [ OK ] threshold_stable_count: 28 (expected 28)

Track A verification (run_evaluate)
  [ OK ] Level 1 MAVIS_full TPR @ t=1.0: 0.913 (expected 0.913)
  [ OK ] Level 1 monomer_only TPR @ t=1.0: 0.391 (expected 0.391)
  [ OK ] Level 3 HBB Pearson r: 0.894 (expected 0.890)
  [ OK ] Level 2 pathogenic detected @ t=1.0: 24 (expected 24)
  [ OK ] Level 2 pathogenic_gof detected @ t=1.0: 5 (expected 5)

Verification summary: 16/16 passed
[PASS] All v5 framework checks passed.
```

**If numbers diverge after a fresh FoldX rerun:** this is expected if Luke updates any of (a) annotation CSV, (b) structure files, (c) FoldX repair settings. The v5 numbers reproduce against the *current cached intermediate*; a full rerun on the same inputs should land at the same headlines within ±0.01 due to FoldX stochasticity.

---

## What might shift after fresh FoldX rerun

The cached intermediate `mavis_v7_results_with_nbhd.csv` reflects a prior FoldX run. A fresh rerun may produce slightly different ΔΔG values due to:

1. **FoldX stochasticity in BuildModel.** Despite fixed seed, the 5-replicate mean and SD can vary by ~0.05–0.1 kcal/mol between runs.
2. **RepairPDB differences** if any structure files changed.
3. **AlphaFold updates** if any AF DB structures were re-downloaded.

If headlines shift by more than ±0.02, investigate before publishing. The verify script's tolerance is ±0.005 for SA/MC; larger drift signals an input change worth understanding.

---

## Open items / future work

**Pending decisions and work for future sessions:**

1. **KRAS G12D/G12V annotation review.** Hunter 2015 Figure 3 explicitly classifies these as "low RAF affinity" — contradicting current "WT-like" annotations. Pending Luke pulling the figure with institutional access.

2. **Phase 2 literature verification batches 4-12.** Three of 12 paper batches were verified (Kosinski 2010, Kiger 1998, Hunter 2015). Nine batches remain.

3. **CHD validation.** Once methods paper is complete, deploy v5 framework on the 163-variant CHD candidate gene set.

4. **Benchmark expansion.** Target ~60 variants before publication to strengthen per-class statistics.

5. **v6 code unification.** Track A and Track B remain in separate scripts in v5; v6 can dissolve them into one unified evaluation framework. Pure refactor, no methodology change.

6. **Methods paper drafting.** Source: `methods_metrics_sketch_v4.md`. Headline framing locked: paired claim (SA 0.71 / MC 0.75 at Sapozhnikov).

7. **mild_stab annotation category.** Mentioned in earlier sessions as a possible methodological addition; no current annotations affected; deferred.

---

## What this hand-off does NOT cover

- **Structure curation** — outside this scope; assumed inputs are valid AF/PDB structures
- **Variant annotation accuracy** — assumed correct per Phase 2 verification (3/12 batches done; 9 remaining)
- **Manuscript writing** — separate effort; this hand-off only covers the computational framework

---

## Files presented to Luke this session

- `documentation/methods_metrics_sketch_v4.md` — methods paper outline with verified Rubric B numbers
- `documentation/design_decisions.md` — comprehensive provenance documentation
- `documentation/PHASE4_HANDOFF.md` — this file
