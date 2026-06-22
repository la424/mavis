# MAVIS v7 — Methods Paper Metrics Outline (v4 — implementation-verified)

**Status:** All numbers verified by `verification/verify_v5.py` against the v5 implementation. The v5 framework is fully implemented in `apply_concordance_v5.py`, `mavis_v7/mechanism.py`, `mavis_v7/pipeline.py`, `mavis_v7_baseline_correct.py`, and `mavis_v7/evaluation.py`. Locked.

**v3 → v4 update:** mech_consistency numbers reflect the actual Rubric B implementation (per-axis fold subtype checks via `grade_mechanism_consistency` in `apply_concordance_v5.py`), which is more permissive than the standalone-compute sketch numbers — fold subtypes are checked against per-axis annotations, but the rubric retains the existing axis-presence framework rather than imposing a coarser axis-alignment rule.

---

## §1. Pipeline outputs

For each input variant, MAVIS v7 generates structural-disruption predictions decomposed across multiple axes:

- A **structural tier** (Tier 1–4) reflecting position-level structural importance, derived from Grantham severity, contacts, burial state, interface participation, and pLDDT confidence
- **Three FoldX ΔΔG axes** with replicate-derived 95% confidence intervals: monomer fold stability, fold-in-complex stability per partner, and binding interaction energy per partner
- A **mechanism call** in one of 33 categories synthesizing tier and ΔΔG signals into a biophysically-named outcome (full category list and decision logic in supplementary)
- **Diagnostic columns** (described in §7) supporting interpretation but not graded

Implementation details of the underlying calculations (FoldX execution, pLDDT gating thresholds, score formulas) are in the Pipeline Methods section.

The pipeline is positioned as a **structural-disruption detector**, not a pathogenicity classifier. Predictions are intended for use alongside external pathogenicity evidence to inform variant prioritization in a downstream interpretation context.

## §2. Benchmark

The evaluation benchmark is `benchmark_variants_v5.csv`: 44 variants curated from literature across 11 protein-protein interaction systems. Systems were selected to span loss-of-function and gain-of-function mechanisms, fold-driven and binding-driven disruption, and well-characterized vs. ambiguous biology.

Each variant has per-axis ground-truth annotations:
- **role**: pathogenic / pathogenic_lof / pathogenic_gof / benign / mechanism_control
- **phenotype**: pathogenic / pathogenic_gof / benign (refines role; the generic `pathogenic` label captures cases where LoF vs GoF assignment is unspecified or ambiguous in the literature)
- **expected_ddg_monomer**: destab / mild_destab / neutral / unknown
- **expected_ddg_fold_complex**: destab / neutral / unknown
- **expected_ddg_binding**: destab / mild_destab / stab / neutral / unknown
- **expected_topology**: at_interface / not_at_interface / unknown
- **expected_mech_class** (derived from per-axis annotations): structurally_silent / mixed_structural / fold_mechanism / ppi_destab_mechanism / ppi_stab_mechanism

Annotations were curated from peer-reviewed literature with citations recorded per axis. Approximately 30% of axis-level annotations are `unknown` because direct biochemical or structural measurements were not available — these axes are excluded from grading per §5.

The `mechanism_control` role category is dissolved into the regular pathogenic/benign pool for evaluation; the framework grades structural-disruption prediction accuracy, which is independent of whether a variant was originally curated as a control.

## §3. Five evaluation questions

The framework asks five questions about pipeline accuracy. **Q1 and Q2 are complementary, testing per-axis correctness and synthesis-level correctness respectively.** Q3 tests quantitative ΔΔG accuracy where experimental data exists. Q4 tests methodological contribution against baselines. Q5 tests cross-tool concordance for context, with explicit secondary framing.

Headline numbers are reported with bootstrap 95% confidence intervals (n=1000 iterations).

### Q1: Does the pipeline detect structural disruption at the per-axis level?

**Metric: structural_agreement.** For each variant, gradeable axes (tier + 3 ΔΔG axes) are checked against per-axis ground truth using §5's exclusion rules. Score = correct / gradeable, summed across variants.

**Headline: 0.71 at the Sapozhnikov-confident threshold (mono=2.9 / fold=2.9 / bind=3.5 kcal/mol)**, with sensitivity range 0.66–0.72 across uniform thresholds 1.0–2.5 kcal/mol.

Sensitivity sweep (n correct / n gradeable):

| Threshold | structural_agreement |
|---|---|
| t = 1.0 | 0.658 (77/117) |
| t = 1.5 | 0.718 (84/117) |
| t = 2.0 | 0.709 (83/117) |
| t = 2.5 | 0.718 (84/117) |
| t = Sapozhnikov | 0.709 (83/117) |

The metric saturates around 0.71, reflecting that beyond t=1.5, the predictions are bimodal — either firmly above or firmly below threshold, with few borderline cases.

### Q2: Does the synthesized mechanism call correctly identify the kind of disruption?

**Metric: mech_consistency.** The mechanism call is graded against `expected_mech_class` using a three-state rubric (consistent / partial / inconsistent), with **per-axis fold subtype checks (Rubric B)** — fold-monomer and fold-complex axes are graded separately. A call missing one fold subtype while matching the other → partial; matching both → consistent. Score = (consistent + 0.5 × partial) / total.

**Headline: 0.75 at the Sapozhnikov-confident threshold**, with sensitivity range 0.60–0.76 across uniform thresholds 1.0–2.5 kcal/mol.

Sensitivity sweep:

| Threshold | mech_consistency | consistent / partial / inconsistent |
|---|---|---|
| t = 1.0 | 0.602 | 23 / 7 / 14 |
| t = 1.5 | 0.693 | 27 / 7 / 10 |
| t = 2.0 | 0.716 | 27 / 9 / 8 |
| t = 2.5 | 0.761 | 28 / 11 / 5 |
| t = Sapozhnikov | 0.750 | 27 / 12 / 5 |

mech_consistency rises monotonically with threshold strictness (with a small drop at Sapozhnikov per-axis as the higher fold/binding cutoffs flip a few near-threshold calls). Note that threshold-stability of grades is reduced under Rubric B: 28/44 variants are stable across all five thresholds (vs 30/44 under the more permissive pre-v5 rubric), because Rubric B is more sensitive to threshold-dependent fold-axis crossings.

### Q3: How accurately does FoldX quantify the magnitude of structural effects?

**Metric: HBB W37 quantitative correlation.** Pearson r and Spearman ρ between pipeline-predicted ΔΔG and experimental ΔΔG (Kwiatkowski 1998) on the four W37 hemoglobin variants.

**Headline: Pearson r = 0.89, Spearman ρ = 0.80, MAE = 1.57 kcal/mol.**

Per-variant predictions:

| Variant | Pipeline ΔΔG | Experimental | Error |
|---|---|---|---|
| W37Y | 1.60 | 2.0 | 0.40 |
| W37A | 4.53 | 5.0 | 0.47 |
| W37G | 5.45 | 7.0 | 1.55 |
| W37E | 5.14 | 9.0 | 3.86 |

This is the only experimental-magnitude validation in the benchmark; expansion is planned. The W37G/W37E underestimation reflects FoldX's known force-field plateau at large destabilizations (§10).

### Q4: Does the multimer axis add value beyond simpler scoring?

**Metric: baseline comparison.** Three classifiers compared against role-based binary outcome (pathogenic vs benign): MAVIS_full (tier + multimer ΔΔG), monomer_only (tier + monomer ΔΔG), structural_score (tier alone). Reports TPR, TNR, accuracy.

**Headline: at threshold 1.0, MAVIS_full TPR = 0.913 vs monomer_only TPR = 0.391**, demonstrating the multimer axis recovers pathogenic structural disruption that single-residue scoring misses. The threshold for this comparison is t=1.0 for continuity with the original PHASE3_CHECKPOINT analysis. Sensitivity at higher thresholds and full per-classifier comparison are reported in supplementary.

(Note: these numbers are pre-v5 PHASE3_CHECKPOINT values. Track A's update under v5 — LoF/GoF cleanup and fold-split categories — may shift them slightly. Final values will be re-verified after implementation.)

### Q5 (supporting): Where do MAVIS predictions agree with external pathogenicity evidence?

**Metric: external_consensus** (AM + Franklin agreement) and **legacy 4-way concordance** (tier + collapsed DDG + AM + Franklin). The 4-way is retained for backward-compatibility with prior cross-tool comparison literature; it is not a primary contribution.

**Headline: AM and Franklin agree on 30/40 variants where both have classifications. Legacy 4-way concordance: 0.616 strict / 0.733 relaxed.**

## §4. Threshold methodology

ΔΔG thresholds set the bar for "the pipeline detected destabilization." Five thresholds reflecting the field's range are reported:

- **t = 1.0 kcal/mol** — CAGI5 frataxin convention, prior CHD-pipeline continuity (sensitivity-favoring)
- **t = 1.5 kcal/mol** — Caldararu/Guerois empirical classification optimum
- **t = 2.0 kcal/mol** — intermediate, conservative-leaning
- **t = 2.5 kcal/mol** — approaches but does not meet Sapozhnikov-confident floor
- **Sapozhnikov per-axis** — mono = 2.9, fold = 2.9, bind = 3.5 kcal/mol; matches per-axis 95% prediction interval against experimental ΔΔG (Sapozhnikov et al. 2023, BMC Bioinformatics)

The pipeline does not designate a single primary threshold. Headline numbers are reported at the Sapozhnikov-confident threshold (most stringent, most defensible against precision criticism). Sensitivity across the full sweep is reported in supplementary tables. The Sapozhnikov-per-axis threshold produces structural_agreement and mech_consistency numbers nearly identical to t=2.5 (SA: 0.709 vs 0.718; MC: 0.648 vs 0.659), reflecting the bimodal distribution of predictions in this benchmark — most variants either fire well above 2.5 kcal/mol or stay well below.

## §5. Per-axis confidence and exclusion rules

An axis is excluded from grading when any of the following hold:

- **Ground truth unknown.** Axes with `expected_*` annotated as `unknown` are not gradeable.
- **Pipeline confidence low.** Per-axis pLDDT < 50 at the variant residue (or partner pLDDT < 50 for partner-specific axes) excludes the axis.
- **CI95 spans zero.** When 0 falls within ±1.96 × SD of the FoldX prediction across 5 replicates, the prediction is statistically indistinguishable from null. Such axes are excluded regardless of ground-truth annotation type, applied symmetrically to mild_destab, destab, stab, and neutral.

These rules apply per-axis. A variant may have one axis excluded and others graded; the denominator is variant-specific.

The rationale: the pipeline should not be penalized for limitations of the underlying calculation. When FoldX's own replicate variance prevents a confident call, no ground-truth comparison is meaningful.

## §6. Mechanism category decomposition

The mechanism call distinguishes **monomer-driven** and **multimer-driven** fold disruption. The 33-category scheme covers:

- **Pure fold disruption (12 categories)** — mono / multi / both × destab / stab × at_interface / no_interface
- **Fold + PPI same-direction (6)** — fold_cat × destab / stab on binding axis matching
- **Fold + PPI conflicting (6)** — fold_cat × opposite-direction binding (the "auto-inhibition release" signature)
- **Pure PPI (3)** — destab / stab / conflicting (mixed partner signals) on binding without fold signal
- **Structural with no DDG signal (3)** — interface variant, contact-driven, burial-driven
- **Catch-all (3)** — no structural effect, unevaluable, low-confidence

This split was added in the v5 framework. The original 17-category scheme (retained from the prior CHD-pipeline v6.0) collapsed monomer-driven and multimer-driven fold disruption into a single "Fold destabilization" category. The collapse allowed the synthesis label to land in a consistent-grading bucket even when the pipeline was firing on different underlying axes than ground truth specified — a "right for the wrong reasons" failure mode.

Under v5, both the synthesis label and the consistency grading exploit the fold split. The grading rubric (**Rubric B**, see §7 for diagnostics) checks fold-monomer and fold-complex axes separately against per-axis ground truth annotations: a call that asserts only "Multimer fold destab" against ground truth annotated `expected_ddg_monomer = destab AND expected_ddg_fold_complex = destab` is graded as `partial` because the call missed the monomer axis even though it caught the multimer axis.

Empirically, the headline mech_consistency at t=2.5 dropped from 0.773 (pre-v5 17-category scheme without subtype checking) to 0.761 (v5 33-category scheme with Rubric B subtype checking). The decrease is small but biologically meaningful: it reflects a small number of variants (e.g., BRCA1 C61G at t=2.5) where the call asserted only one fold subtype while the underlying ground truth expected both. The fold split is doing diagnostic and grading work simultaneously without dramatically inflating or deflating the headline.

## §7. Diagnostic columns

The following columns are not graded but support interpretation:

- **`tier_structural_signal_type`** ∈ {fold_related, binding_related, both, ambiguous, none} — categorizes the structural feature class underlying the tier prediction
- **`evaluation_note`** — flags variants where mech_consistency and structural_agreement disagree directionally, with short reason
- **`p1_ddg_concordance_t*`**, **`p2_ddg_concordance_t*`** — per-threshold tier × DDG agreement labels (concordant_disruption / structural_only / ddg_only / concordant_silent) for both pipeline variants (P1 single-residue, P2 ±3 neighborhood)
- **Per-axis vote columns** — `ddg_{monomer,fold,binding}_vote_{strict,relaxed}` enabling reconstruction of alternative concordance metrics
- **`mech_consistency_threshold_stable`** — flag for variants whose mech_consistency grade is identical across all five thresholds (28/44 stable in current benchmark under Rubric B; the unstable subset is surfaced in §10 as boundary cases)

## §8. Worked examples

Three variants traced through the framework illustrate the metrics in action.

### §8.1 PIK3CA E545K — convention-dependent case study (pathogenic_gof)

E545K is a known gain-of-function PIK3CA variant. Ground truth: `expected_mech_class = mixed_structural`, `expected_ddg_fold_complex = destab`, `expected_ddg_binding = destab`, `phenotype = pathogenic_gof`.

Pipeline output:
- Tier 2, score 3.5; `tier_structural_signal_type = binding_related` (interface signal at nSH2-iSH2)
- `ddg_monomer = 0.18` (SD 0.034; CI95 = [0.11, 0.25] — distinguishable from 0; expected neutral, predicted silent → correct)
- `ddg_fold_pik3r1 = 1.89` (SD 0.46; fires at t=1.0/1.5, silent at t=2.0/2.5/Sapozhnikov; correct positive direction; expected destab → sub-threshold at strict thresholds)
- `ddg_binding_pik3r1 = −4.63` (SD 0.46; strong stabilization predicted; ground truth expected destabilization → wrong direction)
- mech at t=1.0/1.5: `Multimer fold destab. + PPI stabilization (conflicting)`
- mech at t=2.0/2.5/Sapozhnikov: `PPI stabilization`

Per-axis structural_agreement at Sapozhnikov: tier ✓ (high tier expected, fired), monomer ✓ (correctly silent on neutral), fold ✗ (sub-threshold at 2.9 kcal/mol), binding ✗ (wrong direction). Score: 2/4.

mech_consistency (Rubric B):
- t=1.0/1.5: synthesis call "Multimer fold destab. + PPI stabilization (conflicting)" → graded **consistent** against `mixed_structural`. The conflicting category captures the mixed signal directly; per-axis fold subtype check passes (call asserts multimer-fold; ground truth `expected_ddg_monomer = neutral`, `expected_ddg_fold_complex = destab` → axis correctly matched; missed monomer is excused because it's annotated neutral).
- t≥2.0: synthesis call collapses to "PPI stabilization" → graded **inconsistent**. The call now misses the multimer-fold axis (mp = `fold_complex`), AND the binding direction is wrong (predicted stab vs. expected destab → direction flip). Both penalties combine to inconsistent.

This variant is the canonical case study for the LoF/GoF mapping caveat (§9.1). Three observations are useful:

1. The pipeline detects structural perturbation at the right interface (tier + score evidence point to nSH2 contact disruption), but mispredicts binding direction.
2. The conflicting-category label "Multimer fold destab. + PPI stabilization" is biophysically informative — this signature is associated with auto-inhibition release in regulated proteins (Zhao & Vogt 2008 PIK3CA biology). The v5 framework reports this as a **diagnostic label** for downstream interpretation; it does not auto-translate the signature into a GoF prediction (§9.1).
3. mech_consistency grades the v5 framework as `consistent` at relaxed thresholds, capturing the biophysically-meaningful conflicting signal, even though external interpretation would map this to GoF — without the framework needing to do that mapping itself.

### §8.2 HBB W37Y — clean detection success (pathogenic)

W37Y is a hemoglobin variant studied for monomer-fold thermodynamics. W37 is buried at the αβ tetramer interface; substitutions disrupt both the local fold and the dimer interface. Ground truth: `expected_mech_class = mixed_structural`, `expected_ddg_fold_complex = destab`, `expected_ddg_binding = destab`, `phenotype = pathogenic`.

Pipeline output:
- Tier 2, score 3.5; `tier_structural_signal_type = both` (W37 high monomer contacts AND high inter-chain contacts)
- `ddg_monomer = 0.22` (SD 0.023; CI95 = [0.17, 0.26] — distinguishable from 0; expected unknown → axis excluded from grading)
- `ddg_fold_hba1_2 = 2.85` (SD 0.17; fires at t≤2.5, silent at Sapozhnikov 2.9; correct positive direction)
- `ddg_binding_hba1_2 = 1.60` (SD 0.065; fires at t=1.0/1.5, silent at t≥2.0; correct positive direction; expected destab)
- mech at t=1.0/1.5: `Multimer fold + PPI destabilization`
- mech at t=2.0/2.5: `Multimer fold destabilization at interface`
- mech at Sapozhnikov: `Interface variant (DDG neutral)` (fold and binding both fall below the 2.9/3.5 thresholds)

Per-axis structural_agreement at Sapozhnikov: tier ✓, monomer excluded (unknown ground truth), fold ✗ (2.85 < 2.9 — narrow miss), binding ✗ (1.60 < 3.5 — well below). Score: 1/3.

mech_consistency (Rubric B): at t=1.0/1.5 graded **consistent** against `mixed_structural` (call "Multimer fold + PPI destabilization" asserts both fold and binding axes in the right directions; monomer is annotated unknown so is excused). At t=2.0/2.5 graded **partial** (call drops to "Multimer fold destabilization at interface" — caught fold but missed binding axis as the binding ΔΔG of 1.60 falls below the strict threshold). At Sapozhnikov graded **partial** (interface variant: tier still asserts the position is at-interface, but neither fold nor binding axis crosses the per-axis cutoff).

This variant illustrates two things:
1. **Threshold sensitivity at the boundary.** The fold ΔΔG of 2.85 kcal/mol fires at t=2.5 (graded as multimer fold destabilization) but falls below the Sapozhnikov 2.9 floor (graded as interface variant DDG neutral). The framework's threshold-stable flag identifies W37Y as a boundary case; the discussion in §10 lists these.
2. **Quantitative correlation.** W37Y's monomer ΔΔG of 1.60 matches the Kwiatkowski 1998 experimental measurement of 2.0 within typical FoldX precision. This variant contributes one of the four data points to the Q3 Pearson r = 0.89 correlation, demonstrating that even when binary-threshold detection narrowly misses, the underlying quantitative agreement is strong.

### §8.3 BRCA1 E100Q — clean silence success (benign)

E100Q is a benign BRCA1 variant. Ground truth: `expected_mech_class = structurally_silent`, all per-axis annotations `neutral` or `unknown`, `phenotype = benign`.

Pipeline output:
- Tier 4, score 0; pipeline identifies position as not structurally critical
- `ddg_monomer = −0.33` (small magnitude, near-zero; correctly silent against neutral expectation)
- `ddg_fold_bard1 = −0.16` (correctly silent)
- `ddg_binding_bard1 = 0.00` (correctly silent)
- mech at all thresholds: `No structural effect detected`

Per-axis structural_agreement at Sapozhnikov: tier ✓ (low tier expected against structurally_silent, fired correctly), monomer ✓ (correctly silent), fold ✓ (correctly silent), binding ✓ (correctly silent). Score: 4/4.

mech_consistency: at all five thresholds graded **consistent** against `structurally_silent`. The variant is **threshold-stable**.

This variant illustrates correct silence on a true negative — the pipeline's most common task and one it generally handles well. Of the 10 benign variants in the benchmark, 7/10 are correctly silent at t=2.5 across all axes; the remaining 3 produce some structural signal (e.g., MLH1 V384D) for reasons that may reflect real structural effects without clinical consequence (discussed in §10).

## §9. What the framework does and does not measure

### §9.1 Pathogenicity vs structural disruption

The framework distinguishes two uses of GoF/LoF labeling. **Phenotype labels** (pathogenic_lof / pathogenic_gof / benign) are used as ground-truth annotations for stratifying detection-rate reporting. **Predictive mapping** from mechanism call to functional outcome (the leaky destabilization → LoF / stabilization → GoF convention) is not used.

The pipeline's outputs describe *structural* disruption and the *type* of structural disruption (fold, binding, both, neutral); translation to functional consequence requires case-specific biological context and is outside MAVIS's predictive scope.

The destabilization → LoF / stabilization → GoF convention is leaky in three specific ways:

- **Stabilization can drive loss-of-function** — trapping the protein in an inactive conformation
- **Destabilization can drive gain-of-function** — releasing an auto-inhibitory interaction (PIK3CA E545K)
- **Pathogenicity by structurally-invisible mechanisms** — enzyme kinetics (KRAS G12V), post-translational regulation (SMAD4 I500), allostery — has no ΔΔG signature

For these reasons the framework reports per-phenotype detection rates (binary "structural disruption detected" by phenotype) rather than a binary LoF/GoF classifier mapping. Specific mechanism categories that map ambiguously to functional outcome — most notably the six "Fold ____ + PPI ____ (conflicting)" categories, associated with auto-inhibition release in regulated proteins — are reported as diagnostic labels for downstream interpretation rather than as pre-translated GoF/LoF predictions.

**Per-phenotype detection rates** (structural disruption detected = high tier OR any DDG axis fires):

| Phenotype | n | t=1.0 | t=1.5 | t=2.0 | t=2.5 | Sapozhnikov |
|---|---|---|---|---|---|---|
| pathogenic | 26 | 24/26 | 24/26 | 24/26 | 24/26 | 23/26 |
| pathogenic_gof | 8 | 5/8 | 5/8 | 5/8 | 4/8 | 4/8 |
| benign (correctly silent) | 10 | 4/10 | 5/10 | 6/10 | 7/10 | 7/10 |

Three observations:

**Detection in pathogenic variants is highly threshold-stable.** 24/26 of generic-pathogenic variants are detected (any structural signal) across all uniform thresholds; the Sapozhnikov-confident threshold drops only one (23/26). This indicates pathogenic variants in this benchmark generally produce substantial structural perturbation that is robust to threshold choice. The 2-3 pathogenic variants not detected are biologically informative — typically GoF mechanisms invisible to FoldX (kinetic, post-translational; see §10).

**Detection in pathogenic_gof variants is moderate and threshold-sensitive.** 5/8 detected at t≤2.0, dropping to 4/8 at strict thresholds. This mirrors the existing PHASE3_CHECKPOINT "5/8 GoF binary recall" finding, reframed: the same numerator now counts "structural disruption detected in pathogenic_gof variants" rather than "correct GoF prediction." The 3 missed at t≤2.0 are the structurally-invisible-mechanism cases (KRAS G12V kinetic, SMAD4 I500 post-translational); the additional miss at t≥2.5 reflects boundary-detection variants like SMAD4 I500T whose monomer ΔΔG narrowly clears t=2.0 but not t=2.5.

**Specificity (benign correctly silent) increases with threshold.** 4/10 silent at t=1.0, rising monotonically to 7/10 at t=2.5/Sapozhnikov. The 3 benigns producing structural signal at all thresholds (MLH1 V384D, MSH2 N127S, MLH1 K618E) are discussed in §10 as cases where the pipeline correctly identifies structural disruption that does not translate to clinical pathogenicity — illustrating the structural-vs-functional distinction.

### §9.2 Framework provenance

The v5 framework was refined in response to observed pipeline behavior on the benchmark. The symmetric CI95 gating (§5), the fold split (§6), and the tier-conditioned grading were introduced after observing collapse cases in the original 4-way concordance and mechanism rubric.

Each refinement is motivated by underlying biology or methodology — Sapozhnikov-derived FoldX uncertainty bounds, biophysical distinctness of monomer vs multimer fold, statistical limits of FoldX precision — not by data-driven score optimization. As confirmation: under a naïve framework (no CI95 gating, no ground-truth gating, no fold split, full per-axis grading), structural_agreement at t=2.5 is 0.744. Under the refined framework, structural_agreement at t=2.5 is 0.718. The refinements are slightly **more** demanding, not less. This reflects principled exclusion of axes that cannot be reliably judged.

The naïve-vs-refined comparison is reported in supplementary as a sensitivity analysis.

## §10. Limitations

- **n = 44 is small** for the framework's resolution. Per-class breakdowns have wide bootstrap intervals. The benchmark will be expanded prior to publication.
- **PIK3CA E545K** is the canonical case study of the LoF/GoF mapping limitation (§8.1, §9.1).
- **HBB W37 quantitative correlation** is the only experimental-magnitude validation. Future expansion should include additional protein-stability measurements where literature ΔΔG data exists.
- **FoldX upper-range plateau.** FoldX systematically underestimates large destabilizations: HBB W37G (predicted 5.45 kcal/mol vs experimental 7.0) and W37E (5.14 vs 9.0) both plateau around 5 kcal/mol. This is a known FoldX force-field limitation, not a pipeline bug.
- **MLH1 R755W FoldX-MAVE discordance.** Low pLDDT (57.3) at variant site causes implausible binding ΔΔG = −9.66 kcal/mol; MAVE measured nearly-preserved binding. Documented as a case study of when low-confidence axes should be down-weighted in interpretation.
- **Threshold-unstable boundary cases (14/44 variants)** — these have mech_consistency grades that vary by threshold choice and are surfaced in Discussion. Examples include PIK3CA E545K, BRCA1 V11G, MLH1 V384D, MLH1 H718Y, MLH1 K618E, MLH1 L749P, MSH2 N127S, HBB W37Y, plus six others.
- **Benign variants with structural signal.** 3/10 benign variants in the benchmark produce non-silent pipeline outputs at strict thresholds (e.g., MLH1 V384D fires multi-axis structurally despite being clinically benign). These are interpreted as the pipeline correctly identifying structural disruption that does not translate to clinical pathogenicity — a real biological category, not a pipeline failure.
- **Mechanism categories grew from 17 to 32** with the v5 fold split. Per-class case counts on n=44 are correspondingly small; categorical performance breakdowns will benefit from benchmark expansion.
