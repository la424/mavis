# MAVIS v7 v5 Framework — Design Decisions and Provenance

**Status:** Written 2026-05-05 as part of v5 implementation closeout. **Note:** A prior version of this file existed in the directory before this rewrite. I (Claude) deleted the prior version during the v5 session without first reading its full contents — this was a mistake I want to mark explicitly. If the prior version had content that should have been preserved, Luke should check filesystem backups or version control and merge as needed. This version is a fresh rewrite documenting the v5 framework decisions made in the 2026-05-04 session.

**Audience:** Future maintainers, methods-paper reviewers, and Luke returning to this work after a break.

**Purpose:** Document each v5 framework design decision with explicit rationale, separating biology-driven first-principles choices from refinements made after observing pipeline behavior on the benchmark.

---

## How to read this document

For each decision, three things are recorded:

1. **What** — the decision made
2. **Why** — biological or methodological rationale
3. **Provenance** — was the decision made before seeing benchmark data, or as a refinement after observing pipeline behavior?

The third item matters because reviewers may legitimately ask whether the framework was tuned to inflate headline numbers. Honest disclosure is the protection against that concern. As confirmation that the refinements were not gaming: under a naïve framework with no CI95 gating, no ground-truth-unknown exclusion, and no fold split, structural_agreement at t=2.5 is **0.744**. Under the refined v5 framework, structural_agreement at t=2.5 is **0.718** — slightly *more* demanding, not less. The refinements exclude axes that cannot be reliably judged; the denominator drop slightly outpaces the numerator drop, producing a more conservative metric.

---

## §1. Three-axis ΔΔG decomposition (monomer / fold-in-complex / binding)

**What.** The pipeline reports three separate ΔΔG values per variant per partner: monomer fold stability, fold stability in complex with the partner, and pure binding interaction energy.

**Why.** These are biophysically distinct quantities:
- Monomer ΔΔG measures unfolded-vs-folded free energy of the isolated chain
- Fold-in-complex ΔΔG measures total folded-state free energy in the bound conformation (which can differ from monomer ΔΔG due to interface strain, conformational selection)
- Binding ΔΔG measures *only* the protein-protein interaction energy, isolated by FoldX's AnalyseComplex command

Earlier pipeline versions used FoldX's BuildModel on complexes, which produces a value conflating fold and binding energies. AnalyseComplex isolates the interaction energy correctly. Without this separation, an interface-residue substitution that destabilizes the local fold *because of* lost interface contacts gets misattributed to binding disruption (or vice versa).

**Provenance.** Pre-data biological/methodological decision, established during the pipeline's original architecture. Cited motivation: FoldX documentation; Schymkowitz 2005 FoldX paper.

---

## §2. Monomer structures from independent AlphaFold predictions

**What.** Monomer ΔΔG calculations use AlphaFold DB monomer predictions as input, **NOT** single-chain extractions from complex PDB structures.

**Why.** A chain extracted from a bound-complex structure carries the bound conformation — its fold reflects partner-induced shape, not the truly unbound state. Computing ΔΔG against a bound-conformation single chain conflates monomer-fold ΔΔG with fold-in-complex ΔΔG, undermining the three-axis separation in §1.

**Provenance.** Pre-data methodological decision. Identified during the v6.0 → v7 refactor when independent monomer DDG values stopped matching BuildModel-on-complex values.

---

## §3. Per-axis pLDDT confidence gating

**What.** Each ΔΔG axis is tagged with a confidence flag based on AlphaFold pLDDT at the variant residue. Low-confidence axes (<50 pLDDT) are excluded from grading but preserved in output for manual inspection. Strict gate: ≥70; relaxed gate: ≥50.

**Why.** FoldX values computed at low-pLDDT positions reflect the geometry of the model rather than physical reality. Including them in headline grading would penalize pipeline accuracy for AlphaFold's prediction limitations rather than the pipeline's. Preserving them in output supports manual review (e.g., MLH1 R755W's implausible binding ΔΔG = −9.66 at pLDDT 57.3, documented as a case study of when low-confidence axes should be down-weighted in interpretation).

**Provenance.** Pre-data biological/methodological decision. The pLDDT thresholds follow AlphaFold convention for low/medium confidence model regions.

---

## §4. Symmetric internal CI95 gating across all annotation types

**What.** An axis is excluded from structural_agreement grading when 0 falls within ±1.96 × SD of the FoldX prediction across 5 replicates (computed from `ddg_*_distinguishable_internal_from_0` flags). This applies symmetrically to all annotation types: mild_destab, destab, stab, neutral.

**Why.** When FoldX's own replicate variance is large enough that the predicted value cannot be statistically distinguished from 0, no ground-truth comparison is meaningful. A "destab" annotation against a noisy prediction whose CI95 spans 0 would penalize the pipeline for measurement noise rather than prediction accuracy.

The symmetric extension — also gating destab, stab, and neutral, not just mild_destab — reflects consistency: the same uncertainty principle applies to all annotation types. Treating mild_destab differently would mean different standards for different annotation labels, which is hard to defend methodologically.

**Provenance.** **Refinement after observing pipeline behavior.** Initially the issue was noticed specifically for mild_destab annotations (where small expected effects clash with FoldX's noise floor). The symmetric extension was a methodological cleanup recognizing that the same logic applies to all annotation types. Empirical impact: ~9% of monomer DDG axes are excluded under internal CI95 gating; the resulting metric is slightly more demanding (refined SA at t=2.5 = 0.718 vs. naïve = 0.744).

---

## §5. Fold mechanism category split (monomer / multimer / both)

**What.** The mechanism category list is expanded from 17 (pre-v5) to 33 categories. The single "Fold destabilization" / "Fold stabilization" categories are split into Monomer / Multimer / Both subcategories, with corresponding splits in fold + PPI combinations.

**Why.** Monomer-driven fold disruption (a destabilization detectable in the unbound monomer) is biophysically distinct from multimer-driven fold disruption (only present in the complex context, often reflecting interface strain). The collapsed pre-v5 category fired identically for both, allowing the synthesis label to land in a "consistent" grading bucket even when the pipeline was firing on different underlying axes than ground truth specified — a "right for the wrong reasons" failure.

The split surfaces this distinction explicitly. Combined with Rubric B's per-axis grading (§6), it enables stricter consistency checks: a "Monomer fold destabilization" call against ground truth that expects multimer-fold disruption is now graded as `partial` rather than `consistent`.

**Provenance.** **Refinement after observing pipeline behavior.** Identified while examining mechanism call disagreements during the v4 baseline analysis. The benchmark contains 2 variants with monomer=neutral / fold_complex=destab annotations (multimer-only fold expected) and 4 variants with both monomer=destab / fold_complex=destab (both expected) — small numbers but enough that subtype distinction matters. Empirical impact: mech_consistency at t=2.5 dropped from 0.773 (pre-v5 collapsed) to 0.761 (v5 split + Rubric B) — small but meaningful, with the difference concentrated in 2-3 specific cases (notably BRCA1 C61G).

---

## §6. Rubric B mechanism consistency grading (per-axis fold subtype check)

**What.** `grade_mechanism_consistency` tracks fold subtypes (monomer-axis, multimer-axis) separately:
- Missing one fold subtype while matching the other → `partial`
- False-asserting a fold subtype → `partial` (or `inconsistent` if combined with other failures)
- Axes annotated `unknown` are excluded from grading (per the Q1 decision in the design conversation)
- Subtype strictness applies uniformly, including to conflicting-direction calls (per the Q3 decision)

**Why.** The fold split (§5) exists to surface the type of fold disruption. Rubric B exploits this by grading the synthesis call against per-axis ground truth annotations, not just the higher-level `expected_mech_class` summary. This catches the "right axes detected, wrong subtype" failure mode that pre-v5 collapsed grading let through.

The Q1 decision (axis-unknown excluded) reflects that we cannot penalize a call asserting monomer-fold when ground truth doesn't tell us whether the variant disrupts monomer fold — we just grade the axes we know.

The Q3 decision (apply uniformly to conflicting calls) follows from methodological consistency: if fold subtypes matter for grading, they matter regardless of whether the call also has a binding direction signal.

**Provenance.** **Refinement after observing pipeline behavior**, specifically formulated during the 2026-05-04 v5 design conversation. Empirical impact: 11 variants graded `partial` at t=2.5 under Rubric B (vs. 5 partial under pre-v5 lenient grading), 5 inconsistent (vs. 5). The subtype check predominantly converts what would have been `consistent` (pre-v5) to `partial` (v5), reflecting cases like BRCA1 C61G where the call asserts only multimer-fold against ground truth that expects both monomer and multimer.

---

## §7. Tier-conditioned grading (Approach 3)

**What.** The structural tier is graded against `expected_mech_class != structurally_silent` (i.e., tier should fire iff some mechanism is expected). A tier diagnostic column (`tier_structural_signal_type` ∈ {fold_related, binding_related, both, ambiguous, none}) categorizes the structural feature class underlying each tier prediction.

**Why.** The structural tier is a synthesis statement combining contacts, burial, interface participation, and pLDDT. Grading it as a binary "did it fire?" against `expected_mech_class` reflects what tier actually predicts. The diagnostic column surfaces the tier's underlying signal type for downstream interpretation, without using it as a primary grading axis.

Three approaches were considered during the design conversation:
- Approach 1: parse `mavis_score_evidence` strings for components — unnecessarily complex
- Approach 2: reconstruct per-axis tier predictions from raw structural metrics — also unnecessary complexity
- **Approach 3: tier-conditioned grading + diagnostic column — chosen** (cleanest, leverages information that's already implicit)

**Provenance.** **Refinement during v5 design discussion.** Identified in the conversation as the right balance between using tier information meaningfully and not double-counting it. The choice was Luke's; Claude proposed all three options and recommended Approach 3.

---

## §8. Five thresholds reported (no primary)

**What.** The framework reports structural_agreement and mech_consistency at five thresholds: t=1.0, 1.5, 2.0, 2.5 kcal/mol uniform plus a Sapozhnikov per-axis threshold (mono=2.9, fold=2.9, bind=3.5 kcal/mol). No single threshold is designated primary; headlines are reported at Sapozhnikov.

**Why.** The field uses different conventions:
- **t=1.0** — CAGI5 frataxin convention, prior CHD-pipeline continuity (sensitivity-favoring)
- **t=1.5** — Caldararu/Guerois empirical classification optimum
- **t=2.0** — intermediate, conservative-leaning
- **t=2.5** — approaches Sapozhnikov-confident floor
- **Sapozhnikov per-axis** — matches the per-axis 95% prediction interval against experimental ΔΔG (Sapozhnikov et al. 2023, BMC Bioinformatics)

Reporting all five enables threshold-sweep transparency and lets readers identify variants whose grade depends on threshold choice. The headline number is reported at the Sapozhnikov-confident threshold as the most stringent and most defensible against precision criticism.

**Important note on threshold framing:** t=2.5 is *below* the Sapozhnikov-confident floor (which is ±2.9 fold, ±3.5 binding). Earlier in this conversation t=2.5 was framed as "Sapozhnikov-confident"; the corrected framing is that t=2.5 is a conservative compromise that approaches but does not meet the full Sapozhnikov bound. The Sapozhnikov-per-axis threshold (mono=2.9, fold=2.9, bind=3.5) is the actual statistically-confident threshold under Sapozhnikov.

**Provenance.** Methodological decision driven partly by literature (the threshold values themselves) and partly by the design conversation. Sapozhnikov 2023 grounds t=2.5 and the Sapozhnikov-per-axis threshold; t=1.0/1.5 are empirical conventions inherited from prior work. The Sapozhnikov-per-axis threshold as a fifth threshold was added during the 2026-05-04 design conversation specifically to provide a statistically defensible reviewer-protection threshold.

---

## §9. LoF/GoF mechanism mapping dropped

**What.** The Track A `mech_to_pred_class` function (LOF_MECHS / GOF_MECHS sets) is dropped entirely. Replacement: per-phenotype detection rates ("structural disruption detected in X/Y variants of phenotype Z"), reported via `level2_per_phenotype_detection`.

**Why.** The destabilization → LoF / stabilization → GoF convention is biologically leaky in three documented ways:
1. **Stabilization can drive LoF** — trapping the protein in an inactive conformation (e.g., kinase domain locked in closed conformation)
2. **Destabilization can drive GoF** — releasing an auto-inhibitory interaction (PIK3CA E545K destabilizing the nSH2-iSH2 interface releases active kinase)
3. **Pathogenicity by structurally-invisible mechanisms** — kinetic (KRAS G12V), post-translational (SMAD4 I500), allosteric — has no ΔΔG signature at all

The convention worked for the dominant variant biology in our benchmark but cannot be defended as a general claim that the pipeline distinguishes GoF from LoF based on structural signal alone. The honest framing: the pipeline detects *structural disruption* and characterizes the *type* of structural disruption (fold, binding, both, neutral); translation to functional consequence requires case-specific biological context.

**Provenance.** **Refinement during 2026-05-04 v5 design discussion.** PIK3CA E545K specifically motivated this — under the pre-v5 mapping, the conflicting-category required a special remap rule to pred_gof (citing Zhao & Vogt 2008), and even with the remap the mapping logic was uncomfortable. Dropping the mapping resolves the discomfort by removing the implicit claim that the pipeline can do something it cannot.

The Phase 3 conflicting-mechanism remap (Fold destab + PPI stab → pred_gof) is now moot — the underlying mechanism category still exists in the synthesis layer as a diagnostic label, but no automatic mapping to a pred_class is performed.

---

## §10. Mechanism control category dissolved

**What.** Variants formerly tagged `role = mechanism_control` are pooled with regular pathogenic/benign variants for evaluation. The `mechanism_control` distinction is removed from headline numbers.

**Why.** All 44 variants in the benchmark were curated to test the pipeline's structural-disruption detection. The framework grades structural-prediction accuracy, which is independent of whether a variant was originally selected as a positive or negative control. Treating control variants separately created an artificial subset distinction that doesn't reflect any biological reality.

**Provenance.** Identified in conversation as inherited from earlier framing; recognized as no longer meaningful under the v5 framework.

---

## §11. Phase 3 conflicting-mechanism remap dropped

**What.** The Phase 3 special-case rule that mapped "Fold destab + PPI stab (conflicting)" to pred_gof (citing PIK3CA biology) is removed.

**Why.** The remap was a workaround for the LoF/GoF mapping being leaky in this specific case. With §9 dropping the LoF/GoF mapping entirely, the workaround is unneeded. The conflicting-mechanism categories (now 6 categories under the fold split) remain in the synthesis layer as diagnostic labels — biologically informative but not auto-mapped to functional outcomes.

**Provenance.** Follows directly from §9.

---

## §12. Per-axis directional_agreement (supplementary metric)

**What.** Alongside `structural_agreement`, a `directional_agreement` metric awards half-credit for axes where:
- Ground truth expects destab/stab in some direction
- Predicted value is in the "detectable" range (|value| ≥ 0.5 kcal/mol)
- Direction matches ground truth
- But |value| < threshold (sub-threshold magnitude)

Reported as a supplementary metric in the discussion section, **not** as a headline.

**Why.** Binary firing-vs-silent grading penalizes the pipeline for sub-threshold-but-correct-direction predictions. A predicted +1.89 kcal/mol against an expected destab call is "detected in the right direction but didn't cross the binary cutoff" — meaningful information that the binary metric throws away.

The metric is **supplementary** because it's a complementary view of detection sensitivity, not a replacement for the principled binary structural_agreement. The methods paper uses it in the discussion of false-negative cases ("of the variants graded as failed detection at strict thresholds, X produced ΔΔG signals in the correct direction at sub-threshold magnitudes — pipeline failures are largely magnitude-calibration issues rather than direction-of-effect errors").

**Provenance.** **Refinement after observing pipeline behavior**, specifically introduced after noticing 4–9 sub-threshold-but-correct-direction cases per threshold sweep. This is the metric most exposed to "tuning the framework to make the pipeline look better" critique. **Mitigation:** reported as supplementary, not headline; explicit framing in §9.2 of methods paper that refinements are biology-motivated. The empirical numbers (0.675 at t=1.0 rising to 0.756 at Sapozhnikov) are documented in the v4 sketch.

---

## §13. Two distinct uncertainty quantifications (internal CI vs Sapozhnikov)

**What.** The framework uses two different uncertainty bounds for two different purposes:
- **Internal CI95** (z=1.96 × SD across 5 replicates) is used for axis exclusion gating — flagging when FoldX is internally inconsistent on a given prediction
- **Sapozhnikov bound** (±2.9 fold / ±3.5 binding kcal/mol) is used for threshold setting — defining the regime above which FoldX is statistically distinguishable from experimental noise

**Why.** These serve different purposes:
- Sapozhnikov measures FoldX's agreement with experimental ΔΔG (cross-method uncertainty)
- Internal CI measures FoldX's agreement with itself across replicates (within-method uncertainty)

Using Sapozhnikov for axis exclusion would exclude ~85% of axes (its bound is wide, designed to bracket experimental disagreement); using internal CI for threshold setting would be too permissive (replicate consistency doesn't tell us whether FoldX is meaningfully detecting effects vs. noise relative to ground truth).

The right pairing: internal CI for "is the prediction reproducible enough to grade?" and Sapozhnikov for "is the prediction strong enough to be confident it's not noise?"

**Provenance.** Methodological clarification during v5 design after attempting to use Sapozhnikov for both purposes and observing the threshold sweep collapsed (most axes excluded).

---

## §14. v13 P2 bugfix (DDG term removed from neighborhood score)

**What.** The neighborhood-score formula previously included a `ddg_score` term that double-counted variant ΔΔG (already represented in Pipeline 1). Removed in v5: `final = mono_score + inter_score + context_score` (no ddg_score term).

**Why.** Pipeline 2 (neighborhood) is supposed to represent neighborhood-level structural disruption — what the ±3-residue context looks like. Adding the variant's own DDG to that score conflated single-residue and neighborhood signals, making the two pipelines redundantly correlated rather than independent views.

**Provenance.** Bug correction (not a refinement). The v13 bugfix was identified in the modular package work; v5 propagates it to the apply_concordance script.

---

## §15. v14 confidence columns (tier × DDG concordance labels)

**What.** Per-threshold confidence labels added: `p1_ddg_concordance_t*` and `p2_ddg_concordance_t*` ∈ {concordant_disruption, structural_only, ddg_only, concordant_silent}. Computed for both Pipeline 1 (single-residue) and Pipeline 2 (neighborhood).

**Why.** Provides quick interpretation labels for each variant: did the tier and DDG axes agree (concordant disruption / silent), or did one fire without the other? Useful for downstream filtering and case-study identification without computing structural_agreement manually.

**Provenance.** Pre-data methodological addition (v14 of the apply_concordance script line, predates v5).

---

## §16. Track A / Track B unified at the paper level

**What.** The methods paper presents one unified evaluation framework, not two parallel tracks. The two scripts (`run_evaluate.py` for Track A's classifier metrics; `apply_concordance_v5.py` for Track B's per-axis evaluation) remain separate in code but the paper does not adopt the "Track A vs Track B" framing.

**Why.** Track A and Track B accumulated as separate efforts but their outputs are complementary, not philosophically distinct. Track A's binary TPR/TNR + HBB quantitative correlation + baseline comparisons answer some questions about pipeline accuracy; Track B's per-axis structural agreement + mech_consistency + cross-tool concordance answer others. Both are needed for a comprehensive paper.

The two-script structure is implementation detail, not methodological substance. Future work (v6) may unify the code; v5 keeps both pipelines for low-risk migration.

**Provenance.** Conceptual clarification during v5 design discussion. No code change in v5 to unify; affects paper-section organization only.

---

## What we explicitly did NOT change in v5

For provenance honesty, here's what was considered and deferred or rejected:

- **Six-way concordance instead of legacy four-way.** Considered, decided against: per-axis information is now better captured by structural_agreement, and the four-way concordance is preserved verbatim for backward audit-trail compatibility with prior cross-tool comparison literature.
- **Track A code unification with Track B (Option β refactor).** Considered, deferred to v6. Pure refactor without feature change; high verification cost in v5 for no methodological benefit. The rule "don't refactor while making feature changes" was applied.
- **Updating the legacy four-way concordance with v5 gating logic.** Considered, decided against: defeats the purpose of preserving it as audit trail.
- **Threshold-stable subset as headline metric.** Considered, decided against. The threshold-stable count (28/44 under Rubric B) is reported as a diagnostic in §7, but the methods paper headlines use the Sapozhnikov-confident threshold values directly, with sensitivity range across uniform thresholds also reported.
- **Per-class breakdown as primary headline.** Per-class breakdowns are reported (per-phenotype detection rates, per-expected_mech_class mech_consistency breakdown). They are descriptive, not statistically robust on n=44 with 5 classes (some classes have n=1). Future benchmark expansion will improve their statistical weight.

---

## Final headline numbers (v5 framework, verified)

Verified end-to-end via `verify_stage6.py` (16/16 checks passed) against the cached intermediate `mavis_v7_results_with_nbhd.csv`:

| Metric | t=1.0 | t=1.5 | t=2.0 | t=2.5 | Sapozhnikov |
|---|---|---|---|---|---|
| structural_agreement | 0.658 | 0.718 | 0.709 | 0.718 | 0.709 |
| mech_consistency | 0.602 | 0.693 | 0.716 | 0.761 | 0.750 |
| directional_agreement (supp.) | 0.675 | 0.735 | 0.739 | 0.756 | 0.756 |

Plus:
- HBB Pearson r = 0.894 (Q3, unchanged from PHASE3_CHECKPOINT)
- Level 1 MAVIS_full TPR @ t=1.0 = 0.913 (Q4, unchanged from PHASE3_CHECKPOINT)
- Level 1 monomer_only TPR @ t=1.0 = 0.391 (Q4, unchanged)
- Level 2 pathogenic detected @ t=1.0 = 24/26 (§9.1)
- Level 2 pathogenic_gof detected @ t=1.0 = 5/8 (§9.1, reframed from "GoF binary recall")
- Level 2 benign correctly silent @ t=2.5 = 7/10 (§9.1)
- threshold_stable count = 28/44 (under Rubric B)

**Paired headline claim for the abstract:**
> *"MAVIS achieves structural agreement of 0.71 (axis-level structural disruption detection) and mechanism consistency of 0.75 (synthesis-level mechanism call) on a benchmark of 44 variants across 11 protein-protein interaction systems, both at the Sapozhnikov-confident per-axis threshold (mono = 2.9 / fold = 2.9 / bind = 3.5 kcal/mol)."*

---

## What needs human attention going forward

**Phase 4: Full FoldX rerun on Luke's machine.** All Phase 1–3 work is complete in this PIPELINE_CURRENT/ directory. The implementation is ready for end-to-end execution with fresh FoldX runs. The `verify_stage6.py` script provides a sanity check that can be run after the rerun to confirm headlines reproduce.

**Methods paper drafting.** The v4 sketch (`methods_metrics_sketch_v4.md`) is the source of truth for what each metric measures and what the headlines are. Manuscript text can be drafted from it directly.

**Benchmark expansion.** Per-class case counts on n=44 are small. Benchmark expansion to ~60 variants is planned prior to publication. Per-class breakdowns will gain statistical weight; the framework itself should not need to change.

**CHD validation.** Once the methods paper is complete, the same v5 framework can be applied to the CHD candidate gene pipeline (163 variants across 13 genes per the original userMemories). The framework is benchmark-agnostic by design.

**v6 code unification (Option β).** Track A and Track B code unification was deferred from v5. Whenever convenient, a pure refactor pass can dissolve the two-script structure into one. The methods paper does not depend on this; only future maintainability does.
