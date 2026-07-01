# MAVIS — Results Synthesis for the Two Papers

*Prepared 2026-06-15. This document was written **after** a verification pass against the raw data files, not from memory. Every quantitative claim is tagged:*

- **[V]** = recomputed from the data files in this session and matches.
- **[L]** = attested in the canonical ledger but **not** independently recomputed here.
- **[F]** = **flagged** — could not be reproduced from the released data, or was found to be wrong; must be resolved before it appears in a paper.

*Data files used: `mavis_v7_concordance_annotated.csv` (benchmark, 44 variants × 1255 cols — appears to be the canonical reconciled file under a different name than the ledger's `mavis_v7_concordance_v5_reconciled.csv`); `chd_concordance_results_FIXED.csv` (CHD, 384 rows / 144 variants); `fold_zic3_model_0.{cif,pdb}`; `MAVIS_v7_canonical_benchmark_ledger.md`.*

---

## 1. Origin and rationale

The project grew from a practical coverage problem, not from a thesis. Most proteins do their job in complex, but structural and pathogenicity prediction is overwhelmingly monomer-based — so anything a variant does at an interface, or to the stability of the assembled complex rather than the isolated chain, is invisible to the standard toolkit. The goal was to **fold multimer analysis into structural/pathogenicity prediction**, and the natural testbed was the congenital-heart-disease (CHD) gene set (SHROOM3, ZIC3, and partners), because those genes act through interaction networks where the monomer story is incomplete.

Two capabilities emerged, in sequence. First the intended one — **multimer-aware disruption detection**, scoring a variant inside its complexes (interface energy + complex-fold stability), not just the isolated chain. Then, while building it, the realization that the same FoldX decomposition could attribute disruption to a **specific axis** — monomer fold, complex fold, or binding. That mechanism-of-disruption capability was the extension, and it is **why the benchmark exists**: the 44-variant / 11-system set was assembled specifically to test whether mechanism prediction is real and reliable on variants with characterizable mechanisms.

The theme **"structural disruption ≠ pathogenicity"** is therefore a *finding*, not a founding premise. The four-way concordance frame places the two structural channels (tier + FoldX ΔΔG) beside the two pathogenicity channels (AlphaMissense + Franklin/ClinVar); the concordance landscape is the empirical map of where structural disruption and pathogenicity coincide versus diverge, and the divergences are what make their separability an earned conclusion.

For publication, the order inverts the chronology: the methods/benchmark paper is the prerequisite (validate mechanism prediction, then apply it), so it is the priority despite being the later-built piece.

---

## 2. Benchmark pipeline (methods paper)

### Setup
44 variants across 11 PPI systems, every one a mechanism control. There is no `mechanism_control` subclass; 11 rows that had mislabeled roles were corrected. The tier system is **Grantham severity × contact count, with no ΔΔG term**, so structural-evidence strength is not circularly defined by the ΔΔG it complements. Each variant carries primary-literature-grounded ground-truth tokens per axis (`expected_ddg_{monomer, fold_complex, binding}` ∈ {stab, neutral, destab, unknown}), under a strict standard (Definition B: direct biophysical measurement required; FoldX self-prediction and bare structural-position inference are inadmissible).

### Results
- **[V] Composition:** 44 variants — benign 10, pathogenic 26, pathogenic_gof 6, pathogenic_lof 2.
- **[V] Mislabel corrections (all 11 confirmed present):** KRAS G12D/G12V/Q61H + PIK3CA H1047R → `pathogenic_gof`; BRCA1 C61G, MLH1 R755W, MSH2 G674R/C697F, VHL W117R/Y98H, HBB E6V → `pathogenic`.
- **[V] Structural-evidence tier gradient (reproduces exactly):** T1 12/12 = **100%**, T2 13/16 = **81%**, T3 7/10 = **70%**, T4 2/6 = **33%** pathogenic. This is the cleanest single demonstration that structural-evidence strength correlates with — but does not equal — pathogenicity (T4 is still 33% pathogenic).
- **[V] mech_consistency (raw) = 0.70** ((26 consistent + ½·8 partial) / 43 evaluable = 0.698). **[L] pLDDT-reconciled = 0.73** (applies the four documented interface-exclusions; the +0.03 was not independently recomputed here).
- **[V] structural_agreement (thresholded, t=2.5) = 0.764** (81/106 evaluable axes) — reproducible and consistent with the ledger headline 0.77 (sweep range 0.76–0.80). **This is the number to report as the primary `structural_agreement`.** **[F → diagnosed] The directional variant (ledger 0.773 = 51/66 strict / 0.757 = 53/70 relaxed) is not reproducible from this file, but the reason is now understood:** the ledger's denominator counts axes *per partner chain* for multi-chain systems — the five HBB W37/N102 variants each carry three partner chains, with FoldX false-negatives against the T-state — and applies the four §5 pLDDT interface-exclusions, an enumeration the collapsed vote columns discard (per-partner graded axes total 100, collapsed 64; the script's convention sits between, at 66). Reconstructions over the collapsed columns bracket the value at **0.79–0.82**; the file's stored `directional_agreement_*_t25` column gives 0.792 but **counts `unknown`-token axes and must not be cited as-is**; the relaxed 53/70 cannot be reconstructed here at all (no `expected_ddg_*_relaxed` columns — those live in `mavis_v7_results_relaxed.csv`). **Report the thresholded 0.764 / 0.77 as primary; regenerate the directional from the grading script with its axis convention documented, or omit it.**
- **[V] Claim-arm counts:** PPI-disruption (binding-axis destabilizers) = **13**; complex-fold destabilizers = **11**; **monomer-fold destabilizers = 2** — the undersupported arm. **[F] "Silence-not-benign" is 19, not 29:** the verifiable count of non-benign variants with no destabilizer on any structural axis is 19 (the "29" in prior notes was the neutral-binding-axis count, which includes benign variants — a different and incorrect quantity).
- **[V] Gain-of-function reframing holds, with a presentation caveat:** all 6 GoF variants have `expected_mech_class = structurally_silent` and all 6 score `mech_consistency = consistent`. **However**, the raw `mavis_mechanism` *display label* reads "fold destabilization" for PIK3CA H1047R and SMAD4 I500T (a threshold-based label), while the consistency *metric* — which gates on statistical distinguishability from zero — correctly scores them as silent. **Figures and slides must use the consistency verdict, not the raw label, for GoF variants.**
- **[V] Pipeline 2 (neighborhood tier) — the substantive claims verify.** P1 and P2 have **identical** mech_consistency (both consistent 26 / partial 8 / inconsistent 9 → 0.698 raw), so the neighborhood pipeline adds **no** mechanism resolution. The neighborhood tier **degrades the pathogenicity gradient**: the P1 tier carries an odds ratio of 6.48 (T1/T2 vs T3/T4 for pathogenic vs benign) versus 4.00 for the neighborhood tier, and neighborhood-*elevated* variants are anti-enriched for pathogenicity (58% pathogenic among the 12 elevated vs 84% among the 32 concordant; OR 0.26). **[F → note] The ledger's specific "OR 0.48" does not reproduce under either natural 2×2** (the tier OR is 4.00, the elevated-subset OR is 0.26) — that exact figure needs its contingency table pinned — but the conclusion (neighborhood is a tested-and-rejected alternative) is fully supported. The `pipeline_agreement` split is Concordant-high 28 / Concordant-low 4 / Neighborhood-elevated 12.

### Interpretation
The machinery and the central theme both validate. mech_consistency 0.70 (raw) / 0.73 (reconciled), structural_agreement ~0.76, and a monotonic tier gradient establish that mechanism prediction is reliable for the **binding and complex-fold** axes (13 and 11 graded destabilizers respectively). The **monomer-fold axis is the weak arm** — only 2 graded destabilizers — so the honest framing is that the benchmark validates roughly **two and a half of three axes**. The GoF variants reinforce disruption ≠ pathogenicity by being pathogenic yet scored structurally silent. The "right-for-wrong-reasons" failure mode (a correct call from incorrect axis attribution, or a threshold label that disagrees with the CI-aware metric) is treated as first-class throughout.

---

## 3. CHD pipeline (application paper)

### Setup
144 variants — **11 ZIC3 zinc-finger known-pathogenic controls + 133 patient-population variants** — across 10 per-pairing systems (SHROOM3 × {actin, dvl2, ctnnb1, rock2, cdh2}; ZIC3 × {gli3, kpna1, kpna6, mdfi, tcf7l1}), bidirectionally scored, four-way concordance collapsed to one row per variant. The 11 controls are fully citation-matched.

### Results
- **[V] Collapsed concordance (after the AlphaMissense parse correction):** 45 evaluable (4/4 = 6, 3/4 = 10, 2/4 = 4, 1/4 = 10, 0/4 = 15), 99 unevaluable. The correction moved KPNA6 I498T (3/4 → 4/4) and K424N (2/4 → 3/4); all other collapsed hit counts are unchanged from the pre-correction ledger.
- **[V] Control validation:** structure axis fires on **8/9** evaluable controls, ddG on **7/9**; H286R is the only full structural miss (NLS/trafficking mechanism — correctly silent, 2/4 external-only). Concordance: 4/4 — C253S, C297F, H318N, T323M, K405E; 3/4 — H281Y, R350G, S402P; 2/4 — H286R; out-of-scope (disordered, pLDDT < 30) — P217A (0/2), A447G (1/2).
- **[V] Patient candidates:** eight, with **KPNA6 I498T now at 4/4** (all four channels agree) and seven at 3/4 that each miss a *different* channel:

  | candidate | structure | ddG | AM | Franklin | missing | reading |
  |---|---|---|---|---|---|---|
  | KPNA6 I498T | ✓ | ✓ | ✓ | ✓ | — (4/4) | all four channels agree; AM restored by the parse fix |
  | DVL2 D441Y | ✓ | ✓ | ✓ | — | Franklin | clean novel candidate (clinically unannotated) |
  | SHROOM3 G35V | ✓ | ✓ | ✓ | — | Franklin | clean novel candidate |
  | SHROOM3 G60V | ✓ | ✓ | ✓ | — | Franklin | clean novel candidate |
  | ROCK2 T367M | ✓ | ✓ | — | ✓ | AM | structurally + clinically flagged, AM-discordant |
  | KPNA6 K424N | ✓ | — | ✓ | ✓ | ddG | AM-corrected; structure + AM + Franklin agree, ΔΔG sub-threshold |
  | CTNNB1 R151C | ✓ | — | ✓ | ✓ | ddG | caught by the structure tier, ΔΔG sub-threshold |
  | CTNNB1 T297M | ✓ | — | ✓ | ✓ | ddG | caught by the structure tier, ΔΔG sub-threshold |

- **[V] Evaluability is gated on `best_plddt = max(monomer_plddt, multimer_plddt_max) ≥ 50`** — **not** monomer pLDDT alone. Monomer pLDDT does not separate evaluable from unevaluable (evaluable down to 41.4, unevaluable up to 49.8); `best_plddt` separates cleanly at 50. The gate is system-specific and already credits fold-upon-binding: SHROOM3 N1972D and P1971L (C-terminal SD2 region) are unevaluable in four partner systems but become evaluable in `shroom3_rock2`, where that residue gains order (position pLDDT 51–53) while the monomer stays disordered (41–44).
- **[V] The 99 unevaluable variants are disordered in *every* context** (best_plddt < 50; the ZIC3 N-terminal/importin set tops out at multimer pLDDT 45), so no max-rule change reaches them.
- **[V] Zinc gap:** the ZIC3 monomer is an AlphaFold-3 (Server) **zinc-free** prediction, so FoldX fold-ΔΔG on the C2H2 fingers reflects packing only. The two coordinating cysteines split accordingly: C297F fires at +12.1 kcal/mol via a coincidental Cys→Phe steric clash, while C253S registers only +1.64 (sub-threshold; the coordination-loss term is unmodeled) and is recovered by the structure axis, not ΔΔG. (C253/C297 ligand assignments are literature-grounded; H318-as-ligand is inferred, not confirmed — the automated C2H2 register only cleanly resolved ZF5.) The large `zic3_gli3` binding values are real interface contacts for R350G/K405E/S402P/H318N but artifacts (0 interface contacts) for C297F and T323M — neither artifact changes any collapsed hit count.

### Interpretation
The multimer integration delivers where it should. K405E and S402P are real zic3–gli3 interface variants the pipeline catches *at* the interface (K405E is monomer-neutral, ddg_monomer ≈ 0, yet a genuine 4-contact interface disruptor) — exactly the effect a monomer-only analysis misses. Evaluability is genuinely max-based, and it credits fold-upon-binding (the two SHROOM3 rescues). The control recovery (8/9) leans on the coordination-agnostic structure axis rather than ΔΔG, which is a stronger argument for the multi-channel design than a clean ΔΔG result would be.

The standout candidate: **KPNA6 I498T reaches 4/4** — structure, ΔΔG, AM, and Franklin all converge — once the AlphaMissense parse error is corrected (the AM score had been mis-loaded into the class column across the entire gene, silently zeroing the AM hit). Recovering it is itself a small methods-rigor result: an annotation-pipeline bug was masking a fully-concordant candidate. The remaining seven each miss a *different* channel, which is the more informative pattern — the three Franklin-missing variants (DVL2 D441Y, SHROOM3 G35V/G60V) are the cleanest genuinely-novel candidates (clinically unannotated, so no Franklin hit yet); ROCK2 T367M is structurally and clinically flagged but AM-discordant (AM calls it benign at 0.07 against Tier-1 structure+ΔΔG agreement), the single clearest MAVIS-vs-AM disagreement; and the ddG-missing variants (KPNA6 K424N, CTNNB1 R151C/T297M) were caught by the structure tier where the ΔΔG fell short. Note the threshold dependence: I498T’s 4/4 counts its VUS(high) Franklin as a hit; under strict pathogenic-only Franklin the novel candidates stay at 3/4, since a VUS is not yet a clinical classification. Roughly 69% of the patient cohort (the 99 unevaluable) is structurally inaccessible, so for those variants interpretation rests entirely on AM/Franklin; the structural arm's power lives in the ordered zinc-finger and interface regions.

---

## 4. Did we achieve our goals?

**Original goal — integrate multimer analysis into structural/pathogenicity prediction — yes, with a precise caveat.** The pipeline is multimer-aware, catches interface variants monomer tools miss (K405E/S402P; benchmark binding-disruption n = 13), and the four-way frame operationalizes the integration with pathogenicity. Evaluability is genuinely max(monomer, multimer)-gated and credits fold-upon-binding. The caveat, learned empirically: the integration adds interface-specific calls and a *small* amount of coverage (the two fold-on-binding rescues), but it cannot reach variants that are disordered in every context, which is most of the CHD patient cohort.

**Extension goal — predict the mechanism of disruption — validated at ~2.5 of 3 axes.** mech_consistency 0.70 (raw) / 0.73 (reconciled) and structural_agreement ~0.76 support the binding and complex-fold axes; the monomer-fold axis (2 graded destabilizers) is the unproven arm and is the methods paper's one true pre-submission task.

**Emergent payoff — structural disruption ≠ pathogenicity — earned, not assumed.** Supported by the tier gradient, the GoF variants (mech-consistent with silence), and the 19 structurally-silent non-benign variants.

---

## 5. Limitations (state plainly in both papers)

1. **Metal-coordination energetics are unmodeled** (apo AlphaFold + FoldX), visible in the ZIC3 zinc fingers.
2. **Intrinsically disordered regions are structurally inaccessible**, though the max-pLDDT gate does credit fold-upon-binding where a complex orders the position.
3. **The monomer-fold axis is undersupported** in the benchmark (2 graded destabilizers).
4. **Display label vs. metric can disagree** (GoF): the threshold-based mechanism label is not the CI-aware consistency verdict.

---

## 6. Open items before publication

1. **[Resolved → action] The directional `structural_agreement` figure is not reproducible from the released data.** Its denominator (66) reflects the grading script's per-partner multi-chain axis enumeration plus the §5 pLDDT exclusions (per-partner graded axes = 100, collapsed = 64, script convention = 66); reconstructions bracket 0.79–0.82. **Report the reproducible thresholded value (0.764 / 0.77 @ t=2.5) as primary**, and either regenerate the directional from the grading script with the axis convention documented, or omit it. The relaxed 53/70 additionally requires `mavis_v7_results_relaxed.csv` (no relaxed tokens in the annotated file).
2. **[F] Fix "silence-not-benign."** The verifiable count is 19, not 29; pin the exact definition.
3. **[Mostly resolved] Pipeline 2.** P1 == P2 mech_consistency (both 0.698 raw) and the gradient-degradation are now verified from the data (tier OR 6.48 → 4.00; neighborhood-elevated subset anti-enriched, OR 0.26). Only the ledger's specific "OR 0.48" is unmatched — pin its 2×2 definition, or report the gradient-degradation figures above instead.
4. **[F] Resolve the canonical-filename discrepancy** (`mavis_v7_concordance_annotated.csv` vs the ledger's `mavis_v7_concordance_v5_reconciled.csv`) — confirm they are the same artifact.
5. **Expand the monomer-fold destabilizer cohort** (2 → 6–8 graded; leads: BRCA1 V11G, MLH1 Q542L, SMAD4 R361C, TNNI3 R162W) — this is the methods paper's gating task and will reopen the recompute.
6. **(CHD) Per-axis evaluability gating refinement:** gate each axis on its own context's pLDDT (monomer-fold ← monomer; complex/binding ← complex-position), so a fold-on-binding variant reports its binding/complex calls and suppresses an untrustworthy monomer-fold call. No effect on the current 45/99 split, but it matters on the full production run.

---

## Appendix A — Verification ledger

| Claim | Status | Evidence |
|---|---|---|
| Benchmark composition 10/26/6/2 (n=44) | **[V]** | role value_counts |
| 11 mislabel corrections present | **[V]** | per-variant role lookup |
| Tier gradient 100/81/70/33 | **[V]** | mavis_tier × pathogenicity: 12/12, 13/16, 7/10, 2/6 |
| mech_consistency raw 0.70 | **[V]** | (26 + ½·8)/43 = 0.698 |
| mech_consistency reconciled 0.73 | **[L]** | ledger; 4 pLDDT exclusions not recomputed |
| structural_agreement thresholded ~0.77 | **[V]** | 81/106 = 0.764 (in ledger range 0.76–0.80) |
| structural_agreement directional 0.773 (51/66) | **[F → diagnosed]** | not reproducible from collapsed columns; ledger counts per-partner multi-chain axes + §5 pLDDT exclusions (100 per-partner / 64 collapsed / 66 script); reconstructions 0.79–0.82; report thresholded 0.764 instead |
| PPI-disruption n=13 | **[V]** | binding-axis destab tokens = 13 |
| fold-detection n=2 | **[V]** | monomer-fold destab tokens = 2 (complex-fold = 11) |
| silence-not-benign n=29 | **[F]** | corrected: 19 non-benign with no destab token |
| GoF reframed as structurally silent | **[V]** | 6/6 GoF mech_consistency = consistent vs expected silent |
| GoF raw label vs metric caveat | **[V]** | H1047R/I500T labeled "destab" but scored consistent-silent |
| Pipeline 2 == P1 on mechanism (both 0.70 raw) | **[V]** | mech_consistency identical: consistent 26 / partial 8 / inconsistent 9 → 0.698 for both |
| Neighborhood tier degrades the gradient | **[V]** | tier OR 6.48 (P1) → 4.00 (nbhd); elevated subset anti-enriched (58% vs 84% path, OR 0.26) |
| Neighborhood specific OR = 0.48 | **[F]** | unmatched (tier OR 4.00, elevated-subset OR 0.26); exact 2×2 needs pinning |
| CHD 45 evaluable, 5/10/5/10/15; 99 unevaluable | **[V]** | collapse from FIXED CSV |
| CHD controls 8/9 structure, 7/9 ddG | **[V]** | per-control vote collapse |
| CHD 8 candidates: KPNA6 I498T at 4/4 (AM-corrected), seven at 3/4; channel split | **[V]** | per-candidate vote breakdown |
| Evaluability = max(monomer, multimer) ≥ 50 | **[V]** | best_plddt = max(...); clean separation at 50; 2 fold-on-binding rescues |
| 99 unevaluable disordered in all contexts | **[V]** | best multimer pLDDT ≤ 45 |
| ZIC3 monomer is zinc-free AlphaFold-3 | **[V]** | no ZN/HETATM in structure; AF3 Server provenance |
| C253S +1.64 vs C297F +12.1 (coordination split) | **[V]** | ddg_monomer per control |
| zic3_gli3 binding artifacts C297F/T323M (0 contacts) | **[V]** | multi_gli3_is_interface / inter_contacts |
