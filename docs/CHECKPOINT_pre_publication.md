# MAVIS — Pre-Publication Checkpoint

**Date:** 2026-07-01  **Scope:** both papers (methods/benchmark + CHD application), evaluated in parallel.

**Verdict:** the pipeline is in strong shape and reproducible at the concordance level. One CHD takeaway
moved with the AlphaMissense parse fix (§5); benchmark takeaways are unchanged. Two items remain before
submission — one substantive (the monomer-fold cohort), one production-only.

## 1. Inputs — complete
- **Benchmark:** `benchmark_variants_v5.csv` (44 variants / 11 PPI systems).
- **CHD:** `chd_input_final.csv`, `chd_input_per_system.csv` (bidirectional per-system),
  `data/MAVIS_CHD_variant_domain_interface_map.csv`.
- **Annotations:** `variants_with_alphamissense_and_franklin_expanded.csv` (AlphaMissense + Franklin,
  transposition-corrected), `inputs/AM_variants_mavis_mechanism_test.xlsx`.
- **Structures:** AlphaFold-3 Server (multimers); experimental PDBs where AF2 confidence is insufficient
  (2HHB / HBB, 6XI7 / KRAS-RAF1, 1JM7 / BRCA1-BARD1). `structures/` is gitignored (obtain per docs).

## 2. Pipeline — works as intended, with known limits
The `mavis_v7` engine is byte-identical across the working and release trees; the concordance engine
round-trips byte-for-byte on the canonical CHD file; `scripts/apply_chd_concordance.py` is idempotent on
already-correct input.

Remaining limits:
- **Monomer-fold axis undersupported — n=2 graded destabilizers.** The methods paper's one substantive
  gating task (see §7).
- **Full FoldX-level reproducibility gap.** The pLDDT-reconciled structural layer was not persisted, so
  outputs regenerate from the locked structural CSV (concordance level), not from raw FoldX. State in methods.
- **Directional `structural_agreement` is not reproducible from released columns** → report thresholded 0.77.

## 3. Outputs — inventory
| Artifact | Path | Notes |
|---|---|---|
| Benchmark comprehensive | `reference_outputs/mavis_v7_concordance_annotated.csv` | 44 variants, per-axis, 1255 cols (canonical; formerly `v5_reconciled`) |
| CHD comprehensive | `reference_outputs/chd_concordance_results_FIXED.csv` | 384 rows × 213 cols |
| CHD collapsed | `reference_outputs/chd_concordance_collapsed.csv` + `MAVIS_CHD_concordance_collapsed.xlsx` | 144 variants, one row each |
| Per-variant overview | `reference_outputs/MAVIS_results_summary.xlsx` | 11 sheets (benchmark + CHD per-variant, distributions, candidates, control recall, orthogonal cases) |
| Narrative summary | `reference_outputs/MAVIS_results_summary.docx` | prose writeup + matching tables |
| Benchmark takeaways | `docs/MAVIS_v7_canonical_benchmark_ledger.md` | locked metrics, per-system notes, principles |
| Two-paper synthesis | `docs/MAVIS_results_synthesis.md` | benchmark + CHD results, interpretation, verification ledger |
| Methods / decisions | `docs/methods_metrics_sketch.md`, `docs/design_decisions.md` | |

## 4. Locked headline results
**Benchmark** (P1, t=2.5, pLDDT-reconciled): structural_agreement **0.77** (thresholded — primary);
mech_consistency **0.73** reconciled / 0.70 raw; tier gradient **100 / 81 / 70 / 33** (T1–T4);
AlphaMissense accuracy on confident calls 37/41 = 90%. Claim arms: PPI-disruption **13**, complex-fold
**11**, monomer-fold **2**, silence-not-benign **19**. All 6 GoF variants reframed as structurally silent.
Pipeline 2 (neighborhood) equals P1 on mechanism and degrades the pathogenicity gradient → tested-and-rejected.

**CHD:** 144 variants (133 patient + 11 ZIC3 controls); 45 evaluable, 99 unevaluable (disordered in every
context). Controls: 8/9 structure, 7/9 ddG. **Patient candidates: 8** — KPNA6 I498T at 4/4 (all four
channels), seven at 3/4. ROCK2 T367M is the single clean MAVIS-vs-AM disagreement; DVL2 D441Y the cleanest
novel candidate; H286R correctly scored silent (trafficking / NLS mechanism).

## 5. Did the recent AlphaMissense fix change any takeaways?
- **Benchmark: no.** The fix was CHD-only (kpna1 / kpna6 / tcf7l1); no benchmark gene, metric, or claim moved.
- **CHD: one.** Correcting the whole-gene AlphaMissense score→class transposition restored KPNA6 I498T's
  AM hit, lifting it 3/4 → **4/4** (`concordance` metric) and adding KPNA6 K424N as an 8th candidate
  (2/4 → 3/4). I498T is now AM-concordant (previously logged as AM-discordant). The central CHD thesis is
  intact (ROCK2 T367M, DVL2 D441Y, structural disruption ≠ pathogenicity); under strict pathogenic-only
  Franklin the novel candidates still sit at 3/4. Reflected across the CSVs, the summary workbook, and both
  narrative docs.

## 6. Resolved at this checkpoint (consistency pass)
1. Report thresholded `structural_agreement` **0.77** as primary; directional 0.773 superseded (not
   reproducible from released columns; ledger §6 updated).
2. Silence-not-benign = **19** (verified): non-benign variants with no `destab` token on any structural
   axis; the "29" was the all-silent count, which also counted the 10 benign variants.
3. Pipeline 2 = tested-and-rejected; report tier OR 6.48 → 4.00 and elevated-subset OR 0.26; the ledger's
   "OR 0.48" is superseded (does not reproduce under either natural 2×2).
4. Canonical benchmark file = `mavis_v7_concordance_annotated.csv` (formerly `v5_reconciled`, absent from
   disk; confirmed the same artifact — carries the full raw + pLDDT-reconciled mech_consistency columns).

## 7. Open before submission
- **[Substantive] Expand the monomer-fold destabilizer cohort** (2 → 6–8 graded). Leads BRCA1 V11G,
  MLH1 Q542L, SMAD4 R361C, TNNI3 R162W were confirmed non-gradeable; productive path = unstable-hemoglobin
  variants (HBB already on 2HHB) and/or a new monomer-fold system. Reopens the recompute; gating for the
  methods paper.
- **[Production] CHD per-axis evaluability gating refinement** — gate each axis on its own context's pLDDT
  (monomer-fold ← monomer; complex/binding ← complex-position). No effect on the current 45/99 split;
  matters on the full production run.

## 8. Repository state
`la424/mavis`, branch `main`, in sync with GitHub as of this checkpoint. The FoldX binary and AlphaFold
structures are gitignored — users supply their own (see `README.md` and `docs/`).
