# MAVIS v7 — Canonical Benchmark Ledger (post-verification, post-relaxed-tier)

**Purpose.** The single source of truth for what the 44-variant benchmark *is* and why each call was
made. Supersedes the scattered append blocks for reference purposes (the append blocks remain the
detailed derivation record). For a methods paper where the benchmark is the contribution, this
ledger is a core deliverable, not hygiene.

## 1. Composition
- **44 variants / 11 PPI systems.** Roles (post mech-control retirement): benign 10, pathogenic 26,
  pathogenic_gof 6, pathogenic_lof 2.
- **No `mechanism_control` category** — the 11 originally-mislabeled rows are *mislabel corrected*
  (not reassigned): KRAS G12D/G12V/Q61H + PIK3CA H1047R → `pathogenic_gof`; BRCA1 C61G, MLH1 R755W,
  MSH2 G674R/C697F, VHL W117R/Y98H, HBB E6V → `pathogenic`.

## 2. Three ground-truth tiers (per axis: monomer, complex-fold, binding)
- **STRICT** (`expected_ddg_{axis}`): primary-literature-grounded. Token ∈ {stab, neutral, destab,
  unknown}. `unknown` = no groundable direct measurement (Definition B; FoldX self-prediction and
  bare structural-position inference are NOT admissible grounding).
- **RELAXED** (`expected_ddg_{axis}_relaxed`): strict + directional promotions of contracted axes
  that have indirect-but-specific support (coupled readout | sibling analogy | off-axis measurement
  that speaks to THIS axis). Structural-position inference excluded. One tier.
- **FUNCTIONAL** (not yet a column): off-axis quantitative data that doesn't speak to any structural
  axis (e.g. CaM Ca²⁺-affinity) — noted in evidence, reserved for a possible future functional axis.

## 3. Verification impact (the before/after)
- Graded structural axes: **96 → 80** (net −16; gross 18 contracted to unknown, 2 recovered to neutral).
- Disruptor tokens (destab/stab/mild_destab): **39 → 26** (a third removed as ungroundable or
  wrong-direction).
- Roles reassigned: **11** (mech-control retirement).
- **structural_agreement essentially unchanged (0.775 → 0.773 @ t2.5):** the benchmark was corrected
  without moving the headline — the prior number was right for partly-wrong reasons, now right for
  right reasons.

## 4. The 18 contracted axes — relaxed disposition (6 promoted / 12 stay-unknown)
| variant | axis | strict | relaxed | basis |
|---------|------|--------|---------|-------|
| MSH2 A636P | mono | unknown | **destab** | COUPLED: instability/reduced-expression (Ollila) |
| MSH2 A636P | fold | unknown | unknown | no independent directional readout |
| MSH2 C697F | mono | unknown | unknown | no directional readout |
| MSH2 C697F | fold | unknown | unknown | no directional readout |
| MSH2 N127S | mono | unknown | **neutral** | COUPLED: MMR-proficient functional control |
| MSH2 N127S | fold | unknown | **neutral** | COUPLED: preserved assembly *(unevaluable: no partner)* |
| MSH2 G322D | mono | unknown | **neutral** | COUPLED: near-neutral polymorphism |
| MSH2 G322D | fold | unknown | **neutral** | COUPLED: preserved assembly *(unevaluable: no partner)* |
| VHL W117R | mono | unknown | unknown | structural-position only (excluded) |
| VHL W117R | fold | unknown | unknown | structural-position only (excluded) |
| VHL W117R | bind | unknown | unknown | not in Kishida panel; HIF-axis lesion |
| VHL Y98H | mono | unknown | **neutral** | COUPLED: retained VBC assembly proxy |
| TNNI3 R162W | bind | unknown | unknown | inference-only (Zhou 2013 not in hand) |
| CaM D96V | mono | unknown | unknown | off-axis Ca²⁺ ≠ fold (NMR: fold intact) |
| CaM N98S | mono | unknown | unknown | off-axis Ca²⁺ ≠ fold |
| CaM F142L | mono | unknown | unknown | off-axis Ca²⁺ ≠ fold |
| SMAD4 D351H | fold | unknown | unknown | Shi = oligomerization (binding axis), not fold |
| SMAD4 R361C | fold | unknown | unknown | Shi = oligomerization, not fold |

Promotions: 6 (1 destab, 5 neutral), all COUPLED-class. Conservative by design.

## 5. pLDDT interface-exclusions (4 variants — apply to BOTH metrics, BOTH pipelines)
| variant | partner | pLDDT | excluded axis | batch |
|---------|---------|-------|---------------|-------|
| TNNI3 R145G | tnnc1 | 54.93 | fold+binding | B9 |
| TNNI3 R145Q | tnnc1 | 54.93 | fold+binding | B9 |
| CaM D96V | cacna1c | 69.45 | fold+binding | B10 |
| SMAD4 D351H | smad3 | 57.37 | binding | B11 |

## 6. Locked headline metrics (Pipeline 1, t=2.5, pLDDT-reconciled)
- **structural_agreement (strict, thresholded sweep):** 0.77 @ t2.5 (range 0.76–0.80).
- **structural_agreement (directional, 1.0 floor):** strict 0.757... → reported as the common ruler
  for the strict/relaxed comparison: **strict 0.773 (51/66) / relaxed 0.757 (53/70).**
- **mech_consistency:** 0.70 raw → **0.73 pLDDT-reconciled** @ t2.5.
- **Pipeline 2:** pending local run (see runbook).

## 7. Per-system one-liners (mechanism character)
- **MLH1-PMS2** (Kosinski): mixed; L749P fold→unknown corrected.
- **HBB tetramer/dimer** (Kiger): W37 series R-state contacts — FoldX FN against T-state 2HHB.
- **KRAS-CRAF** (Hunter, closed): G12D/G12V GoF, low-RAF-affinity; Q61H direction was uncertain.
- **BRCA1-BARD1**: C61G (RING/Zn) mislabel-corrected to pathogenic.
- **PIK3CA-PIK3R1**: E545K/E542K released-autoinhibition GoF → FoldX reads as stabilizing.
- **VHL-ElonginC**: L158Q true interface disruptor; W117R fully contracted (HIF-axis); Y98H clean
  ElonginC-silent control.
- **MSH2-MSH6**: regulatory/DNA-binding system, zero true partner disruptors; A636P Walker-A.
- **TNNI3-TNNC1**: regulatory system, zero partner disruptors; R145G/Q binding preserved.
- **CaM-Cav1.2**: Ca²⁺-sensing system; CDI loss with interface preserved-to-strengthened; FoldX
  correctly "stabilizing/neutral" yet blind to the sensing lesion.
- **SMAD4-SMAD3**: R361C clean interface TP; D351H interface FN (AF2 mis-modeled, pLDDT-excluded);
  I500T/V structurally-silent post-translational GoF.

## 8. Cross-cutting principles confirmed by the campaign
- Three-axis separation is load-bearing; "right for the wrong reasons" is a first-class failure mode.
- `unknown` is a valid value; primary literature required; FoldX self-prediction and structural-
  position inference are not admissible grounding (strict tier).
- Off-axis evidence promotes a relaxed token ONLY if it speaks to that axis's direction.
- pLDDT gating applied uniformly across structural_agreement and mech_consistency, both pipelines.

## 9. Artifacts (this session)
- `mavis_v7_results_corrected_v2.csv` — strict tokens, post-recompute.
- `mavis_v7_results_relaxed.csv` — adds parallel `expected_ddg_*_relaxed` columns.
- `mavis_v7_concordance_v5_reconciled.csv` — P1 mech_consistency raw + pLDDT-reconciled, all thresholds.
- `batch12_recompute_bundle.py`, `mech_consistency_plddt_patch.py`, `relaxed_regrounding_walk.py` —
  the auditable edit/grade/promotion logic.
- Reports: Batch-12 cross-system, Track-B mech_consistency, reconciliation summary, P2 runbook,
  this ledger.
