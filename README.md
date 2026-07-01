# MAVIS — Multimer-Aware Variant Impact Scoring

MAVIS is a FoldX-based structural variant-interpretation pipeline. Most structural
variant-effect tools score a mutation against a single protein in isolation. But
disease genes act in **protein complexes**, and a variant's real structural
consequence often only appears in that multimeric context — at an interface, or in
the fold of a subunit as it sits within its complex. MAVIS's central premise is that
variant disruptiveness should be evaluated **in the appropriate multimer**, and that
doing so resolves the **specific mechanism** of disruption rather than emitting a
single undifferentiated score.

**Try it in your browser (no install):** a guided Colab notebook fetches or predicts the
structures, runs MAVIS, and shows per-variant mechanism cards.

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/la424/mavis/blob/main/notebooks/MAVIS_colab.ipynb)

## What it does

For each missense variant, MAVIS computes a **three-axis ΔΔG** profile against
AlphaFold structures using FoldX, decomposing the structural effect by mechanism:

- `ddg_monomer` — destabilization of the isolated subunit's fold
- `ddg_fold_{partner}` — destabilization of that subunit's fold *within the complex*
- `ddg_binding_{partner}` — disruption of the protein–protein interface itself

This decomposition is what monomer-based tools miss: a variant that looks innocuous
on the lone subunit can be a clear interface disruptor in the assembled complex, and
MAVIS tells you which axis is hit. ΔΔG concordance is **pLDDT-gated** (≥70 strict, ≥50
relaxed) to suppress FoldX artifacts at low-confidence positions, and interface calls
are gated by interface-position pLDDT rather than raw contact count. A **four-way
concordance framework** then integrates the structural tier, FoldX ΔΔG, AlphaMissense,
and Franklin/ClinVar annotations.

**Scope (important, but not the headline):** MAVIS predicts *structural disruption and
its mechanism* — not pathogenicity. Structural disruption overlaps with, but is not
identical to, pathogenicity, so MAVIS reports mechanism and evidence strength
**separately** from phenotype, and "no structural effect detected" is never silently
read as "benign." The benchmark below measures how well its structural calls agree
with literature-grounded structural expectations.

## Repository layout

```
scripts/                 core engine + drivers
  mavis_v7/              the MAVIS package (config, foldx_runner, mechanism,
                          evaluation, metrics, concordance, pipeline, ...)
  run_live.py            benchmark driver (live FoldX run)
  apply_concordance_v5.py   four-way concordance (external tools)
  build_report.py        spreadsheet report
  mavis_v7_baseline_correct.py, *_patch.py, relaxed_regrounding_walk.py
                          post-processing / grading steps (see docs/)
run.py                   generic runner — score ANY genes from a YAML config
new_system.py            scaffold a YAML systems block for new genes
run_chd.py               CHD pipeline driver
prepare_chd_input.py     builds CHD variant input
configs/                 YAML systems configs (chd + benchmark worked examples)
benchmark_variants_v5.csv, chd_input_final.csv   variant inputs
inputs/                  cached intermediates for the no-FoldX self-test + AM table
reference_outputs/       canonical result files: comprehensive CSVs, collapsed CHD,
                          plus MAVIS_results_summary.xlsx (11-sheet overview) + .docx narrative
data/                    reference inputs (UniProt domain ranges, variant–domain map)
docs/                    benchmark ledger, results synthesis, methods, design notes,
                          CHECKPOINT_pre_publication.md
verification/            self-test that reproduces the headline metrics
examples/                getting started
```

## Installation

```bash
pip install -r requirements.txt        # pandas, numpy, biopython, openpyxl, pyyaml
```

Two external dependencies are **not** bundled (see "Inputs"):

1. **FoldX 5.x** — proprietary, free for academics from https://foldxsuite.crg.eu/.
   Download your own copy and point MAVIS at it:
   ```bash
   export FOLDX_BINARY=/path/to/foldx
   ```
2. **AlphaFold structures** — you supply monomer + multimer predictions for your
   proteins (the pipeline consumes structures; it does not predict them).

## The three pipelines

All three share one engine; they differ only in input and in whether external
tools are layered on.

**1. Benchmark** (44 variants / 11 PPI systems)
```bash
export FOLDX_BINARY=/path/to/foldx
python scripts/run_live.py          # -> results/mavis_v7_results.csv
# then concordance + evaluation (see scripts/apply_concordance_v5.py --help and docs/)
```

**2. CHD — full (with external concordance)**
```bash
export FOLDX_BINARY=/path/to/foldx
python run_chd.py                                   # structural results
python scripts/apply_concordance_v5.py --help       # then fold in AlphaMissense + Franklin
```

**3. CHD — structural only (multimer structural results, no external tools)**

This is simply the structural stage on its own — run `run_chd.py` and stop. The
output `results/chd_rerun/chd_structural_results.csv` contains the three ΔΔG axes,
the structural tier, and pLDDT gating, with **no** AlphaMissense / Franklin / ClinVar
columns. Do not run `apply_concordance_v5.py`.
```bash
export FOLDX_BINARY=/path/to/foldx
python run_chd.py
```

## Run it on your own genes

MAVIS isn't limited to the genes above — it's config-driven. To score variants in your own
proteins, describe your complexes in a YAML file and supply AlphaFold structures; no Python
editing required.

```bash
# scaffold a config block (prints the YAML + the structure files you need)
python new_system.py --hub MYGENE --partner PARTNERA --partner PARTNERB

# then run — auto-expands a simple gene,ref_aa,position,alt_aa CSV across systems
export FOLDX_BINARY=/path/to/foldx
python run.py --config my_systems.yaml --variants my_variants.csv \
              --structures ./structures --out results/my_run --dry-run
```

`configs/chd_systems.yaml` is a worked example; `configs/benchmark_systems.yaml` covers the
harder cases (x-ray/NMR, position offsets, non-standard chains, multi-chain complexes). Full
walkthrough: **docs/adding_your_own_genes.md**.

## Quick self-test (no FoldX required)

The framework's headline metrics can be reproduced from cached intermediates
without running FoldX. See **`verification/README.md`** for the exact command; it
runs `verification/verify_stage6.py` against `inputs/intermediate/`.

## Benchmark results

On the 44-variant / 11-PPI-system benchmark (Pipeline 1, t = 2.5, pLDDT-reconciled):

| Metric | Value |
|---|---|
| structural_agreement | **0.77** (threshold sweep 0.76–0.80) |
| mech_consistency | **0.73** (pLDDT-reconciled; 0.70 raw) |

Pipeline 1 (Grantham severity × contact count, **no ΔΔG term**) is a calibrated
structural-evidence-**strength** gradient that complements — rather than competes
with — the ΔΔG mechanism call. The neighborhood/Pipeline-2 variant was tested and
**rejected** (it degraded the gradient); it is retained only as a tested
alternative. The canonical derivation lives in `docs/MAVIS_v7_canonical_benchmark_ledger.md`.

## Inputs

**Variant CSV** (one row per variant). Minimal columns for structural scoring:

| column | meaning |
|---|---|
| `gene` | gene symbol (matches the system config) |
| `ref_aa`, `position`, `alt_aa` | reference AA, 1-based residue, alternate AA |
| `system` | which PPI system / partner set (defined in the config) |

The full concordance step additionally uses `AlphaMissense`, `AlphaMissense_pathogenicity`,
and `franklin` columns. System → partner/structure mappings are defined in
`scripts/mavis_v7/build_chd_config.py` (CHD) and `mavis_v7/config.py` (benchmark);
adapt these for your own proteins.

**Structures.** Place AlphaFold monomer + multimer PDBs under `./structures`,
named per the system config. Monomers are downloadable from the AlphaFold DB by
UniProt ID; multimers must be predicted (AlphaFold Server / ColabFold / AlphaFold-Multimer).

## CHD data provenance & required acknowledgment

The CHD variants analyzed here derive from the **Gabriella Miller Kids First
Pediatric Research Program (Kids First)**, supported by the Common Fund of the
Office of the Director, National Institutes of Health. The data were obtained
through the Kids First Data Resource Center and dbGaP, and their use is governed
by a Kids First / dbGaP Data Use Certification. Only distinct-variant-level,
de-identified data are included in this repository - no per-individual genotypes,
allele counts, or sample identifiers.

> dbGaP study accession: _to be added_ - insert the phs###### for the specific
> Kids First CHD dataset, and verify the exact acknowledgment wording your Data
> Use Certification requires.

## Reproducing the full study

Code, the worked example, and the result CSVs live here. The full AlphaFold
structure set (large) is best archived separately on **Zenodo/Figshare** with a
DOI linked from this README, which also keeps a citable record of the exact inputs.
(GitHub's Zenodo integration can additionally mint a DOI for this code on release.)

## Caveats

- "Bring your own structures + FoldX": MAVIS consumes AlphaFold structures and a
  user-supplied FoldX binary; it does not generate structures.
- Variants in disordered / low-pLDDT regions are structurally **unevaluable**
  regardless of substitution severity — this is reported, not silently dropped.

## Citation

See `CITATION.cff`.

## License

MIT — see `LICENSE`.
