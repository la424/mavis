# Adding your own genes

MAVIS is config-driven: you describe your complexes in a YAML file and provide the
AlphaFold structures. No Python editing is required.

## What you need

- **AlphaFold structures** for each gene and each hub-partner complex (you supply these;
  MAVIS consumes structures, it does not predict them).
- **FoldX 5.x** (free for academics, https://foldxsuite.crg.eu/), exposed as
  `export FOLDX_BINARY=/path/to/foldx`.
- A **variant CSV** with columns `gene,ref_aa,position,alt_aa`.

## Step 1 — structures

For every complex you want to score, obtain the **multimer** prediction (hub + partner
together) and a **monomer** prediction for each gene. Put them all in one directory
(passed as `--structures`). The default naming the loader and scaffold use is
`fold_<gene>_model_0.pdb` for monomers and `fold_<hub>_<partner>_model_0.pdb` for
complexes, but files can be named anything as long as the YAML points at the right names.
If your PDB B-factor column is zeroed (common for some AlphaFold exports), drop the
matching `.cif` next to each `.pdb` and reference it via `cif_file` / `monomer_cif` —
MAVIS reads pLDDT from the CIF in that case.

## Step 2 — the systems config (YAML)

Scaffold a starting point:

    python new_system.py --hub MYGENE --partner PARTNERA --partner PARTNERB

That prints a ready-to-fill YAML block and lists the structure filenames it expects. Save
it as e.g. `my_systems.yaml`. A minimal complex:

    systems:
      mygene_partnera:
        structure_type: AF
        pdb_file: fold_mygene_partnera_model_0.pdb
        cif_file: fold_mygene_partnera_model_0.cif   # optional
        genes:
          mygene:
            chain: A
            monomer_pdb: fold_mygene_model_0.pdb
          partnera:
            chain: B
            monomer_pdb: fold_partnera_model_0.pdb

`configs/chd_systems.yaml` is a full real-world example; `configs/benchmark_systems.yaml`
shows the harder cases (x-ray/NMR complexes with HETATM handling, position offsets,
non-standard chain order, and a 4-chain homo-tetramer). The complete field reference is
the docstring at the top of `scripts/mavis_v7/config_io.py`.

Fields worth knowing:
- `chain:` is the gene's chain **in the complex** (drives the interface analysis).
- `monomer_offset` / `multimer_offset` shift CSV residue numbers to match a construct that
  doesn't start at residue 1 (the value is subtracted from the CSV position).
- For x-ray/NMR complexes set `structure_type: xray|nmr`, `plddt_gate: false`, and the
  HETATM keep/strip lists (see the benchmark example).
- For complexes with repeated chains, give `chain_map` and `pairwise_partners` explicitly
  (see the hemoglobin tetramer in the benchmark example).

**Verify chain assignments against the actual PDB** — don't trust the file's chain labels
blindly. (Two benchmark complexes have a non-standard chain order; getting this wrong
silently scores the wrong interface.)

## Step 3 — your variants

A plain CSV is enough:

    gene,ref_aa,position,alt_aa
    MYGENE,G,35,V
    PARTNERA,R,180,Q

You do **not** tag systems yourself. `run.py` auto-expands each variant into one row per
system its gene appears in (hub or partner side). Keep any extra columns (e.g.
`AlphaMissense`, `franklin`) if you plan to run the four-way concordance — they're carried
through.

## Step 4 — run

    export FOLDX_BINARY=/path/to/foldx
    python run.py --config my_systems.yaml --variants my_variants.csv \
                  --structures /path/to/structures --out results/my_run

Add `--dry-run` first to validate the config and confirm every referenced structure is
present (and preview the fan-out) without launching FoldX. Structural results land in
`results/my_run/structural_results.csv` — the **structural-only** output (three ΔΔG axes,
structural tier, pLDDT gating; no external tools).

## Step 5 — four-way concordance (optional)

To layer AlphaMissense + Franklin/ClinVar on top, include `AlphaMissense`,
`AlphaMissense_pathogenicity`, and `franklin` columns in your variant CSV, then run
`scripts/apply_concordance_v5.py` on the structural output (see its `--help`).

## What MAVIS still can't do for you

- Predict structures (bring your own AlphaFold models).
- Guess chain assignments (verify them in the PDB).
- Supply FoldX (proprietary; you install it).
