"""
build_chd_config.py — CHD STRUCTURE_CONFIG for the mavis_v7 package.

DROP-IN: add `build_chd_config()` to config.py (alongside `build_benchmark_config()`),
or import from here. Uses the EXACT same SystemConfig / MonomerSpec / MultimerSpec schema.

DESIGN — Option 1 (per-pairing systems), chosen this session:
  CHD's structure is "one hub gene (shroom3 / zic3), many separate pairwise complex PDBs."
  The package MultimerSpec models ONE pdb per system. So each (hub, partner) pairing becomes
  its OWN SystemConfig — structurally identical to a benchmark two-chain system. The hub's
  variants are run once per partner-system. Zero engine changes.

SCOPE (locked this session):
  - DROPPED: shroom3-actb (short ACTB piece — user wants the longer actin chain),
             shroom3-actb_no_bind, shroom3-cdh2_truncated.
  - KEPT: shroom3 × {actin, dvl2, ctnnb1, rock2, cdh2_cyto};  zic3 × {gli3, kpna1, kpna6, mdfi, tcf7l1}.
  - => 10 systems.
  - BIDIRECTIONAL scoring: every gene's variants are scored against the complex(es) they appear in
    (e.g. GLI3 variants vs zic3_gli3 on the GLI3 side). pairwise_partners is bidirectional for this.
  - ZIC5 (14) + ZIC2 (4) FULLY EXCLUDED (no complexes) → CHD N = 163 - 18 = 145 analyzed variants.

INPUT FILE: variants_with_alphamissense_and_franklin_expanded.csv (163 rows). The ONLY candidate with
AlphaMissense (163/163) + franklin (150/163) columns, which the four-way concordance requires.
(variant_comprehensive_v5_2.csv is the 165-row PRE-dedup structural-metrics OUTPUT, no AM/franklin.)
Apply: gene NOT in {zic5, zic2} → 145 rows. (Dedup already reflected in this file: 165→163.)

OFFSETS: the ONLY offset among KEPT systems is CDH2 cytoplasmic (-745) on shroom3_cdh2_cyto, applied
to CDH2 PARTNER positions only. SHROOM3/ZIC3 hub positions use 0 (as-is).

=============================== VERIFY BEFORE RUNNING ===============================
The notebook (Cell 1) used a find_file() search across dirs and CIF-or-PDB per structure.
The package MonomerSpec/MultimerSpec take a single pdb_file resolved under structure_dir.
You MUST confirm, on your Mac, that these files exist under ~/mavis_v8/structures/ (copy
them there from PyCharmMiscProject/structures if needed) and that chains/types are right:

  [V1] Partner monomer files: the notebook listed monomer structures only for some partners.
       For partners WITHOUT a monomer file (rock2/gli3/kpna1/kpna6/mdfi/tcf7l1 are pdb-only;
       actin has NO monomer entry at all), the partner-monomer ΔΔG axis can't run. Since
       variants are in HUBS only, this is fine — but the schema wants a MonomerSpec per gene.
       Approach below: give partners a MonomerSpec pointing at their pdb-only file where one
       exists; for actin (no monomer file) the partner monomer axis is effectively skipped.
       VERIFY each partner monomer filename actually exists, or set to a known-present file.
  [V2] actin structure: kept `fold_shroom3_actin_chain_model_0.pdb` as the COMPLEX. Confirm
       this is the long/biological actin chain (user's intent), not the short ACTB piece.
  [V3] Chain assignments: Cell 1 had all kept complexes as hub=A, partner=B. Confirm against
       the actual PDBs (the SMAD4 lesson — don't trust the label, check the file). Especially
       confirm none of the CHD complexes has a non-standard (swapped) chain order.
  [V4] Structure type: all kept complexes are AlphaFold (fold_*_model_0) → plddt_gate=True,
       no preprocessing. Confirm none needs crystal/NMR handling.
  [V5] CIF for pLDDT: notebook used CIF for pLDDT on shroom3/zic3/dvl2/ctnnb1 (PDB B-factors=0)
       and PDB B-factors for rock2/gli3/kpna*/mdfi/tcf7l1. The package MonomerSpec has a
       `cif_file` field for exactly this — set below where the notebook used CIF. VERIFY the
       package's structure_loading honors cif_file for pLDDT (it has the field; confirm usage).
====================================================================================
"""

from .config import SystemConfig, MonomerSpec, MultimerSpec


def build_chd_config():
    """
    CHD STRUCTURE_CONFIG — 9 per-pairing systems (Option 1).
    Hubs: shroom3, zic3. Variants are in hubs only.
    """
    configs = {}

    # ---- shared hub monomer specs (CIF for pLDDT — PDB B-factors are 0 for these) ----
    shroom3_mono = MonomerSpec("shroom3", "fold_shroom3_model_0.pdb",
                               cif_file="fold_shroom3_model_0.cif")  # [V5] CIF pLDDT
    zic3_mono = MonomerSpec("zic3", "fold_zic3_model_0.pdb",
                            cif_file="fold_zic3_model_0.cif")        # [V5] CIF pLDDT

    # ---- partner monomer specs (best-effort; partner monomer axis is non-critical) ----
    # Partners with their own monomer file (CIF where notebook used it):
    dvl2_mono   = MonomerSpec("dvl2",   "fold_dvl2_model_0.pdb",   cif_file="fold_dvl2_model_0.cif")
    ctnnb1_mono = MonomerSpec("ctnnb1", "fold_ctnnb1_model_0.pdb", cif_file="fold_ctnnb1_model_0.cif")
    # pdb-only partners (PDB B-factors valid):
    rock2_mono  = MonomerSpec("rock2",  "rock2.pdb")
    gli3_mono   = MonomerSpec("gli3",   "gli3.pdb")
    kpna1_mono  = MonomerSpec("kpna1",  "kpna1.pdb")
    kpna6_mono  = MonomerSpec("kpna6",  "kpna6.pdb")
    mdfi_mono   = MonomerSpec("mdfi",   "mdfi.pdb")
    tcf7l1_mono = MonomerSpec("tcf7l1", "tcf7l1.pdb")
    # [V1] actin: NO standalone monomer file exists, and NO actin variants exist, so the actin
    # partner-monomer axis is never used. Point at the complex file to satisfy the schema; it is
    # inert (no variant is in actin). Confirmed with user: no actin monomer needed.
    actin_mono  = MonomerSpec("actin",  "fold_shroom3_actin_chain_model_0.pdb")
    # cdh2 monomer (CIF for pLDDT — PDB B-factors zeroed). CDH2 variants exist (12), scored
    # partner-side against shroom3_cdh2_cyto; monomer axis uses full-length CDH2 numbering.
    cdh2_mono   = MonomerSpec("cdh2",   "fold_cdh2_model_0.pdb",   cif_file="fold_cdh2_model_0.cif")

    def pair(system, hub_gene, hub_mono, partner_gene, partner_mono, pdb_file,
             hub_chain="A", partner_chain="B", cif_file=None):
        return SystemConfig(
            system=system,
            monomers={hub_gene: hub_mono, partner_gene: partner_mono},
            multimer=MultimerSpec(
                system=system,
                pdb_file=pdb_file,
                cif_file=cif_file,
                structure_type="AF",            # [V4] all kept complexes are AlphaFold
                chain_map={hub_gene: hub_chain, partner_gene: partner_chain},  # [V3] verify chains
                pairwise_partners={hub_gene: [partner_chain], partner_gene: [hub_chain]},
                plddt_gate=True,
            ),
            notes=f"CHD per-pairing system: {hub_gene} vs {partner_gene}. Variants in {hub_gene}.",
        )

    # -------------------- SHROOM3 systems (4) --------------------
    configs["shroom3_actin"] = pair(
        "shroom3_actin", "shroom3", shroom3_mono, "actin", actin_mono,
        "fold_shroom3_actin_chain_model_0.pdb", cif_file="fold_shroom3_actin_chain_model_0.cif")            # [V2] confirm = long actin chain

    configs["shroom3_dvl2"] = pair(
        "shroom3_dvl2", "shroom3", shroom3_mono, "dvl2", dvl2_mono,
        "fold_shroom3_dvl2_model_0.pdb", cif_file="fold_shroom3_dvl2_model_0.cif")

    configs["shroom3_ctnnb1"] = pair(
        "shroom3_ctnnb1", "shroom3", shroom3_mono, "ctnnb1", ctnnb1_mono,
        "fold_shroom3_ctnnb1_model_0.pdb", cif_file="fold_shroom3_ctnnb1_model_0.cif")

    configs["shroom3_rock2"] = pair(
        "shroom3_rock2", "shroom3", shroom3_mono, "rock2", rock2_mono,
        "fold_shroom3_rock2_model_0.pdb")

    # cdh2_cyto: CDH2 cytoplasmic-domain construct. CDH2 (the PARTNER here) variant positions use
    # full-length numbering; the construct starts at full-length residue 746, so a -745 offset maps
    # full-length CDH2 position -> construct-local position. Only matters for CDH2-side (reverse)
    # residue queries; SHROOM3-side (hub) is unaffected. position_offsets is SUBTRACTED from CSV pos.
    configs["shroom3_cdh2_cyto"] = SystemConfig(
        system="shroom3_cdh2_cyto",
        monomers={"shroom3": shroom3_mono, "cdh2": cdh2_mono},
        multimer=MultimerSpec(
            system="shroom3_cdh2_cyto",
            pdb_file="fold_shroom3_cdh2_cytoplasmic_domain.pdb",   # [V] confirm exact filename on disk
            structure_type="AF",
            chain_map={"shroom3": "A", "cdh2": "B"},               # [V3] verify chains
            pairwise_partners={"shroom3": ["B"], "cdh2": ["A"]},
            plddt_gate=True,
            position_offsets={"cdh2": 745},   # full-length CDH2 pos - 745 = cyto-construct pos
        ),
        notes="CDH2 cytoplasmic-domain construct (full-length res 746-906). CDH2 partner positions "
              "offset -745. SHROOM3 variants (hub) unaffected.",
    )

    # -------------------- ZIC3 systems (5) --------------------
    configs["zic3_gli3"] = pair(
        "zic3_gli3", "zic3", zic3_mono, "gli3", gli3_mono,
        "fold_zic3_gli3_model_0.pdb", cif_file="fold_zic3_gli3_model_0.cif")

    configs["zic3_kpna1"] = pair(
        "zic3_kpna1", "zic3", zic3_mono, "kpna1", kpna1_mono,
        "fold_zic3_kpna1_model_0.pdb", cif_file="fold_zic3_kpna1_model_0.cif")

    configs["zic3_kpna6"] = pair(
        "zic3_kpna6", "zic3", zic3_mono, "kpna6", kpna6_mono,
        "fold_zic3_kpna6_model_0.pdb", cif_file="fold_zic3_kpna6_model_0.cif")

    configs["zic3_mdfi"] = pair(
        "zic3_mdfi", "zic3", zic3_mono, "mdfi", mdfi_mono,
        "fold_zic3_mdfi_model_0.pdb", cif_file="fold_zic3_mdfi_model_0.cif")

    configs["zic3_tcf7l1"] = pair(
        "zic3_tcf7l1", "zic3", zic3_mono, "tcf7l1", tcf7l1_mono,
        "fold_zic3_tcf7l1_model_0.pdb", cif_file="fold_zic3_tcf7l1_model_0.cif")

    return configs


# ============================================================================
# Notes on the CHD VARIANT INPUT (for when you run it)
# ============================================================================
# The notebook's variant input was `variants_with_alphamissense_and_franklin_expanded.csv`
# (NOT variant_comprehensive_v5_2.csv — that appears to be an intermediate/older file).
# CONFIRM which CSV is the authoritative 163-variant input WITH AlphaMissense + franklin columns,
# since the four-way concordance layer needs those two columns present.
#
# Variants must be tagged with `system` matching the keys above. Because SHROOM3 variants are
# evaluated against MULTIPLE partner-systems, each SHROOM3 variant needs a row per partner-system
# (shroom3_actin, shroom3_dvl2, shroom3_ctnnb1, shroom3_rock2) — i.e. the same fan-out the
# notebook did via its MULTIMER_STRUCTURES loop. Same for ZIC3 across its 5 partners.
# This is the one structural difference from the benchmark input (where each variant = 1 system).
#
# Known CHD quirks to reproduce (from prior results):
#   - DVL2 R367G/R367Q dedup (165 -> 163 variants)
#   - ZIC5 monomer PDB gap (14 variants skipped) -- NOTE: ZIC5 not in this config; confirm whether
#     ZIC5 was a separate hub that needs adding, or was correctly excluded.
#   - GLI3 / TCF7L1 structure-unevaluable (pipeline should call these "Structure unevaluable")
