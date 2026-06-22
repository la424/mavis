"""
STRUCTURE_CONFIG schema and validation.

Each system entry defines:
  - Structure type (AF / xray / nmr) which controls pLDDT gating
  - Monomer structure(s) per gene
  - Multimer structure with chain map + pairwise partner list
  - HETATM policy for preprocessing (crystal/NMR only)
  - Optional position offsets per construct (e.g. CDH2 numbering)

The config is consumed by:
  - preprocessing.py   — uses hetatm_* and nmr_model fields
  - structure_loading.py — uses chain_map, plddt_gate
  - foldx_runner.py    — uses chain_map, pairwise partners, position_offsets
  - mechanism.py       — uses plddt_gate (crystal/NMR bypass)
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple


# ============================================================================
# Global thresholds (unchanged from v6)
# ============================================================================
TIER_THRESHOLDS = [(5.0, "Tier 1"), (3.0, "Tier 2"), (1.5, "Tier 3")]
DEFAULT_TIER = "Tier 4"

# DDG boundaries (kcal/mol)
DDG_NEUTRAL = 1.0       # |DDG| < 1.0 → neutral
DDG_DESTAB = 1.0        # DDG > 1.0 → destabilizing
DDG_HIGHLY = 2.0        # DDG > 2.0 → highly destabilizing
# mild_destab window for evaluation: 0.5 - 2.0 (overlaps neutral at 0.5-1.0)
DDG_MILD_LOW = 0.5
DDG_MILD_HIGH = 2.0

# pLDDT thresholds
PLDDT_STRICT = 70
PLDDT_RELAXED = 50

# Multi-metric structural score thresholds
DISRUPTION_POINTS = [(20, 4.0), (10, 3.0), (4, 2.0), (1, 1.0)]
CONTACT_DRIVEN_THRESHOLD = 6  # 75th percentile


# ============================================================================
# Per-system config schema
# ============================================================================
@dataclass
class MonomerSpec:
    """Monomer structure for one gene."""
    gene: str
    pdb_file: str                  # AF: fold_{gene}_model_0.pdb
    cif_file: Optional[str] = None  # optional CIF for CIF-only pLDDT
    structure_type: str = "AF"      # "AF" | "xray" | "nmr"
    plddt_gate: bool = True         # False for crystal/NMR
    chain_id: str = "A"             # default AF chain
    position_offset: int = 0        # SUBTRACT from CSV position to get structure position
                                    #   CDH2 truncated: +159  (CSV 200 → struct 41)
                                    #   HBB AF monomer: -1    (CSV 6 mature → struct 7 UniProt)


@dataclass
class MultimerSpec:
    """Multimer complex config for one system."""
    system: str
    pdb_file: str                                    # PDB filename (may be _processed)
    structure_type: str                              # "AF" | "xray" | "nmr"
    chain_map: Dict[str, str]                        # gene -> chain_id
    pairwise_partners: Dict[str, List[str]]          # gene -> [partner_chains]
    plddt_gate: bool                                 # False for crystal/NMR
    position_offsets: Dict[str, int] = field(default_factory=dict)  # gene -> offset (subtract from CSV position)
    cif_file: Optional[str] = None                   # optional CIF for CIF-only pLDDT (AF complexes)

    # Preprocessing fields (used only if structure_type in {"xray","nmr"})
    source_file: Optional[str] = None                # original file before preprocessing
    nmr_model: Optional[int] = None
    hetatm_keep_all: List[str] = field(default_factory=list)
    hetatm_keep_per_chain: Dict[str, List[str]] = field(default_factory=dict)
    hetatm_strip_all: List[str] = field(default_factory=list)
    preprocessing_notes: str = ""


@dataclass
class SystemConfig:
    """Complete config for one benchmark system."""
    system: str
    monomers: Dict[str, MonomerSpec]   # gene -> MonomerSpec
    multimer: MultimerSpec
    notes: str = ""


# ============================================================================
# Validation
# ============================================================================
def validate_config(configs: Dict[str, SystemConfig], structure_dir: Path,
                    preprocessed_dir: Path) -> Tuple[bool, List[str]]:
    """
    Verify all referenced files exist and chain maps are self-consistent.
    Returns (all_ok, list_of_issues).
    """
    structure_dir = Path(structure_dir)
    preprocessed_dir = Path(preprocessed_dir)
    issues = []

    for sys_name, cfg in configs.items():
        # Monomer files
        for gene, mspec in cfg.monomers.items():
            pdb = structure_dir / mspec.pdb_file
            if not pdb.exists():
                issues.append(f"[{sys_name}] monomer {gene}: file not found: {pdb}")

        # Multimer file — check preprocessed first for crystal/NMR, then raw
        multi = cfg.multimer
        if multi.structure_type in ("xray", "nmr"):
            expected_processed = preprocessed_dir / multi.pdb_file
            expected_source = structure_dir / (multi.source_file or multi.pdb_file)
            if not expected_processed.exists() and not expected_source.exists():
                issues.append(f"[{sys_name}] multimer: neither {expected_processed} nor {expected_source} exists")
        else:
            pdb = structure_dir / multi.pdb_file
            if not pdb.exists():
                issues.append(f"[{sys_name}] multimer: file not found: {pdb}")

        # Chain map consistency with pairwise partners
        for gene, partners in multi.pairwise_partners.items():
            if gene not in multi.chain_map:
                issues.append(f"[{sys_name}] gene '{gene}' in pairwise_partners but not in chain_map")
            for p_chain in partners:
                if p_chain not in multi.chain_map.values():
                    issues.append(f"[{sys_name}] pairwise partner chain '{p_chain}' not in chain_map values")

    return len(issues) == 0, issues


# ============================================================================
# Benchmark instance — the 11 systems, 44 variants
# ============================================================================
def build_benchmark_config() -> Dict[str, SystemConfig]:
    """
    Build the MAVIS benchmark STRUCTURE_CONFIG (11 systems).
    Preprocessed crystal/NMR filenames have _processed suffix.
    """
    configs = {}

    # -------- brca1_bard1 (NMR, preprocessed) --------
    configs["brca1_bard1"] = SystemConfig(
        system="brca1_bard1",
        monomers={
            "brca1": MonomerSpec("brca1", "fold_brca1_model_0.pdb"),
            "bard1": MonomerSpec("bard1", "fold_bard1_model_0.pdb"),
        },
        multimer=MultimerSpec(
            system="brca1_bard1",
            pdb_file="1JM7_processed.pdb",
            source_file="1JM7.pdb",
            structure_type="nmr",
            chain_map={"brca1": "A", "bard1": "B"},
            pairwise_partners={"brca1": ["B"], "bard1": ["A"]},
            plddt_gate=False,
            nmr_model=1,
            hetatm_keep_all=["ZN"],
            hetatm_strip_all=["HOH"],
            preprocessing_notes="NMR MODEL 1 extracted. ZN retained (RING domain C3HC4).",
        ),
        notes="RING domain zinc coordination; NMR geometry may increase FoldX SD."
    )

    # -------- pi3k (AF) --------
    configs["pi3k"] = SystemConfig(
        system="pi3k",
        monomers={
            "pik3ca": MonomerSpec("pik3ca", "fold_pik3ca_model_0.pdb"),
            "pik3r1": MonomerSpec("pik3r1", "fold_pik3r1_model_0.pdb"),
        },
        multimer=MultimerSpec(
            system="pi3k",
            pdb_file="fold_pik3ca_pik3r1_model_0.pdb",
            structure_type="AF",
            chain_map={"pik3ca": "A", "pik3r1": "B"},
            pairwise_partners={"pik3ca": ["B"], "pik3r1": ["A"]},
            plddt_gate=True,
        ),
    )

    # -------- mlh1_pms2 (AF) --------
    configs["mlh1_pms2"] = SystemConfig(
        system="mlh1_pms2",
        monomers={
            "mlh1": MonomerSpec("mlh1", "fold_mlh1_model_0.pdb"),
            "pms2": MonomerSpec("pms2", "fold_pms2_model_0.pdb"),
        },
        multimer=MultimerSpec(
            system="mlh1_pms2",
            pdb_file="fold_mlh1_pms2_model_0.pdb",
            structure_type="AF",
            chain_map={"mlh1": "A", "pms2": "B"},
            pairwise_partners={"mlh1": ["B"], "pms2": ["A"]},
            plddt_gate=True,
        ),
    )

    # -------- msh2_msh6 (AF) --------
    configs["msh2_msh6"] = SystemConfig(
        system="msh2_msh6",
        monomers={
            "msh2": MonomerSpec("msh2", "fold_msh2_model_0.pdb"),
            "msh6": MonomerSpec("msh6", "fold_msh6_model_0.pdb"),
        },
        multimer=MultimerSpec(
            system="msh2_msh6",
            pdb_file="fold_msh2_msh6_model_0.pdb",
            structure_type="AF",
            chain_map={"msh2": "A", "msh6": "B"},
            pairwise_partners={"msh2": ["B"], "msh6": ["A"]},
            plddt_gate=True,
        ),
    )

    # -------- vhl_elonginc (AF, non-standard chain order) --------
    configs["vhl_elonginc"] = SystemConfig(
        system="vhl_elonginc",
        monomers={
            "vhl":   MonomerSpec("vhl", "fold_vhl_model_0.pdb"),
            "tceb1": MonomerSpec("tceb1", "fold_tceb1_model_0.pdb"),
        },
        multimer=MultimerSpec(
            system="vhl_elonginc",
            pdb_file="fold_vhl_tceb1_model_0.pdb",
            structure_type="AF",
            chain_map={"vhl": "B", "tceb1": "A"},   # NON-STANDARD
            pairwise_partners={"vhl": ["A"], "tceb1": ["B"]},
            plddt_gate=True,
        ),
        notes="Non-standard chain order: VHL=B, ElonginC=A."
    )

    # -------- kras_craf (X-ray, preprocessed) --------
    configs["kras_craf"] = SystemConfig(
        system="kras_craf",
        monomers={
            "kras": MonomerSpec("kras", "fold_kras_model_0.pdb"),
            "raf1": MonomerSpec("raf1", "fold_raf1_model_0.pdb"),
        },
        multimer=MultimerSpec(
            system="kras_craf",
            pdb_file="6XI7_processed.pdb",
            source_file="6XI7.pdb",
            structure_type="xray",
            chain_map={"kras": "A", "raf1": "B"},
            pairwise_partners={"kras": ["B"], "raf1": ["A"]},
            plddt_gate=False,
            hetatm_keep_per_chain={"B": ["ZN"]},
            hetatm_strip_all=["GNP", "MG", "SO4", "CL", "HOH"],
            preprocessing_notes="GNP removed (FoldX cannot handle). KRAS G12 monomer DDG low-reliability.",
        ),
    )

    # -------- hemoglobin_dimer (X-ray, extracted+preprocessed) --------
    configs["hemoglobin_dimer"] = SystemConfig(
        system="hemoglobin_dimer",
        monomers={
            # HBB/HBA1: AF uses UniProt numbering (incl. initiator Met), CSV uses
            # mature numbering (Met cleaved). Offset = -1 means struct_pos = csv_pos + 1.
            "hbb":  MonomerSpec("hbb", "fold_hbb_model_0.pdb", position_offset=-1),
            "hba1": MonomerSpec("hba1", "fold_hba1_model_0.pdb", position_offset=-1),
        },
        multimer=MultimerSpec(
            system="hemoglobin_dimer",
            pdb_file="2HHB_dimer_AB_processed.pdb",
            source_file="2HHB_dimer_AB.pdb",
            structure_type="xray",
            chain_map={"hbb": "B", "hba1": "A"},
            pairwise_partners={"hbb": ["A"], "hba1": ["B"]},
            plddt_gate=False,
            hetatm_strip_all=["HEM", "HOH", "PO4"],
            preprocessing_notes="HEM removed. E6V >8A from HEM — geometry unaffected.",
        ),
    )

    # -------- hemoglobin_tetramer (X-ray, 4 chains, PAIRWISE AnalyseComplex) --------
    configs["hemoglobin_tetramer"] = SystemConfig(
        system="hemoglobin_tetramer",
        monomers={
            "hbb":  MonomerSpec("hbb", "fold_hbb_model_0.pdb", position_offset=-1),
            "hba1": MonomerSpec("hba1", "fold_hba1_model_0.pdb", position_offset=-1),
        },
        multimer=MultimerSpec(
            system="hemoglobin_tetramer",
            pdb_file="2HHB_processed.pdb",
            source_file="2HHB.pdb",
            structure_type="xray",
            chain_map={"hba1_1": "A", "hbb": "B", "hba1_2": "C", "hbb_2": "D"},
            # For HBB variants (chain B), run AnalyseComplex against A, C, AND D separately
            # All keys must match chain_map keys (gene or instance label)
            pairwise_partners={
                "hbb":    ["A", "C", "D"],   # B vs A (α1), B vs C (α2), B vs D (β2)
                "hba1_1": ["B", "D"],        # A vs B, A vs D  (if α1 had variants)
                "hba1_2": ["B", "D"],        # C vs B, C vs D
                "hbb_2":  ["A", "B", "C"],   # D vs A, D vs B, D vs C
            },
            plddt_gate=False,
            hetatm_strip_all=["HEM", "HOH", "PO4"],
            preprocessing_notes="HEM removed. N102T excluded from Level 3 correlation (4.4A from HEM).",
        ),
        notes="4-chain tetramer. Pairwise AnalyseComplex per partner chain."
    )

    # -------- troponin_ic (AF) --------
    configs["troponin_ic"] = SystemConfig(
        system="troponin_ic",
        monomers={
            "tnni3": MonomerSpec("tnni3", "fold_tnni3_model_0.pdb"),
            "tnnc1": MonomerSpec("tnnc1", "fold_tnnc1_model_0.pdb"),
        },
        multimer=MultimerSpec(
            system="troponin_ic",
            pdb_file="fold_tnni3_tnnc1_model_0.pdb",
            structure_type="AF",
            chain_map={"tnni3": "A", "tnnc1": "B"},
            pairwise_partners={"tnni3": ["B"], "tnnc1": ["A"]},
            plddt_gate=True,
        ),
    )

    # -------- cam_cav12 (AF) --------
    configs["cam_cav12"] = SystemConfig(
        system="cam_cav12",
        monomers={
            "calm1":   MonomerSpec("calm1", "fold_calm1_model_0.pdb"),
            "cacna1c": MonomerSpec("cacna1c", "fold_cacna1c_model_0.pdb"),
        },
        multimer=MultimerSpec(
            system="cam_cav12",
            pdb_file="fold_calm1_cacna1c_model_0.pdb",
            structure_type="AF",
            chain_map={"calm1": "A", "cacna1c": "B"},
            pairwise_partners={"calm1": ["B"], "cacna1c": ["A"]},
            plddt_gate=True,
        ),
        notes="CACNA1C IQ motif region; partial pLDDT (~77 mean in region)."
    )

    # -------- smad4_smad3 (AF, non-standard chain order) --------
    configs["smad4_smad3"] = SystemConfig(
        system="smad4_smad3",
        monomers={
            "smad4": MonomerSpec("smad4", "fold_smad4_model_0.pdb"),
            "smad3": MonomerSpec("smad3", "fold_smad3_model_0.pdb"),
        },
        multimer=MultimerSpec(
            system="smad4_smad3",
            pdb_file="fold_smad4_smad3_model_0.pdb",
            structure_type="AF",
            chain_map={"smad4": "B", "smad3": "A"},   # NON-STANDARD
            pairwise_partners={"smad4": ["A"], "smad3": ["B"]},
            plddt_gate=True,
        ),
        notes="Non-standard chain order: SMAD4=B, SMAD3=A."
    )

    return configs
