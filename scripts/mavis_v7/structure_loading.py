"""
Structure loading, pLDDT extraction, contact/interface computation.

Mostly unchanged from v6, but:
  - pLDDT gating is now config-aware: crystal/NMR structures report pLDDT=100
    uniformly, and downstream gates treat them as "confident" regardless.
  - Monomer lookup is per-system via MonomerSpec.
"""

from pathlib import Path
from typing import Dict, Optional, Tuple, Set
import warnings
warnings.filterwarnings("ignore")

from Bio.PDB import PDBParser, MMCIFParser, NeighborSearch, ShrakeRupley

from .constants import THREE_TO_ONE, MAX_SASA
from .config import SystemConfig, MonomerSpec, MultimerSpec

_pdb_parser = PDBParser(QUIET=True)
_cif_parser = MMCIFParser(QUIET=True)


# ============================================================================
# Low-level loaders
# ============================================================================
def load_pdb(path):
    if path and Path(path).exists():
        try:
            return _pdb_parser.get_structure("s", str(path))
        except Exception:
            pass
    return None


def load_cif(path):
    if path and Path(path).exists():
        try:
            return _cif_parser.get_structure("s", str(path))
        except Exception:
            pass
    return None


# ============================================================================
# pLDDT
# ============================================================================
CRYSTAL_PLDDT_SENTINEL = 100.0


def get_plddt(structure, chain_id, plddt_gate: bool = True):
    """
    Per-residue pLDDT from B-factor column.
    If plddt_gate is False (crystal/NMR), return sentinel 100.0 everywhere
    so downstream >=50 and >=70 thresholds always pass.
    """
    plddt = {}
    if structure is None:
        return plddt
    model = structure[0]
    if chain_id not in model:
        # Fall back to first chain if requested chain missing
        for c in model:
            chain_id = c.id
            break
    if chain_id not in model:
        return plddt

    for res in model[chain_id].get_residues():
        if res.id[0] != ' ':
            continue
        if plddt_gate:
            # AF: real pLDDT from B-factor
            p = None
            if 'CA' in res:
                p = res['CA'].bfactor
            else:
                for atom in res:
                    p = atom.bfactor
                    break
            if p is not None and p > 0:
                plddt[res.id[1]] = round(p, 2)
        else:
            # Crystal/NMR: uniform sentinel
            plddt[res.id[1]] = CRYSTAL_PLDDT_SENTINEL
    return plddt


def get_residue_aa(structure, chain_id="A"):
    """Map {residue_number: one-letter AA}."""
    if structure is None:
        return {}
    model = structure[0]
    if chain_id not in model:
        for c in model:
            chain_id = c.id
            break
    if chain_id not in model:
        return {}
    return {r.id[1]: THREE_TO_ONE.get(r.resname, '?')
            for r in model[chain_id].get_residues() if r.id[0] == ' '}


# ============================================================================
# Contact analysis
# ============================================================================
def count_intra_contacts(structure, chain_id="A", distance=5.0, positions=None):
    """
    Unique residue-residue contacts within a chain, sequence separation >= 3.

    If `positions` is given (list of residue numbers), only those positions
    are scored — O(k) atom lookups via NeighborSearch instead of O(n²) all-pairs.
    If None, computes for every residue (fallback, slow on large proteins).
    """
    if structure is None:
        return {}
    model = structure[0]
    if chain_id not in model:
        for c in model:
            chain_id = c.id
            break
    if chain_id not in model:
        return {}

    chain = model[chain_id]
    residues = [r for r in chain.get_residues() if r.id[0] == ' ']
    res_by_id = {r.id[1]: r for r in residues}
    res_index = {r.id[1]: i for i, r in enumerate(residues)}

    # Build NeighborSearch over all chain atoms
    all_atoms = [a for r in residues for a in r.get_atoms()]
    if not all_atoms:
        return {}
    ns = NeighborSearch(all_atoms)

    target_positions = positions if positions is not None else [r.id[1] for r in residues]
    contacts = {}

    for pos in target_positions:
        res = res_by_id.get(pos)
        if res is None:
            continue
        my_idx = res_index[pos]
        neighbor_set = set()
        for atom in res.get_atoms():
            for nb_res in ns.search(atom.coord, distance, 'R'):
                if nb_res.id[0] != ' ':
                    continue
                other_idx = res_index.get(nb_res.id[1])
                if other_idx is None:
                    continue
                if abs(my_idx - other_idx) < 3:
                    continue
                neighbor_set.add(nb_res.id[1])
        contacts[pos] = len(neighbor_set)
    return contacts


def count_interface(structure, my_chain, partner_chain, distance=5.0):
    """Inter-chain contacts: unique partner residues within distance per residue."""
    if structure is None:
        return {}, set()
    model = structure[0]
    if my_chain not in model or partner_chain not in model:
        return {}, set()

    partner_atoms = list(model[partner_chain].get_atoms())
    if not partner_atoms:
        return {}, set()
    ns = NeighborSearch(partner_atoms)

    inter = {}
    iface = set()
    for res in model[my_chain].get_residues():
        if res.id[0] != ' ':
            continue
        partner_residues = set()
        for atom in res.get_atoms():
            for nb in ns.search(atom.coord, distance, 'R'):
                if nb.id[0] == ' ':
                    partner_residues.add(nb.id[1])
        cnt = len(partner_residues)
        if cnt > 0:
            inter[res.id[1]] = cnt
            iface.add(res.id[1])
    return inter, iface


# ============================================================================
# SASA / burial
# ============================================================================
def compute_sasa(structure, chain_id="A"):
    """Per-residue SASA via Shrake-Rupley."""
    if structure is None:
        return {}
    try:
        sr = ShrakeRupley()
        sr.compute(structure, level="R")
    except Exception:
        return {}
    model = structure[0]
    if chain_id not in model:
        for c in model:
            chain_id = c.id
            break
    if chain_id not in model:
        return {}
    out = {}
    for res in model[chain_id].get_residues():
        if res.id[0] != ' ':
            continue
        if hasattr(res, 'sasa'):
            out[res.id[1]] = round(res.sasa, 2)
    return out


def classify_burial(sasa, aa_one_letter):
    """Classify burial as surface_exposed / partially_buried / buried_core."""
    if sasa is None or aa_one_letter not in MAX_SASA:
        return 'unknown'
    rsa = sasa / MAX_SASA[aa_one_letter]
    if rsa < 0.10:
        return 'buried_core'
    elif rsa < 0.25:
        return 'partially_buried'
    else:
        return 'surface_exposed'


# ============================================================================
# High-level per-system loaders
# ============================================================================
def load_monomer(gene: str, system_cfg: SystemConfig, structure_dir: Path):
    """
    Load monomer structure for one gene. Returns (structure, plddt_dict, monomer_spec).
    Uses CIF first if provided, otherwise PDB. Applies plddt_gate based on monomer type.
    """
    structure_dir = Path(structure_dir)
    mspec = system_cfg.monomers.get(gene.lower())
    if mspec is None:
        return None, {}, None

    # Try CIF first if specified
    if mspec.cif_file:
        path = structure_dir / mspec.cif_file
        struct = load_cif(path)
        if struct is not None:
            plddt = get_plddt(struct, mspec.chain_id, plddt_gate=mspec.plddt_gate)
            if plddt:
                return struct, plddt, mspec

    # PDB fallback (or primary)
    path = structure_dir / mspec.pdb_file
    struct = load_pdb(path)
    if struct is not None:
        plddt = get_plddt(struct, mspec.chain_id, plddt_gate=mspec.plddt_gate)
        return struct, plddt, mspec

    return None, {}, mspec


def load_multimer(system_cfg: SystemConfig, preprocessed_dir: Path, source_dir: Path):
    """
    Load the multimer structure for a system.
    Crystal/NMR systems are loaded from preprocessed_dir; AF from source_dir.
    Returns the Biopython structure object or None.
    """
    multi = system_cfg.multimer
    if multi.structure_type == "AF":
        # AF complex PDBs have zeroed B-factors; real pLDDT lives in the CIF.
        # Prefer CIF when provided (mirrors load_monomer), else fall back to PDB.
        if getattr(multi, "cif_file", None):
            cif_path = Path(source_dir) / multi.cif_file
            struct = load_cif(cif_path)
            if struct is not None:
                return struct
        path = Path(source_dir) / multi.pdb_file
    else:
        path = Path(preprocessed_dir) / multi.pdb_file

    return load_pdb(path)


def get_multimer_chain_plddt(structure, chain_id: str, plddt_gate: bool) -> Dict[int, float]:
    """Per-residue pLDDT for one chain in a multimer. Sentinel-filled if not gated."""
    return get_plddt(structure, chain_id, plddt_gate=plddt_gate)
