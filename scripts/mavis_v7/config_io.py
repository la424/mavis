"""
config_io.py - declarative (YAML) systems config for MAVIS.

Lets users define their OWN genes/complexes in a YAML file instead of editing
Python. `load_config()` parses YAML into the SAME SystemConfig / MonomerSpec /
MultimerSpec dataclasses the engine already uses, so nothing downstream changes.

Also provides:
  - configs_to_yaml(): serialize existing Python configs to YAML (used to generate
    the shipped example configs, and for round-trip testing).
  - scaffold_system_yaml(): emit a ready-to-fill YAML block + the structure
    filenames a new (hub, partners) system needs.
  - autofanout(): expand a simple gene/ref_aa/position/alt_aa variant table into
    the one-row-per-system form the pipeline expects, using the config.

YAML schema (one block per complex under `systems:`):

systems:
  <system_name>:
    structure_type: AF            # AF | xray | nmr   (the COMPLEX)
    plddt_gate: true              # default: true for AF, false for xray/nmr
    pdb_file: <complex>.pdb       # required; looked up under --structures
    cif_file: <complex>.cif       # optional; CIF used for pLDDT (AF B-factors are 0)
    # xray/nmr preprocessing (optional):
    source_file: <raw>.pdb
    nmr_model: 1
    hetatm_keep_all: [ZN]
    hetatm_keep_per_chain: {B: [ZN]}
    hetatm_strip_all: [HOH]
    preprocessing_notes: "..."
    notes: "..."
    genes:                        # one entry per gene/chain in the complex
      <gene>:
        chain: A                  # this gene's chain id in the COMPLEX
        monomer_pdb: fold_<gene>_model_0.pdb   # required
        monomer_cif: fold_<gene>_model_0.cif   # optional (CIF pLDDT)
        monomer_offset: 0         # subtracted from CSV position for the MONOMER axis
        multimer_offset: 0        # subtracted from CSV position for the COMPLEX axes
    # Advanced (multi-chain / instance-labelled complexes, e.g. a homo-tetramer):
    # give chain_map / pairwise_partners explicitly to override what is derived from
    # genes.*.chain. Labels may be gene names or instance labels (e.g. hbb, hbb_2).
    # chain_map: {hba1_1: A, hbb: B, hba1_2: C, hbb_2: D}
    # pairwise_partners: {hbb: [A, C, D], ...}
"""
from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Tuple

from .config import SystemConfig, MonomerSpec, MultimerSpec


def _require_yaml():
    try:
        import yaml
        return yaml
    except ImportError as e:  # pragma: no cover
        raise ImportError(
            "PyYAML is required for YAML configs. Install it with:  pip install pyyaml"
        ) from e


# ---------------------------------------------------------------------------
# Load: YAML file -> {system_name: SystemConfig}
# ---------------------------------------------------------------------------
def load_config(path) -> Dict[str, SystemConfig]:
    yaml = _require_yaml()
    path = Path(path)
    with open(path) as fh:
        data = yaml.safe_load(fh)
    if not data or "systems" not in data or not data["systems"]:
        raise ValueError(f"{path}: top-level 'systems:' mapping not found or empty")
    return {name: _build_system(name, spec) for name, spec in data["systems"].items()}


def _build_system(name: str, s: dict) -> SystemConfig:
    if "pdb_file" not in s:
        raise ValueError(f"[{name}] missing required 'pdb_file'")
    if not s.get("genes"):
        raise ValueError(f"[{name}] missing required 'genes:' block")

    structure_type = s.get("structure_type", "AF")
    plddt_gate = s.get("plddt_gate", structure_type == "AF")
    genes = s["genes"]

    # monomers (one MonomerSpec per gene). NOTE: a gene's 'chain' is its chain in the
    # COMPLEX (-> chain_map); the monomer file's own chain is 'monomer_chain' (AF = A).
    monomers: Dict[str, MonomerSpec] = {}
    for gene, g in genes.items():
        g = g or {}
        if "monomer_pdb" not in g:
            raise ValueError(f"[{name}] gene '{gene}' missing 'monomer_pdb'")
        monomers[gene] = MonomerSpec(
            gene=gene,
            pdb_file=g["monomer_pdb"],
            cif_file=g.get("monomer_cif"),
            structure_type=g.get("monomer_structure_type", "AF"),
            plddt_gate=g.get("monomer_plddt_gate", True),
            chain_id=g.get("monomer_chain", "A"),
            position_offset=int(g.get("monomer_offset", 0)),
        )

    # chain_map: explicit override, else derived from genes.*.chain
    if s.get("chain_map"):
        chain_map = dict(s["chain_map"])
    else:
        chain_map = {}
        for gene, g in genes.items():
            g = g or {}
            if "chain" not in g:
                raise ValueError(
                    f"[{name}] gene '{gene}' needs 'chain:' "
                    f"(or give this system an explicit 'chain_map:')"
                )
            chain_map[gene] = g["chain"]

    # pairwise_partners: explicit override, else each label vs every other chain
    if s.get("pairwise_partners"):
        pairwise = {k: list(v) for k, v in s["pairwise_partners"].items()}
    else:
        all_chains = list(chain_map.values())
        pairwise = {lbl: [c for c in all_chains if c != ch] for lbl, ch in chain_map.items()}

    # multimer position offsets (nonzero only)
    position_offsets: Dict[str, int] = {}
    for gene, g in genes.items():
        g = g or {}
        off = int(g.get("multimer_offset", 0))
        if off:
            position_offsets[gene] = off

    multimer = MultimerSpec(
        system=name,
        pdb_file=s["pdb_file"],
        structure_type=structure_type,
        chain_map=chain_map,
        pairwise_partners=pairwise,
        plddt_gate=plddt_gate,
        position_offsets=position_offsets,
        cif_file=s.get("cif_file"),
        source_file=s.get("source_file"),
        nmr_model=s.get("nmr_model"),
        hetatm_keep_all=list(s.get("hetatm_keep_all") or []),
        hetatm_keep_per_chain={k: list(v) for k, v in (s.get("hetatm_keep_per_chain") or {}).items()},
        hetatm_strip_all=list(s.get("hetatm_strip_all") or []),
        preprocessing_notes=s.get("preprocessing_notes", ""),
    )
    return SystemConfig(system=name, monomers=monomers, multimer=multimer, notes=s.get("notes", ""))


# ---------------------------------------------------------------------------
# Dump: {system_name: SystemConfig} -> YAML string  (lossless round-trip)
# ---------------------------------------------------------------------------
def configs_to_yaml(configs: Dict[str, SystemConfig]) -> str:
    yaml = _require_yaml()
    out = {"systems": {}}
    for name, cfg in configs.items():
        m = cfg.multimer
        genes_block = {}
        for gene, mspec in cfg.monomers.items():
            gb = {"chain": m.chain_map.get(gene, mspec.chain_id), "monomer_pdb": mspec.pdb_file}
            if mspec.cif_file:
                gb["monomer_cif"] = mspec.cif_file
            if mspec.position_offset:
                gb["monomer_offset"] = mspec.position_offset
            if m.position_offsets.get(gene, 0):
                gb["multimer_offset"] = m.position_offsets[gene]
            if mspec.chain_id != "A":
                gb["monomer_chain"] = mspec.chain_id
            genes_block[gene] = gb

        entry = {"structure_type": m.structure_type, "pdb_file": m.pdb_file}
        if m.plddt_gate != (m.structure_type == "AF"):
            entry["plddt_gate"] = m.plddt_gate
        if m.cif_file:
            entry["cif_file"] = m.cif_file
        if m.source_file:
            entry["source_file"] = m.source_file
        if m.nmr_model is not None:
            entry["nmr_model"] = m.nmr_model
        if m.hetatm_keep_all:
            entry["hetatm_keep_all"] = list(m.hetatm_keep_all)
        if m.hetatm_keep_per_chain:
            entry["hetatm_keep_per_chain"] = {k: list(v) for k, v in m.hetatm_keep_per_chain.items()}
        if m.hetatm_strip_all:
            entry["hetatm_strip_all"] = list(m.hetatm_strip_all)
        if m.preprocessing_notes:
            entry["preprocessing_notes"] = m.preprocessing_notes
        if cfg.notes:
            entry["notes"] = cfg.notes
        entry["genes"] = genes_block

        # chain_map: emit only if it has labels beyond the genes, or differs from genes.*.chain
        derived_cm = {g: gb.get("chain") for g, gb in genes_block.items()}
        if set(m.chain_map) != set(genes_block) or m.chain_map != derived_cm:
            entry["chain_map"] = dict(m.chain_map)
        # pairwise: emit only if it differs from the default derivation
        cm = entry.get("chain_map", derived_cm)
        chains = list(cm.values())
        derived_pw = {lbl: [c for c in chains if c != ch] for lbl, ch in cm.items()}
        if m.pairwise_partners != derived_pw:
            entry["pairwise_partners"] = {k: list(v) for k, v in m.pairwise_partners.items()}

        out["systems"][name] = entry
    return yaml.safe_dump(out, sort_keys=False, default_flow_style=False)


# ---------------------------------------------------------------------------
# Scaffold: print a YAML block + needed structure files for a new system
# ---------------------------------------------------------------------------
def scaffold_system_yaml(hub: str, partners: List[str]) -> Tuple[str, List[str]]:
    hub = hub.lower()
    lines = ["systems:"]
    files = {f"fold_{hub}_model_0.pdb"}
    for p in partners:
        p = p.lower()
        sysname = f"{hub}_{p}"
        complex_pdb = f"fold_{hub}_{p}_model_0.pdb"
        files.update({complex_pdb, f"fold_{p}_model_0.pdb"})
        lines += [
            f"  {sysname}:",
            "    structure_type: AF",
            f"    pdb_file: {complex_pdb}",
            "    genes:",
            f"      {hub}:",
            "        chain: A",
            f"        monomer_pdb: fold_{hub}_model_0.pdb",
            f"      {p}:",
            "        chain: B",
            f"        monomer_pdb: fold_{p}_model_0.pdb",
        ]
    return "\n".join(lines) + "\n", sorted(files)


# ---------------------------------------------------------------------------
# Auto fan-out: simple variant table -> one row per (variant x system it appears in)
# ---------------------------------------------------------------------------
def autofanout(df, configs: Dict[str, SystemConfig]):
    import pandas as pd
    gene_to_systems: Dict[str, List[str]] = {}
    for sysname, cfg in configs.items():
        for gene in cfg.monomers:
            gene_to_systems.setdefault(gene.lower(), []).append(sysname)
    rows, unmatched = [], set()
    for _, row in df.iterrows():
        g = str(row["gene"]).strip().lower()
        systems = gene_to_systems.get(g)
        if not systems:
            unmatched.add(g)
            continue
        for s in systems:
            r = row.to_dict()
            r["system"] = s
            rows.append(r)
    return pd.DataFrame(rows), sorted(unmatched)
