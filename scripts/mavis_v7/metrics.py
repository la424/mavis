"""
Per-variant structural metrics: contacts, interface status, burial, pLDDT.

Consumes the variant dataframe (with system/gene/position) and the STRUCTURE_CONFIG,
produces monomer_* and multi_{partner}_* columns used by scoring and the
mechanism classifier.

Handles mixed AF/crystal pLDDT via the sentinel in structure_loading.py.
"""

from pathlib import Path
from typing import Dict, List
import pandas as pd

from .config import SystemConfig
from .constants import (
    THREE_TO_ONE, get_grantham, grantham_severity, get_property_changes,
)
from .structure_loading import (
    load_pdb, load_cif, load_monomer, load_multimer,
    get_plddt, get_residue_aa, count_intra_contacts, count_interface,
    compute_sasa, classify_burial, CRYSTAL_PLDDT_SENTINEL,
)


def _partner_label_for_chain(multi_cfg, chain_id: str) -> str:
    """Given a multimer's chain, return its gene/label (inverse of chain_map)."""
    for gene, ch in multi_cfg.chain_map.items():
        if ch == chain_id:
            return gene
    return f"chain_{chain_id}"


def collect_partner_labels(configs: Dict[str, SystemConfig]) -> List[str]:
    """
    Collect the full set of partner labels across all systems, for column naming.
    Each system's partners become potential multi_{label}_* columns.
    """
    labels = set()
    for sys_cfg in configs.values():
        multi = sys_cfg.multimer
        for gene, partner_chains in multi.pairwise_partners.items():
            for ch in partner_chains:
                pl = _partner_label_for_chain(multi, ch)
                labels.add(pl)
    return sorted(labels)


def compute_monomer_metrics(df: pd.DataFrame, configs: Dict[str, SystemConfig],
                            structure_dir: Path, verbose: bool = True) -> pd.DataFrame:
    """
    Add monomer metrics per variant:
      monomer_plddt, monomer_aa_match, monomer_n_contacts, monomer_burial
    """
    # Gather monomer-space positions of interest per (system, gene)
    positions_needed = {}
    for _, row in df.iterrows():
        key = (row["system"], row["gene"])
        positions_needed.setdefault(key, set()).add(int(row["position_mono"]))

    # Cache: (system, gene) → dict of metrics
    cache = {}

    out_cols = {
        "monomer_plddt": [],
        "monomer_aa_at_pos": [],
        "monomer_aa_match": [],
        "monomer_n_contacts": [],
        "monomer_burial": [],
        "monomer_structure_type": [],
    }

    for _, row in df.iterrows():
        sys_name = row["system"]
        gene = row["gene"]
        pos = int(row["position_mono"])
        ref = row["ref_aa"]

        key = (sys_name, gene)
        if key not in cache:
            cfg = configs[sys_name]
            struct, plddt, mspec = load_monomer(gene, cfg, structure_dir)
            pos_list = sorted(positions_needed.get(key, set()))
            contacts = count_intra_contacts(struct, mspec.chain_id, positions=pos_list) if struct else {}
            aa_map = get_residue_aa(struct, mspec.chain_id) if struct else {}
            # SASA is whole-structure; computed once per (system, gene). Filter to needed positions only.
            sasa_full = compute_sasa(struct, mspec.chain_id) if struct else {}
            sasa = {p: sasa_full[p] for p in pos_list if p in sasa_full}
            cache[key] = {
                "plddt": plddt, "contacts": contacts, "aa_map": aa_map,
                "sasa": sasa, "mspec": mspec,
            }

        c = cache[key]
        mspec = c["mspec"]

        pl_val = c["plddt"].get(pos)
        aa_at = c["aa_map"].get(pos, "")
        aa_match = (aa_at == ref) if aa_at else None
        nc = c["contacts"].get(pos)
        sasa_val = c["sasa"].get(pos)
        burial = classify_burial(sasa_val, ref) if sasa_val is not None else "unknown"

        out_cols["monomer_plddt"].append(pl_val)
        out_cols["monomer_aa_at_pos"].append(aa_at)
        out_cols["monomer_aa_match"].append(aa_match)
        out_cols["monomer_n_contacts"].append(nc)
        out_cols["monomer_burial"].append(burial)
        out_cols["monomer_structure_type"].append(mspec.structure_type if mspec else "")

    for k, v in out_cols.items():
        df[k] = v

    # Substitution severity (Grantham-based)
    df["substitution_severity"] = df.apply(
        lambda r: grantham_severity(get_grantham(r["ref_aa"], r["alt_aa"])), axis=1
    )
    df["grantham_distance"] = df.apply(
        lambda r: get_grantham(r["ref_aa"], r["alt_aa"]), axis=1
    )
    df["property_changes"] = df.apply(
        lambda r: get_property_changes(r["ref_aa"], r["alt_aa"]), axis=1
    )

    if verbose:
        n_with_plddt = df["monomer_plddt"].notna().sum()
        n_mismatch = (df["monomer_aa_match"] == False).sum()
        print(f"  Monomer metrics: {n_with_plddt}/{len(df)} have pLDDT, "
              f"{n_mismatch} AA mismatches")

    return df


def compute_multimer_metrics(df: pd.DataFrame, configs: Dict[str, SystemConfig],
                             preprocessed_dir: Path, source_dir: Path,
                             verbose: bool = True) -> pd.DataFrame:
    """
    Add per-partner multimer metrics:
      multi_{partner}_plddt, multi_{partner}_inter_contacts,
      multi_{partner}_is_interface, multi_{partner}_burial, multi_{partner}_aa_at_pos

    where {partner} is derived from the chain_map (gene/label name, not chain letter).
    """
    all_partners = collect_partner_labels(configs)

    # Initialize columns
    for pl in all_partners:
        df[f"multi_{pl}_plddt"] = pd.NA
        df[f"multi_{pl}_inter_contacts"] = pd.NA
        df[f"multi_{pl}_is_interface"] = False
        df[f"multi_{pl}_burial"] = ""
        df[f"multi_{pl}_aa_at_pos"] = ""

    df["multimer_plddt_max"] = pd.NA
    df["multimer_plddt_min"] = pd.NA

    # Cache per system
    mstruct_cache = {}  # system_name -> loaded multimer structure + derived metrics

    for idx, row in df.iterrows():
        sys_name = row["system"]
        gene = row["gene"]
        pos = int(row["position_multi"])
        ref = row["ref_aa"]
        cfg = configs[sys_name]
        multi = cfg.multimer

        my_chain = multi.chain_map.get(gene)
        if my_chain is None:
            continue

        # Load the multimer structure once per system
        if sys_name not in mstruct_cache:
            struct = load_multimer(cfg, preprocessed_dir, source_dir)
            mstruct_cache[sys_name] = {"struct": struct, "per_partner": {}}

        struct = mstruct_cache[sys_name]["struct"]
        if struct is None:
            continue

        # pLDDT of the variant chain in the multimer
        my_plddt_map = get_plddt(struct, my_chain, plddt_gate=multi.plddt_gate)
        my_aa_map = get_residue_aa(struct, my_chain)

        partner_plddts_at_pos = []

        for p_chain in multi.pairwise_partners.get(gene, []):
            pl = _partner_label_for_chain(multi, p_chain)
            cache_key = (gene, p_chain)
            if cache_key not in mstruct_cache[sys_name]["per_partner"]:
                inter, iface = count_interface(struct, my_chain, p_chain)
                mstruct_cache[sys_name]["per_partner"][cache_key] = (inter, iface)
            inter, iface = mstruct_cache[sys_name]["per_partner"][cache_key]

            # Per-partner pLDDT at the variant site — use partner's chain for gating
            p_plddt_map = get_plddt(struct, p_chain, plddt_gate=multi.plddt_gate)
            # The effective partner pLDDT for this variant site is the partner chain's
            # average nearby residues if at interface, otherwise the variant residue's
            # pLDDT in its own chain. We use variant-chain pLDDT as the gate value here,
            # consistent with v6 semantics.
            variant_site_plddt = my_plddt_map.get(pos)

            df.at[idx, f"multi_{pl}_plddt"] = variant_site_plddt
            df.at[idx, f"multi_{pl}_inter_contacts"] = inter.get(pos, 0)
            df.at[idx, f"multi_{pl}_is_interface"] = pos in iface
            df.at[idx, f"multi_{pl}_aa_at_pos"] = my_aa_map.get(pos, "")

            if variant_site_plddt is not None:
                partner_plddts_at_pos.append(variant_site_plddt)

        if partner_plddts_at_pos:
            df.at[idx, "multimer_plddt_max"] = max(partner_plddts_at_pos)
            df.at[idx, "multimer_plddt_min"] = min(partner_plddts_at_pos)

    # best_plddt = max(monomer_plddt, multimer_plddt_max)
    def best_plddt(r):
        vals = []
        for c in ["monomer_plddt", "multimer_plddt_max"]:
            v = r.get(c)
            if pd.notna(v):
                vals.append(float(v))
        return max(vals) if vals else None
    df["best_plddt"] = df.apply(best_plddt, axis=1)

    # Confidence flag per axis
    df["ddg_monomer_confident"] = df["monomer_plddt"].apply(
        lambda v: pd.notna(v) and float(v) >= 50
    )

    if verbose:
        n_any_iface = sum(
            df.apply(
                lambda r: any(r.get(f"multi_{pl}_is_interface") for pl in all_partners),
                axis=1,
            )
        )
        print(f"  Multimer metrics: {n_any_iface}/{len(df)} at ≥1 interface")

    return df


def populate_ddg_confidence(df: pd.DataFrame, partner_labels: List[str]) -> pd.DataFrame:
    """After FoldX results are merged, mark per-partner confidence flags."""
    for pl in partner_labels:
        plddt_col = f"multi_{pl}_plddt"
        conf_col = f"ddg_{pl}_confident"
        if plddt_col in df.columns:
            df[conf_col] = df[plddt_col].apply(
                lambda v: pd.notna(v) and float(v) >= 50
            )
        else:
            df[conf_col] = False
    return df
