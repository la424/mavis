"""
Variant loading with system threading and deduplication.

Accepts the benchmark v5 CSV schema (system, gene, ref_aa, alt_aa, position, ...)
or the legacy combined-variant-string format. Dedup key is
(system, gene, position, ref_aa, alt_aa) — allows the same variant across
different systems (e.g., same position in different complexes).
"""

from pathlib import Path
from typing import Dict, List, Optional
import pandas as pd

from .config import SystemConfig


def load_variants(
    csv_path: Path,
    configs: Dict[str, SystemConfig],
    verbose: bool = True,
) -> pd.DataFrame:
    """
    Load the variant CSV and validate against configs.
    Applies position_offsets if defined in the system's multimer config.
    """
    csv_path = Path(csv_path)
    # utf-8-sig to handle BOM
    df = pd.read_csv(csv_path, encoding="utf-8-sig")

    # Normalize column names
    df.columns = [c.strip() for c in df.columns]
    required = ["system", "gene", "ref_aa", "alt_aa", "position"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(
            f"Variant CSV missing required columns: {missing}\n"
            f"Available columns: {list(df.columns)}"
        )

    # Normalize values
    df["system"] = df["system"].astype(str).str.strip()
    df["gene"] = df["gene"].astype(str).str.strip().str.lower()
    df["ref_aa"] = df["ref_aa"].astype(str).str.strip().str.upper()
    df["alt_aa"] = df["alt_aa"].astype(str).str.strip().str.upper()
    df["position"] = df["position"].astype(int)

    # Variant string for display
    df["variant"] = df["ref_aa"] + df["position"].astype(str) + df["alt_aa"]

    n_before = len(df)

    # Dedup on (system, gene, position, ref, alt)
    dedup_key = ["system", "gene", "position", "ref_aa", "alt_aa"]
    dupes = df.duplicated(dedup_key, keep=False)
    if dupes.any():
        if verbose:
            print(f"⚠ {dupes.sum()} duplicate variant rows found — keeping first:")
            for _, group in df[dupes].groupby(dedup_key):
                r = group.iloc[0]
                print(f"    {r['system']} / {r['gene']} {r['variant']} × {len(group)}")
        df = df.drop_duplicates(dedup_key, keep="first").reset_index(drop=True)

    # Validate systems exist in configs
    unknown_systems = set(df["system"]) - set(configs.keys())
    if unknown_systems:
        raise ValueError(
            f"Unknown systems in variant CSV (not in STRUCTURE_CONFIG): {unknown_systems}"
        )

    # Validate genes exist in per-system monomer specs
    unknown = []
    for _, row in df.iterrows():
        sys_cfg = configs[row["system"]]
        if row["gene"] not in sys_cfg.monomers:
            unknown.append(f"{row['system']}/{row['gene']}")
    if unknown:
        raise ValueError(
            f"Variants reference genes not in system monomer specs:\n  " +
            "\n  ".join(sorted(set(unknown)))
        )

    # Apply position offsets is now deferred to the consumers (monomer vs multimer),
    # since monomer AF structures and multimer crystal structures can have different
    # numbering conventions. Each consumer uses:
    #   position_mono  = position - monomer_spec.position_offset
    #   position_multi = position - multimer.position_offsets.get(gene, 0)
    # The raw CSV `position` column is preserved unchanged.
    df["position_mono_offset"] = 0
    df["position_multi_offset"] = 0
    for idx, row in df.iterrows():
        sys_cfg = configs[row["system"]]
        mspec = sys_cfg.monomers.get(row["gene"])
        if mspec:
            df.at[idx, "position_mono_offset"] = mspec.position_offset
        df.at[idx, "position_multi_offset"] = sys_cfg.multimer.position_offsets.get(row["gene"], 0)

    df["position_mono"] = df["position"] - df["position_mono_offset"]
    df["position_multi"] = df["position"] - df["position_multi_offset"]

    if verbose:
        print(f"✓ Loaded {len(df)} variants from {csv_path.name}")
        if n_before != len(df):
            print(f"  (deduplicated from {n_before})")
        by_sys = df.groupby("system").size().sort_values(ascending=False)
        for sys, n in by_sys.items():
            print(f"    {sys:<22} {n} variants")

    return df
