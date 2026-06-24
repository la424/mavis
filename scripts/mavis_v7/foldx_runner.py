"""
FoldX wrapper functions.

Three-axis DDG framework:
  Layer 2: BuildModel on monomer                 → ddg_monomer
  Layer 3: BuildModel on complex                 → ddg_fold_{partner}
  Layer 4: AnalyseComplex pairwise per partner   → ddg_binding_{partner}

v7 changes:
  - n_runs=5 with per-run value recording + SD
  - Pairwise AnalyseComplex — iterates over each partner chain individually,
    so multi-chain complexes (e.g. hemoglobin tetramer) produce one binding
    value per (variant, partner_chain) pair
  - RepairPDB cache shared across all calls (keyed by resolved path)
  - DRY_RUN mode for testing pipeline plumbing without FoldX installed
"""

import os
import shutil
import subprocess
from pathlib import Path
from statistics import mean, stdev
from typing import Dict, List, Optional, Tuple
import random

# Module-level cache — persists across calls within one Python session
_REPAIRED_CACHE: Dict[str, Path] = {}


# ============================================================================
# RepairPDB
# ============================================================================
def repair_pdb(
    structure_path: Path,
    work_dir: Path,
    foldx_binary: Path,
    rotabase: Optional[Path] = None,
    timeout: int = 600,
    verbose: bool = False,
    dry_run: bool = False,
) -> Path:
    """
    Run FoldX RepairPDB. Returns path to repaired structure (or original on failure).
    Cached per resolved input path.
    """
    if dry_run:
        return Path(structure_path)
    structure_path = Path(structure_path)
    cache_key = str(structure_path.resolve())
    if cache_key in _REPAIRED_CACHE:
        cached = _REPAIRED_CACHE[cache_key]
        if cached.exists():
            return cached

    work_dir = Path(work_dir).resolve()
    work_dir.mkdir(parents=True, exist_ok=True)
    struct_name = structure_path.name

    shutil.copy2(structure_path, work_dir / struct_name)
    if rotabase and rotabase.exists() and not (work_dir / "rotabase.txt").exists():
        shutil.copy2(rotabase, work_dir / "rotabase.txt")

    cmd = [
        str(foldx_binary),
        "--command=RepairPDB",
        f"--pdb={struct_name}",
        f"--output-dir={work_dir}",
    ]

    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True, timeout=timeout, cwd=str(work_dir)
        )
        if result.returncode != 0:
            print(f"    RepairPDB FAILED (rc={result.returncode}) for {struct_name}")
            print(f"    STDERR tail: {result.stderr[-400:]}")
            _REPAIRED_CACHE[cache_key] = structure_path
            return structure_path

        base = struct_name.replace(".pdb", "")
        repaired = work_dir / f"{base}_Repair.pdb"
        if not repaired.exists():
            for f in work_dir.glob("*_Repair.pdb"):
                repaired = f
                break

        if repaired.exists():
            _REPAIRED_CACHE[cache_key] = repaired
            return repaired
        else:
            if verbose:
                print(f"    RepairPDB: no _Repair.pdb found for {struct_name}, using original")
            _REPAIRED_CACHE[cache_key] = structure_path
            return structure_path

    except subprocess.TimeoutExpired:
        if verbose:
            print(f"    RepairPDB timeout for {struct_name}")
        _REPAIRED_CACHE[cache_key] = structure_path
        return structure_path
    except Exception as e:
        if verbose:
            print(f"    RepairPDB exception for {struct_name}: {e}")
        _REPAIRED_CACHE[cache_key] = structure_path
        return structure_path


# ============================================================================
# BuildModel (monomer DDG and fold-in-complex DDG)
# ============================================================================
def build_model(
    structure_path: Path,
    chain_id: str,
    ref_aa: str,
    position: int,
    alt_aa: str,
    work_dir: Path,
    foldx_binary: Path,
    rotabase: Optional[Path] = None,
    n_runs: int = 5,
    timeout: int = 600,
    dry_run: bool = False,
    verbose: bool = False,
) -> Tuple[Optional[float], Optional[float], List[float], Optional[Path]]:
    """
    Run FoldX BuildModel. Returns (ddg_mean, ddg_sd, all_runs, mutant_pdb_path).

    ddg_mean    — mean DDG across n_runs (kcal/mol), None on failure
    ddg_sd      — stdev across runs, 0.0 for n=1, None on failure
    all_runs    — raw per-run DDG values
    mutant_pdb  — path to first-run mutant PDB (used by multimer flow for AnalyseComplex)
    """
    if dry_run:
        # Produce reproducible mock values so tests are deterministic
        rng = random.Random(f"{structure_path.name}:{chain_id}{position}:{ref_aa}{alt_aa}")
        mean_val = round(rng.gauss(0.5, 1.2), 4)
        runs = [round(mean_val + rng.gauss(0, 0.3), 4) for _ in range(n_runs)]
        return round(mean(runs), 4), round(stdev(runs), 4), runs, None

    work_dir = Path(work_dir).resolve()
    work_dir.mkdir(parents=True, exist_ok=True)
    struct_name = structure_path.name

    if not (work_dir / struct_name).exists():
        shutil.copy2(structure_path, work_dir / struct_name)
    if rotabase and rotabase.exists() and not (work_dir / "rotabase.txt").exists():
        shutil.copy2(rotabase, work_dir / "rotabase.txt")

    mut_str = f"{ref_aa}{chain_id}{position}{alt_aa};"
    mut_file = work_dir / "individual_list.txt"
    mut_file.write_text(mut_str + "\n")

    cmd = [
        str(foldx_binary),
        "--command=BuildModel",
        f"--pdb={struct_name}",
        "--mutant-file=individual_list.txt",
        f"--numberOfRuns={n_runs}",
        f"--output-dir={work_dir}",
    ]

    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True, timeout=timeout, cwd=str(work_dir)
        )
        if result.returncode != 0:
            print(f"    FoldX BuildModel FAILED (rc={result.returncode}) in {work_dir}")
            print(f"    STDERR: {result.stderr[:500]}")
            print(f"    STDOUT tail: {result.stdout[-500:]}")
            return None, None, [], None

        # Locate Dif file
        struct_base = struct_name.replace(".pdb", "")
        dif_file = work_dir / f"Dif_{struct_base}.fxout"
        if not dif_file.exists():
            candidates = list(work_dir.glob(f"Dif_{struct_base}*.fxout")) or list(work_dir.glob("Dif_*.fxout"))
            if candidates:
                dif_file = candidates[0]

        if not dif_file.exists():
            if verbose:
                print(f"    No Dif output in {work_dir}")
            return None, None, [], None

        ddg_values = []
        with open(dif_file) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith(("Pdb", "#")):
                    continue
                parts = line.split("\t")
                if len(parts) >= 2:
                    try:
                        ddg_values.append(float(parts[1]))
                    except ValueError:
                        continue

        if not ddg_values:
            return None, None, [], None

        ddg_mean = round(mean(ddg_values), 4)
        ddg_sd = round(stdev(ddg_values), 4) if len(ddg_values) > 1 else 0.0

        # Mutant structures — FoldX 5.1 naming:
        #   n_runs=1 → {base}_1.pdb
        #   n_runs>1 → {base}_1_0.pdb, {base}_1_1.pdb, ...
        # Return ALL mutant replicates so AnalyseComplex can be averaged.
        mutant_pdbs = sorted([
            f for f in work_dir.glob(f"{struct_base}_1*.pdb")
            if f.name != struct_name
            and not f.name.startswith("WT_")
            and f.name != f"{struct_base}.pdb"
        ])

        return ddg_mean, ddg_sd, ddg_values, mutant_pdbs

    except subprocess.TimeoutExpired:
        if verbose:
            print(f"    FoldX timeout for {ref_aa}{position}{alt_aa}")
        return None, None, [], None
    except Exception as e:
        if verbose:
            print(f"    FoldX exception: {e}")
        return None, None, [], None


# ============================================================================
# AnalyseComplex (pairwise interaction energy)
# ============================================================================
def analyse_complex(
    structure_path: Path,
    chain_a: str,
    chain_b: str,
    work_dir: Path,
    foldx_binary: Path,
    rotabase: Optional[Path] = None,
    timeout: int = 300,
    dry_run: bool = False,
    verbose: bool = False,
) -> Optional[float]:
    """
    Run FoldX AnalyseComplex for one chain pair. Returns interaction energy (kcal/mol).
    """
    if dry_run:
        rng = random.Random(f"IE:{structure_path.name}:{chain_a}:{chain_b}")
        return round(rng.gauss(-15.0, 3.0), 4)

    work_dir = Path(work_dir).resolve()
    work_dir.mkdir(parents=True, exist_ok=True)
    struct_name = structure_path.name

    if not (work_dir / struct_name).exists():
        shutil.copy2(structure_path, work_dir / struct_name)
    if rotabase and rotabase.exists() and not (work_dir / "rotabase.txt").exists():
        shutil.copy2(rotabase, work_dir / "rotabase.txt")

    chains = f"{chain_a},{chain_b}"
    cmd = [
        str(foldx_binary),
        "--command=AnalyseComplex",
        f"--pdb={struct_name}",
        f"--analyseComplexChains={chains}",
        f"--output-dir={work_dir}",
    ]

    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True, timeout=timeout, cwd=str(work_dir)
        )
        if result.returncode != 0:
            print(f"    FoldX AnalyseComplex FAILED (rc={result.returncode}) in {work_dir}")
            print(f"    STDERR: {result.stderr[:400]}")
            print(f"    STDOUT tail: {result.stdout[-400:]}")
            return None

        # Deterministic file ordering
        interaction_files = sorted(work_dir.glob("Interaction_*_AC.fxout"))
        if not interaction_files:
            print(f"    AnalyseComplex wrote no Interaction file in {work_dir}")
            print(f"    STDOUT tail: {result.stdout[-400:]}")
            return None

        for f in interaction_files:
            with open(f) as fh:
                for line in fh:
                    if line.startswith(("Pdb", "#")) or not line.strip():
                        continue
                    parts = line.strip().split("\t")
                    if len(parts) >= 6:
                        try:
                            return float(parts[5])  # Interaction Energy
                        except ValueError:
                            continue

        print(f"    AnalyseComplex produced Interaction file but no parseable IE in {work_dir}")
        return None
    except subprocess.TimeoutExpired:
        print(f"    AnalyseComplex TIMEOUT in {work_dir}")
        return None
    except Exception as e:
        print(f"    AnalyseComplex EXCEPTION in {work_dir}: {e}")
        return None


# ============================================================================
# Three-axis DDG for one variant × one system
# ============================================================================
def compute_three_axis_ddg(
    variant: Dict,
    system_cfg,
    structure_dir: Path,
    preprocessed_dir: Path,
    foldx_work_dir: Path,
    foldx_binary: Path,
    rotabase: Optional[Path] = None,
    n_runs: int = 5,
    dry_run: bool = False,
    verbose: bool = True,
) -> Dict:
    """
    Compute all three DDG axes for one variant in one system.
    Returns a dict with per-axis values, SDs, per-run lists, and per-partner binding DDGs.

    variant is a dict with keys: gene, ref_aa, alt_aa, position, position_struct, variant
    """
    gene = variant["gene"].lower()
    ref = variant["ref_aa"]
    alt = variant["alt_aa"]
    pos_mono = int(variant["position_mono"])
    pos_multi = int(variant["position_multi"])

    out = {
        "ddg_monomer": None, "ddg_monomer_sd": None, "ddg_monomer_runs": [],
        "ddg_fold_by_partner": {},           # partner_chain -> (fold_mean, fold_sd, runs)
        "ddg_binding_by_partner": {},        # partner_chain -> binding mean
        "ddg_binding_sd_by_partner": {},     # partner_chain -> binding SD
        "ddg_binding_runs_by_partner": {},   # partner_chain -> [5 binding DDG replicates]
        "foldx_errors": [],
    }

    # ------- Axis 2: Monomer DDG (BuildModel on AF monomer) -------
    mspec = system_cfg.monomers.get(gene)
    if mspec is None:
        out["foldx_errors"].append(f"No monomer spec for {gene}")
        return out

    mono_path = Path(structure_dir) / mspec.pdb_file
    if not mono_path.exists():
        out["foldx_errors"].append(f"Monomer file missing: {mono_path}")
    else:
        repaired_mono = repair_pdb(
            mono_path, foldx_work_dir / "repaired",
            foldx_binary, rotabase, dry_run=dry_run,
        )
        wd = foldx_work_dir / "monomer" / f"{variant['system']}_{gene}_{ref}{pos_mono}{alt}"
        mean_val, sd_val, runs, _ = build_model(
            repaired_mono, mspec.chain_id, ref, pos_mono, alt, wd,
            foldx_binary, rotabase, n_runs=n_runs,
            dry_run=dry_run, verbose=verbose,
        )
        out["ddg_monomer"] = mean_val
        out["ddg_monomer_sd"] = sd_val
        out["ddg_monomer_runs"] = runs
        if mean_val is None:
            out["foldx_errors"].append(f"monomer BuildModel failed")

    # ------- Axis 3 + 4: Multimer DDG (fold + binding per partner) -------
    multi = system_cfg.multimer
    my_chain = multi.chain_map.get(gene)
    if my_chain is None:
        out["foldx_errors"].append(f"{gene} not in multimer chain_map for {variant['system']}")
        return out

    partner_chains = multi.pairwise_partners.get(gene, [])
    if not partner_chains:
        out["foldx_errors"].append(f"No partners defined for {gene} in {variant['system']}")
        return out

    # Load preprocessed or AF multimer
    if multi.structure_type == "AF":
        multimer_path = Path(structure_dir) / multi.pdb_file
    else:
        multimer_path = Path(preprocessed_dir) / multi.pdb_file

    if not multimer_path.exists():
        out["foldx_errors"].append(f"Multimer missing: {multimer_path}")
        return out

    # Repair multimer once, then one BuildModel (produces mutant complex PDB)
    repaired_multi = repair_pdb(
        multimer_path, foldx_work_dir / "repaired",
        foldx_binary, rotabase, dry_run=dry_run,
    )

    # BuildModel on the complex — fold DDG + mutant structure for AnalyseComplex
    wd_multi = (
        foldx_work_dir / "multimer" / f"{variant['system']}_{gene}_{ref}{pos_multi}{alt}"
    )
    fold_mean, fold_sd, fold_runs, mutant_pdbs = build_model(
        repaired_multi, my_chain, ref, pos_multi, alt, wd_multi,
        foldx_binary, rotabase, n_runs=n_runs,
        dry_run=dry_run, verbose=verbose,
    )

    if fold_mean is None:
        out["foldx_errors"].append("multimer BuildModel failed")
        return out

    # For DRY_RUN, we don't have real mutant PDBs — synthesize binding DDG directly
    if dry_run:
        for p_chain in partner_chains:
            out["ddg_fold_by_partner"][p_chain] = (fold_mean, fold_sd, fold_runs)
            # mock binding delta
            rng = random.Random(f"BIND:{variant['system']}:{gene}{pos_multi}:{p_chain}")
            out["ddg_binding_by_partner"][p_chain] = round(rng.gauss(0.3, 1.5), 4)
            out["ddg_binding_sd_by_partner"][p_chain] = 0.0
            out["ddg_binding_runs_by_partner"][p_chain] = [out["ddg_binding_by_partner"][p_chain]]
        return out

    if not mutant_pdbs:
        out["foldx_errors"].append("mutant PDBs not found after BuildModel")
        # still record fold result
        for p_chain in partner_chains:
            out["ddg_fold_by_partner"][p_chain] = (fold_mean, fold_sd, fold_runs)
        return out

    # Pairwise AnalyseComplex: WT once, mutant for EACH replicate, then average
    for p_chain in partner_chains:
        # WT (repaired multimer) — run once; WT doesn't vary across replicates
        wt_dir = wd_multi / f"wt_{my_chain}{p_chain}"
        ie_wt = analyse_complex(
            repaired_multi, my_chain, p_chain, wt_dir,
            foldx_binary, rotabase, dry_run=False, verbose=verbose,
        )

        # Mutant — AnalyseComplex on each of the 5 replicate mutant PDBs
        ie_mut_values = []
        for i, mpdb in enumerate(mutant_pdbs):
            if not mpdb.exists():
                continue
            mut_dir = wd_multi / f"mut_{my_chain}{p_chain}_{i}"
            ie_mut = analyse_complex(
                mpdb, my_chain, p_chain, mut_dir,
                foldx_binary, rotabase, dry_run=False, verbose=verbose,
            )
            if ie_mut is not None:
                ie_mut_values.append(ie_mut)

        if ie_wt is not None and ie_mut_values:
            ddg_runs = [round(ie_m - ie_wt, 4) for ie_m in ie_mut_values]
            ddg_binding_mean = round(mean(ddg_runs), 4)
            ddg_binding_sd = round(stdev(ddg_runs), 4) if len(ddg_runs) > 1 else 0.0
            out["ddg_binding_by_partner"][p_chain] = ddg_binding_mean
            out["ddg_binding_sd_by_partner"][p_chain] = ddg_binding_sd
            out["ddg_binding_runs_by_partner"][p_chain] = ddg_runs
        else:
            out["foldx_errors"].append(
                f"AnalyseComplex failed for {my_chain} vs {p_chain} "
                f"(ie_wt={ie_wt}, n_mut_ok={len(ie_mut_values)})"
            )

        # Fold axis is chain-independent; replicate per partner for consistency
        out["ddg_fold_by_partner"][p_chain] = (fold_mean, fold_sd, fold_runs)

    return out
