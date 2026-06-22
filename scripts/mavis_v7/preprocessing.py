"""
Structure preprocessing for FoldX input.

Applies HETATM stripping and NMR MODEL extraction to crystal/NMR structures
per the policy encoded in STRUCTURE_CONFIG. AF structures bypass preprocessing
entirely (returned as-is).

Called once at pipeline start. Writes _processed.pdb files to a dedicated
directory; downstream code reads from there.
"""

from pathlib import Path
from typing import Dict, Tuple
import json
from datetime import datetime

from .config import SystemConfig, MultimerSpec


# ============================================================================
# Core transformation
# ============================================================================
def _extract_model(lines, model_num):
    """Extract the specified MODEL's ATOM/HETATM records from NMR ensemble."""
    in_model = False
    out = []
    for line in lines:
        rec = line[:6].strip()
        if rec == "MODEL":
            try:
                m = int(line.split()[1])
            except (ValueError, IndexError):
                continue
            in_model = (m == model_num)
            continue
        elif rec == "ENDMDL":
            if in_model:
                in_model = False
            continue
        if rec in ("ATOM", "HETATM"):
            if in_model:
                out.append(line)
        else:
            # Header/metadata records: pass through once
            out.append(line)
    return out


def _filter_hetatm(lines, keep_all, keep_per_chain, strip_all):
    """Apply HETATM policy."""
    out = []
    for line in lines:
        rec = line[:6].strip()
        if rec != "HETATM":
            out.append(line)
            continue
        res_name = line[17:20].strip()
        chain = line[21]

        if res_name in strip_all:
            continue
        if res_name in keep_all:
            out.append(line)
            continue
        if res_name in keep_per_chain.get(chain, []):
            out.append(line)
            continue
        # Default: strip unknown HETATM
        continue
    return out


def _audit_lines(lines):
    """Count ATOM chains and HETATM residues by chain."""
    chains = {}
    hetatm = {}
    models = 0
    for line in lines:
        rec = line[:6].strip()
        if rec == "MODEL":
            models += 1
        elif rec == "ATOM":
            ch = line[21]
            chains[ch] = chains.get(ch, 0) + 1
        elif rec == "HETATM":
            rn = line[17:20].strip()
            ch = line[21]
            hetatm[f"{ch}:{rn}"] = hetatm.get(f"{ch}:{rn}", 0) + 1
    return chains, hetatm, models


# ============================================================================
# Per-system preprocessing
# ============================================================================
def preprocess_multimer(multi: MultimerSpec, source_dir: Path, output_dir: Path) -> Tuple[Path, Dict]:
    """
    Preprocess one crystal/NMR multimer. Returns (output_path, stats).
    For AF structures, returns (source_path, empty_stats) without modification.
    """
    source_dir = Path(source_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # AF: no preprocessing needed
    if multi.structure_type == "AF":
        src = source_dir / multi.pdb_file
        return src, {"skipped": True, "reason": "AF structure"}

    # Crystal / NMR
    if multi.source_file:
        src = source_dir / multi.source_file
    else:
        src = source_dir / multi.pdb_file

    if not src.exists():
        raise FileNotFoundError(f"Source structure missing: {src}")

    with open(src) as f:
        lines = f.readlines()

    stats = {
        "system": multi.system,
        "source": src.name,
        "structure_type": multi.structure_type,
        "input_lines": len(lines),
    }

    # Input HETATM audit
    _, hetatm_in, models_in = _audit_lines(lines)
    stats["hetatm_input"] = dict(hetatm_in)
    stats["models_input"] = models_in

    # NMR: extract MODEL
    if multi.nmr_model is not None:
        lines = _extract_model(lines, multi.nmr_model)
        stats["nmr_model_extracted"] = f"MODEL {multi.nmr_model} of {models_in}"

    # HETATM filtering
    lines = _filter_hetatm(
        lines,
        keep_all=multi.hetatm_keep_all,
        keep_per_chain=multi.hetatm_keep_per_chain,
        strip_all=multi.hetatm_strip_all,
    )

    # Output audit
    chains_out, hetatm_out, _ = _audit_lines(lines)
    stats["atom_chains"] = dict(chains_out)
    stats["hetatm_output"] = dict(hetatm_out)
    stats["output_lines"] = len(lines)

    out_path = output_dir / multi.pdb_file
    with open(out_path, "w") as f:
        f.writelines(lines)
    stats["output_file"] = out_path.name
    stats["notes"] = multi.preprocessing_notes

    return out_path, stats


def preprocess_all(configs: Dict[str, SystemConfig], source_dir: Path,
                   output_dir: Path, verbose: bool = True) -> Dict:
    """
    Preprocess all multimers referenced in configs. Returns a provenance dict
    that should be written to disk for the audit trail.
    """
    source_dir = Path(source_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    all_stats = []
    if verbose:
        print("=" * 70)
        print(f"MAVIS v7 — preprocessing {len(configs)} systems")
        print("=" * 70)

    for sys_name, cfg in configs.items():
        try:
            out_path, stats = preprocess_multimer(cfg.multimer, source_dir, output_dir)
            all_stats.append({"system": sys_name, "ok": True, "stats": stats})

            if verbose:
                if stats.get("skipped"):
                    print(f"  {sys_name:<22} AF — no preprocessing")
                else:
                    n_het_in = sum(stats["hetatm_input"].values())
                    n_het_out = sum(stats["hetatm_output"].values())
                    model_info = f"[{stats.get('nmr_model_extracted')}] " if "nmr_model_extracted" in stats else ""
                    print(f"  {sys_name:<22} {stats['structure_type']} {model_info}"
                          f"HETATM: {n_het_in} → {n_het_out}   → {stats['output_file']}")
        except Exception as e:
            all_stats.append({"system": sys_name, "ok": False, "error": str(e)})
            if verbose:
                print(f"  {sys_name:<22} ERROR: {e}")

    provenance = {
        "generated": datetime.now().isoformat(),
        "mavis_version": "v7",
        "source_dir": str(source_dir),
        "output_dir": str(output_dir),
        "systems_processed": sum(1 for s in all_stats if s["ok"]),
        "systems_total": len(all_stats),
        "details": all_stats,
    }

    log_path = output_dir / "preprocessing_provenance.json"
    with open(log_path, "w") as f:
        json.dump(provenance, f, indent=2, default=str)

    if verbose:
        print(f"\n  Provenance: {log_path}")
        ok = provenance["systems_processed"]
        total = provenance["systems_total"]
        print(f"  Done: {ok}/{total} systems processed successfully")

    return provenance
