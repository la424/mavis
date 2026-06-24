"""
mavis_helpers.py — pure-Python logic behind the MAVIS guided Colab notebook.

Everything here is importable and unit-testable OUTSIDE Colab (no GPU, no FoldX,
no google.colab). Colab-only steps (pip installs, GPU checks, files.upload /
files.download, ColabFold execution, AlphaFold Server) live in notebook cells,
not here. Keeping the heavy logic in this module is what lets us validate the
notebook locally.
"""
from __future__ import annotations
import os, re, glob, json, subprocess, urllib.request, urllib.parse
from pathlib import Path

AA1 = set("ACDEFGHIKLMNPQRSTVWY")

# --------------------------------------------------------------------------
# 1. Variant parsing
# --------------------------------------------------------------------------
def parse_variants(text):
    """Parse free-text variant lines into a DataFrame [gene, ref_aa, position, alt_aa].

    Accepts one variant per line, e.g.:
        SHROOM3 G1003R
        tnni3, R162W
        ZIC3  C253S
    Blank lines and lines starting with '#' are ignored.
    Returns (df, errors) where errors is a list of (lineno, raw, reason).
    """
    import pandas as pd
    rows, errors = [], []
    for i, raw in enumerate(text.splitlines(), 1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        parts = re.split(r"[,\s]+", line)
        if len(parts) < 2:
            errors.append((i, raw, "need at least GENE and a variant (e.g. 'GENE R162W')"))
            continue
        gene, tok = parts[0], parts[1]
        m = re.fullmatch(r"([A-Za-z])(\d+)([A-Za-z])", tok)
        if m:
            ref, pos, alt = m.group(1).upper(), int(m.group(2)), m.group(3).upper()
        elif len(parts) >= 4 and parts[2].isdigit():
            ref, pos, alt = parts[1].upper(), int(parts[2]), parts[3].upper()
        else:
            errors.append((i, raw, f"could not parse variant token '{tok}' (want like R162W)"))
            continue
        if ref not in AA1 or alt not in AA1:
            errors.append((i, raw, f"non-standard amino acid in '{tok}'"))
            continue
        rows.append({"gene": gene.lower(), "ref_aa": ref, "position": pos, "alt_aa": alt})
    df = pd.DataFrame(rows, columns=["gene", "ref_aa", "position", "alt_aa"])
    return df, errors

def variant_label(ref_aa, position, alt_aa):
    return f"{ref_aa}{position}{alt_aa}"

# --------------------------------------------------------------------------
# 2. UniProt resolution + sequence
# --------------------------------------------------------------------------
def resolve_uniprot(gene, organism_id=9606, reviewed=True, timeout=30):
    """Resolve a gene symbol to a UniProt accession (human reviewed by default).

    Picks the candidate whose PRIMARY gene name actually matches `gene` -- UniProt's
    relevance ordering is not reliable for exact gene lookups (e.g. a gene_exact:TNNC1
    query can return TNNI3 first). Returns (accession, candidates), or (None, []).
    """
    gene = (gene or "").strip()
    if not gene:
        raise ValueError("empty gene symbol")
    if not str(organism_id).strip():
        organism_id = 9606
    q = f"gene_exact:{gene} AND organism_id:{organism_id}"
    if reviewed:
        q += " AND reviewed:true"
    url = ("https://rest.uniprot.org/uniprotkb/search?query="
           + urllib.parse.quote(q)
           + "&fields=accession,gene_primary,protein_name,length&format=tsv&size=5")
    req = urllib.request.Request(
        url, headers={"User-Agent": "MAVIS-colab/1.0 (+https://github.com/la424/mavis)"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        text = r.read().decode()
    lines = [ln for ln in text.strip().splitlines() if ln]
    if len(lines) < 2:
        return None, []
    header = lines[0].split("\t")
    cands = [dict(zip(header, ln.split("\t"))) for ln in lines[1:]]
    gcol = "Gene Names (primary)"
    exact = [c for c in cands if c.get(gcol, "").strip().lower() == gene.lower()]
    chosen = (exact or cands)[0]
    acc = chosen.get("Entry") or list(chosen.values())[0]
    return acc, cands

def fetch_uniprot_sequence(accession, timeout=30):
    """Return the canonical protein sequence (string) for a UniProt accession."""
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.fasta"
    with urllib.request.urlopen(url, timeout=timeout) as r:
        fasta = r.read().decode()
    return "".join(ln.strip() for ln in fasta.splitlines() if ln and not ln.startswith(">"))

# --------------------------------------------------------------------------
# 3. AlphaFold DB monomer fetch
# --------------------------------------------------------------------------
def alphafold_urls(accession, timeout=30):
    """Return (pdbUrl, cifUrl) for the latest AlphaFold DB model via the AF-DB API."""
    api = f"https://alphafold.ebi.ac.uk/api/prediction/{accession}"
    data = json.load(urllib.request.urlopen(api, timeout=timeout))
    if not data:
        raise RuntimeError(f"no AlphaFold DB entry for {accession}")
    e = data[0]
    return e.get("pdbUrl"), e.get("cifUrl")

def fetch_alphafold_monomer(accession, out_dir, fmt="pdb", timeout=120):
    """Download the latest AlphaFold DB monomer model (pdb or cif) via the AF-DB API.

    Uses the API so it tracks the current model version (v4/v5/v6/...) instead of
    a hardcoded URL. Returns local Path.
    """
    out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)
    pdb_url, cif_url = alphafold_urls(accession, timeout=timeout)
    url = cif_url if fmt == "cif" else pdb_url
    if not url:
        raise RuntimeError(f"AlphaFold DB has no {fmt} model for {accession}")
    dest = out_dir / os.path.basename(url)
    urllib.request.urlretrieve(url, dest)
    if dest.stat().st_size == 0:
        raise RuntimeError(f"empty download from {url}")
    return dest

# --------------------------------------------------------------------------
# 4. Residue counting + size routing + CIF->PDB
# --------------------------------------------------------------------------
def count_residues(path, chain=None):
    """Count protein residues via unique (chain, resseq) on CA atoms. PDB or mmCIF."""
    path = Path(path)
    if path.suffix.lower() == ".cif":
        return _count_residues_cif(path, chain)
    seen = set()
    with open(path) as fh:
        for ln in fh:
            if (ln.startswith("ATOM") or ln.startswith("HETATM")) and ln[12:16].strip() == "CA":
                ch = ln[21]
                if chain and ch != chain:
                    continue
                seen.add((ch, ln[22:27].strip()))
    return len(seen)

def _count_residues_cif(path, chain=None):
    """Minimal mmCIF CA counter (no external deps)."""
    seen, cols, idx = set(), [], {}
    with open(path) as fh:
        for ln in fh:
            s = ln.strip()
            if s.startswith("_atom_site."):
                cols.append(s); continue
            if cols and (s.startswith("ATOM") or s.startswith("HETATM")):
                if not idx:
                    idx = {name: k for k, name in enumerate(cols)}
                f = s.split()
                try:
                    if f[idx["_atom_site.label_atom_id"]].strip('"') != "CA":
                        continue
                    ch = f[idx.get("_atom_site.auth_asym_id", idx.get("_atom_site.label_asym_id"))]
                    res = f[idx.get("_atom_site.auth_seq_id", idx.get("_atom_site.label_seq_id"))]
                except Exception:
                    continue
                if chain and ch != chain:
                    continue
                seen.add((ch, res))
            elif cols and idx and s and not s.startswith(("ATOM", "HETATM", "_", "#")):
                break
    return len(seen)

def route_by_size(total_residues, colabfold_limit=1400):
    """'colabfold' for complexes up to ~colabfold_limit residues, else 'alphafold_server'."""
    return "colabfold" if total_residues <= colabfold_limit else "alphafold_server"

def cif_to_pdb(cif_path, pdb_path=None):
    """Convert mmCIF -> PDB. Prefers gemmi (handles huge structures / hybrid-36)."""
    cif_path = Path(cif_path)
    pdb_path = Path(pdb_path) if pdb_path else cif_path.with_suffix(".pdb")
    try:
        import gemmi
        st = gemmi.read_structure(str(cif_path))
        st.setup_entities()
        Path(pdb_path).write_text(st.make_pdb_string())
        return pdb_path
    except ImportError:
        from Bio.PDB import MMCIFParser, PDBIO
        s = MMCIFParser(QUIET=True).get_structure("s", str(cif_path))
        io = PDBIO(); io.set_structure(s); io.save(str(pdb_path))
        return pdb_path

# --------------------------------------------------------------------------
# 5. ColabFold / AlphaFold Server input prep
# --------------------------------------------------------------------------
def slice_sequence(seq, start, end):
    """1-based inclusive slice for domain scoping."""
    return seq[start - 1:end]

def colabfold_query(sequences):
    """Join chain sequences for a ColabFold multimer query: 'SEQ_A:SEQ_B'."""
    return ":".join(sequences)

def write_fasta(name, sequences, out_dir):
    """Write a ColabFold-style multimer FASTA (':'-joined chains). Returns Path."""
    out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)
    p = out_dir / f"{name}.fasta"
    p.write_text(f">{name}\n{colabfold_query(sequences)}\n")
    return p

# --------------------------------------------------------------------------
# 6. Config generation (+ validation via the repo's own loader)
# --------------------------------------------------------------------------
def build_config_yaml(systems, out_path):
    """Write a MAVIS YAML config.

    systems: list of dicts:
      {"name","structure_type"(AF|xray|nmr),"complex_file",
       "genes":[{"gene","chain","monomer_file","monomer_offset"?,"multimer_offset"?,"monomer_cif"?}, ...]}
    Emits the schema load_config expects: system-level 'pdb_file' (required) +
    per-gene 'monomer_pdb' (required), 'chain', optional offsets.
    """
    import yaml
    doc = {"systems": {}}
    for s in systems:
        genes_block = {}
        for g in s["genes"]:
            gb = {"chain": g["chain"], "monomer_pdb": os.path.basename(g["monomer_file"])}
            if g.get("monomer_cif"):
                gb["monomer_cif"] = os.path.basename(g["monomer_cif"])
            if int(g.get("monomer_offset", 0)):
                gb["monomer_offset"] = int(g["monomer_offset"])
            if int(g.get("multimer_offset", 0)):
                gb["multimer_offset"] = int(g["multimer_offset"])
            genes_block[g["gene"].lower()] = gb
        entry = {"structure_type": s.get("structure_type", "AF"),
                 "pdb_file": os.path.basename(s["complex_file"]),
                 "genes": genes_block}
        if s.get("cif_file"):
            entry["cif_file"] = os.path.basename(s["cif_file"])
        doc["systems"][s["name"]] = entry
    Path(out_path).write_text(yaml.safe_dump(doc, sort_keys=False))
    return Path(out_path)

def validate_config(yaml_path, scripts_dir):
    """Load the YAML through the repo's config_io.load_config (raises on malformed)."""
    import sys
    scripts_dir = str(scripts_dir)
    if scripts_dir not in sys.path:
        sys.path.insert(0, scripts_dir)
    from mavis_v7.config_io import load_config
    return load_config(str(yaml_path))

# --------------------------------------------------------------------------
# 7. Results: per-variant mechanism summary + cards
# --------------------------------------------------------------------------
def _first(g, col):
    if col in g and g[col].notna().any():
        return g[col].dropna().iloc[0]
    return float("nan")

def _max(g, col):
    if col in g and g[col].notna().any():
        return float(g[col].dropna().max())
    return float("nan")

def _any_true(g, col):
    return bool(g[col].fillna(False).any()) if col in g else False

def _tier_num(v):
    import pandas as pd
    if pd.isna(v):
        return float("nan")
    m = re.search(r"(\d+)", str(v))
    return int(m.group(1)) if m else float("nan")

def _min_tier(g):
    import numpy as np
    if "mavis_tier" not in g:
        return float("nan")
    nums = [t for t in (_tier_num(v) for v in g["mavis_tier"])
            if not (isinstance(t, float) and np.isnan(t))]
    return min(nums) if nums else float("nan")

_SEVERITY = [(r"interface", 4), (r"binding", 4), (r"fold", 3),
             (r"destab", 3), (r"monomer", 2), (r"no .*effect|silent|none", 0)]
def _mech_severity(s):
    s = str(s).lower()
    for pat, sev in _SEVERITY:
        if re.search(pat, s):
            return sev
    return 1

def partners_present(df):
    return [c[len("ddg_binding_"):] for c in df.columns
            if c.startswith("ddg_binding_") and not c.endswith("_sd")]

def _mech_for_partner(g, partner):
    best, best_sev = None, -1
    for _, row in g.iterrows():
        mech = row.get("mavis_mechanism")
        if mech is None or isinstance(mech, float):
            continue
        rp = row.get("mavis_mechanism_partner")
        if partner is None or rp == partner or isinstance(rp, float):
            sev = _mech_severity(mech)
            if sev > best_sev:
                best, best_sev = mech, sev
    return best

def overall_mechanism(g):
    best, best_sev = "No structural effect detected", -1
    for v in g.get("mavis_mechanism", []):
        if v is None or isinstance(v, float):
            continue
        sev = _mech_severity(v)
        if sev > best_sev:
            best, best_sev = v, sev
    return best

def per_partner_table(df):
    """Long tidy summary: one row per (gene, variant, partner) with the 3 axes + pLDDT.

    This is the human-readable downloadable summary (vs. the full wide results CSV).
    """
    import pandas as pd, numpy as np
    partners = partners_present(df)
    rows = []
    for (gene, variant), g in df.groupby(["gene", "variant"], sort=False):
        mono = dict(ddg_monomer=_first(g, "ddg_monomer"),
                    ddg_monomer_sd=_first(g, "ddg_monomer_sd"),
                    monomer_plddt=_max(g, "monomer_plddt"),
                    monomer_confident=_any_true(g, "ddg_monomer_confident"))
        tier = _min_tier(g)
        had = False
        for p in partners:
            fb, bb = f"ddg_fold_{p}", f"ddg_binding_{p}"
            if (fb in g and g[fb].notna().any()) or (bb in g and g[bb].notna().any()):
                had = True
                rows.append(dict(gene=gene, variant=variant, partner=p, tier=tier, **mono,
                                 ddg_fold=_first(g, fb), ddg_fold_sd=_first(g, f"{fb}_sd"),
                                 ddg_binding=_first(g, bb), ddg_binding_sd=_first(g, f"{bb}_sd"),
                                 interface_plddt=_max(g, f"multi_{p}_plddt"),
                                 mechanism=_mech_for_partner(g, p)))
        if not had:
            rows.append(dict(gene=gene, variant=variant, partner=None, tier=tier, **mono,
                             ddg_fold=np.nan, ddg_fold_sd=np.nan, ddg_binding=np.nan,
                             ddg_binding_sd=np.nan, interface_plddt=np.nan,
                             mechanism=overall_mechanism(g)))
    cols = ["gene", "variant", "partner", "tier", "mechanism",
            "ddg_monomer", "ddg_monomer_sd", "monomer_plddt", "monomer_confident",
            "ddg_fold", "ddg_fold_sd", "ddg_binding", "ddg_binding_sd", "interface_plddt"]
    return pd.DataFrame(rows)[cols]

def _fmt(v, sd=None, nd=2):
    import pandas as pd
    if v is None or (isinstance(v, float) and pd.isna(v)):
        return "\u2014"
    s = f"{v:+.{nd}f}" if isinstance(v, (int, float)) else str(v)
    if sd is not None and not (isinstance(sd, float) and pd.isna(sd)):
        s += f" \u00b1 {sd:.2f}"
    return s

def _plddt_note(p):
    import pandas as pd
    if p is None or (isinstance(p, float) and pd.isna(p)):
        return ""
    if p >= 70:
        return f"pLDDT {p:.0f} (high)"
    if p >= 50:
        return f"pLDDT {p:.0f} (moderate)"
    return f"pLDDT {p:.0f} (low \u2014 gated out)"

def render_cards_html(df):
    """Build readable per-variant HTML mechanism cards (for display in the notebook)."""
    import pandas as pd
    tbl = per_partner_table(df)
    out = ['<div style="font-family:-apple-system,Segoe UI,Roboto,sans-serif">']
    for (gene, variant), g in tbl.groupby(["gene", "variant"], sort=False):
        tier = g["tier"].iloc[0]
        tier_txt = f"Tier {int(tier)}" if pd.notna(tier) else "Tier \u2014"
        overall = max((m for m in g["mechanism"] if isinstance(m, str)),
                      key=_mech_severity, default="No structural effect detected")
        mono = _fmt(g["ddg_monomer"].iloc[0], g["ddg_monomer_sd"].iloc[0])
        mono_p = _plddt_note(g["monomer_plddt"].iloc[0])
        out.append(
            '<div style="border:1px solid #d0d7de;border-radius:10px;padding:12px 14px;margin:10px 0">'
            f'<div style="font-size:16px;font-weight:600">{gene.upper()} {variant} '
            f'<span style="color:#57606a;font-weight:500">\u2014 {tier_txt} \u00b7 {overall}</span></div>'
            f'<div style="margin-top:6px;font-size:13px"><b>Monomer fold</b>: '
            f'\u0394\u0394G {mono} kcal/mol &nbsp;<span style="color:#57606a">{mono_p}</span></div>')
        partnered = g[g["partner"].notna()]
        if len(partnered):
            out.append('<table style="border-collapse:collapse;margin-top:8px;font-size:13px;width:100%">'
                       '<tr style="text-align:left;color:#57606a">'
                       '<th style="padding:3px 8px">Partner</th>'
                       '<th style="padding:3px 8px">Fold-in-complex \u0394\u0394G</th>'
                       '<th style="padding:3px 8px">Binding \u0394\u0394G</th>'
                       '<th style="padding:3px 8px">Interface</th>'
                       '<th style="padding:3px 8px">Call</th></tr>')
            for _, r in partnered.iterrows():
                out.append('<tr style="border-top:1px solid #eaeef2">'
                           f'<td style="padding:3px 8px"><b>{str(r["partner"]).upper()}</b></td>'
                           f'<td style="padding:3px 8px">{_fmt(r["ddg_fold"], r["ddg_fold_sd"])}</td>'
                           f'<td style="padding:3px 8px">{_fmt(r["ddg_binding"], r["ddg_binding_sd"])}</td>'
                           f'<td style="padding:3px 8px">{_plddt_note(r["interface_plddt"])}</td>'
                           f'<td style="padding:3px 8px">{r["mechanism"] or "\u2014"}</td></tr>')
            out.append('</table>')
        out.append('</div>')
    out.append('</div>')
    return "\n".join(out)

# --------------------------------------------------------------------------
# 8. Optional four-way concordance
# --------------------------------------------------------------------------
def concordance_template(variants_df):
    """Editable template for manually-curated external annotations.

    Columns are exactly what apply_concordance_v5.py --external requires:
    gene, variant, AM pathogenicity, AM class, franklin.
    """
    base = variants_df.copy()
    base["variant"] = [variant_label(r.ref_aa, r.position, r.alt_aa)
                       for r in base.itertuples(index=False)]
    out = base[["gene", "variant"]].drop_duplicates().reset_index(drop=True)
    out["AM pathogenicity"] = ""   # 0-1 AlphaMissense score
    out["AM class"] = ""           # likely_benign | ambiguous | likely_pathogenic
    out["franklin"] = ""           # Franklin / ClinVar classification text
    return out

def run_concordance(structural_csv, annotations_df, outdir, scripts_dir, python_exe="python"):
    """Write annotations to xlsx and invoke apply_concordance_v5.py. Returns (master_csv, proc)."""
    outdir = Path(outdir); outdir.mkdir(parents=True, exist_ok=True)
    ext_xlsx = outdir / "external_annotations.xlsx"
    annotations_df.to_excel(ext_xlsx, index=False)
    script = str(Path(scripts_dir) / "apply_concordance_v5.py")
    cmd = [python_exe, script, "--results", str(structural_csv),
           "--external", str(ext_xlsx), "--outdir", str(outdir)]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    master = outdir / "mavis_v7_concordance.csv"
    return (master if master.exists() else None), proc
