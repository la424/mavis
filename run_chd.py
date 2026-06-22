"""CHD pipeline entry point (full concordance).

Runs the multimer-aware structural scoring over the CHD variant set, producing
the per-variant structural results. External-tool concordance (AlphaMissense,
Franklin/ClinVar) is layered on afterwards via scripts/apply_concordance_v5.py.
For the structural-only variant (no external tools), see README "CHD pipelines".

Requires:
  - FOLDX_BINARY env var pointing to your FoldX 5.x binary
  - AlphaFold structures in ./structures (see README "Inputs")
"""
import os
import sys
from pathlib import Path

# Make the bundled mavis_v7 package importable when run from the repo root.
sys.path.insert(0, str(Path(__file__).resolve().parent / "scripts"))

from mavis_v7.build_chd_config import build_chd_config
from mavis_v7.pipeline import run_pipeline

res = run_pipeline(
    configs=build_chd_config(),
    variants_csv=Path('chd_input_final.csv'),
    structure_dir=Path('structures'),
    preprocessed_dir=Path('processed'),
    output_dir=Path('results/chd_rerun'),
    foldx_binary=Path(os.environ.get("FOLDX_BINARY", "foldx")),
    rotabase=None, n_runs=5, dry_run=False, verbose=True,
)
df = res['df']
df.to_csv('results/chd_rerun/chd_structural_results.csv', index=False)
print('WROTE', len(df), 'rows -> results/chd_rerun/chd_structural_results.csv')
