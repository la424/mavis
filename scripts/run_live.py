"""Benchmark pipeline entry point (live FoldX run).

Runs the multimer-aware structural scoring over the 44-variant / 11-PPI-system
benchmark. Produces results/mavis_v7_results.csv, which the concordance and
evaluation steps consume (see README). To reproduce the headline metrics WITHOUT
FoldX, use the cached self-test in README instead.

Requires:
  - FOLDX_BINARY env var pointing to your FoldX 5.x binary
  - AlphaFold structures in ./structures (see README "Inputs")
Run from the repository root.
"""
import os
import sys
import warnings
warnings.filterwarnings('ignore')
from pathlib import Path

# Make the bundled mavis_v7 package importable.
sys.path.insert(0, str(Path(__file__).resolve().parent))

from mavis_v7.config import build_benchmark_config
from mavis_v7.pipeline import run_pipeline

FOLDX = Path(os.environ.get("FOLDX_BINARY", "foldx"))

results = run_pipeline(
    configs=build_benchmark_config(),
    variants_csv=Path('benchmark_variants_v5.csv'),
    structure_dir=Path('structures'),
    preprocessed_dir=Path('processed'),
    output_dir=Path('results'),
    foldx_binary=FOLDX,
    rotabase=None,
    n_runs=5,
    dry_run=False,
    verbose=True,
)
df = results['df']
df.to_csv('results/mavis_v7_results.csv', index=False)
print('Wrote', len(df), 'variants to results/mavis_v7_results.csv')
