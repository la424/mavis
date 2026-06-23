"""
MAVIS v7 — Multimer-Aware Variant Impact Scoring
Generalized, multimer-aware structural variant-interpretation pipeline:
evaluates variant disruptiveness in the appropriate complex and resolves the
specific mechanism of disruption (it is a structural-disruption / mechanism
pipeline, not a pathogenicity classifier).

Accepts arbitrary monomers, multimers, and variants via STRUCTURE_CONFIG.
Three-axis DDG framework (monomer fold, fold-in-complex, binding) with
pairwise AnalyseComplex for multi-chain complexes.
"""

__version__ = "7.0.0"
