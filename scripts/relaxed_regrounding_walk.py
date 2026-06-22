#!/usr/bin/env python3
"""
MAVIS v7 — Relaxed-tier re-grounding walk over the 18 contracted axes.

For each axis that went graded->unknown across Batches 5-11, classify as either:
  - RELAXED-DIRECTIONAL: indirect-but-specific evidence supports a sign
      (coupled/downstream readout | sibling analogy | off-axis quantitative measurement).
      -> writes expected_ddg_{axis}_relaxed = {destab|stab|neutral}, with provenance + evidence class.
  - STAYS-UNKNOWN: only structural-position inference, or no directional evidence at all.
      -> relaxed column stays 'unknown'.

EVIDENCE CLASSES (the only admissible relaxed sources; structural-position inference EXCLUDED):
  COUPLED   - a downstream/coupled readout implies the axis direction
  ANALOGY   - a directly-measured sibling substitution at the same residue implies direction
  OFFAXIS   - a direct quantitative measurement on a different axis implies this axis leans

Each entry cites the batch and the specific observation. NO direction is invented here;
every relaxed call traces to a logged primary-literature finding from the eleven batches.
"""

# (gene, variant, axis) -> dict(decision, token, evclass, batch, evidence)
# decision in {'relaxed','unknown'}; token only meaningful if 'relaxed'
WALK = {

 # ---- MSH2-MSH6 (Batch 7) ----
 # A636P: ATPase/Walker-A region; Ollila/Lutzen show MMR-dead + reduced expression.
 # mono/fold went unknown for lack of a fold-stability assay. Is there a directional lean?
 ('msh2','A636P','mono'): dict(decision='relaxed', token='destab', evclass='COUPLED', batch='B7',
    evidence='Reduced cellular expression / instability phenotype (Ollila 2008) is a downstream '
             'readout implying reduced fold stability; not a fold-stability assay, but directionally destab.'),
 ('msh2','A636P','fold'): dict(decision='unknown', token=None, evclass=None, batch='B7',
    evidence='Same instability readout already used for the monomer lean; the complex-fold axis has '
             'no independent directional readout (would double-count the same observation). Stays unknown.'),
 # C697F: directly tested, binding preserved (-> bind neutral, already recovered). mono/fold:
 ('msh2','C697F','mono'): dict(decision='unknown', token=None, evclass=None, batch='B7',
    evidence='Lutzen reduced-expression is mild and the variant is MMR-functional in some assays; '
             'no specific directional fold-stability readout. Structural-position only -> excluded. Stays unknown.'),
 ('msh2','C697F','fold'): dict(decision='unknown', token=None, evclass=None, batch='B7',
    evidence='No directional complex-fold readout. Stays unknown.'),
 # N127S, G322D: seed neutral on mono/fold; these are MMR-functional / low-impact. Neutral *lean* defensible?
 ('msh2','N127S','mono'): dict(decision='relaxed', token='neutral', evclass='COUPLED', batch='B7',
    evidence='MMR-proficient / near-WT functional behavior (variant treated as a functional control) is a '
             'downstream readout consistent with a preserved (neutral) fold. Directional-neutral lean.'),
 ('msh2','N127S','fold'): dict(decision='relaxed', token='neutral', evclass='COUPLED', batch='B7',
    evidence='Same MMR-proficiency; preserved assembly implies neutral complex-fold. Directional-neutral.'),
 ('msh2','G322D','mono'): dict(decision='relaxed', token='neutral', evclass='COUPLED', batch='B7',
    evidence='G322D is a common polymorphism / functionally near-neutral; downstream functional '
             'normality implies neutral fold. Directional-neutral lean.'),
 ('msh2','G322D','fold'): dict(decision='relaxed', token='neutral', evclass='COUPLED', batch='B7',
    evidence='Same; preserved function implies neutral complex-fold. Directional-neutral.'),

 # ---- VHL-ElonginC (Batch 6) ----
 # W117R: full contraction. mono/fold/bind all unknown. Evidence?
 ('vhl','W117R','mono'): dict(decision='unknown', token=None, evclass=None, batch='B6',
    evidence='"W117 is a beta-core hydrophobic residue" is structural-position inference (FoldX-circular) '
             '-> EXCLUDED by the relaxed rule. No measured fold-stability readout. Stays unknown.'),
 ('vhl','W117R','fold'): dict(decision='unknown', token=None, evclass=None, batch='B6',
    evidence='Same structural-position basis; excluded. No directional readout. Stays unknown.'),
 ('vhl','W117R','bind'): dict(decision='unknown', token=None, evclass=None, batch='B6',
    evidence='W117 not in Kishida ElonginC panel; no measurement either way. The lesion is on the HIF '
             'axis (Ohh 2000), not ElonginC. No directional ElonginC-binding readout. Stays unknown.'),
 # Y98H: directly tested binding-preserved (bind neutral, primary). mono went neutral->unknown.
 ('vhl','Y98H','mono'): dict(decision='relaxed', token='neutral', evclass='COUPLED', batch='B6',
    evidence='Y98H retains VBC assembly (Kishida) and is "no apparent structural role" by retained '
             'complex formation; preserved assembly is a downstream readout implying neutral fold. '
             'Directional-neutral lean (assembly proxy, not a fold ddG).'),

 # ---- TNNI3 (Batch 9) ----
 ('tnni3','R162W','bind'): dict(decision='unknown', token=None, evclass=None, batch='B9',
    evidence='The "stab" token was cross-variant inference transplanted from R145G; explicitly '
             'unsupported. Zhou 2013 (rat R163W) NOT in hand. No directional R162W cTnI-cTnC readout. '
             'Stays unknown (VERIFICATION-PENDING). A future direct readout could promote it.'),

 # ---- CaM-Cav1.2 (Batch 10) ----
 # mono went mild_destab->unknown for all three. The Ca2+-affinity numbers are OFF-AXIS quantitative
 # BUT they speak to a FUNCTIONAL (Ca2+-sensing) axis, NOT the fold-stability axis. Crotti NMR shows
 # the global fold is INTACT. Per the tightened OFFAXIS rule (off-axis evidence promotes ONLY if it
 # implies THIS axis's direction), the Ca2+ data does not speak to fold -> all three STAY UNKNOWN.
 # (DECIDED - Luke, this session.)
 ('calm1','D96V','mono'): dict(decision='unknown', token=None, evclass=None, batch='B10',
    evidence='C-lobe Ca2+ affinity loss 13.6x (Crotti Fig 5C) is OFF-AXIS functional (sensing), and '
             'Crotti NMR shows the fold is INTACT. Off-axis evidence does not imply the fold-axis '
             'direction -> stays unknown in relaxed too.'),
 ('calm1','N98S','mono'): dict(decision='unknown', token=None, evclass=None, batch='B10',
    evidence='C-domain Ca2+ affinity loss (Sondergaard N97S) is off-axis functional; fold intact. '
             'Does not speak to fold axis -> stays unknown.'),
 ('calm1','F142L','mono'): dict(decision='unknown', token=None, evclass=None, batch='B10',
    evidence='C-lobe Ca2+ affinity loss 5.4x (Crotti) is off-axis functional; fold intact. '
             'Does not speak to fold axis -> stays unknown.'),

 # ---- SMAD4-SMAD3 (Batch 11) ----
 # fold went destab->unknown for D351H, R361C. Shi measures oligomerization (binding axis), not fold.
 # Is there a directional FOLD lean distinct from the binding axis?
 ('smad4','D351H','fold'): dict(decision='unknown', token=None, evclass=None, batch='B11',
    evidence='Shi co-IP/gel-filtration is an oligomerization (binding-axis) readout, already grounding '
             'the binding token; it is not a fold-stability readout. No independent directional '
             'fold-stability evidence. The massive FoldX fold ddG is structural-position/self-prediction '
             '-> excluded. Stays unknown.'),
 ('smad4','R361C','fold'): dict(decision='unknown', token=None, evclass=None, batch='B11',
    evidence='Same as D351H: Shi measures oligomerization, not fold stability. No independent fold '
             'directional readout. Stays unknown.'),
}

def summarize():
    from collections import Counter
    relaxed=[(k,v) for k,v in WALK.items() if v['decision']=='relaxed']
    stays=[(k,v) for k,v in WALK.items() if v['decision']=='unknown']
    print(f'=== Relaxed re-grounding walk: {len(WALK)} contracted axes examined ===')
    print(f'  RELAXED-DIRECTIONAL (promoted): {len(relaxed)}')
    print(f'  STAYS-UNKNOWN:                  {len(stays)}')
    print()
    print('  Promoted, by evidence class:')
    ec=Counter(v['evclass'] for k,v in relaxed)
    for c,n in ec.items(): print(f'    {c}: {n}')
    print()
    print('  RELAXED-DIRECTIONAL detail:')
    for (k,v) in relaxed:
        print('    {:7s} {:7s} {:5s} -> {:8s} [{}/{}]'.format(k[0],k[1],k[2],v['token'],v['evclass'],v['batch']))
    print()
    print('  STAYS-UNKNOWN detail:')
    for (k,v) in stays:
        reason = 'structural-position (excluded)' if 'structural-position' in v['evidence'] or 'self-prediction' in v['evidence'] else \
                 ('no directional readout' if 'no ' in v['evidence'].lower() else 'other')
        print(f'    {k[0]:7s} {k[1]:7s} {k[2]:5s}    ({reason})')

if __name__=='__main__':
    summarize()
