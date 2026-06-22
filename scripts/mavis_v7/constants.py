"""
Amino acid data, Grantham distances, SASA, properties.
Unchanged from v6. Separated for importability.
"""

THREE_TO_ONE = {
    'ALA':'A','CYS':'C','ASP':'D','GLU':'E','PHE':'F','GLY':'G','HIS':'H',
    'ILE':'I','LYS':'K','LEU':'L','MET':'M','ASN':'N','PRO':'P','GLN':'Q',
    'ARG':'R','SER':'S','THR':'T','VAL':'V','TRP':'W','TYR':'Y',
}
ONE_TO_THREE = {v: k for k, v in THREE_TO_ONE.items()}

AA_PROPERTIES = {
    'A':{'size':'small','charge':'neutral','hydrophobic':True},
    'R':{'size':'large','charge':'positive','hydrophobic':False},
    'N':{'size':'medium','charge':'neutral','hydrophobic':False},
    'D':{'size':'medium','charge':'negative','hydrophobic':False},
    'C':{'size':'small','charge':'neutral','hydrophobic':True},
    'E':{'size':'medium','charge':'negative','hydrophobic':False},
    'Q':{'size':'medium','charge':'neutral','hydrophobic':False},
    'G':{'size':'small','charge':'neutral','hydrophobic':False},
    'H':{'size':'medium','charge':'positive','hydrophobic':False},
    'I':{'size':'medium','charge':'neutral','hydrophobic':True},
    'L':{'size':'medium','charge':'neutral','hydrophobic':True},
    'K':{'size':'large','charge':'positive','hydrophobic':False},
    'M':{'size':'medium','charge':'neutral','hydrophobic':True},
    'F':{'size':'large','charge':'neutral','hydrophobic':True},
    'P':{'size':'small','charge':'neutral','hydrophobic':False},
    'S':{'size':'small','charge':'neutral','hydrophobic':False},
    'T':{'size':'small','charge':'neutral','hydrophobic':False},
    'W':{'size':'large','charge':'neutral','hydrophobic':True},
    'Y':{'size':'large','charge':'neutral','hydrophobic':False},
    'V':{'size':'small','charge':'neutral','hydrophobic':True},
}

# Max SASA — Tien et al 2013, theoretical Gly-X-Gly
MAX_SASA = {
    'A':129,'R':274,'N':195,'D':193,'C':167,'E':223,'Q':225,'G':104,
    'H':224,'I':197,'L':201,'K':236,'M':224,'F':240,'P':159,'S':155,
    'T':172,'V':174,'W':285,'Y':263,
}

# Grantham distance matrix
GRANTHAM = {
    ('A','R'):112,('A','N'):111,('A','D'):126,('A','C'):195,('A','Q'):91,('A','E'):107,
    ('A','G'):60,('A','H'):86,('A','I'):94,('A','L'):96,('A','K'):106,('A','M'):84,
    ('A','F'):113,('A','P'):27,('A','S'):99,('A','T'):58,('A','W'):148,('A','Y'):112,('A','V'):64,
    ('R','N'):86,('R','D'):96,('R','C'):180,('R','Q'):43,('R','E'):54,('R','G'):125,
    ('R','H'):29,('R','I'):97,('R','L'):102,('R','K'):26,('R','M'):91,('R','F'):97,
    ('R','P'):103,('R','S'):110,('R','T'):71,('R','W'):101,('R','Y'):77,('R','V'):96,
    ('N','D'):23,('N','C'):139,('N','Q'):46,('N','E'):42,('N','G'):80,('N','H'):68,
    ('N','I'):149,('N','L'):153,('N','K'):94,('N','M'):142,('N','F'):158,('N','P'):91,
    ('N','S'):46,('N','T'):65,('N','W'):174,('N','Y'):143,('N','V'):133,
    ('D','C'):154,('D','Q'):61,('D','E'):45,('D','G'):94,('D','H'):81,('D','I'):168,
    ('D','L'):172,('D','K'):101,('D','M'):160,('D','F'):177,('D','P'):108,('D','S'):65,
    ('D','T'):85,('D','W'):181,('D','Y'):160,('D','V'):152,
    ('C','Q'):154,('C','E'):170,('C','G'):159,('C','H'):174,('C','I'):198,('C','L'):198,
    ('C','K'):202,('C','M'):196,('C','F'):205,('C','P'):169,('C','S'):112,('C','T'):149,
    ('C','W'):215,('C','Y'):194,('C','V'):192,
    ('Q','E'):29,('Q','G'):87,('Q','H'):24,('Q','I'):109,('Q','L'):113,('Q','K'):53,
    ('Q','M'):101,('Q','F'):116,('Q','P'):76,('Q','S'):68,('Q','T'):42,('Q','W'):130,
    ('Q','Y'):99,('Q','V'):96,
    ('E','G'):98,('E','H'):40,('E','I'):134,('E','L'):138,('E','K'):56,('E','M'):126,
    ('E','F'):140,('E','P'):93,('E','S'):80,('E','T'):65,('E','W'):152,('E','Y'):122,('E','V'):121,
    ('G','H'):98,('G','I'):135,('G','L'):138,('G','K'):127,('G','M'):127,('G','F'):153,
    ('G','P'):42,('G','S'):56,('G','T'):59,('G','W'):184,('G','Y'):147,('G','V'):109,
    ('H','I'):94,('H','L'):99,('H','K'):32,('H','M'):87,('H','F'):100,('H','P'):77,
    ('H','S'):89,('H','T'):47,('H','W'):115,('H','Y'):83,('H','V'):84,
    ('I','L'):5,('I','K'):102,('I','M'):10,('I','F'):21,('I','P'):95,('I','S'):142,
    ('I','T'):89,('I','W'):61,('I','Y'):33,('I','V'):29,
    ('L','K'):107,('L','M'):15,('L','F'):22,('L','P'):98,('L','S'):145,('L','T'):92,
    ('L','W'):61,('L','Y'):36,('L','V'):32,
    ('K','M'):95,('K','F'):102,('K','P'):103,('K','S'):121,('K','T'):78,('K','W'):110,
    ('K','Y'):85,('K','V'):97,
    ('M','F'):28,('M','P'):87,('M','S'):135,('M','T'):81,('M','W'):67,('M','Y'):36,('M','V'):21,
    ('F','P'):114,('F','S'):155,('F','T'):103,('F','W'):40,('F','Y'):22,('F','V'):50,
    ('P','S'):74,('P','T'):38,('P','W'):147,('P','Y'):110,('P','V'):68,
    ('S','T'):58,('S','W'):177,('S','Y'):144,('S','V'):124,
    ('T','W'):128,('T','Y'):92,('T','V'):69,
    ('W','Y'):37,('W','V'):88,
    ('Y','V'):55,
}

BURIAL_RANK = {'unknown': 0, 'surface_exposed': 1, 'partially_buried': 2, 'buried_core': 3}
RANK_TO_BURIAL = {v: k for k, v in BURIAL_RANK.items()}


# Helpers used throughout the package
def sf(v, d=0.0):
    """Safe float conversion."""
    import pandas as pd
    if pd.isna(v) or v is None:
        return d
    try:
        return float(v)
    except (TypeError, ValueError):
        return d


def si(v, d=0):
    """Safe int conversion."""
    import pandas as pd
    if pd.isna(v) or v is None:
        return d
    try:
        return int(v)
    except (TypeError, ValueError):
        return d


def ss(v):
    """Safe string conversion."""
    import pandas as pd
    return '' if pd.isna(v) or v is None else str(v)


def sb(v, d=False):
    """Safe bool conversion."""
    import pandas as pd
    if pd.isna(v) or v is None:
        return d
    return bool(v)


def get_grantham(a1, a2):
    """Look up Grantham distance; symmetric."""
    import pandas as pd
    if pd.isna(a1) or pd.isna(a2):
        return -1
    a1, a2 = str(a1).upper(), str(a2).upper()
    if a1 == a2:
        return 0
    return GRANTHAM.get((a1, a2), GRANTHAM.get((a2, a1), -1))


def classify_grantham(d):
    """Classify Grantham distance into conservative/moderate/radical."""
    import pandas as pd
    if pd.isna(d) or d is None or d < 0:
        return 'unknown'
    d = int(d)
    if d <= 50:
        return 'conservative'
    elif d <= 100:
        return 'moderately_conservative'
    elif d <= 150:
        return 'moderately_radical'
    else:
        return 'radical'


def grantham_severity(d):
    """Grantham severity score (0-4)."""
    import pandas as pd
    if pd.isna(d) or d is None or d < 0:
        return 0.0
    return min(4.0, float(d) / 53.75)


def get_property_changes(r, a):
    """Summarize size/charge/hydrophobicity changes."""
    import pandas as pd
    if pd.isna(r) or pd.isna(a):
        return 'unknown'
    p1 = AA_PROPERTIES.get(str(r).upper(), {})
    p2 = AA_PROPERTIES.get(str(a).upper(), {})
    if not p1 or not p2:
        return 'unknown'
    ch = []
    for k in ['size', 'charge', 'hydrophobic']:
        if p1.get(k) != p2.get(k):
            ch.append(f"{k}:{p1[k]}->{p2[k]}")
    return ';'.join(ch) if ch else 'none'
