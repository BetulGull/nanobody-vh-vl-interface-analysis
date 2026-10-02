"""
Stage 4 — Centralized Sequence Chemistry Definitions

Single source of truth for the Kyte-Doolittle hydropathy scale and residue
classifications, so no script redefines these independently.

IMPORTANT: values here are INTRINSIC sequence-level chemistry (Kyte-Doolittle
1982 hydropathy index), computed with no reference to solvent accessibility
or structural burial. Never call this "surface hydrophobicity" or "exposed
hydrophobic burden" -- that requires SASA (Stage 5+).
"""

THREE_TO_ONE = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q',
    'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
    'MET': 'M', 'MSE': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T',
    'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
}

# Kyte & Doolittle (1982) hydropathy index, keyed by one-letter amino acid.
KYTE_DOOLITTLE = {
    'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5, 'Q': -3.5,
    'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5, 'L': 3.8, 'K': -3.9,
    'M': 1.9, 'F': 2.8, 'P': -1.6, 'S': -0.8, 'T': -0.7, 'W': -0.9,
    'Y': -1.3, 'V': 4.2,
}

ACIDIC = {'D', 'E'}
BASIC = {'K', 'R'}
AROMATIC = {'F', 'W', 'Y'}
GLYCINE = {'G'}
PROLINE = {'P'}
CYSTEINE = {'C'}
HISTIDINE = {'H'}

CANONICAL_20 = set(KYTE_DOOLITTLE.keys())


def get_kd(one_letter: str) -> float:
    """Kyte-Doolittle intrinsic hydropathy for a one-letter amino acid code."""
    return KYTE_DOOLITTLE.get(one_letter)


def residue_class_flags(one_letter: str) -> dict:
    """Boolean class membership flags for one residue (sequence-level, no SASA)."""
    kd = get_kd(one_letter)
    return {
        'residue_acidic': one_letter in ACIDIC,
        'residue_basic': one_letter in BASIC,
        'residue_aromatic': one_letter in AROMATIC,
        'residue_positive_KD': (kd is not None) and (kd > 0),
        'residue_glycine': one_letter in GLYCINE,
        'residue_proline': one_letter in PROLINE,
        'residue_cysteine': one_letter in CYSTEINE,
        'residue_histidine': one_letter in HISTIDINE,
    }
