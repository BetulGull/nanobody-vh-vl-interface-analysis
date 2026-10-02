"""
Stage 3 — FreeSASA Wrapper

Deterministic, explicitly-configured FreeSASA calculation on variable-domain
atom subsets extracted by stage3_variable_domain.py.

FreeSASA version: freesasa-python 2.2.1 (conda-forge)
Algorithm: Lee & Richards (explicit; not left to library default drift)
Probe radius: 1.4 Å (explicit)
Classifier: FreeSASA default (NACCESS-derived van der Waals radii; natively
recognizes MSE — verified, see stage3_variable_domain_validation.py note)
"""

from typing import Dict, List, Any, Tuple
import string
import freesasa

FREESASA_PACKAGE_VERSION = "2.2.1"
PROBE_RADIUS_ANGSTROM = 1.4
ALGORITHM = freesasa.LeeRichards

freesasa.setVerbosity(freesasa.nowarnings)


def get_parameters() -> freesasa.Parameters:
    params = freesasa.Parameters()
    params.setAlgorithm(ALGORITHM)
    params.setProbeRadius(PROBE_RADIUS_ANGSTROM)
    return params


def _build_chain_surrogate_map(atoms: List[Dict[str, Any]]) -> Dict[str, str]:
    """
    Map each distinct real auth_asym_id present in `atoms` to a guaranteed-
    unique single-character surrogate label.

    REQUIRED WORKAROUND: freesasa.Structure.addAtom silently truncates any
    chainLabel to its first character (confirmed empirically: a real chain
    named 'HHH' comes back keyed as 'H' in result.residueAreas()). Relying on
    freesasa to echo back our real chain ID is therefore unsafe whenever a
    PDB uses extended/multi-character chain naming (encountered in this
    cohort — e.g. PDB entries with chains 'HHH'/'LLL'/'AAA'). Assigning our
    own single-character surrogate per call (at most 2 distinct chains per
    calc_sasa call here: VH alone, or VH+VL) and translating back ourselves
    avoids any dependency on freesasa's truncated echo, and rules out any
    truncation collision (e.g. two different real chains both starting with
    'H') since surrogates are assigned by first-appearance order, not by the
    real label's first character.
    """
    surrogate_map = {}
    available = iter(string.ascii_uppercase + string.digits)
    for atom in atoms:
        real_chain = atom['auth_asym_id']
        if real_chain not in surrogate_map:
            surrogate_map[real_chain] = next(available)
    return surrogate_map


def build_structure(atoms: List[Dict[str, Any]], surrogate_map: Dict[str, str]) -> freesasa.Structure:
    """
    Build a freesasa.Structure from an explicit atom list (as produced by
    stage3_variable_domain.extract_variable_domain_atoms), preserving
    insertion codes by packing them into the residue-number string (verified
    behavior: freesasa keys residueAreas() by this exact string, e.g. '111A'),
    and using single-character surrogate chain labels (see
    _build_chain_surrogate_map) to avoid freesasa's chain-label truncation.
    """
    s = freesasa.Structure()
    for atom in atoms:
        residue_number_str = f"{atom['auth_seq_id']}{atom['ins_code']}"
        s.addAtom(
            atom['atom_name'],
            atom['residue_name'],
            residue_number_str,
            surrogate_map[atom['auth_asym_id']],
            atom['x'], atom['y'], atom['z'],
        )
    return s


def calc_sasa(atoms: List[Dict[str, Any]]) -> Tuple[Dict[Tuple[str, str], float], float]:
    """
    Compute SASA for a given atom set.

    Returns:
        (per_residue_sasa, total_sasa)
        per_residue_sasa: dict mapping (real_chain_id, resnum_str) -> total residue SASA (Å²)
            keyed by the TRUE auth_asym_id (translated back from freesasa's
            surrogate label — see _build_chain_surrogate_map), and the packed
            'auth_seq_id+ins_code' string used in build_structure.
        total_sasa: total SASA (Å²) as reported by FreeSASA directly (result.totalArea())
    """
    if not atoms:
        return {}, 0.0

    surrogate_map = _build_chain_surrogate_map(atoms)
    reverse_map = {v: k for k, v in surrogate_map.items()}

    structure = build_structure(atoms, surrogate_map)
    params = get_parameters()
    result = freesasa.calc(structure, params)

    residue_areas = result.residueAreas()
    per_residue = {}
    for surrogate_chain, resmap in residue_areas.items():
        real_chain = reverse_map[surrogate_chain]
        for resnum_str, area in resmap.items():
            per_residue[(real_chain, resnum_str)] = area.total

    total_sasa = result.totalArea()
    return per_residue, total_sasa


def get_freesasa_config_record() -> Dict[str, Any]:
    """Return the exact numerical configuration used, for audit/reproducibility."""
    return {
        'freesasa_python_package_version': FREESASA_PACKAGE_VERSION,
        'algorithm': 'LeeRichards',
        'probe_radius_angstrom': PROBE_RADIUS_ANGSTROM,
        'classifier': 'default (NACCESS-derived)',
    }
