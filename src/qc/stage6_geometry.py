"""
Stage 6 — CDR3 Geometry: Gordon rho/compactness (resolved) and
Zhang FR2_CDR3_Dis (resolved). Gordon theta/phi are NOT implemented
(see results/audits/stage6_gordon_method_definition.md).
"""

import sys
from pathlib import Path
from typing import Dict, Optional, Tuple, List
import numpy as np
import re

sys.path.insert(0, str(Path(__file__).parent))
from mmcif_qc import MMCIFStructure

CDR3_START, CDR3_END = 105, 117
GORDON_ANCHOR_POSITIONS = [102, 103, 118, 119]  # per corroborating author code (see method-definition doc)
GORDON_QC_ANCHOR_WINDOW = [102, 103, 104, 118, 119, 120]  # 3-residue-each-side QC window per paper text
BACKBONE_ATOMS = {'N', 'CA', 'C'}  # Gordon anchor/CDR3 centroid uses N,CA,C (not O) per author code

_RESIDUE_STR_RE = re.compile(r'^(-?\d+)([A-Za-z]*)$')


def _parse_residue_str(residue_str: str) -> Tuple[int, str]:
    m = _RESIDUE_STR_RE.match(residue_str)
    return int(m.group(1)), m.group(2)


def get_cdr3_real_residues(structure: MMCIFStructure, chain_id: str, model_id: str) -> List[Dict]:
    """Real CDR3 residues (105-117) in true sequence order, from build_residue_scheme."""
    scheme = structure.build_residue_scheme(chain_id, model_id)
    return [r for r in scheme if r['auth_seq_id'] is not None and CDR3_START <= r['auth_seq_id'] <= CDR3_END]


def compute_cdr3_length(structure: MMCIFStructure, chain_id: str, model_id: str) -> Dict:
    cdr3_residues = get_cdr3_real_residues(structure, chain_id, model_id)
    sequence_length = len(cdr3_residues)
    resolved = [r for r in cdr3_residues if r['observed']]
    resolved_length = len(resolved)
    coverage = resolved_length / sequence_length if sequence_length else None
    return {
        'cdr3_sequence_length': sequence_length,
        'cdr3_resolved_length': resolved_length,
        'cdr3_coordinate_coverage': coverage,
    }


def _get_backbone_centroid(structure: MMCIFStructure, chain_id: str, model_id: str, positions: List[int]) -> Optional[np.ndarray]:
    """Centroid of N,CA,C atoms across a list of base IMGT positions (no insertion codes considered
    -- used only for the fixed 4-point Gordon anchor set, none of which carry insertions in this cohort)."""
    coords = []
    for pos in positions:
        atoms = structure.get_atoms_for_chain_residue(chain_id, pos, model_id, '')
        if not atoms:
            return None
        atoms = structure.select_altloc(atoms)
        backbone = {a['name']: a for a in atoms if a['name'] in BACKBONE_ATOMS}
        if len(backbone) < 3:
            return None
        for name in BACKBONE_ATOMS:
            a = backbone[name]
            coords.append([a['x'], a['y'], a['z']])
    if not coords:
        return None
    return np.mean(coords, axis=0)


def _get_cdr3_backbone_centroid(structure: MMCIFStructure, chain_id: str, model_id: str) -> Optional[np.ndarray]:
    cdr3_residues = get_cdr3_real_residues(structure, chain_id, model_id)
    coords = []
    for r in cdr3_residues:
        if not r['observed']:
            return None  # any unresolved real CDR3 residue blocks rho (needs full CDR3 backbone)
        atoms = structure.get_atoms_for_chain_residue(chain_id, r['auth_seq_id'], model_id, r['ins_code'])
        if not atoms:
            return None
        atoms = structure.select_altloc(atoms)
        backbone = {a['name']: a for a in atoms if a['name'] in BACKBONE_ATOMS}
        if len(backbone) < 3:
            return None
        for name in BACKBONE_ATOMS:
            a = backbone[name]
            coords.append([a['x'], a['y'], a['z']])
    if not coords:
        return None
    return np.mean(coords, axis=0)


def _check_qc_anchor_window(structure: MMCIFStructure, chain_id: str, model_id: str) -> bool:
    """Paper's 3-residue-each-side backbone-completeness QC criterion (102-104, 118-120)."""
    for pos in GORDON_QC_ANCHOR_WINDOW:
        atoms = structure.get_atoms_for_chain_residue(chain_id, pos, model_id, '')
        if not atoms:
            return False
        atoms = structure.select_altloc(atoms)
        backbone = {a['name'] for a in atoms if a['name'] in BACKBONE_ATOMS}
        if len(backbone) < 3:
            return False
    return True


def compute_gordon_rho_compactness(structure: MMCIFStructure, chain_id: str, model_id: str,
                                     cdr3_sequence_length: int) -> Dict:
    result = {
        'cdr3_rho': None, 'cdr3_compactness': None, 'cdr3_theta': None, 'cdr3_phi': None,
        'gordon_geometry_valid': False, 'geometry_missing_reason': None,
    }

    qc_ok = _check_qc_anchor_window(structure, chain_id, model_id)
    if not qc_ok:
        result['geometry_missing_reason'] = 'GORDON_ANCHOR_QC_WINDOW_INCOMPLETE_BACKBONE'
        return result

    anchor_centroid = _get_backbone_centroid(structure, chain_id, model_id, GORDON_ANCHOR_POSITIONS)
    if anchor_centroid is None:
        result['geometry_missing_reason'] = 'GORDON_ANCHOR_ATOMS_MISSING'
        return result

    cdr3_centroid = _get_cdr3_backbone_centroid(structure, chain_id, model_id)
    if cdr3_centroid is None:
        result['geometry_missing_reason'] = 'GORDON_CDR3_BACKBONE_INCOMPLETE'
        return result

    rho = float(np.linalg.norm(cdr3_centroid - anchor_centroid))
    if rho <= 0 or cdr3_sequence_length == 0:
        result['geometry_missing_reason'] = 'GORDON_RHO_ZERO_OR_NO_CDR3'
        return result

    compactness = cdr3_sequence_length / rho

    result['cdr3_rho'] = rho
    result['cdr3_compactness'] = compactness
    result['gordon_geometry_valid'] = True
    result['geometry_missing_reason'] = 'THETA_PHI_NOT_IMPLEMENTED_NO_SUPPLEMENT'  # rho/compactness valid; theta/phi never computed
    return result


def _get_residue_centroid_all_atoms(structure: MMCIFStructure, chain_id: str, model_id: str,
                                      auth_seq_id: int, ins_code: str) -> Optional[np.ndarray]:
    atoms = structure.get_atoms_for_chain_residue(chain_id, auth_seq_id, model_id, ins_code)
    if not atoms:
        return None
    atoms = structure.select_altloc(atoms)
    atoms = [a for a in atoms if a['element'].strip().upper() not in ('H', 'D')]
    if not atoms:
        return None
    coords = np.array([[a['x'], a['y'], a['z']] for a in atoms])
    return coords.mean(axis=0)


def compute_zhang_fr2_cdr3_dis(structure: MMCIFStructure, chain_id: str, model_id: str) -> Dict:
    result = {'FR2_CDR3_Dis': None, 'zhang_metric_valid': False, 'zhang_missing_reason': None}

    pos42_centroid = _get_residue_centroid_all_atoms(structure, chain_id, model_id, 42, '')
    if pos42_centroid is None:
        result['zhang_missing_reason'] = 'POSITION_42_UNRESOLVED'
        return result

    cdr3_residues = [r for r in get_cdr3_real_residues(structure, chain_id, model_id) if r['observed']]
    if len(cdr3_residues) <= 4:
        result['zhang_missing_reason'] = 'CDR3_TOO_SHORT_FOR_INTERNAL_SET'
        return result

    # True sequence order via label_seq_id (already the order returned by get_cdr3_real_residues);
    # exclude first two and last two REAL residues in that order.
    internal_residues = cdr3_residues[2:-2]
    if not internal_residues:
        result['zhang_missing_reason'] = 'CDR3_TOO_SHORT_FOR_INTERNAL_SET'
        return result

    distances = []
    for r in internal_residues:
        centroid = _get_residue_centroid_all_atoms(structure, chain_id, model_id, r['auth_seq_id'], r['ins_code'])
        if centroid is not None:
            distances.append(float(np.linalg.norm(pos42_centroid - centroid)))

    if not distances:
        result['zhang_missing_reason'] = 'INTERNAL_CDR3_RESIDUES_UNRESOLVED'
        return result

    result['FR2_CDR3_Dis'] = min(distances)
    result['zhang_metric_valid'] = True
    return result
