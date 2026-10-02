"""
Stage 6 — CDR3 Removal Counterfactual and Shieldable-Interface Target

CRITICAL: the no-CDR3 state is a geometric SASA counterfactual only — same
native coordinates, CDR3 atoms simply omitted from the SASA occlusion
calculation. Never described as a physical deletion mutant.

CDR3 positions cannot shield themselves (once removed they have no
before/after comparison), so the shielding target set
(SHIELDABLE_EMPIRICAL_INTERFACE) is the Stage-3 28-position primary
interface MINUS all CDR3 (IMGT 105-117) positions.
"""

import sys
from pathlib import Path
from typing import List, Dict, Set, Tuple
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from mmcif_qc import MMCIFStructure
from stage3_variable_domain import extract_variable_domain_atoms
from stage3_sasa import calc_sasa

CDR3_START, CDR3_END = 105, 117


def get_shieldable_interface_positions(primary_map: pd.DataFrame) -> pd.DataFrame:
    """FULL_EMPIRICAL_INTERFACE (28 positions) minus CDR3 (105-117) positions."""
    full_interface = primary_map[primary_map['primary_binary_interface'] == True].copy()
    is_cdr3 = (full_interface['imgt_position'] >= CDR3_START) & (full_interface['imgt_position'] <= CDR3_END)
    shieldable = full_interface[~is_cdr3].copy()
    cdr3_interface = full_interface[is_cdr3].copy()
    return full_interface, cdr3_interface, shieldable


def extract_no_cdr3_atoms(vhh_atoms: List[Dict]) -> List[Dict]:
    """Same atoms/coordinates as vhh_atoms, with CDR3 (105-117) atoms omitted."""
    return [a for a in vhh_atoms if not (CDR3_START <= a['auth_seq_id'] <= CDR3_END)]


def compute_native_and_nocdr3_sasa(structure: MMCIFStructure, hchain: str, model_id: str) -> Tuple[dict, dict, list, list]:
    """
    Returns (per_res_native, per_res_nocdr3, native_atoms, nocdr3_atoms).
    per_res_* keyed by (chain, resnum_str) as returned by calc_sasa.
    """
    native_atoms = extract_variable_domain_atoms(structure, hchain, model_id)
    nocdr3_atoms = extract_no_cdr3_atoms(native_atoms)

    per_res_native, _ = calc_sasa(native_atoms)
    per_res_nocdr3, _ = calc_sasa(nocdr3_atoms)

    return per_res_native, per_res_nocdr3, native_atoms, nocdr3_atoms
