"""
Stage 7 — Pairwise Mutual Buried Surface Area (BSA) between CDR3 and the
former-VH-VL-interface packing target set.

Terminology (mandated): "pairwise mutual buried surface area", never just
"the BSA". Distinct from Stage 6 "direct shielding" (VH-vs-VH-without-CDR3
counterfactual) and distinct from Arpeggio "specific interactions" (atom-
pair contact typing) -- this module computes purely geometric SASA burial
between two native-coordinate atom subsets, same FreeSASA engine validated
in Stages 3, 5, 6.

Convention:
  SASA_A_alone   = total SASA of CDR3 atoms alone (isolated calc_sasa call)
  SASA_B_alone   = total SASA of target-set atoms alone (isolated calc_sasa call)
  SASA_AB_total  = total SASA of CDR3+target atoms together (native coordinates)
  BSA_pair_total = SASA_A_alone + SASA_B_alone - SASA_AB_total
  BSA_CDR3_side  = SASA_A_alone - sum(per-residue AB-complex SASA over CDR3 residues)
  BSA_target_side= SASA_B_alone - sum(per-residue AB-complex SASA over target residues)
Consistency check (required): BSA_CDR3_side + BSA_target_side == BSA_pair_total
(exact by construction -- both are the same additive decomposition of the
same total-area difference; a mismatch would indicate a residue-membership
bug, not a real physical discrepancy).
"""

import sys
from pathlib import Path
from typing import List, Dict, Set, Tuple

sys.path.insert(0, str(Path(__file__).parent))
from stage3_sasa import calc_sasa

CDR3_START, CDR3_END = 105, 117

NONADJACENT_PACKING_TARGET: Set[int] = {40, 42, 44, 49, 50, 52, 55, 66, 67, 68, 69, 103, 119, 120}
FULL_SHIELDABLE_TARGET: Set[int] = {40, 42, 44, 49, 50, 52, 55, 66, 67, 68, 69, 103, 118, 119, 120}


def select_cdr3_atoms(vhh_atoms: List[Dict]) -> List[Dict]:
    return [a for a in vhh_atoms if CDR3_START <= a['auth_seq_id'] <= CDR3_END]


def select_target_atoms(vhh_atoms: List[Dict], target_set: Set[int]) -> List[Dict]:
    return [a for a in vhh_atoms if a['auth_seq_id'] in target_set]


def _resnum_key(a: Dict) -> str:
    icode = a['ins_code'] if a['ins_code'] else ''
    return f"{a['auth_seq_id']}{icode}"


def compute_pairwise_bsa(vhh_atoms: List[Dict], target_set: Set[int]) -> Dict:
    cdr3_atoms = select_cdr3_atoms(vhh_atoms)
    target_atoms = select_target_atoms(vhh_atoms, target_set)

    if not cdr3_atoms or not target_atoms:
        return {
            'n_cdr3_atoms': len(cdr3_atoms), 'n_target_atoms': len(target_atoms),
            'SASA_A_alone': None, 'SASA_B_alone': None, 'SASA_AB_total': None,
            'BSA_CDR3_side': None, 'BSA_target_side': None, 'BSA_pair_total': None,
            'bsa_consistency_ok': None,
            'status': 'TARGET_NOT_EVALUABLE',
        }

    per_res_a, sasa_a_alone = calc_sasa(cdr3_atoms)
    per_res_b, sasa_b_alone = calc_sasa(target_atoms)
    per_res_ab, sasa_ab_total = calc_sasa(cdr3_atoms + target_atoms)

    cdr3_resnum_keys = {(a['auth_asym_id'], _resnum_key(a)) for a in cdr3_atoms}
    target_resnum_keys = {(a['auth_asym_id'], _resnum_key(a)) for a in target_atoms}

    sum_a_in_complex = sum(area for key, area in per_res_ab.items() if key in cdr3_resnum_keys)
    sum_b_in_complex = sum(area for key, area in per_res_ab.items() if key in target_resnum_keys)

    bsa_cdr3_side = sasa_a_alone - sum_a_in_complex
    bsa_target_side = sasa_b_alone - sum_b_in_complex
    bsa_pair_total = sasa_a_alone + sasa_b_alone - sasa_ab_total

    consistency_ok = abs((bsa_cdr3_side + bsa_target_side) - bsa_pair_total) < 1e-6

    return {
        'n_cdr3_atoms': len(cdr3_atoms), 'n_target_atoms': len(target_atoms),
        'SASA_A_alone': sasa_a_alone, 'SASA_B_alone': sasa_b_alone, 'SASA_AB_total': sasa_ab_total,
        'BSA_CDR3_side': bsa_cdr3_side, 'BSA_target_side': bsa_target_side, 'BSA_pair_total': bsa_pair_total,
        'bsa_consistency_ok': consistency_ok,
        'status': 'OK',
    }
