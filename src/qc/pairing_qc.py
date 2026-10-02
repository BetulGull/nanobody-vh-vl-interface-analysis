"""
Pairing QC — VH-VL physical-proximity sanity check.

Permissive validation that paired structures have both chains and they are
in physical proximity (within 5.0 Å). No SASA/BSA/interface-residue derivation.
"""

import numpy as np
from typing import List, Tuple
import logging
import re

logger = logging.getLogger(__name__)


_RESIDUE_STR_RE = re.compile(r'^(-?\d+)([A-Za-z]*)$')


def _parse_residue_str(residue_str: str) -> Tuple[int, str]:
    """Split '112A' -> (112, 'A'); '112' -> (112, '')."""
    m = _RESIDUE_STR_RE.match(residue_str)
    if not m:
        raise ValueError(f"Cannot parse residue string: {residue_str!r}")
    return int(m.group(1)), m.group(2)


class PairingQuality:
    """
    Validate VH-VL physical pairing.
    """

    def __init__(self, threshold_angstrom: float = 5.0):
        self.threshold = threshold_angstrom

    def check_physical_pair(
        self,
        vh_atoms: List[Tuple[float, float, float]],
        vl_atoms: List[Tuple[float, float, float]],
    ) -> Tuple[bool, float]:
        """
        Check if VH and VL are physically proximal.

        Args:
            vh_atoms: List of (x, y, z) coordinates for VH atoms
            vl_atoms: List of (x, y, z) coordinates for VL atoms

        Returns:
            (is_paired, min_distance)
        """
        if not vh_atoms or not vl_atoms:
            return False, float('inf')

        vh_coords = np.array(vh_atoms)
        vl_coords = np.array(vl_atoms)

        distances = np.linalg.norm(vh_coords[:, np.newaxis, :] - vl_coords[np.newaxis, :, :], axis=2)
        min_dist = np.min(distances)

        is_paired = bool(min_dist <= self.threshold)

        return is_paired, float(min_dist)

    def get_variable_domain_atoms(
        self,
        mmcif_structure,
        chain_id: str,
        model_id: str,
    ) -> List[Tuple[float, float, float]]:
        """
        Extract variable-domain atom coordinates (IMGT 1-128) for a chain,
        including every insertion-coded residue observed.

        Args:
            mmcif_structure: MMCIFStructure instance
            chain_id: Chain ID
            model_id: Model number

        Returns:
            List of (x, y, z) tuples for all variable-domain atoms
        """
        full_mapping = mmcif_structure.build_full_imgt_mapping_with_insertions(chain_id, model_id)

        atoms_out = []
        for imgt_pos in range(1, 129):
            for residue_str in full_mapping.get(imgt_pos, []):
                base_pos, ins_code = _parse_residue_str(residue_str)
                residue_atoms = mmcif_structure.get_atoms_for_chain_residue(
                    chain_id, base_pos, model_id, ins_code
                )
                if residue_atoms:
                    residue_atoms = mmcif_structure.select_altloc(residue_atoms)
                    for atom in residue_atoms:
                        atoms_out.append((atom['x'], atom['y'], atom['z']))

        return atoms_out
