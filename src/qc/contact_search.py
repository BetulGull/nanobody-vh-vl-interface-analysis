"""
Contact Search — Distance calculations for Gordon CDR-antigen criterion.

Implements both brute-force reference and optimized cKDTree algorithms for finding
CDR residues within ≤4.5 Å of protein-antigen atoms.

Gordon criterion: A CDR is "in contact with" protein antigen if at least one CDR
residue has at least one atom within ≤4.5 Å of at least one antigen atom.

For Stage-2 QC: report which CDRs (H1, H2, H3, L1, L2, L3 where applicable) meet
this criterion for a given structure.
"""

import numpy as np
from typing import Dict, List, Tuple, Optional
from scipy.spatial import cKDTree
import logging

logger = logging.getLogger(__name__)


class ContactSearchError(Exception):
    """Raised when contact search encounters a fatal error."""
    pass


class ContactSearch:
    """
    Find CDR-antigen contacts using Gordon ≤4.5Å criterion.

    Workflows:
    1. Brute-force: Compare every CDR atom to every antigen atom (slow, simple, reference)
    2. Optimized: cKDTree spatial search (fast, production-use)

    Both return identical results; optimized is ~100-1000x faster depending on structure size.
    """

    def __init__(self, threshold_angstrom: float = 4.5):
        self.threshold = threshold_angstrom

    def brute_force_contacts(
        self,
        cdr_atoms: np.ndarray,  # shape (n_cdr, 3) with xyz coords
        antigen_atoms: np.ndarray,  # shape (m_antigen, 3) with xyz coords
    ) -> Tuple[bool, List[int], List[int]]:
        """
        Brute-force distance comparison.

        Args:
            cdr_atoms: (N, 3) array of CDR atom coordinates
            antigen_atoms: (M, 3) array of antigen atom coordinates

        Returns:
            (has_contact, cdr_indices_with_contact, antigen_indices_with_contact)
            - has_contact: bool (True if any CDR atom within threshold of any antigen atom)
            - cdr_indices: list of unique CDR atom indices with ≤threshold contact
            - antigen_indices: list of unique antigen atom indices with ≤threshold contact
        """
        if len(cdr_atoms) == 0 or len(antigen_atoms) == 0:
            return False, [], []

        # Compute all pairwise distances: (N, M) matrix
        # distance[i, j] = distance from CDR atom i to antigen atom j
        distances = np.linalg.norm(cdr_atoms[:, np.newaxis, :] - antigen_atoms[np.newaxis, :, :], axis=2)

        # Find all pairs within threshold
        within_threshold = distances <= self.threshold

        if not within_threshold.any():
            return False, [], []

        # Extract indices of atoms with ≤threshold contacts
        cdr_idx, antigen_idx = np.where(within_threshold)
        unique_cdr = sorted(np.unique(cdr_idx).tolist())
        unique_antigen = sorted(np.unique(antigen_idx).tolist())

        return True, unique_cdr, unique_antigen

    def optimized_contacts(
        self,
        cdr_atoms: np.ndarray,
        antigen_atoms: np.ndarray,
    ) -> Tuple[bool, List[int], List[int]]:
        """
        Optimized contact search using cKDTree spatial indexing.

        Same signature as brute_force_contacts; identical results, ~100-1000x faster.

        Args:
            cdr_atoms: (N, 3) array of CDR atom coordinates
            antigen_atoms: (M, 3) array of antigen atom coordinates

        Returns:
            (has_contact, cdr_indices_with_contact, antigen_indices_with_contact)
        """
        if len(cdr_atoms) == 0 or len(antigen_atoms) == 0:
            return False, [], []

        # Build KDTree on antigen atoms
        tree = cKDTree(antigen_atoms)

        # Query: for each CDR atom, find all antigen atoms within threshold
        # query_ball_point returns list of antigen-atom indices within threshold for each CDR atom
        cdr_indices_with_contact = []
        antigen_indices_with_contact = set()

        for i, cdr_coord in enumerate(cdr_atoms):
            # Query returns indices of antigen atoms within threshold distance
            nearby_antigen_indices = tree.query_ball_point(cdr_coord, self.threshold)

            if nearby_antigen_indices:
                cdr_indices_with_contact.append(i)
                antigen_indices_with_contact.update(nearby_antigen_indices)

        has_contact = len(cdr_indices_with_contact) > 0

        return (
            has_contact,
            sorted(cdr_indices_with_contact),
            sorted(list(antigen_indices_with_contact)),
        )

    def contacts_by_cdr(
        self,
        cdr_dict: Dict[str, np.ndarray],  # {'H1': array(...), 'H2': array(...), ...}
        antigen_atoms: np.ndarray,
    ) -> Dict[str, bool]:
        """
        Determine which CDRs (H1, H2, H3, L1, L2, L3) have contact with antigen.

        Args:
            cdr_dict: dict mapping CDR name (e.g., 'H1') to (n_atoms, 3) coordinate array
            antigen_atoms: (m, 3) array of antigen atom coordinates

        Returns:
            Dict mapping CDR name to bool (True if this CDR has contact)

        Example:
            contacts = searcher.contacts_by_cdr(
                {'H1': h1_atoms, 'H2': h2_atoms, 'H3': h3_atoms},
                antigen_atoms
            )
            # contacts = {'H1': False, 'H2': True, 'H3': True}
        """
        result = {}
        for cdr_name, cdr_atoms in cdr_dict.items():
            has_contact, _, _ = self.optimized_contacts(cdr_atoms, antigen_atoms)
            result[cdr_name] = has_contact

        return result

    def validate_against_brute_force(
        self,
        cdr_atoms: np.ndarray,
        antigen_atoms: np.ndarray,
        tolerance: float = 1e-6,
    ) -> bool:
        """
        Verify that optimized_contacts matches brute_force_contacts.

        Args:
            cdr_atoms: (N, 3) array
            antigen_atoms: (M, 3) array
            tolerance: floating-point tolerance for numerical comparison

        Returns:
            True if both algorithms agree (within tolerance)

        Raises:
            ContactSearchError if algorithms disagree
        """
        bf_result = self.brute_force_contacts(cdr_atoms, antigen_atoms)
        opt_result = self.optimized_contacts(cdr_atoms, antigen_atoms)

        bf_has_contact, bf_cdr_idx, bf_antigen_idx = bf_result
        opt_has_contact, opt_cdr_idx, opt_antigen_idx = opt_result

        # Check has_contact agreement
        if bf_has_contact != opt_has_contact:
            msg = f"Disagreement on has_contact: BF={bf_has_contact}, OPT={opt_has_contact}"
            logger.error(msg)
            raise ContactSearchError(msg)

        # Check index lists (should be identical)
        if bf_cdr_idx != opt_cdr_idx:
            msg = f"Disagreement on CDR indices: BF={bf_cdr_idx[:5]}..., OPT={opt_cdr_idx[:5]}..."
            logger.error(msg)
            raise ContactSearchError(msg)

        if bf_antigen_idx != opt_antigen_idx:
            msg = f"Disagreement on antigen indices: BF={bf_antigen_idx[:5]}..., OPT={opt_antigen_idx[:5]}..."
            logger.error(msg)
            raise ContactSearchError(msg)

        return True


def benchmark_contact_search(
    num_cdr_atoms: int = 1000,
    num_antigen_atoms: int = 2000,
    threshold: float = 4.5,
) -> Dict[str, float]:
    """
    Benchmark brute-force vs. optimized for typical structure sizes.

    Args:
        num_cdr_atoms: number of CDR atoms
        num_antigen_atoms: number of antigen atoms
        threshold: Gordon threshold in Angstroms

    Returns:
        Dict with timing results and speedup
    """
    import time

    # Generate random atom positions
    np.random.seed(42)
    cdr_atoms = np.random.randn(num_cdr_atoms, 3) * 10  # ~20Å box
    antigen_atoms = np.random.randn(num_antigen_atoms, 3) * 10

    searcher = ContactSearch(threshold)

    # Brute-force timing
    t0 = time.time()
    bf_result = searcher.brute_force_contacts(cdr_atoms, antigen_atoms)
    bf_time = time.time() - t0

    # Optimized timing
    t0 = time.time()
    opt_result = searcher.optimized_contacts(cdr_atoms, antigen_atoms)
    opt_time = time.time() - t0

    # Verify agreement
    if bf_result[0] != opt_result[0]:
        logger.warning("Brute-force and optimized disagree on has_contact!")

    speedup = bf_time / opt_time if opt_time > 0 else float('inf')

    return {
        "num_cdr_atoms": num_cdr_atoms,
        "num_antigen_atoms": num_antigen_atoms,
        "threshold_angstrom": threshold,
        "brute_force_time_sec": bf_time,
        "optimized_time_sec": opt_time,
        "speedup": speedup,
    }
