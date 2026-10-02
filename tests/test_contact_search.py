"""Unit tests for contact_search.py"""

import numpy as np
import pytest
from pathlib import Path

# Add src to path
import sys
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from qc.contact_search import ContactSearch, ContactSearchError


class TestContactSearch:
    """Test contact detection algorithms."""

    def test_empty_cdr_atoms(self):
        """Empty CDR array should return no contact."""
        searcher = ContactSearch(threshold_angstrom=4.5)
        cdr_atoms = np.array([]).reshape(0, 3)
        antigen_atoms = np.array([[0, 0, 0], [1, 1, 1]])

        has_contact, cdr_idx, antigen_idx = searcher.brute_force_contacts(cdr_atoms, antigen_atoms)
        assert has_contact is False
        assert len(cdr_idx) == 0
        assert len(antigen_idx) == 0

    def test_empty_antigen_atoms(self):
        """Empty antigen array should return no contact."""
        searcher = ContactSearch(threshold_angstrom=4.5)
        cdr_atoms = np.array([[0, 0, 0], [1, 1, 1]])
        antigen_atoms = np.array([]).reshape(0, 3)

        has_contact, cdr_idx, antigen_idx = searcher.brute_force_contacts(cdr_atoms, antigen_atoms)
        assert has_contact is False

    def test_single_close_contact(self):
        """Single CDR atom within threshold should detect contact."""
        searcher = ContactSearch(threshold_angstrom=5.0)
        cdr_atoms = np.array([[0, 0, 0]])
        antigen_atoms = np.array([[2, 0, 0]])  # distance = 2.0

        bf_has, bf_cdr, bf_antigen = searcher.brute_force_contacts(cdr_atoms, antigen_atoms)
        opt_has, opt_cdr, opt_antigen = searcher.optimized_contacts(cdr_atoms, antigen_atoms)

        assert bf_has is True
        assert opt_has is True
        assert bf_cdr == [0]
        assert opt_cdr == [0]
        assert bf_antigen == [0]
        assert opt_antigen == [0]

    def test_single_far_contact(self):
        """Single CDR atom beyond threshold should return no contact."""
        searcher = ContactSearch(threshold_angstrom=4.5)
        cdr_atoms = np.array([[0, 0, 0]])
        antigen_atoms = np.array([[10, 0, 0]])  # distance = 10.0

        bf_has, bf_cdr, bf_antigen = searcher.brute_force_contacts(cdr_atoms, antigen_atoms)
        opt_has, opt_cdr, opt_antigen = searcher.optimized_contacts(cdr_atoms, antigen_atoms)

        assert bf_has is False
        assert opt_has is False

    def test_multiple_atoms_partial_contact(self):
        """Multiple atoms with some within threshold."""
        searcher = ContactSearch(threshold_angstrom=4.5)
        cdr_atoms = np.array([
            [0, 0, 0],      # 0: will contact antigen 0 (dist=2.0)
            [10, 0, 0],     # 1: will not contact anything (dist≥10)
            [1, 1, 1],      # 2: will contact antigen 0 (dist≈1.7)
        ])
        antigen_atoms = np.array([
            [2, 0, 0],      # 0: contacted by CDR 0, 2
            [20, 0, 0],     # 1: too far
        ])

        bf_has, bf_cdr, bf_antigen = searcher.brute_force_contacts(cdr_atoms, antigen_atoms)
        opt_has, opt_cdr, opt_antigen = searcher.optimized_contacts(cdr_atoms, antigen_atoms)

        assert bf_has is True
        assert opt_has is True
        assert set(bf_cdr) == {0, 2}
        assert set(opt_cdr) == {0, 2}
        assert set(bf_antigen) == {0}
        assert set(opt_antigen) == {0}

    def test_brute_force_vs_optimized_small(self):
        """Brute-force and optimized must agree on small random data."""
        searcher = ContactSearch(threshold_angstrom=4.5)
        np.random.seed(42)
        cdr_atoms = np.random.randn(10, 3) * 5
        antigen_atoms = np.random.randn(15, 3) * 5

        bf_result = searcher.brute_force_contacts(cdr_atoms, antigen_atoms)
        opt_result = searcher.optimized_contacts(cdr_atoms, antigen_atoms)

        assert bf_result[0] == opt_result[0]  # has_contact
        assert bf_result[1] == opt_result[1]  # cdr_idx
        assert bf_result[2] == opt_result[2]  # antigen_idx

    def test_brute_force_vs_optimized_large(self):
        """Brute-force and optimized must agree on larger random data."""
        searcher = ContactSearch(threshold_angstrom=4.5)
        np.random.seed(123)
        cdr_atoms = np.random.randn(100, 3) * 10
        antigen_atoms = np.random.randn(200, 3) * 10

        bf_result = searcher.brute_force_contacts(cdr_atoms, antigen_atoms)
        opt_result = searcher.optimized_contacts(cdr_atoms, antigen_atoms)

        assert bf_result[0] == opt_result[0]
        assert bf_result[1] == opt_result[1]
        assert bf_result[2] == opt_result[2]

    def test_contacts_by_cdr(self):
        """Test CDR-level contact detection."""
        searcher = ContactSearch(threshold_angstrom=4.5)

        cdr_dict = {
            'H1': np.array([[0, 0, 0]]),           # within threshold of antigen
            'H2': np.array([[10, 0, 0]]),          # far from antigen
            'H3': np.array([[1, 0, 0]]),           # within threshold of antigen
        }
        antigen_atoms = np.array([[2, 0, 0], [20, 0, 0]])

        contacts = searcher.contacts_by_cdr(cdr_dict, antigen_atoms)

        assert contacts['H1'] is True   # CDR H1 within 4.5Å
        assert contacts['H2'] is False  # CDR H2 too far
        assert contacts['H3'] is True   # CDR H3 within 4.5Å

    def test_validate_against_brute_force_agreement(self):
        """Validation should succeed when algorithms agree."""
        searcher = ContactSearch(threshold_angstrom=4.5)
        np.random.seed(456)
        cdr_atoms = np.random.randn(50, 3) * 8
        antigen_atoms = np.random.randn(100, 3) * 8

        # Should not raise
        result = searcher.validate_against_brute_force(cdr_atoms, antigen_atoms)
        assert result is True

    def test_validate_threshold_exact_boundary(self):
        """Test atoms exactly at threshold distance."""
        searcher = ContactSearch(threshold_angstrom=4.5)
        cdr_atoms = np.array([[0, 0, 0]])
        # Antigen atom exactly 4.5 Angstroms away
        antigen_atoms = np.array([[4.5, 0, 0]])

        bf_has, _, _ = searcher.brute_force_contacts(cdr_atoms, antigen_atoms)
        opt_has, _, _ = searcher.optimized_contacts(cdr_atoms, antigen_atoms)

        # At exact boundary, should be included (≤ threshold)
        assert bf_has is True
        assert opt_has is True


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
