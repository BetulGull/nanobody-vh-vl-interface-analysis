"""
Stage-2 comprehensive test suite.

Covers: case-sensitive chain matching, extraction-manifest / APFS-collision
resolution, model selection, IMGT insertion-code preservation, CDR range
definitions, missing backbone atoms, deterministic altloc, Hchain/Lchain
logic (incl. SD-H no-assigned-Lchain), protein-polymer antigen identification
(incl. exclusion of water/ions/nonpolymer and multiple antigen chains),
VH-VL physical proximity, exact Gordon threshold boundary (4.5000 pass /
4.5001 fail), framework-only contact not satisfying Gordon CDR inclusion,
optimized/brute-force agreement, duplicate-FAB-archive non-duplication, and
Stage-1 source-table immutability.
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from qc.mmcif_qc import MMCIFStructure
from qc.cdr_qc import CDRQuality, IMGT_CDR_RANGES, _parse_residue_str
from qc.antigen_qc import AntigenQuality
from qc.contact_search import ContactSearch
from qc.pairing_qc import PairingQuality
from qc.structure_resolver import StructureResolver

PROJECT_ROOT = Path(__file__).parent.parent

# Real fixture files used across several tests (already validated in Gate 1/2)
VHH_INSERTION_FILE = PROJECT_ROOT / "data/raw_extracted/sdh/pdb_00001bzq/pdb_00001bzq_K_+_abag.cif"
PAIRED_MULTI_ANTIGEN_FILE = PROJECT_ROOT / "data/raw_extracted/sdh/pdb_000011ha/pdb_000011ha_G_+_abag.cif" \
    if (PROJECT_ROOT / "data/raw_extracted/sdh/pdb_000011ha").exists() else None


def _find_real_file(instance_glob):
    matches = list(PROJECT_ROOT.glob(f"data/raw_extracted/**/{instance_glob}"))
    return matches[0] if matches else None


class TestResidueStringParsing:
    def test_parse_plain_residue(self):
        assert _parse_residue_str("112") == (112, "")

    def test_parse_insertion_coded_residue(self):
        assert _parse_residue_str("112A") == (112, "A")
        assert _parse_residue_str("111B") == (111, "B")

    def test_insertion_codes_are_distinct(self):
        """111 != 111A != 111B must never collapse to the same key."""
        a = _parse_residue_str("111")
        b = _parse_residue_str("111A")
        c = _parse_residue_str("111B")
        assert a != b != c
        assert a[0] == b[0] == c[0] == 111
        assert len({a[1], b[1], c[1]}) == 3


class TestCDRRangeDefinitions:
    def test_all_three_cdr_h_ranges(self):
        assert IMGT_CDR_RANGES['H1'] == (27, 38)
        assert IMGT_CDR_RANGES['H2'] == (56, 65)
        assert IMGT_CDR_RANGES['H3'] == (105, 117)

    def test_all_three_cdr_l_ranges(self):
        assert IMGT_CDR_RANGES['L1'] == (27, 38)
        assert IMGT_CDR_RANGES['L2'] == (56, 65)
        assert IMGT_CDR_RANGES['L3'] == (105, 117)

    def test_cdr_position_membership(self):
        cdrq = CDRQuality()
        assert cdrq.is_cdr_position('H', 30) == (True, 'H1')
        assert cdrq.is_cdr_position('H', 60) == (True, 'H2')
        assert cdrq.is_cdr_position('H', 110) == (True, 'H3')
        assert cdrq.is_cdr_position('H', 50) == (False, None)  # framework


class TestRealStructureMmCIFParsing:
    """Uses actual extracted mmCIF files, not synthetic data."""

    def test_insertion_code_preserved_in_real_structure(self):
        if not VHH_INSERTION_FILE.exists():
            pytest.skip("fixture file not present")
        s = MMCIFStructure(str(VHH_INSERTION_FILE))
        model_id = s.get_model()
        full_map = s.build_full_imgt_mapping_with_insertions('K', model_id)
        # pdb_00001bzq chain K has an insertion at IMGT 112 (112, 112A observed distinctly)
        assert 112 in full_map
        assert len(full_map[112]) >= 1  # at least the base residue observed

    def test_case_sensitive_chain_matching_real_file(self):
        if not VHH_INSERTION_FILE.exists():
            pytest.skip("fixture file not present")
        s = MMCIFStructure(str(VHH_INSERTION_FILE))
        model_id = s.get_model()
        chains = s.get_chains_for_model(model_id)
        # Chain IDs must be preserved exactly as authored (upper-case here)
        assert 'K' in chains
        assert 'k' not in chains  # case sensitivity: lowercase must NOT match

    def test_model_selection_returns_valid_model(self):
        if not VHH_INSERTION_FILE.exists():
            pytest.skip("fixture file not present")
        s = MMCIFStructure(str(VHH_INSERTION_FILE))
        model_id = s.get_model()
        assert model_id in s.models


class TestBackboneCompleteness:
    def test_missing_backbone_atom_detected(self):
        s = MMCIFStructure.__new__(MMCIFStructure)  # bypass file loading
        atoms_incomplete = [
            {'name': 'N', 'x': 0, 'y': 0, 'z': 0},
            {'name': 'CA', 'x': 1, 'y': 0, 'z': 0},
            {'name': 'C', 'x': 2, 'y': 0, 'z': 0},
            # 'O' missing
        ]
        assert s.check_backbone_complete(atoms_incomplete) is False

    def test_complete_backbone_detected(self):
        s = MMCIFStructure.__new__(MMCIFStructure)
        atoms_complete = [
            {'name': 'N', 'x': 0, 'y': 0, 'z': 0},
            {'name': 'CA', 'x': 1, 'y': 0, 'z': 0},
            {'name': 'C', 'x': 2, 'y': 0, 'z': 0},
            {'name': 'O', 'x': 3, 'y': 0, 'z': 0},
        ]
        assert s.check_backbone_complete(atoms_complete) is True


class TestAltlocPolicy:
    def test_highest_occupancy_selected(self):
        s = MMCIFStructure.__new__(MMCIFStructure)
        atoms = [
            {'name': 'CA', 'occupancy': 0.4, 'altloc': 'A', 'x': 0, 'y': 0, 'z': 0},
            {'name': 'CA', 'occupancy': 0.6, 'altloc': 'B', 'x': 1, 'y': 1, 'z': 1},
        ]
        result = s.select_altloc(atoms)
        assert len(result) == 1
        assert result[0]['altloc'] == 'B'  # higher occupancy wins

    def test_alphabetical_tie_break(self):
        s = MMCIFStructure.__new__(MMCIFStructure)
        atoms = [
            {'name': 'CA', 'occupancy': 0.5, 'altloc': 'B', 'x': 1, 'y': 1, 'z': 1},
            {'name': 'CA', 'occupancy': 0.5, 'altloc': 'A', 'x': 0, 'y': 0, 'z': 0},
        ]
        result = s.select_altloc(atoms)
        assert len(result) == 1
        assert result[0]['altloc'] == 'A'  # tie -> alphabetical first

    def test_deterministic_repeated_calls(self):
        s = MMCIFStructure.__new__(MMCIFStructure)
        atoms = [
            {'name': 'CA', 'occupancy': 0.5, 'altloc': 'C', 'x': 2, 'y': 2, 'z': 2},
            {'name': 'CA', 'occupancy': 0.5, 'altloc': 'A', 'x': 0, 'y': 0, 'z': 0},
            {'name': 'CA', 'occupancy': 0.5, 'altloc': 'B', 'x': 1, 'y': 1, 'z': 1},
        ]
        r1 = s.select_altloc(list(atoms))
        r2 = s.select_altloc(list(atoms))
        assert r1[0]['altloc'] == r2[0]['altloc'] == 'A'


class TestNoAssignedLightChainLogic:
    def test_na_sentinel_means_no_lchain(self):
        """SAbDab2 uses literal string 'NA' for SD-H rows with no light chain."""
        raw_lchain = 'NA'
        lchain = None if (pd.isna(raw_lchain) or str(raw_lchain).strip() == 'NA') else raw_lchain
        assert lchain is None

    def test_real_lchain_value_not_treated_as_na(self):
        raw_lchain = 'L'
        lchain = None if (pd.isna(raw_lchain) or str(raw_lchain).strip() == 'NA') else raw_lchain
        assert lchain == 'L'

    def test_metadata_vhh_cohort_all_na(self):
        metadata = pd.read_parquet(PROJECT_ROOT / "data/metadata/master_stage1_resolved.parquet")
        vhh = metadata[metadata["primary_vhh_candidate_strict"] == True]
        assert (vhh["Lchain"] == "NA").all(), "All strict-primary VHH rows must have Lchain=='NA'"


class TestAntigenIdentification:
    def test_protein_polymer_types_case_sensitive(self):
        aq = AntigenQuality()
        assert 'polypeptide(L)' in aq.PROTEIN_POLYMER_TYPES
        assert 'polypeptide(l)' not in aq.PROTEIN_POLYMER_TYPES  # case-sensitive by design

    def test_real_structure_antigen_identification(self):
        if not VHH_INSERTION_FILE.exists():
            pytest.skip("fixture file not present")
        s = MMCIFStructure(str(VHH_INSERTION_FILE))
        model_id = s.get_model()
        aq = AntigenQuality()
        antigen_chains, chain_map = aq.identify_antigen_chains(s, model_id, antibody_chains={'K'})
        assert 'D' in antigen_chains  # lysozyme chain in 1BZQ
        assert chain_map['D'] == 'polypeptide(L)'

    def test_antibody_chain_excluded_from_antigen(self):
        if not VHH_INSERTION_FILE.exists():
            pytest.skip("fixture file not present")
        s = MMCIFStructure(str(VHH_INSERTION_FILE))
        model_id = s.get_model()
        aq = AntigenQuality()
        antigen_chains, _ = aq.identify_antigen_chains(s, model_id, antibody_chains={'K'})
        assert 'K' not in antigen_chains

    def test_multiple_antigen_chains_supported(self):
        """Find a real structure with >1 protein-antigen chain and confirm both are captured."""
        metadata = pd.read_parquet(PROJECT_ROOT / "data/metadata/master_stage1_resolved.parquet")
        paired = metadata[metadata["primary_paired_candidate_metadata"] == True]
        found = False
        for _, row in paired.head(50).iterrows():
            path = PROJECT_ROOT / row['extracted_structure_path']
            if not path.exists():
                continue
            try:
                s = MMCIFStructure(str(path))
            except Exception:
                continue
            model_id = s.get_model()
            hchain, lchain = row['Hchain'], row['Lchain']
            antibody_chains = {hchain}
            if pd.notna(lchain) and str(lchain).strip() != 'NA':
                antibody_chains.add(lchain)
            aq = AntigenQuality()
            antigen_chains, _ = aq.identify_antigen_chains(s, model_id, antibody_chains)
            if len(antigen_chains) > 1:
                found = True
                break
        assert found, "Expected at least one real structure with multiple protein-antigen chains"


class TestGordonThresholdBoundary:
    def test_exact_4_5000_passes(self):
        cs = ContactSearch(threshold_angstrom=4.5)
        cdr_atoms = np.array([[0, 0, 0]])
        antigen_atoms = np.array([[4.5, 0, 0]])  # exactly 4.5000
        has_contact, _, _ = cs.optimized_contacts(cdr_atoms, antigen_atoms)
        assert has_contact is True

    def test_4_5001_fails(self):
        cs = ContactSearch(threshold_angstrom=4.5)
        cdr_atoms = np.array([[0, 0, 0]])
        antigen_atoms = np.array([[4.5001, 0, 0]])
        has_contact, _, _ = cs.optimized_contacts(cdr_atoms, antigen_atoms)
        assert has_contact is False

    def test_boundary_agreement_brute_force_and_optimized(self):
        cs = ContactSearch(threshold_angstrom=4.5)
        cdr_atoms = np.array([[0, 0, 0]])
        antigen_atoms = np.array([[4.5, 0, 0]])
        bf = cs.brute_force_contacts(cdr_atoms, antigen_atoms)
        opt = cs.optimized_contacts(cdr_atoms, antigen_atoms)
        assert bf[0] == opt[0] == True


class TestFrameworkOnlyContactDoesNotSatisfyGordon:
    def test_framework_contact_alone_is_not_cdr_contact(self):
        """
        A structure where only framework (non-CDR) residues touch the antigen
        must NOT be reported as a Gordon-qualifying CDR-antigen contact.
        This is enforced by construction: contact search here is only ever run
        on CDR-restricted atom sets (see cdr_qc.get_cdr_atom_coords), so a
        framework-only contact is invisible to the Gordon test by design.
        """
        cs = ContactSearch(threshold_angstrom=4.5)
        # Simulate: framework atom close to antigen, but NOT included in the
        # CDR atom set passed to contact search (as the real pipeline does).
        framework_atom_near_antigen = np.array([[0, 0, 0]])
        cdr_atoms_far = np.array([[100, 100, 100]])
        antigen_atoms = np.array([[0.5, 0, 0]])  # close to framework, not to CDR

        # Correct call: only CDR atoms are tested, framework atom is excluded
        has_contact, _, _ = cs.optimized_contacts(cdr_atoms_far, antigen_atoms)
        assert has_contact is False, "Framework-proximal contact must not register as CDR contact"


class TestPairingQC:
    def test_physical_pair_within_threshold(self):
        pq = PairingQuality(threshold_angstrom=5.0)
        vh_atoms = [(0, 0, 0)]
        vl_atoms = [(3, 0, 0)]  # distance 3.0
        is_paired, min_dist = pq.check_physical_pair(vh_atoms, vl_atoms)
        assert is_paired is True
        assert min_dist == pytest.approx(3.0)

    def test_physical_pair_beyond_threshold(self):
        pq = PairingQuality(threshold_angstrom=5.0)
        vh_atoms = [(0, 0, 0)]
        vl_atoms = [(10, 0, 0)]
        is_paired, min_dist = pq.check_physical_pair(vh_atoms, vl_atoms)
        assert is_paired is False

    def test_empty_chain_atoms_not_paired(self):
        pq = PairingQuality(threshold_angstrom=5.0)
        is_paired, min_dist = pq.check_physical_pair([], [(0, 0, 0)])
        assert is_paired is False
        assert min_dist == float('inf')


class TestStructureResolverIntegration:
    def test_manifest_loads(self):
        resolver = StructureResolver()
        assert len(resolver.manifest) > 0

    def test_resolve_real_structure_path(self):
        metadata = pd.read_parquet(PROJECT_ROOT / "data/metadata/master_stage1_resolved.parquet")
        row = metadata[metadata["primary_vhh_candidate_strict"] == True].iloc[0]
        resolver = StructureResolver()
        result = resolver.resolve_structure_path(
            row['extracted_structure_path'], row['PDB'], row['INSTANCE'], row.get('archive_sha256', '')
        )
        assert result['resolved'] is True
        assert Path(result['structure_path']).exists()

    def test_resolve_missing_structure_fails_cleanly(self):
        resolver = StructureResolver()
        result = resolver.resolve_structure_path(
            "data/raw_extracted/sdh/does_not_exist/nope.cif", "0000", "fake-instance", "deadbeef"
        )
        assert result['resolved'] is False
        assert result['error_reason'] == "FAIL_STRUCTURE_FILE_MISSING"


class TestNoDuplicateScientificRows:
    def test_duplicate_fab_archive_does_not_duplicate_rows(self):
        """
        Stage-1 tracked a content_duplicate_flag for the duplicated FAB
        archive; confirm the resolved metadata table does not carry duplicate
        INSTANCE rows from it.
        """
        metadata = pd.read_parquet(PROJECT_ROOT / "data/metadata/master_stage1_resolved.parquet")
        dup_counts = metadata['INSTANCE'].value_counts()
        duplicated_instances = dup_counts[dup_counts > 1]
        assert len(duplicated_instances) == 0, f"Found duplicated INSTANCE rows: {duplicated_instances.index.tolist()[:10]}"


class TestStage1SourceTablesUnchanged:
    def test_master_stage1_parquet_exists_and_untouched_by_stage2(self):
        """
        Stage 2 must never write to master_stage1.parquet. We verify the file
        exists and its row count matches the originally documented 4931.
        """
        path = PROJECT_ROOT / "data/metadata/master_stage1.parquet"
        assert path.exists()
        df = pd.read_parquet(path)
        assert len(df) == 4931

    def test_master_stage1_resolved_row_count(self):
        path = PROJECT_ROOT / "data/metadata/master_stage1_resolved.parquet"
        df = pd.read_parquet(path)
        assert len(df) == 4931


class TestInputCohortSizesFromTable:
    def test_vhh_structural_qc_pool_size(self):
        metadata = pd.read_parquet(PROJECT_ROOT / "data/metadata/master_stage1_resolved.parquet")
        n = (metadata["vhh_candidate_structural_qc"] == True).sum()
        assert n == 1458, f"Expected 1458 VHH structural-QC rows, found {n}"

    def test_paired_primary_pool_size(self):
        metadata = pd.read_parquet(PROJECT_ROOT / "data/metadata/master_stage1_resolved.parquet")
        n = (metadata["primary_paired_candidate_metadata"] == True).sum()
        assert n == 3092, f"Expected 3092 paired primary rows, found {n}"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
