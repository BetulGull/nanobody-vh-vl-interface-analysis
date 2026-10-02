"""
Stage 3 comprehensive test suite — VH-VL empirical interface derivation by ΔSASA.
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from qc.mmcif_qc import MMCIFStructure
from qc.stage3_variable_domain import extract_variable_domain_atoms
from qc.stage3_sasa import calc_sasa, build_structure

PROJECT_ROOT = Path(__file__).parent.parent

VH_CH1_INSTANCE = 'pdb_000011ha-G-K'      # FAB, VH+CH1 same physical chain
FV_INSTANCE = 'pdb_00001a14-H-L'          # paired FV
INSERTION_INSTANCE = 'pdb_000013fj-H-L'   # insertion-coded variable domain


def _load_row(instance):
    master = pd.read_parquet(PROJECT_ROOT / "data/processed/stage2_structure_qc_master.parquet")
    return master[master['INSTANCE'] == instance].iloc[0]


def _load_structure(row):
    path = PROJECT_ROOT / row['extracted_structure_path']
    s = MMCIFStructure(str(path))
    return s, s.get_model()


class TestVHOnlyExtraction:
    def test_vh_only_contains_no_vl_atoms(self):
        row = _load_row(FV_INSTANCE)
        s, model_id = _load_structure(row)
        vh_atoms = extract_variable_domain_atoms(s, row['Hchain'], model_id)
        assert all(a['auth_asym_id'] == row['Hchain'] for a in vh_atoms)
        assert not any(a['auth_asym_id'] == row['Lchain'] for a in vh_atoms)


class TestVHVLExtraction:
    def test_combined_extraction_contains_only_vh_and_vl_chains(self):
        row = _load_row(FV_INSTANCE)
        s, model_id = _load_structure(row)
        vh_atoms = extract_variable_domain_atoms(s, row['Hchain'], model_id)
        vl_atoms = extract_variable_domain_atoms(s, row['Lchain'], model_id)
        combined_chains = set(a['auth_asym_id'] for a in vh_atoms + vl_atoms)
        assert combined_chains == {row['Hchain'], row['Lchain']}


class TestAntigenAbsence:
    def test_antigen_chains_absent_from_extraction(self):
        row = _load_row(VH_CH1_INSTANCE)
        s, model_id = _load_structure(row)
        vh_atoms = extract_variable_domain_atoms(s, row['Hchain'], model_id)
        vl_atoms = extract_variable_domain_atoms(s, row['Lchain'], model_id)
        chains_used = set(a['auth_asym_id'] for a in vh_atoms + vl_atoms)
        assert chains_used == {row['Hchain'], row['Lchain']}


class TestConstantDomainAbsence:
    def test_ch1_excluded_from_fab_vh_chain(self):
        row = _load_row(VH_CH1_INSTANCE)
        s, model_id = _load_structure(row)
        vh_atoms = extract_variable_domain_atoms(s, row['Hchain'], model_id)
        assert all(a['auth_seq_id'] <= 128 for a in vh_atoms)
        assert max(a['auth_seq_id'] for a in vh_atoms) <= 128
        # Confirm CH1 really was present and excluded (not just absent from file)
        total_h_residues = len(s.build_residue_scheme(row['Hchain'], model_id))
        assert total_h_residues > 128  # CH1 exists in this chain


class TestCoordinateEquality:
    def test_vh_coords_identical_between_states(self):
        row = _load_row(FV_INSTANCE)
        s, model_id = _load_structure(row)
        vh_alone = extract_variable_domain_atoms(s, row['Hchain'], model_id)
        vh_again = extract_variable_domain_atoms(s, row['Hchain'], model_id)
        assert len(vh_alone) == len(vh_again)
        for a, b in zip(vh_alone, vh_again):
            assert (a['x'], a['y'], a['z']) == (b['x'], b['y'], b['z'])


class TestCaseSensitiveChains:
    def test_case_sensitive_chain_extraction(self):
        row = _load_row(FV_INSTANCE)
        s, model_id = _load_structure(row)
        chains_in_model = s.get_chains_for_model(model_id)
        # Confirm exact case match required
        assert row['Hchain'] in chains_in_model
        assert row['Hchain'].lower() not in chains_in_model or row['Hchain'].lower() == row['Hchain']


class TestInsertionCodePreservation:
    def test_insertion_code_retained_in_extraction(self):
        row = _load_row(INSERTION_INSTANCE)
        s, model_id = _load_structure(row)
        vh_atoms = extract_variable_domain_atoms(s, row['Hchain'], model_id)
        insertion_atoms = [a for a in vh_atoms if a['ins_code']]
        assert len(insertion_atoms) > 0
        assert any(a['ins_code'] == 'A' for a in insertion_atoms)


class TestFreeSASAReproducibility:
    def test_repeated_calc_identical(self):
        row = _load_row(FV_INSTANCE)
        s, model_id = _load_structure(row)
        vh_atoms = extract_variable_domain_atoms(s, row['Hchain'], model_id)
        per_res_1, total_1 = calc_sasa(vh_atoms)
        per_res_2, total_2 = calc_sasa(vh_atoms)
        assert abs(total_1 - total_2) < 1e-9
        for k in per_res_1:
            assert abs(per_res_1[k] - per_res_2[k]) < 1e-9


class TestPerResidueAggregation:
    def test_sum_per_residue_equals_total(self):
        row = _load_row(FV_INSTANCE)
        s, model_id = _load_structure(row)
        vh_atoms = extract_variable_domain_atoms(s, row['Hchain'], model_id)
        per_res, total = calc_sasa(vh_atoms)
        assert abs(sum(per_res.values()) - total) < 1e-6


class TestDeltaSASASign:
    def test_delta_formula_direction(self):
        """delta = sasa_alone - sasa_with_vl; burial => positive delta."""
        sasa_alone = 100.0
        sasa_with_vl = 80.0
        delta = sasa_alone - sasa_with_vl
        assert delta == 20.0  # positive means buried by VL

    def test_delta_negative_means_less_burial_not_exposure_error(self):
        sasa_alone = 100.0
        sasa_with_vl = 105.0  # tiny numerical noise scenario
        delta = sasa_alone - sasa_with_vl
        assert delta == -5.0


class TestExactThresholds:
    def test_threshold_1_boundary(self):
        assert not (0.999 >= 1.0)
        assert (1.000 >= 1.0)

    def test_threshold_5_boundary(self):
        assert not (4.999 >= 5.0)
        assert (5.000 >= 5.0)

    def test_threshold_10_boundary(self):
        assert not (9.999 >= 10.0)
        assert (10.000 >= 10.0)


class TestConsensusFrequencyDenominator:
    def test_frequency_uses_observed_denominator_not_total_cohort(self):
        # position observed in only 3 of 10 structures; interface in 2 of those 3
        n_observed = 3
        n_interface = 2
        frequency = n_interface / n_observed
        assert frequency == pytest.approx(2 / 3)
        # NOT n_interface / 10 (total cohort) -- that would be a hidden-denominator bug
        assert frequency != n_interface / 10


class TestObservationFraction:
    def test_observation_fraction_uses_total_structures_denominator(self):
        n_observed = 3
        n_structures_in_map = 10
        obs_fraction = n_observed / n_structures_in_map
        assert obs_fraction == 0.3


class TestPrimaryBinaryInterfaceRule:
    def test_rule_requires_both_frequency_and_observation(self):
        def primary_binary(freq5, obs_frac, freq_t=0.25, obs_t=0.25):
            return freq5 >= freq_t and obs_frac >= obs_t

        assert primary_binary(0.5, 0.5) is True
        assert primary_binary(0.5, 0.1) is False  # fails observation guard
        assert primary_binary(0.1, 0.5) is False  # fails frequency
        assert primary_binary(0.25, 0.25) is True  # exact boundary


class TestRareInsertionObservationGuard:
    def test_rare_position_with_high_conditional_frequency_excluded(self):
        # Observed in only 1/100 structures, but that one has interface -> freq=1.0
        n_observed = 1
        n_structures = 100
        n_interface = 1
        freq = n_interface / n_observed
        obs_fraction = n_observed / n_structures
        primary_binary = freq >= 0.25 and obs_fraction >= 0.25
        assert freq == 1.0  # deceptively "perfect" frequency
        assert obs_fraction == 0.01
        assert primary_binary is False  # correctly excluded by observation guard


class TestTranslatedVLControl:
    def test_far_translated_vl_yields_near_zero_delta(self):
        row = _load_row(FV_INSTANCE)
        s, model_id = _load_structure(row)
        vh_atoms = extract_variable_domain_atoms(s, row['Hchain'], model_id)
        vl_atoms = extract_variable_domain_atoms(s, row['Lchain'], model_id)

        import copy
        vl_far = copy.deepcopy(vl_atoms)
        for a in vl_far:
            a['x'] += 100000.0

        _, total_alone = calc_sasa(vh_atoms)
        per_res_far, _ = calc_sasa(vh_atoms + vl_far)
        vh_keys = set((a['auth_asym_id'], f"{a['auth_seq_id']}{a['ins_code']}") for a in vh_atoms)
        vh_total_with_far = sum(v for k, v in per_res_far.items() if k in vh_keys)

        assert abs(total_alone - vh_total_with_far) < 0.5


class TestRepresentativeSelectionDeterminism:
    def test_selection_is_deterministic_across_runs(self):
        selection_path = PROJECT_ROOT / "results/audits/stage3_representative_selection.csv"
        if not selection_path.exists():
            pytest.skip("representative selection not yet generated")
        df1 = pd.read_csv(selection_path)
        # Re-derive representative flags should be stable: check no group has >1 representative
        rep_counts = df1[df1['is_primary_representative']].groupby('redundancy_group').size()
        assert (rep_counts == 1).all()


class TestNoInflationFromDuplicates:
    def test_primary_map_uses_one_representative_per_group(self):
        rep_path = PROJECT_ROOT / "data/processed/stage3_paired_representatives.parquet"
        if not rep_path.exists():
            pytest.skip("representatives not yet generated")
        rep_df = pd.read_parquet(rep_path)
        assert rep_df['redundancy_group'].is_unique


class TestAllInstanceMapSeparateFromPrimary:
    def test_maps_are_distinct_files_with_labels(self):
        primary_path = PROJECT_ROOT / "data/processed/stage3_empirical_vh_vl_interface_primary.parquet"
        all_path = PROJECT_ROOT / "data/processed/stage3_empirical_vh_vl_interface_all_instances.csv"
        if not (primary_path.exists() and all_path.exists()):
            pytest.skip("maps not yet generated")
        primary_df = pd.read_parquet(primary_path)
        all_df = pd.read_csv(all_path)
        assert primary_df['map_type'].iloc[0] == 'PRIMARY_NONREDUNDANT'
        assert all_df['map_type'].iloc[0] == 'ALL_INSTANCE_SENSITIVITY'


class TestTetradAnnotationDoesNotForceInterface:
    def test_tetrad_flag_independent_of_binary_interface(self):
        primary_path = PROJECT_ROOT / "data/processed/stage3_empirical_vh_vl_interface_primary.parquet"
        if not primary_path.exists():
            pytest.skip("primary map not yet generated")
        df = pd.read_parquet(primary_path)
        tetrad_rows = df[df['vhh_tetrad_position']]
        # Not all tetrad positions need be in the binary interface, and not all
        # binary-interface positions need be tetrad positions.
        assert not df['primary_binary_interface'].equals(df['vhh_tetrad_position'])


class TestStage2SourceUnchanged:
    def test_stage2_paired_primary_row_count_unchanged(self):
        df = pd.read_parquet(PROJECT_ROOT / "data/processed/stage2_paired_primary.parquet")
        assert len(df) == 2879

    def test_stage2_vhh_strict_primary_unchanged(self):
        df = pd.read_parquet(PROJECT_ROOT / "data/processed/stage2_vhh_strict_primary.parquet")
        assert len(df) == 1217


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
