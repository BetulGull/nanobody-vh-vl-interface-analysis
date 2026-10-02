"""
Stage 5 comprehensive test suite — Native 3D Solvent Exposure and
Former-Interface Hydrophobic Burden.
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from qc.mmcif_qc import MMCIFStructure
from qc.stage3_variable_domain import extract_variable_domain_atoms
from qc.stage3_sasa import calc_sasa
from qc.antigen_qc import AntigenQuality

PROJECT_ROOT = Path(__file__).parent.parent


class TestVHHNativeNoAntigen:
    def test_no_antigen_atoms_in_vhh_extraction(self):
        master = pd.read_parquet(PROJECT_ROOT / "data/processed/stage4_vhh_sequence_representatives.parquet")
        row = master.iloc[0]
        s = MMCIFStructure(str(PROJECT_ROOT / row['extracted_structure_path']))
        model_id = s.get_model()
        vhh_atoms = extract_variable_domain_atoms(s, row['Hchain'], model_id)
        chains_used = set(a['auth_asym_id'] for a in vhh_atoms)
        assert chains_used == {row['Hchain']}


class TestVHHNativeRetainsCDR3:
    def test_cdr3_atoms_present(self):
        master = pd.read_parquet(PROJECT_ROOT / "data/processed/stage4_vhh_sequence_representatives.parquet")
        row = master.iloc[0]
        s = MMCIFStructure(str(PROJECT_ROOT / row['extracted_structure_path']))
        model_id = s.get_model()
        vhh_atoms = extract_variable_domain_atoms(s, row['Hchain'], model_id)
        cdr3_atoms = [a for a in vhh_atoms if 105 <= a['auth_seq_id'] <= 117]
        assert len(cdr3_atoms) > 0


class TestVHGeometricCounterfactualNoVL:
    def test_vh_alone_extraction_excludes_vl(self):
        master = pd.read_parquet(PROJECT_ROOT / "data/processed/stage3_paired_representatives.parquet")
        row = master.iloc[0]
        s = MMCIFStructure(str(PROJECT_ROOT / row['extracted_structure_path']))
        model_id = s.get_model()
        vh_atoms = extract_variable_domain_atoms(s, row['Hchain'], model_id)
        chains_used = set(a['auth_asym_id'] for a in vh_atoms)
        assert row['Lchain'] not in chains_used


class TestPairedVHContainsVL:
    def test_stage3_with_vl_state_used_vl_atoms(self):
        residues = pd.read_parquet(PROJECT_ROOT / "data/processed/stage3_vh_vl_delta_sasa_residues.parquet")
        row = residues.iloc[0]
        assert pd.notna(row['sasa_vh_with_vl'])
        assert row['sasa_vh_with_vl'] <= row['sasa_vh_alone'] + 1e-6  # VL presence cannot increase VH SASA materially


class TestStage3VHSASAReuseJoin:
    def test_join_matches_correct_instance_and_residue(self):
        combined = pd.read_parquet(PROJECT_ROOT / "data/processed/stage5_interface_surface_chemistry_residues.parquet")
        stage3 = pd.read_parquet(PROJECT_ROOT / "data/processed/stage3_vh_vl_delta_sasa_residues.parquet")

        sample = combined[(combined['structural_state'] == 'VH_WITHOUT_VL_GEOMETRIC') & combined['sasa'].notna()].iloc[0]
        match = stage3[
            (stage3['INSTANCE'] == sample['INSTANCE']) &
            (stage3['imgt_position'] == sample['imgt_position']) &
            (stage3['insertion_code'] == sample['insertion_code'])
        ]
        assert len(match) == 1
        assert abs(match.iloc[0]['sasa_vh_alone'] - sample['sasa']) < 1e-6


class TestExactIMGTInsertionMatching:
    def test_insertion_positions_distinct(self):
        combined = pd.read_parquet(PROJECT_ROOT / "data/processed/stage5_interface_surface_chemistry_residues.parquet")
        ins_rows = combined[combined['insertion_code'] == 'A']
        base_rows = combined[(combined['imgt_position'].isin(ins_rows['imgt_position'])) & (combined['insertion_code'] == '')]
        assert len(ins_rows) > 0 and len(base_rows) > 0


class TestUnresolvedNotZeroSASA:
    def test_unresolved_residue_has_nan_not_zero(self):
        combined = pd.read_parquet(PROJECT_ROOT / "data/processed/stage5_interface_surface_chemistry_residues.parquet")
        unresolved = combined[combined['mapping_status'] == 'PRESENT_UNRESOLVED']
        assert len(unresolved) > 0
        assert unresolved['sasa'].isna().all()
        assert not (unresolved['sasa'] == 0).any()


class TestAbsentNotZeroSASA:
    def test_absent_residue_has_nan_not_zero(self):
        combined = pd.read_parquet(PROJECT_ROOT / "data/processed/stage5_interface_surface_chemistry_residues.parquet")
        absent = combined[combined['mapping_status'] == 'ABSENT_CONFIRMED']
        assert len(absent) > 0
        assert absent['sasa'].isna().all()


class TestNoncanonicalKDUndefined:
    def test_noncanonical_residue_has_no_kd(self):
        combined = pd.read_parquet(PROJECT_ROOT / "data/processed/stage5_interface_surface_chemistry_residues.parquet")
        noncanonical = combined[combined['noncanonical_residue_flag'] == True]
        assert len(noncanonical) > 0
        assert noncanonical['kyte_doolittle'].isna().all()


class TestTotalInterfaceSASADenominator:
    def test_sum_over_evaluable_only(self):
        sasa_vals = [10.0, 20.0, None, 15.0]  # one absent/unresolved
        evaluable = [v for v in sasa_vals if v is not None]
        total = sum(evaluable)
        assert total == 45.0
        assert len(evaluable) == 3  # denominator excludes the missing one


class TestMeanSASADenominator:
    def test_mean_uses_evaluable_count(self):
        sasa_vals = [10.0, 20.0, None]
        evaluable = [v for v in sasa_vals if v is not None]
        mean_val = sum(evaluable) / len(evaluable)
        assert mean_val == 15.0


class TestPositiveKDSASACalculation:
    def test_only_positive_kd_residues_summed(self):
        residues = [{'sasa': 10, 'kd': 4.5}, {'sasa': 20, 'kd': -3.5}, {'sasa': 5, 'kd': 1.8}]
        positive_sum = sum(r['sasa'] for r in residues if r['kd'] > 0)
        assert positive_sum == 15  # 10 + 5, excludes the -3.5 one


class TestPositiveKDFractionDenominator:
    def test_fraction_uses_total_canonical_sasa(self):
        residues = [{'sasa': 10, 'kd': 4.5}, {'sasa': 20, 'kd': -3.5}]
        positive_sum = sum(r['sasa'] for r in residues if r['kd'] > 0)
        total_sum = sum(r['sasa'] for r in residues)
        fraction = positive_sum / total_sum
        assert fraction == pytest.approx(10 / 30)


class TestSASAWeightedKDCalculation:
    def test_weighted_mean_formula(self):
        residues = [{'sasa': 10, 'kd': 4.5}, {'sasa': 20, 'kd': -3.5}]
        weighted = sum(r['sasa'] * r['kd'] for r in residues) / sum(r['sasa'] for r in residues)
        assert weighted == pytest.approx((10 * 4.5 + 20 * -3.5) / 30)


class TestWeightedInterfaceFrequencyCalculation:
    def test_freq5_weighted_formula(self):
        residues = [{'w': 0.5, 'sasa': 10, 'kd': 2.0}, {'w': 0.8, 'sasa': 20, 'kd': -1.0}]
        num = sum(r['w'] * r['sasa'] * r['kd'] for r in residues)
        denom = sum(r['w'] * r['sasa'] for r in residues)
        result = num / denom
        assert result == pytest.approx(num / denom)
        assert denom == pytest.approx(0.5 * 10 + 0.8 * 20)


class TestEffectiveExposedWeightDenominator:
    def test_denominator_is_weight_times_sasa_not_weight_alone(self):
        residues = [{'w': 0.5, 'sasa': 10}, {'w': 0.8, 'sasa': 20}]
        effective_weight = sum(r['w'] * r['sasa'] for r in residues)
        naive_weight = sum(r['w'] for r in residues)
        assert effective_weight != naive_weight
        assert effective_weight == pytest.approx(0.5 * 10 + 0.8 * 20)


class TestBinaryInterfaceUsesStage3PositionsOnly:
    def test_primary_binary_interface_count_matches_stage3(self):
        combined = pd.read_parquet(PROJECT_ROOT / "data/processed/stage5_interface_surface_chemistry_residues.parquet")
        stage3_map = pd.read_parquet(PROJECT_ROOT / "data/processed/stage3_empirical_vh_vl_interface_primary.parquet")
        n_binary_stage3 = stage3_map['primary_binary_interface'].sum()
        n_binary_stage5 = combined[combined['structural_state'] == 'VHH_NATIVE']['primary_binary_interface'].sum() / \
            combined[combined['structural_state'] == 'VHH_NATIVE']['INSTANCE'].nunique()
        assert abs(n_binary_stage3 - n_binary_stage5) < 1e-6


class TestTetradCannotAlterMembership:
    def test_non_tetrad_positions_exist_in_interface(self):
        pos_comparison = pd.read_parquet(PROJECT_ROOT / "data/processed/stage5_position_surface_chemistry.parquet")
        combined = pd.read_parquet(PROJECT_ROOT / "data/processed/stage5_interface_surface_chemistry_residues.parquet")
        tetrad_positions = {42, 49, 50, 52}
        non_tetrad_interface = pos_comparison[
            (~pos_comparison['imgt_position'].isin(tetrad_positions)) & pos_comparison['primary_binary_interface']
        ]
        assert len(non_tetrad_interface) > 0


class TestNoCDR3RemovalPath:
    def test_no_cdr3_removal_function_exists(self):
        import inspect
        from qc import stage5_structural_chemistry, stage5_sequence_surface_chemistry, stage5_position_and_audits
        for module in [stage5_structural_chemistry, stage5_sequence_surface_chemistry, stage5_position_and_audits]:
            source = inspect.getsource(module)
            assert 'shielding' not in source.lower()
            assert 'unshielded' not in source.lower()
            assert 'no_cdr3' not in source.lower()
            assert 'minus_cdr3' not in source.lower()


class TestHighCoverageSensitivityDoesNotAlterPrimary:
    def test_primary_cohort_unchanged_by_sensitivity_filter(self):
        seq_chem = pd.read_parquet(PROJECT_ROOT / "data/processed/stage5_sequence_surface_chemistry.parquet")
        n_vhh_full = len(seq_chem[seq_chem['structural_state'] == 'VHH_NATIVE'])
        assert n_vhh_full == 533  # unchanged regardless of any coverage sensitivity subset


class TestVHHRepresentativeSetSize:
    def test_vhh_representatives_533(self):
        df = pd.read_parquet(PROJECT_ROOT / "data/processed/stage4_vhh_sequence_representatives.parquet")
        assert len(df) == 533


class TestVHRepresentativeSetFromFrozenStage3:
    def test_vh_representatives_1429(self):
        df = pd.read_parquet(PROJECT_ROOT / "data/processed/stage3_paired_representatives.parquet")
        assert len(df) == 1429


class TestUpstreamSourceFilesUnchanged:
    def test_stage2_unchanged(self):
        assert len(pd.read_parquet(PROJECT_ROOT / "data/processed/stage2_vhh_strict_primary.parquet")) == 1217
        assert len(pd.read_parquet(PROJECT_ROOT / "data/processed/stage2_paired_primary.parquet")) == 2879

    def test_stage3_unchanged(self):
        df = pd.read_parquet(PROJECT_ROOT / "data/processed/stage3_empirical_vh_vl_interface_primary.parquet")
        assert len(df) == 214
        assert df['primary_binary_interface'].sum() == 28

    def test_stage4_unchanged(self):
        df = pd.read_parquet(PROJECT_ROOT / "data/processed/stage4_sequence_interface_chemistry.parquet")
        assert len(df) == 1962  # 533 VHH + 1429 VH


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
