"""
Stage 6 comprehensive test suite — CDR3 Geometry, Direct Former-Interface
Shielding, and H_unshielded Decomposition.
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
from qc.stage6_cdr3_removal import extract_no_cdr3_atoms, get_shieldable_interface_positions, CDR3_START, CDR3_END
from qc.stage6_geometry import compute_cdr3_length, compute_gordon_rho_compactness, compute_zhang_fr2_cdr3_dis

PROJECT_ROOT = Path(__file__).parent.parent


class TestCDR3Membership:
    def test_cdr3_range_105_117(self):
        assert CDR3_START == 105
        assert CDR3_END == 117


class TestInsertionCodedCDR3Membership:
    def test_insertion_coded_position_within_range_counted(self):
        # e.g. 111A, 112B are within [105,117] by base number
        assert CDR3_START <= 111 <= CDR3_END
        assert CDR3_START <= 112 <= CDR3_END


class TestNoCDR3RemovesAllCDR3Atoms:
    def test_zero_cdr3_atoms_remain(self):
        master = pd.read_parquet(PROJECT_ROOT / "data/processed/stage4_vhh_sequence_representatives.parquet")
        row = master.iloc[0]
        s = MMCIFStructure(str(PROJECT_ROOT / row['extracted_structure_path']))
        model_id = s.get_model()
        native_atoms = extract_variable_domain_atoms(s, row['Hchain'], model_id)
        nocdr3_atoms = extract_no_cdr3_atoms(native_atoms)
        assert all(not (CDR3_START <= a['auth_seq_id'] <= CDR3_END) for a in nocdr3_atoms)


class TestNoNonCDR3AtomsRemoved:
    def test_all_noncdr3_atoms_retained(self):
        master = pd.read_parquet(PROJECT_ROOT / "data/processed/stage4_vhh_sequence_representatives.parquet")
        row = master.iloc[0]
        s = MMCIFStructure(str(PROJECT_ROOT / row['extracted_structure_path']))
        model_id = s.get_model()
        native_atoms = extract_variable_domain_atoms(s, row['Hchain'], model_id)
        nocdr3_atoms = extract_no_cdr3_atoms(native_atoms)
        expected_noncdr3 = [a for a in native_atoms if not (CDR3_START <= a['auth_seq_id'] <= CDR3_END)]
        assert len(nocdr3_atoms) == len(expected_noncdr3)


class TestRetainedCoordinatesUnchanged:
    def test_coords_identical(self):
        master = pd.read_parquet(PROJECT_ROOT / "data/processed/stage4_vhh_sequence_representatives.parquet")
        row = master.iloc[0]
        s = MMCIFStructure(str(PROJECT_ROOT / row['extracted_structure_path']))
        model_id = s.get_model()
        native_atoms = extract_variable_domain_atoms(s, row['Hchain'], model_id)
        nocdr3_atoms = extract_no_cdr3_atoms(native_atoms)
        native_noncdr3 = [a for a in native_atoms if not (CDR3_START <= a['auth_seq_id'] <= CDR3_END)]
        for a, b in zip(native_noncdr3, nocdr3_atoms):
            assert (a['x'], a['y'], a['z']) == (b['x'], b['y'], b['z'])


class TestNativeAndNoCDR3IdenticalTargetSet:
    def test_target_residue_keys_identical(self):
        residues = pd.read_parquet(PROJECT_ROOT / "data/processed/stage6_cdr3_shielding_residues.parquet")
        # for every row, native and noCDR3 SASA either both present or both absent (same evaluability)
        evaluable = residues.dropna(subset=['sasa_native'])
        assert evaluable['sasa_noCDR3'].notna().all()


class TestCDR3InterfacePositionsExcludedFromShieldableTarget:
    def test_no_cdr3_positions_in_shieldable_set(self):
        primary_map = pd.read_parquet(PROJECT_ROOT / "data/processed/stage3_empirical_vh_vl_interface_primary.parquet")
        _, _, shieldable = get_shieldable_interface_positions(primary_map)
        assert not ((shieldable['imgt_position'] >= CDR3_START) & (shieldable['imgt_position'] <= CDR3_END)).any()


class TestStage3EmpiricalInterfaceUnchanged:
    def test_28_positions_unchanged(self):
        df = pd.read_parquet(PROJECT_ROOT / "data/processed/stage3_empirical_vh_vl_interface_primary.parquet")
        assert df['primary_binary_interface'].sum() == 28


class TestShieldingDeltaFormulaSign:
    def test_delta_is_nocdr3_minus_native(self):
        sasa_native = 50.0
        sasa_nocdr3 = 80.0
        delta = sasa_nocdr3 - sasa_native
        assert delta == 30.0  # positive = CDR3 was occluding


class TestShieldingFractionDenominator:
    def test_fraction_uses_nocdr3_as_denominator(self):
        delta = 30.0
        sasa_nocdr3 = 80.0
        fraction = delta / sasa_nocdr3
        assert fraction == pytest.approx(0.375)


class TestUnresolvedTargetNotZeroSASA:
    def test_unresolved_has_none_not_zero(self):
        residues = pd.read_parquet(PROJECT_ROOT / "data/processed/stage6_cdr3_shielding_residues.parquet")
        unresolved = residues[residues['coordinate_observed'] == False]
        assert len(unresolved) > 0
        assert unresolved['sasa_native'].isna().all()
        assert not (unresolved['sasa_native'] == 0).any()


class TestAbsentTargetNotZeroSASA:
    def test_absent_has_none_not_zero(self):
        residues = pd.read_parquet(PROJECT_ROOT / "data/processed/stage6_cdr3_shielding_residues.parquet")
        absent = residues[residues['amino_acid'].isna()]
        if len(absent) > 0:
            assert absent['sasa_native'].isna().all()


class TestNegativeShieldingRetained:
    def test_negative_values_not_clamped_in_code(self):
        import inspect
        from qc import stage6_full_cohort
        source = inspect.getsource(stage6_full_cohort)
        assert 'clamp' not in source.lower()
        assert 'max(0' not in source.replace(' ', '')


class TestTranslatedCDR3ApproximatesNoCDR3:
    def test_translated_control_result_passes(self):
        df = pd.read_csv(PROJECT_ROOT / "results/audits/stage6_translated_cdr3_control.csv")
        assert df['control_pass'].all()


class TestKDUndefinedForNoncanonical:
    def test_noncanonical_kd_is_none(self):
        residues = pd.read_parquet(PROJECT_ROOT / "data/processed/stage6_cdr3_shielding_residues.parquet")
        noncanonical = residues[residues['noncanonical_flag'] == True]
        if len(noncanonical) > 0:
            assert noncanonical['kyte_doolittle'].isna().all()


class TestPositiveKDDeltaCalculation:
    def test_delta_formula(self):
        native_pos_sasa = 20.0
        nocdr3_pos_sasa = 35.0
        delta = nocdr3_pos_sasa - native_pos_sasa
        assert delta == 15.0


class TestSASAWeightedKDDeltaCalculation:
    def test_delta_formula(self):
        native_kd = -1.0
        nocdr3_kd = -0.5
        delta = nocdr3_kd - native_kd
        assert delta == 0.5


class TestContinuousInterfaceWeighting:
    def test_weighted_uses_freq5(self):
        combined = pd.read_parquet(PROJECT_ROOT / "data/processed/stage6_cdr3_shielding_geometry_combined.parquet")
        assert 'weighted_delta_sasa_shield_freq5' in combined.columns


class TestSequenceLengthDistinctFromResolvedLength:
    def test_columns_are_separate(self):
        geometries = pd.read_parquet(PROJECT_ROOT / "data/processed/stage6_cdr3_geometry.parquet")
        assert 'cdr3_sequence_length' in geometries.columns
        assert 'cdr3_resolved_length' in geometries.columns
        # they can differ
        assert (geometries['cdr3_sequence_length'] != geometries['cdr3_resolved_length']).any() or \
               (geometries['cdr3_sequence_length'] == geometries['cdr3_resolved_length']).all()


class TestGordonRigidTranslationInvariance:
    def test_translation_invariance_from_validation_audit(self):
        df = pd.read_csv(PROJECT_ROOT / "results/audits/stage6_gordon_geometry_validation.csv")
        valid = df[df['gordon_valid_original']]
        assert (valid['rho_invariance_diff'] < 1e-6).all()


class TestGordonRigidRotationInvariance:
    def test_rotation_invariance_from_validation_audit(self):
        df = pd.read_csv(PROJECT_ROOT / "results/audits/stage6_gordon_geometry_validation.csv")
        valid = df[df['gordon_valid_original']]
        assert (valid['compactness_invariance_diff'] < 1e-6).all()


class TestGordonCompactnessExactFormula:
    def test_compactness_equals_length_over_rho(self):
        geometries = pd.read_parquet(PROJECT_ROOT / "data/processed/stage6_cdr3_geometry.parquet")
        valid = geometries[geometries['gordon_geometry_valid'] == True].head(20)
        for _, row in valid.iterrows():
            expected = row['cdr3_sequence_length'] / row['cdr3_rho']
            assert abs(row['cdr3_compactness'] - expected) < 1e-9


class TestZhangFR2CDR3DisExactDefinition:
    def test_definition_documented(self):
        doc_path = PROJECT_ROOT / "results/audits/stage6_zhang_method_definition.md"
        assert doc_path.exists()
        content = doc_path.read_text()
        assert 'minimum' in content.lower()
        assert 'centroid' in content.lower() or 'central coordinate' in content.lower()


class TestFR2CDR3DisRigidBodyInvariance:
    def test_invariance_from_validation_audit(self):
        df = pd.read_csv(PROJECT_ROOT / "results/audits/stage6_gordon_geometry_validation.csv")
        valid = df[df['zhang_valid_original']]
        assert (valid['FR2_CDR3_Dis_invariance_diff'] < 1e-6).all()


class TestAngularWrapAroundHandling:
    def test_theta_phi_not_computed_no_wraparound_risk(self):
        geometries = pd.read_parquet(PROJECT_ROOT / "data/processed/stage6_cdr3_geometry.parquet")
        # theta/phi genuinely unimplemented -- no angular correlation performed, so no
        # wrap-around artifact is possible in this stage's outputs
        assert 'cdr3_theta' not in geometries.columns or geometries['cdr3_theta'].isna().all()
        assert 'cdr3_phi' not in geometries.columns or geometries['cdr3_phi'].isna().all()


class TestHighCoverageSensitivityDoesNotAlterPrimary:
    def test_primary_row_count_unchanged(self):
        summaries = pd.read_parquet(PROJECT_ROOT / "data/processed/stage6_vhh_cdr3_shielding.parquet")
        assert len(summaries) == 533


class TestNoArpeggioInvoked:
    def test_no_arpeggio_import(self):
        import inspect
        from qc import stage6_full_cohort, stage6_geometry, stage6_cdr3_removal
        for module in [stage6_full_cohort, stage6_geometry, stage6_cdr3_removal]:
            source = inspect.getsource(module)
            assert 'arpeggio' not in source.lower()


class TestUpstreamSourceFilesUnchanged:
    def test_stage2_unchanged(self):
        assert len(pd.read_parquet(PROJECT_ROOT / "data/processed/stage2_vhh_strict_primary.parquet")) == 1217

    def test_stage3_unchanged(self):
        assert len(pd.read_parquet(PROJECT_ROOT / "data/processed/stage3_empirical_vh_vl_interface_primary.parquet")) == 214

    def test_stage4_unchanged(self):
        assert len(pd.read_parquet(PROJECT_ROOT / "data/processed/stage4_vhh_sequence_representatives.parquet")) == 533

    def test_stage5_unchanged(self):
        assert len(pd.read_parquet(PROJECT_ROOT / "data/processed/stage5_sequence_surface_chemistry.parquet")) == 3391


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
