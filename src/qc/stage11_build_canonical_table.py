"""Stage 11 — Canonical per-representative model input table."""
import pandas as pd
from pathlib import Path

ROOT = Path("/Users/betulgul/Desktop/staj")


def load_vhh_meta():
    df = pd.read_parquet(ROOT / "data/processed/stage4_vhh_sequence_representatives.parquet")
    return df[['INSTANCE', 'PDB', 'SABDAB_ID', 'resolution', 'antigen_name']].assign(group='VHH')


def load_vh_meta():
    df = pd.read_parquet(ROOT / "data/processed/stage3_paired_representatives.parquet")
    df = df[df['is_primary_representative'] == True]
    return df[['INSTANCE', 'PDB', 'SABDAB_ID', 'resolution', 'antigen_name']].assign(group='VH')


def main():
    meta = pd.concat([load_vhh_meta(), load_vh_meta()], ignore_index=True)
    assert meta['INSTANCE'].is_unique, "duplicate INSTANCE in base metadata"
    n0 = len(meta)

    # Stage 4: H_intrinsic (per representative_INSTANCE)
    s4 = pd.read_csv(ROOT / "data/processed/stage4_sequence_interface_chemistry.csv")
    s4 = s4.rename(columns={'representative_INSTANCE': 'INSTANCE'})
    s4_cols = ['INSTANCE', 'H_intrinsic_binary_mean_KD', 'primary_interface_position_coverage']
    df = meta.merge(s4[s4_cols], on='INSTANCE', how='left')
    assert len(df) == n0, "Stage-4 join altered row count"

    # Stage 6: shielding + geometry (VHH_NATIVE and VH_WITHOUT_VL_GEOMETRIC both carry
    # H_native_full_interface_sasa_weighted_mean_KD; shielding/geometry VHH-only)
    s6 = pd.read_parquet(ROOT / "data/processed/stage6_cdr3_shielding_geometry_combined.parquet")
    s6_cols = ['INSTANCE', 'H_native_full_interface_sasa_weighted_mean_KD', 'shielding_fraction',
               'delta_sasa_shield', 'cdr3_sequence_length', 'cdr3_rho', 'cdr3_compactness', 'FR2_CDR3_Dis']
    df = df.merge(s6[s6_cols], on='INSTANCE', how='left')
    assert len(df) == n0, "Stage-6 join altered row count"

    # Need VH-side H_native too -- Stage 6 combined table is VHH-only (533 rows).
    # Reload Stage 5 directly, filtered to VH_WITHOUT_VL_GEOMETRIC, for the VH rows.
    s5 = pd.read_parquet(ROOT / "data/processed/stage5_sequence_surface_chemistry.parquet")
    s5_vh = s5[(s5['group'] == 'VH') & (s5['structural_state'] == 'VH_WITHOUT_VL_GEOMETRIC')]
    s5_vh = s5_vh[['INSTANCE', 'interface_sasa_weighted_mean_KD']].rename(
        columns={'interface_sasa_weighted_mean_KD': 'H_native_full_interface_sasa_weighted_mean_KD_vh'})
    df = df.merge(s5_vh, on='INSTANCE', how='left')
    # fill VH rows' H_native column from the Stage-5-direct value (Stage-6 combined has none for VH)
    df['H_native_full_interface_sasa_weighted_mean_KD'] = df['H_native_full_interface_sasa_weighted_mean_KD'].fillna(
        df['H_native_full_interface_sasa_weighted_mean_KD_vh'])
    df = df.drop(columns=['H_native_full_interface_sasa_weighted_mean_KD_vh'])
    assert len(df) == n0, "Stage-5 join altered row count"

    # Stage 7: BSA (VHH only)
    s7 = pd.read_csv(ROOT / "results/tables/stage7_per_vhh_packing_summary.csv")
    df = df.merge(s7[['INSTANCE', 'nonadj_BSA_pair_total']], on='INSTANCE', how='left')
    assert len(df) == n0, "Stage-7 join altered row count"

    # Stage 8: patch + topology + APBS
    s8p = pd.read_parquet(ROOT / "data/processed/stage8_surface_patch_descriptors.parquet")
    s8p_cols = ['INSTANCE', 'former_interface_core_PSH', 'former_interface_core_PPC', 'former_interface_core_PNC']
    df = df.merge(s8p[s8p_cols], on='INSTANCE', how='left')
    assert len(df) == n0, "Stage-8 patch join altered row count"

    s8t = pd.read_parquet(ROOT / "data/processed/stage8_patch_topology.parquet")
    df = df.merge(s8t[['INSTANCE', 'local_thr4.0_largest_patch_residue_count']], on='INSTANCE', how='left')
    assert len(df) == n0, "Stage-8 topology join altered row count"

    s8a = pd.read_parquet(ROOT / "data/processed/stage8_apbs_interface_electrostatics.parquet")
    s8a_cols = ['INSTANCE', 'interface_potential_positive_mean', 'interface_potential_negative_mean', 'interface_potential_absolute_mean']
    df = df.merge(s8a[s8a_cols], on='INSTANCE', how='left')
    assert len(df) == n0, "Stage-8 APBS join altered row count"

    # Stage 9: engagement summary + 15-position and 28-position engagement fractions
    s9 = pd.read_parquet(ROOT / "data/processed/stage9_antigen_engagement_summary.parquet")
    s9_cols = ['INSTANCE', 'framework_contact_fraction', 'cdr3_contact_fraction',
               'total_unique_antibody_contact_residues', 'total_unique_antigen_contact_residues']
    df = df.merge(s9[s9_cols], on='INSTANCE', how='left')
    assert len(df) == n0, "Stage-9 summary join altered row count"

    vhh_shield = pd.read_csv(ROOT / "data/interim/stage9/vhh_shieldable_interface_engagement.csv")
    vh_shield = pd.read_csv(ROOT / "data/interim/stage9/vh_shieldable_interface_engagement.csv")
    shield_eng = pd.concat([vhh_shield, vh_shield], ignore_index=True).rename(
        columns={'fraction_target_contacting': 'shieldable_engagement_fraction',
                 'n_target_positions_contacting': 'shieldable_n_contacting',
                 'n_target_positions_total': 'shieldable_n_total'})
    df = df.merge(shield_eng[['INSTANCE', 'shieldable_engagement_fraction', 'shieldable_n_contacting', 'shieldable_n_total']],
                  on='INSTANCE', how='left')
    assert len(df) == n0, "Stage-9 shieldable-engagement join altered row count"

    vhh_full = pd.read_csv(ROOT / "data/interim/stage9/vhh_full_interface_engagement.csv")
    vh_full = pd.read_csv(ROOT / "data/interim/stage9/vh_full_interface_engagement.csv")
    full_eng = pd.concat([vhh_full, vh_full], ignore_index=True).rename(
        columns={'fraction_target_contacting': 'full_interface_engagement_fraction',
                 'n_target_positions_contacting': 'full_n_contacting', 'n_target_positions_total': 'full_n_total'})
    df = df.merge(full_eng[['INSTANCE', 'full_interface_engagement_fraction', 'full_n_contacting', 'full_n_total']],
                  on='INSTANCE', how='left')
    assert len(df) == n0, "Stage-9 full-engagement join altered row count"

    assert df['INSTANCE'].is_unique, "duplicate INSTANCE after all joins"
    assert len(df) == n0 == 1962, f"final row count {len(df)} != expected 1962"

    df.to_parquet(ROOT / "data/processed/stage11_model_input.parquet", index=False)
    df.to_csv(ROOT / "data/processed/stage11_model_input.csv", index=False)
    print(f"Canonical table: {len(df)} rows, {len(df.columns)} columns")
    print(df.groupby('group').size())
    return df


if __name__ == '__main__':
    main()
