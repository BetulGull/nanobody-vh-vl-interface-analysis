"""Stage 11 — Authoritative Input Variable Dictionary."""
import pandas as pd
from pathlib import Path

ROOT = Path("/Users/betulgul/Desktop/staj")

VARS = [
    dict(variable='H_intrinsic_binary_mean_KD', source_stage=4, source_file='data/processed/stage4_sequence_interface_chemistry.csv',
         exact_column='H_intrinsic_binary_mean_KD', group_applicability='VHH,VH', definition='Intrinsic (sequence-only) KD hydrophobicity mean over the frozen 28-position former-interface, binary-membership weighted', units='KD units', primary_or_secondary='primary'),
    dict(variable='H_native_full_interface_sasa_weighted_mean_KD', source_stage=6, source_file='data/processed/stage6_cdr3_shielding_geometry_combined.parquet',
         exact_column='H_native_full_interface_sasa_weighted_mean_KD', group_applicability='VHH,VH', definition='Native structural exposure: SASA-weighted mean KD over the former interface. VHH=VHH_NATIVE state; VH=VH_WITHOUT_VL_GEOMETRIC state (renamed from stage5 interface_sasa_weighted_mean_KD)', units='KD units', primary_or_secondary='primary'),
    dict(variable='shielding_fraction', source_stage=6, source_file='data/processed/stage6_cdr3_shielding_geometry_combined.parquet',
         exact_column='shielding_fraction', group_applicability='VHH only', definition='Fraction of native SASA at the 15-position shieldable interface attributable to CDR3 (native vs no-CDR3 geometric counterfactual)', units='fraction [0,1]', primary_or_secondary='primary'),
    dict(variable='delta_sasa_shield', source_stage=6, source_file='data/processed/stage6_cdr3_shielding_geometry_combined.parquet',
         exact_column='delta_sasa_shield', group_applicability='VHH only', definition='Absolute SASA shielded by CDR3 at the 15-position shieldable interface', units='Å²', primary_or_secondary='secondary'),
    dict(variable='cdr3_sequence_length', source_stage=6, source_file='data/processed/stage6_cdr3_shielding_geometry_combined.parquet',
         exact_column='cdr3_sequence_length', group_applicability='VHH only', definition='Real-residue CDR3 sequence length (IMGT 105-117, insertion-aware)', units='residues', primary_or_secondary='primary'),
    dict(variable='cdr3_rho', source_stage=6, source_file='data/processed/stage6_cdr3_shielding_geometry_combined.parquet',
         exact_column='cdr3_rho', group_applicability='VHH only', definition='Gordon reach: centroid distance CDR3 backbone to 102/103/118/119 anchor centroid', units='Å', primary_or_secondary='primary'),
    dict(variable='cdr3_compactness', source_stage=6, source_file='data/processed/stage6_cdr3_shielding_geometry_combined.parquet',
         exact_column='cdr3_compactness', group_applicability='VHH only', definition='Gordon compactness = cdr3_sequence_length / rho', units='dimensionless', primary_or_secondary='primary'),
    dict(variable='FR2_CDR3_Dis', source_stage=6, source_file='data/processed/stage6_cdr3_shielding_geometry_combined.parquet',
         exact_column='FR2_CDR3_Dis', group_applicability='VHH only', definition='Zhang: min distance from IMGT 42 centroid to internal CDR3 residue centroids', units='Å', primary_or_secondary='primary'),
    dict(variable='nonadj_BSA_pair_total', source_stage=7, source_file='results/tables/stage7_per_vhh_packing_summary.csv',
         exact_column='nonadj_BSA_pair_total', group_applicability='VHH only', definition='Pairwise mutual buried surface area, CDR3 vs 14-position NONADJACENT_PACKING_TARGET', units='Å²', primary_or_secondary='secondary (geometric corroboration for shielding)'),
    dict(variable='former_interface_core_PSH', source_stage=8, source_file='data/processed/stage8_surface_patch_descriptors.parquet',
         exact_column='former_interface_core_PSH', group_applicability='VHH,VH', definition='TAP-style adapted hydrophobic patch score, 28-position former-interface core', units='dimensionless (KD-normalized patch score)', primary_or_secondary='primary'),
    dict(variable='former_interface_core_PPC', source_stage=8, source_file='data/processed/stage8_surface_patch_descriptors.parquet',
         exact_column='former_interface_core_PPC', group_applicability='VHH,VH', definition='Adapted positive discrete charge-patch score, former-interface core', units='dimensionless', primary_or_secondary='secondary'),
    dict(variable='former_interface_core_PNC', source_stage=8, source_file='data/processed/stage8_surface_patch_descriptors.parquet',
         exact_column='former_interface_core_PNC', group_applicability='VHH,VH', definition='Adapted negative discrete charge-patch score, former-interface core', units='dimensionless', primary_or_secondary='secondary'),
    dict(variable='local_thr4.0_largest_patch_residue_count', source_stage=8, source_file='data/processed/stage8_patch_topology.parquet',
         exact_column='local_thr4.0_largest_patch_residue_count', group_applicability='VHH,VH', definition='PROJECT-SPECIFIC largest connected exposed-hydrophobic-residue patch, local (former-interface vicinity) scope, primary 4.0Å threshold', units='residue count', primary_or_secondary='primary'),
    dict(variable='interface_potential_positive_mean', source_stage=8, source_file='data/processed/stage8_apbs_interface_electrostatics.parquet',
         exact_column='interface_potential_positive_mean', group_applicability='VHH,VH', definition='APBS potential (kT/e) sampled at exposed former-interface atoms, positive-component mean', units='kT/e', primary_or_secondary='secondary'),
    dict(variable='interface_potential_negative_mean', source_stage=8, source_file='data/processed/stage8_apbs_interface_electrostatics.parquet',
         exact_column='interface_potential_negative_mean', group_applicability='VHH,VH', definition='APBS potential (kT/e) sampled at exposed former-interface atoms, negative-component mean', units='kT/e', primary_or_secondary='secondary'),
    dict(variable='framework_contact_fraction', source_stage=9, source_file='data/processed/stage9_antigen_engagement_summary.parquet',
         exact_column='framework_contact_fraction', group_applicability='VHH,VH(VH-only)', definition='Framework-region fraction of heavy-domain antigen-contact residues (4.5Å)', units='fraction [0,1]', primary_or_secondary='primary'),
    dict(variable='cdr3_contact_fraction', source_stage=9, source_file='data/processed/stage9_antigen_engagement_summary.parquet',
         exact_column='cdr3_contact_fraction', group_applicability='VHH,VH(VH-only)', definition='CDR3 fraction of heavy-domain antigen-contact residues (4.5Å)', units='fraction [0,1]', primary_or_secondary='primary'),
    dict(variable='shieldable_engagement_fraction', source_stage=9, source_file='data/interim/stage9/{vhh,vh}_shieldable_interface_engagement.csv',
         exact_column='fraction_target_contacting (target=15-position SHIELDABLE_NONCDR3)', group_applicability='VHH,VH(VH-only)', definition='Fraction of the 15-position non-CDR3 shieldable former interface contacting antigen (4.5Å)', units='fraction [0,1], denominator=15', primary_or_secondary='primary'),
    dict(variable='full_interface_engagement_fraction', source_stage=9, source_file='data/interim/stage9/{vhh,vh}_full_interface_engagement.csv',
         exact_column='fraction_target_contacting (target=28-position FULL_EMPIRICAL_INTERFACE)', group_applicability='VHH,VH(VH-only)', definition='Fraction of the 28-position full former interface contacting antigen (4.5Å)', units='fraction [0,1], denominator=25 unique positions', primary_or_secondary='secondary'),
    dict(variable='resolution', source_stage=2, source_file='data/processed/stage4_vhh_sequence_representatives.parquet / stage3_paired_representatives.parquet',
         exact_column='resolution', group_applicability='VHH,VH', definition='Crystallographic resolution', units='Å', primary_or_secondary='sensitivity covariate'),
    dict(variable='primary_interface_position_coverage', source_stage=4, source_file='data/processed/stage4_sequence_interface_chemistry.csv',
         exact_column='primary_interface_position_coverage', group_applicability='VHH,VH', definition='Fraction of 28-position former interface structurally evaluable', units='fraction [0,1]', primary_or_secondary='high-coverage sensitivity covariate'),
    dict(variable='antigen_name', source_stage=2, source_file='data/processed/stage4_vhh_sequence_representatives.parquet / stage3_paired_representatives.parquet',
         exact_column='antigen_name', group_applicability='VHH,VH', definition='Free-text antigen name from SAbDab2 metadata (exact identity only, not a validated family/class label)', units='text', primary_or_secondary='antigen-composition sensitivity grouping variable'),
]


def build():
    df = pd.DataFrame(VARS)
    df['missing_n'] = None
    df['missing_fraction'] = None
    df['notes'] = ''
    out = ROOT / "results/audits/stage11_input_variable_dictionary.csv"
    df.to_csv(out, index=False)
    print(f"Wrote {out}, {len(df)} variables")
    return df


if __name__ == '__main__':
    build()
