"""
Stage 5 — Position-Level Structural Comparison, Coverage/Noncanonical Audits,
High-Coverage Sensitivity, Interface-vs-Control, Three-State Comparison,
Intrinsic-vs-Structural Correlation.
"""

import sys
from pathlib import Path
import pandas as pd
import numpy as np
from scipy.stats import pearsonr, spearmanr
import logging

PROJECT_ROOT = Path(__file__).parent.parent.parent
PROCESSED_DIR = PROJECT_ROOT / "data" / "processed"
TABLES_DIR = PROJECT_ROOT / "results" / "tables"
AUDITS_DIR = PROJECT_ROOT / "results" / "audits"

logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger(__name__)

CANONICAL_20 = set('ARNDCQEGHILKMFPSTWYV')


def position_stats(combined: pd.DataFrame, structural_state: str) -> pd.DataFrame:
    sub = combined[(combined['structural_state'] == structural_state) & (combined['primary_binary_interface'] == True)]
    rows = []
    for (pos, ins), grp in sub.groupby(['imgt_position', 'insertion_code']):
        n_present = (grp['mapping_status'].isin(['PRESENT_OBSERVED', 'PRESENT_UNRESOLVED'])).sum()
        n_coord = (grp['coordinate_observed'] == True).sum()
        evaluable = grp[grp['sasa'].notna()]
        evaluable_kd = evaluable.dropna(subset=['kyte_doolittle'])

        pos_kd_sasa = evaluable_kd[evaluable_kd['positive_KD'] == True]['sasa'].sum()
        total_kd_sasa = evaluable_kd['sasa'].sum()

        rows.append({
            'imgt_position': pos, 'insertion_code': ins,
            'n_sequence_present': n_present,
            'n_coordinate_observed': n_coord,
            'coordinate_coverage': n_coord / n_present if n_present else None,
            'mean_sasa': evaluable['sasa'].mean() if len(evaluable) else None,
            'median_sasa': evaluable['sasa'].median() if len(evaluable) else None,
            'mean_KD_evaluable': evaluable_kd['kyte_doolittle'].mean() if len(evaluable_kd) else None,
            'mean_sasa_times_KD': (evaluable_kd['sasa'] * evaluable_kd['kyte_doolittle']).mean() if len(evaluable_kd) else None,
            'median_sasa_times_KD': (evaluable_kd['sasa'] * evaluable_kd['kyte_doolittle']).median() if len(evaluable_kd) else None,
            'positiveKD_sasa_fraction': pos_kd_sasa / total_kd_sasa if total_kd_sasa and total_kd_sasa > 0 else None,
        })
    df = pd.DataFrame(rows)
    df.columns = ['imgt_position', 'insertion_code'] + [f'{c}_{structural_state}' for c in df.columns[2:]]
    return df


def main():
    combined = pd.read_parquet(PROCESSED_DIR / "stage5_interface_surface_chemistry_residues.parquet")

    # ============ Position-level structural comparison ============
    vhh_pos = position_stats(combined, 'VHH_NATIVE')
    vh_without = position_stats(combined, 'VH_WITHOUT_VL_GEOMETRIC')
    vh_with = position_stats(combined, 'VH_WITH_VL_PAIRED')

    meta = combined.drop_duplicates(subset=['imgt_position', 'insertion_code'])[
        ['imgt_position', 'insertion_code', 'region', 'primary_binary_interface', 'interface_frequency_5']
    ]

    pos_comparison = meta.merge(vhh_pos, on=['imgt_position', 'insertion_code'], how='left')
    pos_comparison = pos_comparison.merge(vh_without, on=['imgt_position', 'insertion_code'], how='left')
    pos_comparison = pos_comparison.merge(vh_with, on=['imgt_position', 'insertion_code'], how='left')

    pos_comparison['delta_mean_sasa'] = pos_comparison['mean_sasa_VHH_NATIVE'] - pos_comparison['mean_sasa_VH_WITHOUT_VL_GEOMETRIC']
    pos_comparison['delta_mean_sasa_KD'] = pos_comparison['mean_sasa_times_KD_VHH_NATIVE'] - pos_comparison['mean_sasa_times_KD_VH_WITHOUT_VL_GEOMETRIC']
    pos_comparison['delta_positiveKD_sasa_fraction'] = pos_comparison['positiveKD_sasa_fraction_VHH_NATIVE'] - pos_comparison['positiveKD_sasa_fraction_VH_WITHOUT_VL_GEOMETRIC']

    pos_comparison = pos_comparison.sort_values(['imgt_position', 'insertion_code']).reset_index(drop=True)
    pos_out_csv = PROCESSED_DIR / "stage5_position_surface_chemistry.csv"
    pos_out_parquet = PROCESSED_DIR / "stage5_position_surface_chemistry.parquet"
    pos_comparison.to_csv(pos_out_csv, index=False)
    pos_comparison.to_parquet(pos_out_parquet, index=False)
    logger.info(f"Wrote {pos_out_csv} / .parquet: {len(pos_comparison)} positions")

    interface_only = pos_comparison[pos_comparison['primary_binary_interface'] == True]
    logger.info("\nLargest |delta_mean_sasa_KD| at primary interface positions:")
    logger.info(interface_only.reindex(interface_only['delta_mean_sasa_KD'].abs().sort_values(ascending=False).index)
                [['imgt_position', 'insertion_code', 'region', 'delta_mean_sasa', 'delta_mean_sasa_KD']].head(5).to_string(index=False))

    # ============ Coverage output ============
    seq_chem = pd.read_parquet(PROCESSED_DIR / "stage5_sequence_surface_chemistry.parquet")
    coverage_out = TABLES_DIR / "stage5_interface_structural_coverage.csv"
    coverage_summary = seq_chem.groupby(['group', 'structural_state'])[
        ['interface_sequence_coverage', 'interface_coordinate_coverage_given_present', 'interface_chemistry_coverage_given_present']
    ].describe()
    coverage_summary.to_csv(coverage_out)
    logger.info(f"\nWrote {coverage_out}")

    low_coverage = seq_chem[seq_chem['interface_coordinate_coverage_given_present'] < 0.5]
    logger.info(f"Sequences with interface_coordinate_coverage_given_present < 0.5: {len(low_coverage)} (retained, not excluded)")

    # ============ High-coverage sensitivity (>=0.90) ============
    high_cov = seq_chem[seq_chem['interface_coordinate_coverage_given_present'] >= 0.90]
    sensitivity_rows = []
    for state in ['VHH_NATIVE', 'VH_WITHOUT_VL_GEOMETRIC', 'VH_WITH_VL_PAIRED']:
        full = seq_chem[seq_chem['structural_state'] == state]
        hc = high_cov[high_cov['structural_state'] == state]
        sensitivity_rows.append({
            'structural_state': state,
            'n_full': len(full), 'n_high_coverage': len(hc),
            'interface_sasa_weighted_mean_KD_full': full['interface_sasa_weighted_mean_KD'].mean(),
            'interface_sasa_weighted_mean_KD_high_coverage': hc['interface_sasa_weighted_mean_KD'].mean(),
            'interface_positiveKD_sasa_fraction_full': full['interface_positiveKD_sasa_fraction'].mean(),
            'interface_positiveKD_sasa_fraction_high_coverage': hc['interface_positiveKD_sasa_fraction'].mean(),
        })
    sensitivity_df = pd.DataFrame(sensitivity_rows)
    sensitivity_out = TABLES_DIR / "stage5_high_coverage_sensitivity.csv"
    sensitivity_df.to_csv(sensitivity_out, index=False)
    logger.info(f"\nWrote {sensitivity_out}")
    logger.info(sensitivity_df.to_string(index=False))

    # ============ Noncanonical residue audit ============
    interface_present = combined[(combined['primary_binary_interface'] == True) &
                                   (combined['mapping_status'].isin(['PRESENT_OBSERVED', 'PRESENT_UNRESOLVED']))]
    noncanonical = interface_present[
        interface_present['amino_acid'].notna() & (~interface_present['amino_acid'].isin(CANONICAL_20))
    ]
    noncanonical_out = AUDITS_DIR / "stage5_noncanonical_residues.csv"
    noncanonical.to_csv(noncanonical_out, index=False)
    logger.info(f"\nWrote {noncanonical_out}: {len(noncanonical)} noncanonical observations within primary interface")

    # ============ Three-state comparison ============
    three_state_rows = []
    for state in ['VHH_NATIVE', 'VH_WITHOUT_VL_GEOMETRIC', 'VH_WITH_VL_PAIRED']:
        sub = seq_chem[seq_chem['structural_state'] == state]
        for col in ['interface_total_sasa', 'interface_mean_sasa_per_residue', 'interface_positiveKD_sasa',
                    'interface_positiveKD_sasa_fraction', 'interface_sasa_weighted_mean_KD',
                    'weighted_interface_sasa_mean_KD_freq5']:
            vals = sub[col].dropna()
            three_state_rows.append({
                'structural_state': state, 'descriptor': col, 'N': len(vals),
                'mean': vals.mean(), 'median': vals.median(), 'SD': vals.std(),
                'Q25': vals.quantile(0.25), 'Q75': vals.quantile(0.75),
            })
    three_state_df = pd.DataFrame(three_state_rows)
    three_state_out = TABLES_DIR / "stage5_three_state_surface_comparison.csv"
    three_state_df.to_csv(three_state_out, index=False)
    logger.info(f"\nWrote {three_state_out}")

    # ============ Interface vs non-interface framework control ============
    stage4_control = pd.read_csv(PROJECT_ROOT / "results/tables/stage4_interface_vs_noninterface_framework.csv")
    control_positions = set(zip(
        stage4_control[stage4_control['category'] == 'CONTROL_NONINTERFACE_FRAMEWORK']['imgt_position'],
        stage4_control[stage4_control['category'] == 'CONTROL_NONINTERFACE_FRAMEWORK']['insertion_code'].fillna(''),
    ))

    def structural_summary_for_positions(structural_state, pos_pairs, label):
        sub = combined[
            (combined['structural_state'] == structural_state) &
            (combined.set_index(['imgt_position', 'insertion_code']).index.isin(pos_pairs)) &
            (combined['sasa'].notna())
        ]
        sub_kd = sub.dropna(subset=['kyte_doolittle'])
        pos_kd_sasa = sub_kd[sub_kd['positive_KD'] == True]['sasa'].sum()
        total_kd_sasa = sub_kd['sasa'].sum()
        weighted_kd = (sub_kd['sasa'] * sub_kd['kyte_doolittle']).sum() / total_kd_sasa if total_kd_sasa else None
        return {
            'category': label, 'structural_state': structural_state,
            'n_observations': len(sub),
            'mean_exposed_sasa': sub['sasa'].mean(),
            'sasa_weighted_mean_KD': weighted_kd,
            'positiveKD_sasa_fraction': pos_kd_sasa / total_kd_sasa if total_kd_sasa else None,
        }

    interface_pos_pairs = set(zip(
        pos_comparison[pos_comparison['primary_binary_interface']]['imgt_position'],
        pos_comparison[pos_comparison['primary_binary_interface']]['insertion_code'],
    ))

    control_rows = []
    for state in ['VHH_NATIVE', 'VH_WITHOUT_VL_GEOMETRIC']:
        control_rows.append(structural_summary_for_positions(state, interface_pos_pairs, 'EMPIRICAL_INTERFACE'))
        control_rows.append(structural_summary_for_positions(state, control_positions, 'CONTROL_NONINTERFACE_FRAMEWORK'))
    control_df = pd.DataFrame(control_rows)
    control_out = TABLES_DIR / "stage5_interface_vs_noninterface_framework.csv"
    control_df.to_csv(control_out, index=False)
    logger.info(f"\nWrote {control_out} (descriptive only, no inferential significance test)")
    logger.info(control_df.to_string(index=False))

    # ============ Intrinsic vs structural correlation ============
    corr_rows = []
    for group, state in [('VHH', 'VHH_NATIVE'), ('VH', 'VH_WITHOUT_VL_GEOMETRIC')]:
        sub = seq_chem[(seq_chem['group'] == group) & (seq_chem['structural_state'] == state)].dropna(
            subset=['H_intrinsic_binary_mean_KD', 'interface_sasa_weighted_mean_KD']
        )
        if len(sub) > 2:
            pear = pearsonr(sub['H_intrinsic_binary_mean_KD'], sub['interface_sasa_weighted_mean_KD'])
            spear = spearmanr(sub['H_intrinsic_binary_mean_KD'], sub['interface_sasa_weighted_mean_KD'])
        else:
            pear, spear = (None, None), (None, None)

        sub2 = seq_chem[(seq_chem['group'] == group) & (seq_chem['structural_state'] == state)].dropna(
            subset=['H_intrinsic_weighted_freq5', 'weighted_interface_sasa_mean_KD_freq5']
        )
        if len(sub2) > 2:
            pear2 = pearsonr(sub2['H_intrinsic_weighted_freq5'], sub2['weighted_interface_sasa_mean_KD_freq5'])
            spear2 = spearmanr(sub2['H_intrinsic_weighted_freq5'], sub2['weighted_interface_sasa_mean_KD_freq5'])
        else:
            pear2, spear2 = (None, None), (None, None)

        corr_rows.append({
            'group': group, 'structural_state': state, 'N_binary': len(sub),
            'pearson_binary_intrinsic_vs_structural': pear[0], 'spearman_binary_intrinsic_vs_structural': spear[0],
            'N_weighted': len(sub2),
            'pearson_weighted_intrinsic_vs_structural': pear2[0], 'spearman_weighted_intrinsic_vs_structural': spear2[0],
            'intrinsic_mean': sub['H_intrinsic_binary_mean_KD'].mean(),
            'structural_mean': sub['interface_sasa_weighted_mean_KD'].mean(),
        })
    corr_df = pd.DataFrame(corr_rows)
    corr_out = TABLES_DIR / "stage5_intrinsic_vs_structural_chemistry.csv"
    corr_df.to_csv(corr_out, index=False)
    logger.info(f"\nWrote {corr_out}")
    logger.info(corr_df.to_string(index=False))

    return pos_comparison, sensitivity_df, control_df, corr_df, noncanonical


if __name__ == "__main__":
    main()
