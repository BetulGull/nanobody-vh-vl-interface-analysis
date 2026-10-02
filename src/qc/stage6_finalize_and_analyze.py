"""
Stage 6 — Finalize Outputs, Negative-Shielding Audit, Correlations,
Position-Level Map, Sensitivity, Noncanonical Audit, Combined Table.
"""

import sys
from pathlib import Path
import pandas as pd
import numpy as np
from scipy.stats import pearsonr, spearmanr
import logging

sys.path.insert(0, str(Path(__file__).parent.parent))

PROJECT_ROOT = Path(__file__).parent.parent.parent
CHECKPOINT_DIR = PROJECT_ROOT / "data" / "interim" / "stage6" / "checkpoints"
PROCESSED_DIR = PROJECT_ROOT / "data" / "processed"
TABLES_DIR = PROJECT_ROOT / "results" / "tables"
AUDITS_DIR = PROJECT_ROOT / "results" / "audits"

logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger(__name__)


def load_final(prefix):
    files = sorted(CHECKPOINT_DIR.glob(f"{prefix}_*.parquet"))
    return pd.read_parquet(files[-1])


def main():
    residues = load_final("residues")
    summaries = load_final("summaries")
    geometries = load_final("geometries")

    # --- Core outputs ---
    residues.to_parquet(PROCESSED_DIR / "stage6_cdr3_shielding_residues.parquet", index=False)
    summaries.to_parquet(PROCESSED_DIR / "stage6_vhh_cdr3_shielding.parquet", index=False)
    summaries.to_csv(PROCESSED_DIR / "stage6_vhh_cdr3_shielding.csv", index=False)
    geometries.to_parquet(PROCESSED_DIR / "stage6_cdr3_geometry.parquet", index=False)
    geometries.to_csv(PROCESSED_DIR / "stage6_cdr3_geometry.csv", index=False)
    logger.info(f"Wrote core Stage-6 outputs: {len(residues)} residue rows, {len(summaries)} VHH summaries, {len(geometries)} geometries")

    # --- Negative shielding audit ---
    residue_deltas = residues['delta_sasa_shield'].dropna()
    neg_residue = {
        'n_total_residue_observations': len(residue_deltas),
        'n_delta_negative': (residue_deltas < 0).sum(),
        'min_delta': residue_deltas.min(),
        'n_lt_minus_0.1': (residue_deltas < -0.1).sum(),
        'n_lt_minus_0.5': (residue_deltas < -0.5).sum(),
        'n_lt_minus_1.0': (residue_deltas < -1.0).sum(),
    }
    struct_deltas = summaries['delta_sasa_shield'].dropna()
    neg_struct = {
        'n_total_structures': len(struct_deltas),
        'n_delta_negative': (struct_deltas < 0).sum(),
        'min_delta': struct_deltas.min(),
        'n_lt_minus_0.1': (struct_deltas < -0.1).sum(),
        'n_lt_minus_0.5': (struct_deltas < -0.5).sum(),
        'n_lt_minus_1.0': (struct_deltas < -1.0).sum(),
    }
    neg_df = pd.DataFrame([{'level': 'residue', **neg_residue}, {'level': 'structure', **neg_struct}])
    neg_out = AUDITS_DIR / "stage6_negative_shielding.csv"
    neg_df.to_csv(neg_out, index=False)
    logger.info(f"\nWrote {neg_out}")
    logger.info(neg_df.to_string(index=False))

    # --- Noncanonical residue audit ---
    noncanonical = residues[residues['noncanonical_flag'] == True]
    noncanonical.to_csv(AUDITS_DIR / "stage6_noncanonical_residues.csv", index=False)
    logger.info(f"\nNoncanonical residues in shieldable interface: {len(noncanonical)}")

    # --- Shielding distribution table ---
    dist_rows = []
    for col in ['delta_sasa_shield', 'shielding_fraction']:
        vals = summaries[col].dropna()
        dist_rows.append({
            'descriptor': col, 'N': len(vals), 'mean': vals.mean(), 'median': vals.median(),
            'SD': vals.std(), 'Q25': vals.quantile(0.25), 'Q75': vals.quantile(0.75),
            'min': vals.min(), 'max': vals.max(),
        })
    thresholds = [1, 5, 10, 25, 50]
    threshold_fractions = {}
    for t in thresholds:
        threshold_fractions[f'fraction_delta_shield_lt_{t}'] = (summaries['delta_sasa_shield'] < t).mean()
        threshold_fractions[f'fraction_delta_shield_ge_{t}'] = (summaries['delta_sasa_shield'] >= t).mean()
    dist_df = pd.DataFrame(dist_rows)
    dist_out = TABLES_DIR / "stage6_shielding_distribution.csv"
    dist_df.to_csv(dist_out, index=False)
    pd.DataFrame([threshold_fractions]).to_csv(TABLES_DIR / "stage6_shielding_threshold_fractions.csv", index=False)
    logger.info(f"\nWrote {dist_out}")
    logger.info(dist_df.to_string(index=False))
    logger.info(f"\nThreshold fractions: {threshold_fractions}")

    # --- High-coverage sensitivity ---
    high_cov = summaries[
        (summaries['shieldable_interface_coordinate_coverage'] >= 0.90)
    ].merge(geometries[['INSTANCE', 'cdr3_coordinate_coverage']], on='INSTANCE', how='left')
    high_cov = high_cov[high_cov['cdr3_coordinate_coverage'] >= 0.90]

    sensitivity_rows = [{
        'n_full': len(summaries), 'n_high_coverage': len(high_cov),
        'delta_sasa_shield_mean_full': summaries['delta_sasa_shield'].mean(),
        'delta_sasa_shield_mean_high_coverage': high_cov['delta_sasa_shield'].mean(),
        'shielding_fraction_mean_full': summaries['shielding_fraction'].mean(),
        'shielding_fraction_mean_high_coverage': high_cov['shielding_fraction'].mean(),
    }]
    sens_df = pd.DataFrame(sensitivity_rows)
    sens_out = TABLES_DIR / "stage6_high_coverage_sensitivity.csv"
    sens_df.to_csv(sens_out, index=False)
    logger.info(f"\nWrote {sens_out}")
    logger.info(sens_df.to_string(index=False))

    # --- Position-level shielding map ---
    pos_rows = []
    for (pos, ins), grp in residues.groupby(['imgt_position', 'insertion_code']):
        evaluable = grp.dropna(subset=['delta_sasa_shield'])
        n_present = (grp['coordinate_observed'] | grp['amino_acid'].notna()).sum()
        n_coord_observed = grp['coordinate_observed'].sum()
        deltas = evaluable['delta_sasa_shield']
        pos_rows.append({
            'imgt_position': pos, 'insertion_code': ins, 'region': grp['region'].iloc[0],
            'n_sequence_present': int(n_present), 'n_coordinate_observed': int(n_coord_observed),
            'n_shielding_evaluable': len(evaluable),
            'mean_sasa_native': evaluable['sasa_native'].mean() if len(evaluable) else None,
            'mean_sasa_noCDR3': evaluable['sasa_noCDR3'].mean() if len(evaluable) else None,
            'mean_delta_sasa_shield': deltas.mean() if len(deltas) else None,
            'median_delta_sasa_shield': deltas.median() if len(deltas) else None,
            'q25': deltas.quantile(0.25) if len(deltas) else None,
            'q75': deltas.quantile(0.75) if len(deltas) else None,
            'fraction_delta_shield_gt_1': (deltas > 1).mean() if len(deltas) else None,
            'fraction_delta_shield_gt_5': (deltas > 5).mean() if len(deltas) else None,
            'fraction_delta_shield_gt_10': (deltas > 10).mean() if len(deltas) else None,
            'mean_interface_frequency_5': grp['interface_frequency_5'].iloc[0],
        })
    pos_df = pd.DataFrame(pos_rows).sort_values(['imgt_position', 'insertion_code']).reset_index(drop=True)
    pos_df.to_csv(PROCESSED_DIR / "stage6_position_shielding_summary.csv", index=False)
    pos_df.to_parquet(PROCESSED_DIR / "stage6_position_shielding_summary.parquet", index=False)
    logger.info(f"\nPosition-level shielding map: {len(pos_df)} positions")
    logger.info(pos_df[['imgt_position', 'insertion_code', 'region', 'mean_delta_sasa_shield']].to_string(index=False))

    # --- Combined mechanism table ---
    stage4_chem = pd.read_parquet(PROCESSED_DIR / "stage4_sequence_interface_chemistry.parquet")
    stage4_vhh = stage4_chem[stage4_chem['group'] == 'VHH'].rename(columns={'representative_INSTANCE': 'INSTANCE'})[
        ['INSTANCE', 'H_intrinsic_binary_mean_KD', 'H_intrinsic_weighted_freq5']
    ]
    stage5_seq = pd.read_parquet(PROCESSED_DIR / "stage5_sequence_surface_chemistry.parquet")
    stage5_vhh = stage5_seq[stage5_seq['structural_state'] == 'VHH_NATIVE'][
        ['INSTANCE', 'interface_sasa_weighted_mean_KD', 'interface_positiveKD_sasa_fraction', 'interface_total_sasa']
    ].rename(columns={
        'interface_sasa_weighted_mean_KD': 'H_native_full_interface_sasa_weighted_mean_KD',
        'interface_positiveKD_sasa_fraction': 'H_native_full_interface_positiveKD_fraction',
        'interface_total_sasa': 'H_native_full_interface_total_sasa',
    })

    combined = summaries.merge(geometries.drop(columns=['sequence_group_id', 'SABDAB_ID', 'PDB']), on='INSTANCE', how='left')
    combined = combined.merge(stage4_vhh, on='INSTANCE', how='left')
    combined = combined.merge(stage5_vhh, on='INSTANCE', how='left')

    combined.to_parquet(PROCESSED_DIR / "stage6_cdr3_shielding_geometry_combined.parquet", index=False)
    combined.to_csv(PROCESSED_DIR / "stage6_cdr3_shielding_geometry_combined.csv", index=False)
    logger.info(f"\nWrote combined mechanism table: {len(combined)} rows, {len(combined.columns)} columns")

    # --- Relationship: CDR3 length vs shielding ---
    def corr_report(x, y, label):
        sub = pd.DataFrame({'x': x, 'y': y}).dropna()
        if len(sub) > 2:
            pear = pearsonr(sub['x'], sub['y'])
            spear = spearmanr(sub['x'], sub['y'])
            return {'relationship': label, 'N': len(sub), 'pearson_r': pear[0], 'pearson_p_exploratory': pear[1],
                    'spearman_r': spear[0], 'spearman_p_exploratory': spear[1]}
        return {'relationship': label, 'N': len(sub), 'pearson_r': None, 'pearson_p_exploratory': None,
                'spearman_r': None, 'spearman_p_exploratory': None}

    geom_corr_rows = [
        corr_report(combined['cdr3_sequence_length'], combined['delta_sasa_shield'], 'cdr3_length_vs_delta_sasa_shield'),
        corr_report(combined['cdr3_sequence_length'], combined['shielding_fraction'], 'cdr3_length_vs_shielding_fraction'),
        corr_report(combined['cdr3_rho'], combined['delta_sasa_shield'], 'rho_vs_delta_sasa_shield'),
        corr_report(combined['cdr3_rho'], combined['shielding_fraction'], 'rho_vs_shielding_fraction'),
        corr_report(combined['cdr3_compactness'], combined['delta_sasa_shield'], 'compactness_vs_delta_sasa_shield'),
        corr_report(combined['cdr3_compactness'], combined['shielding_fraction'], 'compactness_vs_shielding_fraction'),
        corr_report(combined['FR2_CDR3_Dis'], combined['delta_sasa_shield'], 'FR2_CDR3_Dis_vs_delta_sasa_shield'),
        corr_report(combined['FR2_CDR3_Dis'], combined['shielding_fraction'], 'FR2_CDR3_Dis_vs_shielding_fraction'),
        corr_report(combined['cdr3_rho'], combined['weighted_delta_sasa_shield_freq5'], 'rho_vs_weighted_delta_shield_freq5'),
    ]
    geom_corr_df = pd.DataFrame(geom_corr_rows)
    geom_corr_out = TABLES_DIR / "stage6_geometry_shielding_correlations.csv"
    geom_corr_df.to_csv(geom_corr_out, index=False)
    logger.info(f"\nWrote {geom_corr_out} (exploratory correlations, no formal inference)")
    logger.info(geom_corr_df.to_string(index=False))

    # --- Relationship: H_intrinsic / H_native vs shielding ---
    intrinsic_corr_rows = [
        corr_report(combined['H_intrinsic_binary_mean_KD'], combined['delta_sasa_shield'], 'H_intrinsic_binary_mean_KD_vs_delta_sasa_shield'),
        corr_report(combined['H_intrinsic_binary_mean_KD'], combined['shielding_fraction'], 'H_intrinsic_binary_mean_KD_vs_shielding_fraction'),
        corr_report(combined['H_intrinsic_weighted_freq5'], combined['weighted_delta_sasa_shield_freq5'], 'H_intrinsic_weighted_freq5_vs_weighted_delta_shield_freq5'),
        corr_report(combined['native_sasa_weighted_mean_KD'], combined['delta_sasa_shield'], 'H_native_shieldable_vs_delta_sasa_shield'),
        corr_report(combined['native_sasa_weighted_mean_KD'], combined['shielding_fraction'], 'H_native_shieldable_vs_shielding_fraction'),
    ]
    intrinsic_corr_df = pd.DataFrame(intrinsic_corr_rows)
    intrinsic_corr_out = TABLES_DIR / "stage6_intrinsic_shielding_correlations.csv"
    intrinsic_corr_df.to_csv(intrinsic_corr_out, index=False)
    logger.info(f"\nWrote {intrinsic_corr_out} (exploratory correlations, no formal inference)")
    logger.info(intrinsic_corr_df.to_string(index=False))

    return combined, dist_df, geom_corr_df, intrinsic_corr_df


if __name__ == "__main__":
    main()
