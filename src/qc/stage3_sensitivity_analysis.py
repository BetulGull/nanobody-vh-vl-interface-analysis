"""
Stage 3 — Sensitivity Analyses: threshold combinations, redundancy comparison,
region-level summary, tetrad annotation summary, negative-ΔSASA audit.
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

DELTA_THRESHOLDS = [1, 5, 10]
FREQUENCY_THRESHOLDS = [0.10, 0.25, 0.50]
OBSERVATION_THRESHOLDS = [0.10, 0.25, 0.50]
TETRAD_POSITIONS = {42, 49, 50, 52}
PRIMARY_DELTA, PRIMARY_FREQ, PRIMARY_OBS = 5, 0.25, 0.25

IMGT_REGIONS = {
    'FR1': (1, 26), 'CDR1': (27, 38), 'FR2': (39, 55), 'CDR2': (56, 65),
    'FR3': (66, 104), 'CDR3': (105, 117), 'FR4': (118, 128),
}


def threshold_sensitivity(residue_df: pd.DataFrame, n_structures: int, primary_positions: set) -> pd.DataFrame:
    """For each (delta_threshold, freq_threshold, obs_threshold) combo, recompute consensus set."""
    rows = []

    # Precompute per-position n_observed and per-threshold interface counts once
    grouped = residue_df.groupby(['imgt_position', 'insertion_code'])
    base = grouped.agg(
        n_observed=('delta_sasa_vl_raw', 'count'),
        n_int_1=('interface_dSASA_1', 'sum'),
        n_int_5=('interface_dSASA_5', 'sum'),
        n_int_10=('interface_dSASA_10', 'sum'),
    ).reset_index()
    base['observation_fraction'] = base['n_observed'] / n_structures
    base['freq_1'] = base['n_int_1'] / base['n_observed']
    base['freq_5'] = base['n_int_5'] / base['n_observed']
    base['freq_10'] = base['n_int_10'] / base['n_observed']
    base['position_key'] = list(zip(base['imgt_position'], base['insertion_code']))

    freq_col_map = {1: 'freq_1', 5: 'freq_5', 10: 'freq_10'}

    for delta_t in DELTA_THRESHOLDS:
        freq_col = freq_col_map[delta_t]
        for freq_t in FREQUENCY_THRESHOLDS:
            for obs_t in OBSERVATION_THRESHOLDS:
                mask = (base[freq_col] >= freq_t) & (base['observation_fraction'] >= obs_t)
                positions_included = set(base.loc[mask, 'position_key'])
                n_positions = len(positions_included)

                overlap = len(positions_included & primary_positions)
                union = len(positions_included | primary_positions)
                jaccard = overlap / union if union else 1.0

                n_tetrad = sum(1 for (pos, ins) in positions_included if pos in TETRAD_POSITIONS)

                rows.append({
                    'delta_sasa_threshold': delta_t,
                    'frequency_threshold': freq_t,
                    'observation_fraction_threshold': obs_t,
                    'n_consensus_positions': n_positions,
                    'overlap_with_primary': overlap,
                    'jaccard_with_primary': jaccard,
                    'n_tetrad_positions_included': n_tetrad,
                    'is_primary_definition': (delta_t == PRIMARY_DELTA and freq_t == PRIMARY_FREQ and obs_t == PRIMARY_OBS),
                })

    return pd.DataFrame(rows)


def redundancy_sensitivity(primary_map: pd.DataFrame, all_map: pd.DataFrame) -> dict:
    """Compare primary nonredundant map vs all-instance map."""
    merged = primary_map.merge(
        all_map, on=['imgt_position', 'insertion_code'], suffixes=('_primary', '_all')
    )

    pearson_freq = pearsonr(merged['interface_frequency_5_primary'], merged['interface_frequency_5_all'])
    spearman_freq = spearmanr(merged['interface_frequency_5_primary'], merged['interface_frequency_5_all'])

    mean_abs_diff = (merged['interface_frequency_5_primary'] - merged['interface_frequency_5_all']).abs().mean()
    max_abs_diff = (merged['interface_frequency_5_primary'] - merged['interface_frequency_5_all']).abs().max()

    primary_binary = set(zip(
        merged.loc[merged['primary_binary_interface_primary'], 'imgt_position'],
        merged.loc[merged['primary_binary_interface_primary'], 'insertion_code'],
    ))
    all_binary = set(zip(
        merged.loc[merged['primary_binary_interface_all'], 'imgt_position'],
        merged.loc[merged['primary_binary_interface_all'], 'insertion_code'],
    ))
    jaccard_binary = len(primary_binary & all_binary) / len(primary_binary | all_binary) if (primary_binary | all_binary) else 1.0

    return {
        'pearson_r_interface_frequency_5': pearson_freq[0],
        'pearson_p': pearson_freq[1],
        'spearman_r_interface_frequency_5': spearman_freq[0],
        'spearman_p': spearman_freq[1],
        'mean_abs_diff_interface_frequency_5': mean_abs_diff,
        'max_abs_diff_interface_frequency_5': max_abs_diff,
        'binary_interface_jaccard': jaccard_binary,
        'n_positions_compared': len(merged),
        'n_primary_binary_interface': len(primary_binary),
        'n_all_instance_binary_interface': len(all_binary),
    }


def region_summary(primary_map: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for region in IMGT_REGIONS:
        sub = primary_map[primary_map['region'] == region]
        rows.append({
            'region': region,
            'n_positions_total': len(sub),
            'n_primary_consensus_positions': sub['primary_binary_interface'].sum(),
            'mean_interface_frequency_5': sub['interface_frequency_5'].mean(),
            'median_interface_frequency_5': sub['interface_frequency_5'].median(),
            'mean_delta_sasa': sub['mean_delta_sasa'].mean(),
        })
    return pd.DataFrame(rows)


def tetrad_summary(primary_map: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for pos in sorted(TETRAD_POSITIONS):
        sub = primary_map[primary_map['imgt_position'] == pos]
        if len(sub) == 0:
            rows.append({'imgt_position': pos, 'n_observed': 0})
            continue
        # combine across insertion-code variants at this position (rare) by summing observed n and weighting means
        n_obs_total = sub['n_observed'].sum()
        mean_delta = np.average(sub['mean_delta_sasa'], weights=sub['n_observed']) if n_obs_total else None
        median_delta = sub['median_delta_sasa'].median()
        freq1 = np.average(sub['interface_frequency_1'], weights=sub['n_observed']) if n_obs_total else None
        freq5 = np.average(sub['interface_frequency_5'], weights=sub['n_observed']) if n_obs_total else None
        freq10 = np.average(sub['interface_frequency_10'], weights=sub['n_observed']) if n_obs_total else None
        in_primary = bool(sub['primary_binary_interface'].any())
        rows.append({
            'imgt_position': pos,
            'n_observed': n_obs_total,
            'mean_delta_sasa': mean_delta,
            'median_delta_sasa': median_delta,
            'frequency_ge_1': freq1,
            'frequency_ge_5': freq5,
            'frequency_ge_10': freq10,
            'enters_primary_binary_interface': in_primary,
        })
    return pd.DataFrame(rows)


def negative_delta_audit(residue_df: pd.DataFrame) -> dict:
    deltas = residue_df['delta_sasa_vl_raw']
    negative = deltas[deltas < 0]
    return {
        'n_total_residue_observations': len(deltas),
        'n_delta_negative': len(negative),
        'min_delta': deltas.min(),
        'n_delta_lt_minus_0.1': (deltas < -0.1).sum(),
        'n_delta_lt_minus_0.5': (deltas < -0.5).sum(),
        'n_delta_lt_minus_1.0': (deltas < -1.0).sum(),
        'negative_fraction': len(negative) / len(deltas) if len(deltas) else 0,
        'negative_mean': negative.mean() if len(negative) else None,
        'negative_std': negative.std() if len(negative) else None,
    }


def main():
    residue_df = pd.read_parquet(PROCESSED_DIR / "stage3_vh_vl_delta_sasa_residues.parquet")
    primary_map = pd.read_parquet(PROCESSED_DIR / "stage3_empirical_vh_vl_interface_primary.parquet")
    # keep_default_na=False: prevents pandas from turning empty-string
    # insertion_code ('') into NaN on CSV read, which would silently break
    # the merge join against the parquet-loaded primary_map (where empty
    # insertion codes are the literal string '') in redundancy_sensitivity().
    all_map = pd.read_csv(
        PROCESSED_DIR / "stage3_empirical_vh_vl_interface_all_instances.csv",
        keep_default_na=False,
    )

    primary_positions = set(zip(
        primary_map.loc[primary_map['primary_binary_interface'], 'imgt_position'],
        primary_map.loc[primary_map['primary_binary_interface'], 'insertion_code'],
    ))

    primary_residues = residue_df[residue_df['is_primary_representative'] == True]
    n_primary_structures = primary_residues['INSTANCE'].nunique()

    # --- Threshold sensitivity ---
    threshold_df = threshold_sensitivity(primary_residues, n_primary_structures, primary_positions)
    threshold_out = TABLES_DIR / "stage3_interface_threshold_sensitivity.csv"
    threshold_df.to_csv(threshold_out, index=False)
    logger.info(f"Wrote {threshold_out} ({len(threshold_df)} combinations)")

    # --- Redundancy sensitivity ---
    redundancy_result = redundancy_sensitivity(primary_map, all_map)
    redundancy_df = pd.DataFrame([redundancy_result])
    redundancy_out = TABLES_DIR / "stage3_redundancy_sensitivity.csv"
    redundancy_df.to_csv(redundancy_out, index=False)
    logger.info(f"Wrote {redundancy_out}")
    logger.info(redundancy_df.T.to_string(header=False))

    # --- Region summary ---
    region_df = region_summary(primary_map)
    region_out = TABLES_DIR / "stage3_region_summary.csv"
    region_df.to_csv(region_out, index=False)
    logger.info(f"Wrote {region_out}")

    # --- Tetrad summary ---
    tetrad_df = tetrad_summary(primary_map)
    tetrad_out = TABLES_DIR / "stage3_tetrad_summary.csv"
    tetrad_df.to_csv(tetrad_out, index=False)
    logger.info(f"Wrote {tetrad_out}")
    logger.info(tetrad_df.to_string(index=False))

    # --- Negative delta audit ---
    neg_result = negative_delta_audit(residue_df)
    neg_df = pd.DataFrame([neg_result])
    neg_out = AUDITS_DIR / "stage3_negative_delta_sasa.csv"
    neg_df.to_csv(neg_out, index=False)
    logger.info(f"Wrote {neg_out}")
    logger.info(neg_df.T.to_string(header=False))

    # --- Primary interface position list ---
    primary_interface_list = primary_map[primary_map['primary_binary_interface']][
        ['imgt_position', 'region', 'n_observed', 'observation_fraction',
         'interface_frequency_5', 'mean_delta_sasa', 'median_delta_sasa', 'vhh_tetrad_position']
    ].copy()
    primary_list_out = TABLES_DIR / "stage3_primary_interface_positions.csv"
    primary_interface_list.to_csv(primary_list_out, index=False)
    logger.info(f"Wrote {primary_list_out}: {len(primary_interface_list)} positions")

    return threshold_df, redundancy_result, region_df, tetrad_df, neg_result


if __name__ == "__main__":
    main()
