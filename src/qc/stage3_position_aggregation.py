"""
Stage 3 — Position-Level Aggregation

Builds the PRIMARY (nonredundant representatives) and ALL-INSTANCE
(sensitivity) empirical VH-VL interface maps from the long residue table.
"""

import sys
from pathlib import Path
import pandas as pd
import numpy as np
import logging

PROJECT_ROOT = Path(__file__).parent.parent.parent
PROCESSED_DIR = PROJECT_ROOT / "data" / "processed"

logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger(__name__)

TETRAD_POSITIONS = {42, 49, 50, 52}

PRIMARY_DELTA_THRESHOLD = 5
PRIMARY_FREQUENCY_THRESHOLD = 0.25
PRIMARY_OBSERVATION_FRACTION_THRESHOLD = 0.25


def aggregate_positions(residue_df: pd.DataFrame, n_structures: int, map_type: str) -> pd.DataFrame:
    """
    Aggregate the long residue table by exact IMGT VH position
    (imgt_number + insertion_code), matching spec §15.
    """
    rows = []
    group_cols = ['imgt_position', 'insertion_code']

    for (imgt_pos, ins_code), grp in residue_df.groupby(group_cols):
        n_obs = len(grp)
        deltas = grp['delta_sasa_vl_raw']

        row = {
            'imgt_position': imgt_pos,
            'imgt_number': imgt_pos,
            'insertion_code': ins_code,
            'region': grp['region'].iloc[0],
            'n_observed': n_obs,
            'observation_fraction': n_obs / n_structures if n_structures else None,
            'mean_delta_sasa': deltas.mean(),
            'median_delta_sasa': deltas.median(),
            'sd_delta_sasa': deltas.std(ddof=0) if n_obs > 1 else 0.0,
            'q25_delta_sasa': deltas.quantile(0.25),
            'q75_delta_sasa': deltas.quantile(0.75),
            'max_delta_sasa': deltas.max(),
            'n_interface_1': grp['interface_dSASA_1'].sum(),
            'n_interface_5': grp['interface_dSASA_5'].sum(),
            'n_interface_10': grp['interface_dSASA_10'].sum(),
            'interface_frequency_1': grp['interface_dSASA_1'].sum() / n_obs if n_obs else None,
            'interface_frequency_5': grp['interface_dSASA_5'].sum() / n_obs if n_obs else None,
            'interface_frequency_10': grp['interface_dSASA_10'].sum() / n_obs if n_obs else None,
            'vhh_tetrad_position': imgt_pos in TETRAD_POSITIONS,
            'map_type': map_type,
        }
        rows.append(row)

    df = pd.DataFrame(rows).sort_values(['imgt_position', 'insertion_code']).reset_index(drop=True)

    df['primary_binary_interface'] = (
        (df['interface_frequency_5'] >= PRIMARY_FREQUENCY_THRESHOLD) &
        (df['observation_fraction'] >= PRIMARY_OBSERVATION_FRACTION_THRESHOLD)
    )
    df['w_interface_frequency_5'] = df['interface_frequency_5']
    df['w_mean_delta_sasa'] = df['mean_delta_sasa']
    df['w_median_delta_sasa'] = df['median_delta_sasa']

    return df


def main():
    residue_df = pd.read_parquet(PROCESSED_DIR / "stage3_vh_vl_delta_sasa_residues.parquet")

    # --- PRIMARY map: nonredundant representatives only ---
    primary_residues = residue_df[residue_df['is_primary_representative'] == True]
    n_primary_structures = primary_residues['INSTANCE'].nunique()
    logger.info(f"Primary map: {n_primary_structures} representative structures, {len(primary_residues)} residue rows")

    primary_map = aggregate_positions(primary_residues, n_primary_structures, 'PRIMARY_NONREDUNDANT')
    primary_map.to_csv(PROCESSED_DIR / "stage3_empirical_vh_vl_interface_primary.csv", index=False)
    primary_map.to_parquet(PROCESSED_DIR / "stage3_empirical_vh_vl_interface_primary.parquet", index=False)
    logger.info(f"Wrote primary map: {len(primary_map)} VH IMGT positions")

    # --- ALL-INSTANCE sensitivity map ---
    n_all_structures = residue_df['INSTANCE'].nunique()
    logger.info(f"All-instance map: {n_all_structures} structures, {len(residue_df)} residue rows")

    all_map = aggregate_positions(residue_df, n_all_structures, 'ALL_INSTANCE_SENSITIVITY')
    all_map.to_csv(PROCESSED_DIR / "stage3_empirical_vh_vl_interface_all_instances.csv", index=False)
    logger.info(f"Wrote all-instance map: {len(all_map)} VH IMGT positions")

    n_primary_interface = primary_map['primary_binary_interface'].sum()
    logger.info(f"\nPrimary binary interface positions: {n_primary_interface}/{len(primary_map)}")

    return primary_map, all_map


if __name__ == "__main__":
    main()
