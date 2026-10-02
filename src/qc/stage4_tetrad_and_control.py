"""
Stage 4 — Tetrad Subset Analysis and Empirical-Interface-vs-Control Comparison
"""

import sys
from pathlib import Path
import pandas as pd
import numpy as np
import logging

PROJECT_ROOT = Path(__file__).parent.parent.parent
PROCESSED_DIR = PROJECT_ROOT / "data" / "processed"
TABLES_DIR = PROJECT_ROOT / "results" / "tables"

logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger(__name__)

TETRAD_POSITIONS = {42, 49, 50, 52}
CONTROL_MIN_OBSERVATION_FRACTION = 0.25


def main():
    comparison = pd.read_parquet(PROCESSED_DIR / "stage4_position_chemistry_comparison.parquet")
    long_df = pd.read_parquet(PROCESSED_DIR / "stage4_interface_projection_long.parquet")

    # --- Tetrad comparison table ---
    tetrad_rows = comparison[comparison['imgt_position'].isin(TETRAD_POSITIONS)].copy()
    tetrad_out = TABLES_DIR / "stage4_tetrad_comparison.csv"
    tetrad_rows.to_csv(tetrad_out, index=False)
    logger.info(f"Wrote {tetrad_out}: {len(tetrad_rows)} tetrad rows")
    logger.info(tetrad_rows[['imgt_position', 'interface_frequency_5', 'mean_KD_VHH', 'mean_KD_VH',
                              'delta_mean_KD', 'JSD', 'primary_binary_interface']].to_string(index=False))

    # --- Tetrad-only vs full-interface comparison ---
    full_interface = comparison[comparison['primary_binary_interface'] == True]
    tetrad_in_interface = full_interface[full_interface['imgt_position'].isin(TETRAD_POSITIONS)]

    n_tetrad_in_full_interface = len(tetrad_in_interface)
    n_full_interface = len(full_interface)

    total_freq_weight_full = full_interface['interface_frequency_5'].sum()
    total_freq_weight_tetrad = tetrad_in_interface['interface_frequency_5'].sum()
    fraction_weight_carried_by_tetrad = (
        total_freq_weight_tetrad / total_freq_weight_full if total_freq_weight_full else None
    )

    summary_rows = [
        {'subset': 'TETRAD_ONLY (42,49,50,52)', 'n_positions': len(tetrad_rows),
         'n_in_primary_binary_interface': (tetrad_rows['primary_binary_interface']).sum(),
         'mean_interface_frequency_5': tetrad_rows['interface_frequency_5'].mean(),
         'mean_delta_mean_KD': tetrad_rows['delta_mean_KD'].mean(),
         'mean_JSD': tetrad_rows['JSD'].mean()},
        {'subset': 'FULL_PRIMARY_INTERFACE (28 positions)', 'n_positions': n_full_interface,
         'n_in_primary_binary_interface': n_full_interface,
         'mean_interface_frequency_5': full_interface['interface_frequency_5'].mean(),
         'mean_delta_mean_KD': full_interface['delta_mean_KD'].mean(),
         'mean_JSD': full_interface['JSD'].mean()},
    ]
    tetrad_vs_full_df = pd.DataFrame(summary_rows)
    tetrad_vs_full_df['n_tetrad_positions_in_full_interface'] = n_tetrad_in_full_interface
    tetrad_vs_full_df['fraction_of_interface_frequency_weight_from_tetrad'] = fraction_weight_carried_by_tetrad

    tetrad_vs_full_out = TABLES_DIR / "stage4_tetrad_vs_full_interface.csv"
    tetrad_vs_full_df.to_csv(tetrad_vs_full_out, index=False)
    logger.info(f"\nWrote {tetrad_vs_full_out}")
    logger.info(tetrad_vs_full_df.to_string(index=False))
    logger.info(f"\nTetrad positions in full 28-position interface: {n_tetrad_in_full_interface}/4")
    logger.info(f"Fraction of interface-frequency weight from tetrad: {fraction_weight_carried_by_tetrad:.3f}"
                if fraction_weight_carried_by_tetrad is not None else "N/A")

    # --- Interface-vs-non-interface-framework control ---
    framework_regions = {'FR1', 'FR2', 'FR3', 'FR4'}
    control_candidates = comparison[
        (comparison['region'].isin(framework_regions)) &
        (comparison['primary_binary_interface'] == False) &
        (comparison['observation_fraction'] >= CONTROL_MIN_OBSERVATION_FRACTION)
    ].copy()

    interface_positions = comparison[comparison['primary_binary_interface'] == True].copy()

    interface_positions['category'] = 'EMPIRICAL_INTERFACE'
    control_candidates['category'] = 'CONTROL_NONINTERFACE_FRAMEWORK'

    combined = pd.concat([interface_positions, control_candidates], ignore_index=True)
    control_out = TABLES_DIR / "stage4_interface_vs_noninterface_framework.csv"
    combined.to_csv(control_out, index=False)
    logger.info(f"\nWrote {control_out}: {len(interface_positions)} interface + {len(control_candidates)} control positions")

    interface_summary = {
        'category': 'EMPIRICAL_INTERFACE',
        'n_positions': len(interface_positions),
        'mean_delta_mean_KD': interface_positions['delta_mean_KD'].mean(),
        'mean_abs_delta_mean_KD': interface_positions['delta_mean_KD'].abs().mean(),
        'mean_JSD': interface_positions['JSD'].mean(),
        'mean_delta_acidic_fraction': interface_positions['delta_acidic_fraction'].mean(),
        'mean_delta_basic_fraction': interface_positions['delta_basic_fraction'].mean(),
    }
    control_summary = {
        'category': 'CONTROL_NONINTERFACE_FRAMEWORK',
        'n_positions': len(control_candidates),
        'mean_delta_mean_KD': control_candidates['delta_mean_KD'].mean(),
        'mean_abs_delta_mean_KD': control_candidates['delta_mean_KD'].abs().mean(),
        'mean_JSD': control_candidates['JSD'].mean(),
        'mean_delta_acidic_fraction': control_candidates['delta_acidic_fraction'].mean(),
        'mean_delta_basic_fraction': control_candidates['delta_basic_fraction'].mean(),
    }
    summary_df = pd.DataFrame([interface_summary, control_summary])
    logger.info("\nInterface vs control descriptive summary:")
    logger.info(summary_df.to_string(index=False))

    return tetrad_rows, tetrad_vs_full_df, combined, summary_df


if __name__ == "__main__":
    main()
