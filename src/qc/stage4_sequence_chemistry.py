"""
Stage 4 — Per-Sequence H_intrinsic Chemistry Descriptors

Computed strictly from the primary binary empirical interface (28 positions)
and from the continuous weighted map. Absent/ambiguous positions never
contribute an amino acid or chemistry value; denominators reflect only
positions that actually contributed (never hidden).
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


def compute_sequence_chemistry(long_df: pd.DataFrame) -> pd.DataFrame:
    """One row per (group, sequence_group_id, representative_INSTANCE)."""
    rows = []

    n_primary_total = long_df[long_df['group'] == long_df['group'].iloc[0]]['primary_binary_interface'].sum() \
        if len(long_df) else 0
    # n_primary_total should be computed per full map, not per sequence subset;
    # recompute properly from a single sequence's full position set below.

    for (group, seq_group_id, rep_instance), grp in long_df.groupby(
        ['group', 'sequence_group_id', 'representative_INSTANCE']
    ):
        n_primary_positions_total = grp['primary_binary_interface'].sum()  # constant per map, same for every sequence
        binary_grp = grp[grp['primary_binary_interface'] == True]
        present_binary = binary_grp[binary_grp['residue_present_in_sequence'] == True]
        n_present_binary = len(present_binary)

        primary_interface_position_coverage = (
            n_present_binary / n_primary_positions_total if n_primary_positions_total else None
        )

        # Coordinate coverage among interface residues genuinely PRESENT in the sequence
        if n_present_binary > 0:
            n_coord_observed = present_binary['coordinate_observed_if_available'].sum()
            primary_interface_coordinate_coverage_given_present = n_coord_observed / n_present_binary
        else:
            primary_interface_coordinate_coverage_given_present = None

        kd_vals = present_binary['kyte_doolittle'].dropna()
        h_binary_mean_kd = kd_vals.mean() if len(kd_vals) else None
        h_binary_sum_kd = kd_vals.sum() if len(kd_vals) else None

        # Weighted (continuous) descriptors over ALL present positions in the full map (not just binary)
        present_all = grp[grp['residue_present_in_sequence'] == True]

        def weighted_kd(weight_col):
            sub = present_all.dropna(subset=[weight_col, 'kyte_doolittle'])
            w_sum = sub[weight_col].sum()
            if w_sum == 0 or len(sub) == 0:
                return None, 0.0
            return (sub[weight_col] * sub['kyte_doolittle']).sum() / w_sum, w_sum

        h_weighted_freq5, eff_w_freq5 = weighted_kd('interface_frequency_5')
        h_weighted_mean_dsasa, eff_w_mean_dsasa = weighted_kd('mean_delta_sasa')
        h_weighted_median_dsasa, eff_w_median_dsasa = weighted_kd('median_delta_sasa')

        # Residue-class fractions within the primary binary interface (present residues only)
        n_p = n_present_binary
        acidic_fraction = present_binary['residue_acidic'].sum() / n_p if n_p else None
        basic_fraction = present_binary['residue_basic'].sum() / n_p if n_p else None
        aromatic_fraction = present_binary['residue_aromatic'].sum() / n_p if n_p else None
        positive_kd_fraction = present_binary['residue_positive_KD'].sum() / n_p if n_p else None
        glycine_fraction = present_binary['residue_glycine'].sum() / n_p if n_p else None
        proline_fraction = present_binary['residue_proline'].sum() / n_p if n_p else None
        cysteine_fraction = present_binary['residue_cysteine'].sum() / n_p if n_p else None
        histidine_fraction = present_binary['residue_histidine'].sum() / n_p if n_p else None

        n_acidic = present_binary['residue_acidic'].sum() if n_p else 0
        n_basic = present_binary['residue_basic'].sum() if n_p else 0
        nominal_sidechain_charge = (n_basic - n_acidic) if n_p else None

        rows.append({
            'group': group, 'sequence_group_id': seq_group_id, 'representative_INSTANCE': rep_instance,
            'primary_interface_position_coverage': primary_interface_position_coverage,
            'primary_interface_coordinate_coverage_given_present': primary_interface_coordinate_coverage_given_present,
            'n_primary_positions_present': n_present_binary,
            'H_intrinsic_binary_mean_KD': h_binary_mean_kd,
            'H_intrinsic_binary_sum_KD': h_binary_sum_kd,
            'H_intrinsic_weighted_freq5': h_weighted_freq5,
            'effective_weight_sum_freq5': eff_w_freq5,
            'H_intrinsic_weighted_mean_dSASA': h_weighted_mean_dsasa,
            'effective_weight_sum_mean_dSASA': eff_w_mean_dsasa,
            'H_intrinsic_weighted_median_dSASA': h_weighted_median_dsasa,
            'effective_weight_sum_median_dSASA': eff_w_median_dsasa,
            'acidic_fraction': acidic_fraction,
            'basic_fraction': basic_fraction,
            'aromatic_fraction': aromatic_fraction,
            'positive_KD_fraction': positive_kd_fraction,
            'glycine_fraction': glycine_fraction,
            'proline_fraction': proline_fraction,
            'cysteine_fraction': cysteine_fraction,
            'histidine_fraction': histidine_fraction,
            'nominal_sidechain_charge': nominal_sidechain_charge,
        })

    return pd.DataFrame(rows)


def main():
    long_df = pd.read_parquet(PROCESSED_DIR / "stage4_interface_projection_long.parquet")

    # SABDAB_ID column also needed in output per spec
    chem_df = compute_sequence_chemistry(long_df)

    # attach SABDAB_ID (constant within sequence_group_id/representative_INSTANCE)
    sabdab_lookup = long_df.drop_duplicates(subset=['representative_INSTANCE'])[
        ['representative_INSTANCE', 'SABDAB_ID']
    ]
    chem_df = chem_df.merge(sabdab_lookup, on='representative_INSTANCE', how='left')

    cols_order = [
        'group', 'sequence_group_id', 'SABDAB_ID', 'representative_INSTANCE',
        'primary_interface_position_coverage', 'primary_interface_coordinate_coverage_given_present',
        'n_primary_positions_present',
        'H_intrinsic_binary_mean_KD', 'H_intrinsic_binary_sum_KD',
        'H_intrinsic_weighted_freq5', 'effective_weight_sum_freq5',
        'H_intrinsic_weighted_mean_dSASA', 'effective_weight_sum_mean_dSASA',
        'H_intrinsic_weighted_median_dSASA', 'effective_weight_sum_median_dSASA',
        'acidic_fraction', 'basic_fraction', 'aromatic_fraction', 'positive_KD_fraction',
        'glycine_fraction', 'proline_fraction', 'cysteine_fraction', 'histidine_fraction',
        'nominal_sidechain_charge',
    ]
    chem_df = chem_df[cols_order]

    out_parquet = PROCESSED_DIR / "stage4_sequence_interface_chemistry.parquet"
    out_csv = PROCESSED_DIR / "stage4_sequence_interface_chemistry.csv"
    chem_df.to_parquet(out_parquet, index=False)
    chem_df.to_csv(out_csv, index=False)
    logger.info(f"Wrote {out_parquet} / .csv: {len(chem_df)} rows")

    logger.info("\nH_intrinsic_binary_mean_KD by group:")
    logger.info(chem_df.groupby('group')['H_intrinsic_binary_mean_KD'].describe().to_string())

    # --- Coverage output ---
    coverage_out = TABLES_DIR / "stage4_interface_projection_coverage.csv"
    coverage_summary = chem_df.groupby('group')[
        ['primary_interface_position_coverage', 'primary_interface_coordinate_coverage_given_present']
    ].describe()
    coverage_summary.to_csv(coverage_out)
    logger.info(f"Wrote {coverage_out}")

    low_coverage = chem_df[chem_df['primary_interface_position_coverage'] < 0.5]
    logger.info(f"\nSequences with primary_interface_position_coverage < 0.5: {len(low_coverage)} (retained, not excluded)")

    return chem_df


if __name__ == "__main__":
    main()
