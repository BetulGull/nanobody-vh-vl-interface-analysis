"""
Stage 5 — Per-Sequence Structural Surface Chemistry Descriptors
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


def compute_sequence_surface_chemistry(combined: pd.DataFrame) -> pd.DataFrame:
    rows = []

    for (group, structural_state, seq_group_id, instance), grp in combined.groupby(
        ['group', 'structural_state', 'sequence_group_id', 'INSTANCE']
    ):
        n_total = grp['primary_binary_interface'].sum()  # constant per full map, same for every sequence

        binary = grp[grp['primary_binary_interface'] == True]
        n_present = (binary['mapping_status'].isin(['PRESENT_OBSERVED', 'PRESENT_UNRESOLVED'])).sum()
        n_coord_observed = (binary['coordinate_observed'] == True).sum()
        # chemistry_evaluable == has a real SASA value (canonical/noncanonical alike -- SASA
        # itself doesn't require canonical identity, only coordinates)
        n_chemistry_evaluable = binary['sasa'].notna().sum()

        interface_sequence_coverage = n_present / n_total if n_total else None
        interface_coordinate_coverage_given_present = n_coord_observed / n_present if n_present else None
        interface_chemistry_coverage_given_present = n_chemistry_evaluable / n_present if n_present else None

        evaluable = binary[binary['sasa'].notna()]
        # KD-based stats need canonical KD too (noncanonical residues have KD=None)
        evaluable_kd = evaluable.dropna(subset=['kyte_doolittle'])

        interface_total_sasa = evaluable['sasa'].sum() if len(evaluable) else None
        interface_mean_sasa_per_residue = evaluable['sasa'].mean() if len(evaluable) else None

        positive_kd_rows = evaluable_kd[evaluable_kd['positive_KD'] == True]
        interface_positiveKD_sasa = positive_kd_rows['sasa'].sum() if len(evaluable_kd) else None
        canonical_sasa_sum = evaluable_kd['sasa'].sum()
        interface_positiveKD_sasa_fraction = (
            interface_positiveKD_sasa / canonical_sasa_sum if canonical_sasa_sum and canonical_sasa_sum > 0 else None
        )

        if len(evaluable_kd) > 0 and evaluable_kd['sasa'].sum() != 0:
            interface_sasa_weighted_mean_KD = (evaluable_kd['sasa'] * evaluable_kd['kyte_doolittle']).sum() / evaluable_kd['sasa'].sum()
        else:
            interface_sasa_weighted_mean_KD = None
        interface_KD_area_sum = (evaluable_kd['sasa'] * evaluable_kd['kyte_doolittle']).sum() if len(evaluable_kd) else None

        # Continuous weighted (freq5) descriptors -- over ALL present-with-SASA positions in full map
        present_all = grp[(grp['sasa'].notna())].dropna(subset=['kyte_doolittle'])
        w = present_all['interface_frequency_5']
        sasa = present_all['sasa']
        kd = present_all['kyte_doolittle']

        weighted_interface_total_sasa_freq5 = (w * sasa).sum() if len(present_all) else None
        pos_mask = present_all['positive_KD'] == True
        weighted_interface_positiveKD_sasa_freq5 = (w[pos_mask] * sasa[pos_mask]).sum() if pos_mask.any() else 0.0
        effective_exposed_weight_freq5 = (w * sasa).sum() if len(present_all) else None
        weighted_interface_positiveKD_fraction_freq5 = (
            weighted_interface_positiveKD_sasa_freq5 / effective_exposed_weight_freq5
            if effective_exposed_weight_freq5 and effective_exposed_weight_freq5 > 0 else None
        )
        weighted_interface_sasa_mean_KD_freq5 = (
            (w * sasa * kd).sum() / effective_exposed_weight_freq5
            if effective_exposed_weight_freq5 and effective_exposed_weight_freq5 > 0 else None
        )

        rows.append({
            'group': group, 'structural_state': structural_state,
            'sequence_group_id': seq_group_id, 'INSTANCE': instance,
            'SABDAB_ID': grp['SABDAB_ID'].iloc[0], 'PDB': grp['PDB'].iloc[0],
            'n_primary_interface_positions_total': n_total,
            'n_primary_interface_residues_present': n_present,
            'n_primary_interface_residues_coordinate_observed': n_coord_observed,
            'n_primary_interface_residues_chemistry_evaluable': n_chemistry_evaluable,
            'interface_sequence_coverage': interface_sequence_coverage,
            'interface_coordinate_coverage_given_present': interface_coordinate_coverage_given_present,
            'interface_chemistry_coverage_given_present': interface_chemistry_coverage_given_present,
            'interface_total_sasa': interface_total_sasa,
            'interface_mean_sasa_per_residue': interface_mean_sasa_per_residue,
            'interface_positiveKD_sasa': interface_positiveKD_sasa,
            'interface_positiveKD_sasa_fraction': interface_positiveKD_sasa_fraction,
            'interface_sasa_weighted_mean_KD': interface_sasa_weighted_mean_KD,
            'interface_KD_area_sum': interface_KD_area_sum,
            'weighted_interface_total_sasa_freq5': weighted_interface_total_sasa_freq5,
            'weighted_interface_positiveKD_sasa_freq5': weighted_interface_positiveKD_sasa_freq5,
            'weighted_interface_positiveKD_fraction_freq5': weighted_interface_positiveKD_fraction_freq5,
            'weighted_interface_sasa_mean_KD_freq5': weighted_interface_sasa_mean_KD_freq5,
            'effective_exposed_weight_freq5': effective_exposed_weight_freq5,
        })

    return pd.DataFrame(rows)


def main():
    combined = pd.read_parquet(PROCESSED_DIR / "stage5_interface_surface_chemistry_residues.parquet")
    seq_chem = compute_sequence_surface_chemistry(combined)

    # Join Stage-4 intrinsic descriptors where available (VHH and VH primary only)
    stage4_chem = pd.read_parquet(PROCESSED_DIR / "stage4_sequence_interface_chemistry.parquet")
    stage4_chem_renamed = stage4_chem.rename(columns={'representative_INSTANCE': 'INSTANCE'})[
        ['INSTANCE', 'H_intrinsic_binary_mean_KD', 'H_intrinsic_weighted_freq5']
    ]
    seq_chem = seq_chem.merge(stage4_chem_renamed, on='INSTANCE', how='left')

    out_parquet = PROCESSED_DIR / "stage5_sequence_surface_chemistry.parquet"
    out_csv = PROCESSED_DIR / "stage5_sequence_surface_chemistry.csv"
    seq_chem.to_parquet(out_parquet, index=False)
    seq_chem.to_csv(out_csv, index=False)
    logger.info(f"Wrote {out_parquet} / .csv: {len(seq_chem)} rows")

    logger.info("\ninterface_sasa_weighted_mean_KD by structural_state:")
    logger.info(seq_chem.groupby('structural_state')['interface_sasa_weighted_mean_KD'].describe().to_string())

    return seq_chem


if __name__ == "__main__":
    main()
