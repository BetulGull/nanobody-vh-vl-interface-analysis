"""
Stage 4 — Mapping-Missingness Bias Audit

Tests whether the observed VHH-vs-VH chemistry divergence at the empirical
interface could be explained by unequal sequence-mapping confidence between
groups, rather than genuine amino-acid chemistry differences.
"""

import sys
from pathlib import Path
import pandas as pd
import numpy as np
import logging

PROJECT_ROOT = Path(__file__).parent.parent.parent
PROCESSED_DIR = PROJECT_ROOT / "data" / "processed"
TABLES_DIR = PROJECT_ROOT / "results" / "tables"
AUDITS_DIR = PROJECT_ROOT / "results" / "audits"

logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger(__name__)

CONFIDENT_STATES = {'PRESENT_OBSERVED', 'PRESENT_UNRESOLVED'}


def main():
    long_df = pd.read_parquet(PROCESSED_DIR / "stage4_interface_projection_long.parquet")
    interface_long = long_df[long_df['primary_binary_interface'] == True].copy()

    # --- Section 1+4: per-position mapping status by group ---
    rows = []
    for (pos, ins), pos_grp in interface_long.groupby(['imgt_position', 'insertion_code']):
        row = {'imgt_position': pos, 'insertion_code': ins}
        for group in ['VHH', 'VH']:
            g = pos_grp[pos_grp['group'] == group]
            n_total = len(g)
            state_counts = g['sequence_state'].value_counts()
            n_observed = state_counts.get('PRESENT_OBSERVED', 0)
            n_unresolved = state_counts.get('PRESENT_UNRESOLVED', 0)
            n_absent = state_counts.get('ABSENT_CONFIRMED', 0)
            n_ambiguous = state_counts.get('AMBIGUOUS_UNRESOLVED', 0)
            n_confident = n_observed + n_unresolved
            confident_fraction = n_confident / n_total if n_total else None

            row[f'n_total_{group}'] = n_total
            row[f'n_present_observed_{group}'] = n_observed
            row[f'n_present_unresolved_{group}'] = n_unresolved
            row[f'n_absent_confirmed_{group}'] = n_absent
            row[f'n_ambiguous_unresolved_{group}'] = n_ambiguous
            row[f'n_confident_residue_{group}'] = n_confident
            row[f'confident_mapping_fraction_{group}'] = confident_fraction
            row[f'ambiguous_fraction_{group}'] = n_ambiguous / n_total if n_total else None

        row['delta_confident_mapping_fraction'] = (
            row['confident_mapping_fraction_VHH'] - row['confident_mapping_fraction_VH']
        )
        rows.append(row)

    status_df = pd.DataFrame(rows).sort_values(['imgt_position', 'insertion_code']).reset_index(drop=True)
    status_out = AUDITS_DIR / "stage4_interface_mapping_status_by_position.csv"
    status_df.to_csv(status_out, index=False)
    logger.info(f"Wrote {status_out}: {len(status_df)} positions")

    # --- Section 2: differential mapping bias ---
    imbalance_005 = status_df[status_df['delta_confident_mapping_fraction'].abs() >= 0.05]
    imbalance_010 = status_df[status_df['delta_confident_mapping_fraction'].abs() >= 0.10]
    logger.info(f"\nPositions with |delta confident mapping fraction| >= 0.05: {len(imbalance_005)}")
    logger.info(f"Positions with |delta confident mapping fraction| >= 0.10: {len(imbalance_010)}")
    if len(imbalance_005) > 0:
        logger.info(imbalance_005[['imgt_position', 'insertion_code',
                                     'confident_mapping_fraction_VHH', 'confident_mapping_fraction_VH',
                                     'delta_confident_mapping_fraction']].to_string(index=False))

    # Overall ambiguity rate within the 28-position primary interface, per group
    vhh_ambig_rate = status_df['n_ambiguous_unresolved_VHH'].sum() / status_df['n_total_VHH'].sum()
    vh_ambig_rate = status_df['n_ambiguous_unresolved_VH'].sum() / status_df['n_total_VH'].sum()
    logger.info(f"\nOverall VHH ambiguity rate within 28 primary interface positions: {vhh_ambig_rate:.4f}")
    logger.info(f"Overall VH ambiguity rate within 28 primary interface positions: {vh_ambig_rate:.4f}")

    # --- Section 3: high-coverage sensitivity ---
    chem_df = pd.read_parquet(PROCESSED_DIR / "stage4_sequence_interface_chemistry.parquet")
    # primary_interface_position_coverage already = n_present_binary / n_primary_positions_total,
    # exactly the "confident evaluable coverage" definition (n_present_binary only counts
    # residue_present_in_sequence==True, i.e. PRESENT_OBSERVED + PRESENT_UNRESOLVED).
    high_cov = chem_df[chem_df['primary_interface_position_coverage'] >= 0.90]

    sensitivity_rows = []
    for group in ['VHH', 'VH']:
        full = chem_df[chem_df['group'] == group]
        hc = high_cov[high_cov['group'] == group]
        sensitivity_rows.append({
            'group': group,
            'n_full': len(full), 'n_high_coverage_subset': len(hc),
            'fraction_retained': len(hc) / len(full) if len(full) else None,
            'H_intrinsic_binary_mean_KD_full': full['H_intrinsic_binary_mean_KD'].mean(),
            'H_intrinsic_binary_mean_KD_high_coverage': hc['H_intrinsic_binary_mean_KD'].mean(),
            'H_intrinsic_weighted_freq5_full': full['H_intrinsic_weighted_freq5'].mean(),
            'H_intrinsic_weighted_freq5_high_coverage': hc['H_intrinsic_weighted_freq5'].mean(),
            'n_evaluable_residues_median_full': full['n_primary_positions_present'].median(),
            'n_evaluable_residues_median_high_coverage': hc['n_primary_positions_present'].median(),
            'effective_weight_sum_freq5_median_full': full['effective_weight_sum_freq5'].median(),
            'effective_weight_sum_freq5_median_high_coverage': hc['effective_weight_sum_freq5'].median(),
        })
    sensitivity_df = pd.DataFrame(sensitivity_rows)
    sensitivity_out = TABLES_DIR / "stage4_high_coverage_sensitivity.csv"
    sensitivity_df.to_csv(sensitivity_out, index=False)
    logger.info(f"\nWrote {sensitivity_out}")
    logger.info(sensitivity_df.to_string(index=False))

    # --- Section 5: noncanonical residues ---
    CANONICAL_20 = set('ARNDCQEGHILKMFPSTWYV')
    interface_present = interface_long[interface_long['residue_present_in_sequence'] == True]
    noncanonical = interface_present[~interface_present['amino_acid'].isin(CANONICAL_20)]
    n_noncanonical = len(noncanonical)
    logger.info(f"\nNoncanonical residue observations within primary interface: {n_noncanonical}")
    if n_noncanonical > 0:
        logger.info(noncanonical[['group', 'representative_INSTANCE', 'imgt_position',
                                    'insertion_code', 'amino_acid']].to_string(index=False))
        flagged_positions = noncanonical[['imgt_position', 'insertion_code']].drop_duplicates()
        logger.info(f"Primary interface positions flagged for noncanonical residue occurrence: {len(flagged_positions)}")

    # --- Section 6: interface-vs-control coverage check ---
    control_df = pd.read_csv(TABLES_DIR / "stage4_interface_vs_noninterface_framework.csv")

    def confident_fraction_for_positions(pos_ins_pairs, group):
        sub = long_df[
            (long_df['group'] == group) &
            (long_df.set_index(['imgt_position', 'insertion_code']).index.isin(pos_ins_pairs))
        ]
        n_total = len(sub)
        n_confident = sub['sequence_state'].isin(CONFIDENT_STATES).sum()
        return n_confident / n_total if n_total else None

    interface_pairs = list(zip(
        control_df[control_df['category'] == 'EMPIRICAL_INTERFACE']['imgt_position'],
        control_df[control_df['category'] == 'EMPIRICAL_INTERFACE']['insertion_code'].fillna(''),
    ))
    control_pairs = list(zip(
        control_df[control_df['category'] == 'CONTROL_NONINTERFACE_FRAMEWORK']['imgt_position'],
        control_df[control_df['category'] == 'CONTROL_NONINTERFACE_FRAMEWORK']['insertion_code'].fillna(''),
    ))

    coverage_check_rows = []
    for label, pairs in [('EMPIRICAL_INTERFACE', interface_pairs), ('CONTROL_NONINTERFACE_FRAMEWORK', control_pairs)]:
        for group in ['VHH', 'VH']:
            frac = confident_fraction_for_positions(pairs, group)
            coverage_check_rows.append({'category': label, 'group': group, 'mean_confident_mapping_fraction': frac})
    coverage_check_df = pd.DataFrame(coverage_check_rows)
    coverage_check_out = TABLES_DIR / "stage4_interface_vs_control_coverage_check.csv"
    coverage_check_df.to_csv(coverage_check_out, index=False)
    logger.info(f"\nWrote {coverage_check_out}")
    logger.info(coverage_check_df.to_string(index=False))

    return status_df, sensitivity_df, coverage_check_df, n_noncanonical


if __name__ == "__main__":
    main()
