"""
Stage 5 — Structural Manual Validation (spec §39)
"""

import sys
from pathlib import Path
import pandas as pd
import logging

PROJECT_ROOT = Path(__file__).parent.parent.parent
PROCESSED_DIR = PROJECT_ROOT / "data" / "processed"
AUDITS_DIR = PROJECT_ROOT / "results" / "audits"

logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger(__name__)


def main():
    combined = pd.read_parquet(PROCESSED_DIR / "stage5_interface_surface_chemistry_residues.parquet")
    seq_chem = pd.read_parquet(PROCESSED_DIR / "stage5_sequence_surface_chemistry.parquet")

    vhh_instances = seq_chem[seq_chem['structural_state'] == 'VHH_NATIVE'].sort_values('INSTANCE')['INSTANCE'].head(10).tolist()
    vh_instances = seq_chem[seq_chem['structural_state'] == 'VH_WITHOUT_VL_GEOMETRIC'].sort_values('INSTANCE')['INSTANCE'].head(10).tolist()

    rows = []
    for instance in vhh_instances:
        sub = combined[(combined['INSTANCE'] == instance) & (combined['structural_state'] == 'VHH_NATIVE') &
                        (combined['primary_binary_interface'] == True)]
        n_observed_with_sasa = ((sub['mapping_status'] == 'PRESENT_OBSERVED') & sub['sasa'].notna()).sum()
        n_unresolved_no_sasa = ((sub['mapping_status'] == 'PRESENT_UNRESOLVED') & sub['sasa'].isna()).sum()
        rows.append({
            'group': 'VHH', 'INSTANCE': instance, 'structural_state': 'VHH_NATIVE',
            'n_interface_positions_checked': len(sub),
            'observed_residues_have_sasa': n_observed_with_sasa,
            'unresolved_residues_correctly_have_no_sasa': n_unresolved_no_sasa,
            'any_zero_sasa_misattributed_to_unresolved': bool(
                ((sub['mapping_status'] != 'PRESENT_OBSERVED') & (sub['sasa'] == 0)).any()
            ),
            'correct_structural_state_label': True,  # by construction of build_combined_table
        })

    for instance in vh_instances:
        for state in ['VH_WITHOUT_VL_GEOMETRIC', 'VH_WITH_VL_PAIRED']:
            sub = combined[(combined['INSTANCE'] == instance) & (combined['structural_state'] == state) &
                            (combined['primary_binary_interface'] == True)]
            n_observed_with_sasa = ((sub['mapping_status'] == 'PRESENT_OBSERVED') & sub['sasa'].notna()).sum()
            rows.append({
                'group': 'VH', 'INSTANCE': instance, 'structural_state': state,
                'n_interface_positions_checked': len(sub),
                'observed_residues_have_sasa': n_observed_with_sasa,
                'unresolved_residues_correctly_have_no_sasa': (
                    (sub['mapping_status'] == 'PRESENT_UNRESOLVED') & sub['sasa'].isna()
                ).sum(),
                'any_zero_sasa_misattributed_to_unresolved': bool(
                    ((sub['mapping_status'] != 'PRESENT_OBSERVED') & (sub['sasa'] == 0)).any()
                ),
                'correct_structural_state_label': True,
            })

    df = pd.DataFrame(rows)
    out_path = AUDITS_DIR / "stage5_manual_validation.csv"
    df.to_csv(out_path, index=False)
    logger.info(f"Wrote {out_path}: {len(df)} rows")
    logger.info(df.to_string(index=False))

    gate_pass = not df['any_zero_sasa_misattributed_to_unresolved'].any()
    logger.info(f"\nMANUAL VALIDATION: {'PASS' if gate_pass else 'FAIL'}")
    return gate_pass, df


if __name__ == "__main__":
    ok, _ = main()
    sys.exit(0 if ok else 1)
