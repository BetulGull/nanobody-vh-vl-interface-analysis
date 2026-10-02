"""
Stage 4 — Secondary VHH Cohort Sensitivity and All-Instance Redundancy Sensitivity
"""

import sys
from pathlib import Path
import pandas as pd
import numpy as np
import logging

sys.path.insert(0, str(Path(__file__).parent.parent))

from qc.mmcif_qc import MMCIFStructure
from qc.stage4_sequence_source import SequenceSource, compute_ambiguous_base_positions
from qc.stage4_chemistry import THREE_TO_ONE, get_kd, residue_class_flags
from qc.stage4_sequence_chemistry import compute_sequence_chemistry
from qc.stage4_projection import project_group

PROJECT_ROOT = Path(__file__).parent.parent.parent
PROCESSED_DIR = PROJECT_ROOT / "data" / "processed"
TABLES_DIR = PROJECT_ROOT / "results" / "tables"

logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger(__name__)


def main():
    residues = pd.read_parquet(PROCESSED_DIR / "stage3_vh_vl_delta_sasa_residues.parquet")
    ambig = compute_ambiguous_base_positions(residues)
    primary_map = pd.read_parquet(PROCESSED_DIR / "stage3_empirical_vh_vl_interface_primary.parquet")

    # --- Secondary VHH cohort ---
    secondary_master = pd.read_parquet(PROCESSED_DIR / "stage2_vhh_secondary_structural.parquet")
    secondary_master = secondary_master.copy()
    secondary_master['redundancy_group'] = secondary_master['SABDAB_ID']
    secondary_master['is_primary_representative'] = True  # treat every row as its own "representative" for this sensitivity-only table (no dedup claim made)

    secondary_long = project_group(secondary_master, 'VHH_SECONDARY', primary_map, ambig, all_instance=True)
    secondary_chem = compute_sequence_chemistry(secondary_long.rename(columns={'group': 'group'}))

    secondary_out = TABLES_DIR / "stage4_vhh_secondary_sensitivity.csv"
    secondary_chem['provenance_tier'] = 'secondary_structural_sensitivity'
    secondary_chem.to_csv(secondary_out, index=False)
    logger.info(f"Wrote {secondary_out}: {len(secondary_chem)} rows")
    if len(secondary_chem) > 0:
        logger.info(secondary_chem['H_intrinsic_binary_mean_KD'].describe().to_string())

    # --- Primary representative vs all-instance sensitivity ---
    primary_chem = pd.read_parquet(PROCESSED_DIR / "stage4_sequence_interface_chemistry.parquet")

    all_instance_long = pd.read_parquet(PROCESSED_DIR / "stage4_interface_projection_long_all_instances.parquet")
    all_instance_chem = compute_sequence_chemistry(all_instance_long)

    sensitivity_rows = []
    for group in ['VHH', 'VH']:
        prim = primary_chem[primary_chem['group'] == group]
        alln = all_instance_chem[all_instance_chem['group'] == group]
        sensitivity_rows.append({
            'group': group,
            'n_primary_representatives': len(prim),
            'n_all_instance': len(alln),
            'mean_H_intrinsic_binary_mean_KD_primary': prim['H_intrinsic_binary_mean_KD'].mean(),
            'mean_H_intrinsic_binary_mean_KD_all_instance': alln['H_intrinsic_binary_mean_KD'].mean(),
            'mean_H_intrinsic_weighted_freq5_primary': prim['H_intrinsic_weighted_freq5'].mean(),
            'mean_H_intrinsic_weighted_freq5_all_instance': alln['H_intrinsic_weighted_freq5'].mean(),
        })
    sensitivity_df = pd.DataFrame(sensitivity_rows)
    sensitivity_out = TABLES_DIR / "stage4_primary_vs_all_instance_sensitivity.csv"
    sensitivity_df.to_csv(sensitivity_out, index=False)
    logger.info(f"\nWrote {sensitivity_out}")
    logger.info(sensitivity_df.to_string(index=False))

    return secondary_chem, sensitivity_df


if __name__ == "__main__":
    main()
