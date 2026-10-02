"""
Stage 4 — Project the Stage-3 Empirical Interface onto VHH and VH Sequences

Builds the long table: one row per (sequence representative) x (Stage-3
empirical-map IMGT position), for both VHH and conventional-VH groups.
"""

import sys
from pathlib import Path
import pandas as pd
import logging

sys.path.insert(0, str(Path(__file__).parent.parent))

from qc.mmcif_qc import MMCIFStructure
from qc.stage4_sequence_source import SequenceSource, compute_ambiguous_base_positions
from qc.stage4_chemistry import THREE_TO_ONE, get_kd, residue_class_flags

PROJECT_ROOT = Path(__file__).parent.parent.parent
PROCESSED_DIR = PROJECT_ROOT / "data" / "processed"
LOG_DIR = PROJECT_ROOT / "logs"
LOG_DIR.mkdir(exist_ok=True)

logging.basicConfig(
    level=logging.INFO, format='%(asctime)s - %(message)s',
    handlers=[logging.FileHandler(LOG_DIR / "stage4_interface_projection.log"), logging.StreamHandler()],
)
logger = logging.getLogger(__name__)


def project_group(master: pd.DataFrame, group_label: str, primary_map: pd.DataFrame,
                   ambig: set, all_instance: bool = False) -> pd.DataFrame:
    """
    master: dataframe of sequences to project (one row per structure to include)
    group_label: 'VHH' or 'VH'
    all_instance: if True, project every row (sensitivity); else only
                  is_primary_representative==True rows (primary map)
    """
    if not all_instance:
        master = master[master['is_primary_representative'] == True]

    rows = []
    n_errors = 0

    for _, row in master.iterrows():
        instance = row['INSTANCE']
        hchain = row['Hchain']
        sabdab_id = row['SABDAB_ID']
        redundancy_group = row.get('redundancy_group', sabdab_id)
        path = PROJECT_ROOT / row['extracted_structure_path']

        try:
            s = MMCIFStructure(str(path))
            model_id = s.get_model()
            ss = SequenceSource(s, hchain, model_id, ambig)
        except Exception as e:
            n_errors += 1
            logger.error(f"Error on {instance}: {e}")
            continue

        #28-position interface onto VHH/VH sequences; extract amino acid,
        # hydropathy, and residue class for each IMGT position
        for _, prow in primary_map.iterrows():
            pos = int(prow['imgt_position'])
            ins = prow['insertion_code']
            proj = ss.project(pos, ins)

            state = proj['state']
            residue_present = state in ('PRESENT_OBSERVED', 'PRESENT_UNRESOLVED')
            amino_acid = None
            kd = None
            class_flags = {
                'residue_acidic': None, 'residue_basic': None, 'residue_aromatic': None,
                'residue_positive_KD': None, 'residue_glycine': None, 'residue_proline': None,
                'residue_cysteine': None, 'residue_histidine': None,
            }
            if residue_present and proj['amino_acid'] is not None:
                amino_acid = THREE_TO_ONE.get(proj['amino_acid'], 'X')
                kd = get_kd(amino_acid)
                class_flags = residue_class_flags(amino_acid)

            rows.append({
                'group': group_label,
                'sequence_group_id': redundancy_group,
                'SABDAB_ID': sabdab_id,
                'representative_INSTANCE': instance,
                'imgt_position': pos, 'imgt_number': pos, 'insertion_code': ins,
                'region': prow['region'],
                'residue_present_in_sequence': residue_present,
                'sequence_state': state,
                'amino_acid': amino_acid,
                'primary_binary_interface': bool(prow['primary_binary_interface']),
                'interface_frequency_5': prow['interface_frequency_5'],
                'mean_delta_sasa': prow['mean_delta_sasa'],
                'median_delta_sasa': prow['median_delta_sasa'],
                'observation_fraction': prow['observation_fraction'],
                'vhh_tetrad_position': bool(prow['vhh_tetrad_position']),
                'kyte_doolittle': kd,
                **class_flags,
                'coordinate_observed_if_available': proj['coordinate_observed'],
            })

    logger.info(f"{group_label} ({'all-instance' if all_instance else 'primary'}): "
                f"{master['INSTANCE'].nunique()} sequences, {n_errors} errors, {len(rows)} rows")
    return pd.DataFrame(rows)


def main():
    residues = pd.read_parquet(PROCESSED_DIR / "stage3_vh_vl_delta_sasa_residues.parquet")
    ambig = compute_ambiguous_base_positions(residues)
    primary_map = pd.read_parquet(PROCESSED_DIR / "stage3_empirical_vh_vl_interface_primary.parquet")

    vhh_master = pd.read_parquet(PROCESSED_DIR / "stage4_vhh_all_instances_with_redundancy.parquet")
    vh_master = pd.read_parquet(PROCESSED_DIR / "stage3_all_instances_with_redundancy.parquet")

    logger.info("Projecting PRIMARY (nonredundant representative) sets...")
    vhh_primary_long = project_group(vhh_master, 'VHH', primary_map, ambig, all_instance=False)
    vh_primary_long = project_group(vh_master, 'VH', primary_map, ambig, all_instance=False)

    primary_long = pd.concat([vhh_primary_long, vh_primary_long], ignore_index=True)
    out_path = PROCESSED_DIR / "stage4_interface_projection_long.parquet"
    primary_long.to_parquet(out_path, index=False)
    logger.info(f"Wrote {out_path}: {len(primary_long)} rows")

    logger.info("Projecting ALL-INSTANCE sensitivity sets...")
    vhh_all_long = project_group(vhh_master, 'VHH', primary_map, ambig, all_instance=True)
    vh_all_long = project_group(vh_master, 'VH', primary_map, ambig, all_instance=True)
    all_instance_long = pd.concat([vhh_all_long, vh_all_long], ignore_index=True)
    all_out_path = PROCESSED_DIR / "stage4_interface_projection_long_all_instances.parquet"
    all_instance_long.to_parquet(all_out_path, index=False)
    logger.info(f"Wrote {all_out_path}: {len(all_instance_long)} rows")

    return primary_long, all_instance_long


if __name__ == "__main__":
    main()
