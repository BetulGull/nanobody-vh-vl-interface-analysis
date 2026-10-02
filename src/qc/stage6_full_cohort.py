"""
Stage 6 — Full Cohort: CDR3 Shielding + Geometry for all 533 VHH Representatives
"""

import sys
from pathlib import Path
import pandas as pd
import numpy as np
import logging
import traceback

sys.path.insert(0, str(Path(__file__).parent.parent))

from qc.mmcif_qc import MMCIFStructure
from qc.stage3_variable_domain import extract_variable_domain_atoms
from qc.stage3_sasa import calc_sasa
from qc.stage6_cdr3_removal import extract_no_cdr3_atoms, get_shieldable_interface_positions, CDR3_START, CDR3_END
from qc.stage6_geometry import compute_cdr3_length, compute_gordon_rho_compactness, compute_zhang_fr2_cdr3_dis
from qc.stage4_sequence_source import SequenceSource, compute_ambiguous_base_positions
from qc.stage4_chemistry import THREE_TO_ONE, get_kd, residue_class_flags

PROJECT_ROOT = Path(__file__).parent.parent.parent
PROCESSED_DIR = PROJECT_ROOT / "data" / "processed"
CHECKPOINT_DIR = PROJECT_ROOT / "data" / "interim" / "stage6" / "checkpoints"
CHECKPOINT_DIR.mkdir(parents=True, exist_ok=True)
LOG_DIR = PROJECT_ROOT / "logs"
LOG_DIR.mkdir(exist_ok=True)

logging.basicConfig(
    level=logging.INFO, format='%(asctime)s - %(message)s',
    handlers=[logging.FileHandler(LOG_DIR / "stage6_shielding_geometry.log"), logging.StreamHandler()],
)
logger = logging.getLogger(__name__)

CANONICAL_20 = set('ARNDCQEGHILKMFPSTWYV')


def process_one_vhh(row: pd.Series, shieldable_positions_df: pd.DataFrame, ambig: set) -> dict:
    instance = row['INSTANCE']
    hchain = row['Hchain']
    pdb = row['PDB']
    sabdab_id = row['SABDAB_ID']
    seq_group_id = row.get('redundancy_group', sabdab_id)
    path = PROJECT_ROOT / row['extracted_structure_path']

    result = {
        'residue_rows': [],
        'summary': {
            'sequence_group_id': seq_group_id, 'SABDAB_ID': sabdab_id, 'INSTANCE': instance, 'PDB': pdb,
        },
        'geometry': {
            'sequence_group_id': seq_group_id, 'SABDAB_ID': sabdab_id, 'INSTANCE': instance, 'PDB': pdb,
        },
    }

    s = MMCIFStructure(str(path))
    model_id = s.get_model()

    # --- Geometry ---
    length_info = compute_cdr3_length(s, hchain, model_id)
    gordon = compute_gordon_rho_compactness(s, hchain, model_id, length_info['cdr3_sequence_length'])
    zhang = compute_zhang_fr2_cdr3_dis(s, hchain, model_id)
    result['geometry'].update(length_info)
    result['geometry'].update(gordon)
    result['geometry'].update(zhang)

    # --- Shielding ---
    native_atoms = extract_variable_domain_atoms(s, hchain, model_id)
    nocdr3_atoms = extract_no_cdr3_atoms(native_atoms)
    per_res_native, _ = calc_sasa(native_atoms)
    per_res_nocdr3, _ = calc_sasa(nocdr3_atoms)

    ss = SequenceSource(s, hchain, model_id, ambig)

    residue_rows = []
    for _, prow in shieldable_positions_df.iterrows():
        pos, ins = int(prow['imgt_position']), prow['insertion_code']
        proj = ss.project(pos, ins)
        state = proj['state']
        residue_present = state in ('PRESENT_OBSERVED', 'PRESENT_UNRESOLVED')
        coord_observed = state == 'PRESENT_OBSERVED'

        amino_acid = None
        kd = None
        positive_kd = None
        noncanonical = False
        sasa_native = None
        sasa_nocdr3 = None
        delta_shield = None
        sasa_warning = None

        if residue_present and proj['amino_acid'] is not None:
            amino_acid = THREE_TO_ONE.get(proj['amino_acid'], 'X')
            kd = get_kd(amino_acid)
            positive_kd = (kd is not None) and (kd > 0)
            noncanonical = amino_acid not in CANONICAL_20

        if coord_observed:
            key = (hchain, f"{pos}{ins}")
            sasa_native = per_res_native.get(key)
            sasa_nocdr3 = per_res_nocdr3.get(key)
            if sasa_native is not None and sasa_nocdr3 is not None:
                delta_shield = sasa_nocdr3 - sasa_native
            else:
                sasa_warning = "coordinate_observed_but_sasa_lookup_failed"
        else:
            sasa_warning = f"no_coordinates_state={state}"

        residue_rows.append({
            'sequence_group_id': seq_group_id, 'SABDAB_ID': sabdab_id, 'INSTANCE': instance, 'PDB': pdb, 'Hchain': hchain,
            'imgt_position': pos, 'imgt_number': pos, 'insertion_code': ins, 'region': prow['region'],
            'amino_acid': amino_acid, 'kyte_doolittle': kd,
            'interface_frequency_5': prow['interface_frequency_5'],
            'mean_delta_sasa_stage3': prow['mean_delta_sasa'], 'median_delta_sasa_stage3': prow['median_delta_sasa'],
            'sasa_native': sasa_native, 'sasa_noCDR3': sasa_nocdr3, 'delta_sasa_shield': delta_shield,
            'positiveKD': positive_kd,
            'coordinate_observed': coord_observed, 'chemistry_evaluable': (kd is not None) and coord_observed,
            'noncanonical_flag': noncanonical, 'sasa_warning': sasa_warning,
        })

    result['residue_rows'] = residue_rows
    return result


def aggregate_summary(residue_rows: list, seq_group_id, sabdab_id, instance, pdb, n_shieldable_total: int) -> dict:
    df = pd.DataFrame(residue_rows)
    evaluable = df[df['coordinate_observed'] == True].dropna(subset=['sasa_native', 'sasa_noCDR3'])
    n_present = (df['coordinate_observed'] | df['amino_acid'].notna()).sum()
    n_coord_observed = df['coordinate_observed'].sum()

    evaluable_kd = evaluable.dropna(subset=['kyte_doolittle'])

    SASA_interface_native = evaluable['sasa_native'].sum() if len(evaluable) else None
    SASA_interface_noCDR3 = evaluable['sasa_noCDR3'].sum() if len(evaluable) else None
    delta_sasa_shield = (SASA_interface_noCDR3 - SASA_interface_native) if (SASA_interface_native is not None) else None
    shielding_fraction = (delta_sasa_shield / SASA_interface_noCDR3) if (SASA_interface_noCDR3 and SASA_interface_noCDR3 > 0) else None

    w = evaluable['interface_frequency_5']
    weighted_native = (w * evaluable['sasa_native']).sum() if len(evaluable) else None
    weighted_nocdr3 = (w * evaluable['sasa_noCDR3']).sum() if len(evaluable) else None
    weighted_delta_shield_freq5 = (weighted_nocdr3 - weighted_native) if weighted_native is not None else None

    def kd_stats(col):
        pos_rows = evaluable_kd[evaluable_kd['positiveKD'] == True]
        pos_sasa = pos_rows[col].sum()
        total_sasa = evaluable_kd[col].sum()
        frac = pos_sasa / total_sasa if total_sasa and total_sasa > 0 else None
        weighted_kd = (evaluable_kd[col] * evaluable_kd['kyte_doolittle']).sum() / total_sasa if total_sasa and total_sasa > 0 else None
        kd_area_sum = (evaluable_kd[col] * evaluable_kd['kyte_doolittle']).sum() if len(evaluable_kd) else None
        return pos_sasa, frac, weighted_kd, kd_area_sum

    native_pos_sasa, native_frac, native_weighted_kd, native_kd_area = kd_stats('sasa_native')
    nocdr3_pos_sasa, nocdr3_frac, nocdr3_weighted_kd, nocdr3_kd_area = kd_stats('sasa_noCDR3')

    return {
        'sequence_group_id': seq_group_id, 'SABDAB_ID': sabdab_id, 'INSTANCE': instance, 'PDB': pdb,
        'n_shieldable_positions_total': n_shieldable_total,
        'n_shieldable_positions_present': int(n_present),
        'n_shieldable_positions_coordinate_observed': int(n_coord_observed),
        'shieldable_interface_coordinate_coverage': n_coord_observed / n_present if n_present else None,
        'SASA_interface_native': SASA_interface_native,
        'SASA_interface_noCDR3': SASA_interface_noCDR3,
        'delta_sasa_shield': delta_sasa_shield,
        'shielding_fraction': shielding_fraction,
        'weighted_delta_sasa_shield_freq5': weighted_delta_shield_freq5,
        'native_positiveKD_sasa': native_pos_sasa, 'noCDR3_positiveKD_sasa': nocdr3_pos_sasa,
        'delta_positiveKD_sasa_CDR3': (nocdr3_pos_sasa - native_pos_sasa) if (native_pos_sasa is not None and nocdr3_pos_sasa is not None) else None,
        'native_positiveKD_fraction': native_frac, 'noCDR3_positiveKD_fraction': nocdr3_frac,
        'delta_positiveKD_fraction_CDR3': (nocdr3_frac - native_frac) if (native_frac is not None and nocdr3_frac is not None) else None,
        'native_sasa_weighted_mean_KD': native_weighted_kd, 'noCDR3_sasa_weighted_mean_KD': nocdr3_weighted_kd,
        'delta_sasa_weighted_mean_KD_CDR3': (nocdr3_weighted_kd - native_weighted_kd) if (native_weighted_kd is not None and nocdr3_weighted_kd is not None) else None,
        'native_KD_area_sum': native_kd_area, 'noCDR3_KD_area_sum': nocdr3_kd_area,
        'delta_KD_area_sum_CDR3': (nocdr3_kd_area - native_kd_area) if (native_kd_area is not None and nocdr3_kd_area is not None) else None,
        'n_negative_delta_shield': (evaluable['delta_sasa_shield'] < 0).sum() if len(evaluable) else 0,
    }


def main():
    vhh_master = pd.read_parquet(PROCESSED_DIR / "stage4_vhh_all_instances_with_redundancy.parquet")
    vhh_master = vhh_master[vhh_master['is_primary_representative'] == True].reset_index(drop=True)

    primary_map = pd.read_parquet(PROCESSED_DIR / "stage3_empirical_vh_vl_interface_primary.parquet")
    full_interface, cdr3_interface, shieldable = get_shieldable_interface_positions(primary_map)

    residues_stage3 = pd.read_parquet(PROCESSED_DIR / "stage3_vh_vl_delta_sasa_residues.parquet")
    ambig = compute_ambiguous_base_positions(residues_stage3)

    n = len(vhh_master)
    logger.info(f"Processing {n} VHH representatives...")
    logger.info(f"Full interface: {len(full_interface)}, CDR3-in-interface: {len(cdr3_interface)}, shieldable: {len(shieldable)}")

    all_residue_rows = []
    all_summaries = []
    all_geometries = []
    n_errors = 0

    for idx, row in vhh_master.iterrows():
        try:
            result = process_one_vhh(row, shieldable, ambig)
            all_residue_rows.extend(result['residue_rows'])
            summary = aggregate_summary(
                result['residue_rows'], result['summary']['sequence_group_id'], result['summary']['SABDAB_ID'],
                result['summary']['INSTANCE'], result['summary']['PDB'], len(shieldable)
            )
            all_summaries.append(summary)
            all_geometries.append(result['geometry'])
        except Exception as e:
            n_errors += 1
            logger.error(f"Error on {row['INSTANCE']}: {e}")
            logger.error(traceback.format_exc())
            if n_errors >= 20:
                logger.error("STOPPING: too many errors, systematic issue suspected.")
                raise

        if (idx + 1) % 100 == 0 or (idx + 1) == n:
            pd.DataFrame(all_residue_rows).to_parquet(CHECKPOINT_DIR / f"residues_{idx+1:05d}.parquet", index=False)
            pd.DataFrame(all_summaries).to_parquet(CHECKPOINT_DIR / f"summaries_{idx+1:05d}.parquet", index=False)
            pd.DataFrame(all_geometries).to_parquet(CHECKPOINT_DIR / f"geometries_{idx+1:05d}.parquet", index=False)
            logger.info(f"Checkpoint: {idx+1}/{n} processed, {n_errors} errors")

    logger.info(f"Done. {n_errors} errors total.")
    return pd.DataFrame(all_residue_rows), pd.DataFrame(all_summaries), pd.DataFrame(all_geometries)


if __name__ == "__main__":
    main()
