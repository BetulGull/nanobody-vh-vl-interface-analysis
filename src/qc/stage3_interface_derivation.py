"""
Stage 3 — Empirical VH-VL Interface Derivation Orchestrator

For every Stage-2 paired-primary structure, computes per-VH-residue ΔSASA
due to VL presence, checkpointed for resumability.
"""

import sys
from pathlib import Path
import logging
import traceback
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent.parent))

from qc.mmcif_qc import MMCIFStructure, mmCIFParseError
from qc.stage3_variable_domain import extract_variable_domain_atoms
from qc.stage3_sasa import calc_sasa

PROJECT_ROOT = Path(__file__).parent.parent.parent
CHECKPOINT_DIR = PROJECT_ROOT / "data" / "interim" / "stage3" / "checkpoints"
CHECKPOINT_DIR.mkdir(parents=True, exist_ok=True)
LOG_DIR = PROJECT_ROOT / "logs"
LOG_DIR.mkdir(exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(LOG_DIR / "stage3_interface_derivation.log"),
        logging.StreamHandler(),
    ],
)
logger = logging.getLogger(__name__)

IMGT_REGIONS = {
    'FR1': (1, 26), 'CDR1': (27, 38), 'FR2': (39, 55), 'CDR2': (56, 65),
    'FR3': (66, 104), 'CDR3': (105, 117), 'FR4': (118, 128),
}


def get_region(imgt_number: int) -> str:
    for region, (start, end) in IMGT_REGIONS.items():
        if start <= imgt_number <= end:
            return region
    return 'UNKNOWN'


def normalize_lchain(raw):
    if pd.isna(raw) or str(raw).strip() == 'NA':
        return None
    return raw


def process_one_structure(row: pd.Series) -> dict:
    """
    Returns dict with:
      'residue_rows': list of per-VH-residue dicts
      'structure_summary': dict of structure-level metrics
    """
    instance = row['INSTANCE']
    pdb = row['PDB']
    hchain = row['Hchain']
    lchain = normalize_lchain(row['Lchain'])
    sabdab_id = row['SABDAB_ID']
    redundancy_group = row.get('redundancy_group', sabdab_id)
    is_primary = row.get('is_primary_representative', False)

    summary = {
        'PDB': pdb, 'INSTANCE': instance, 'SABDAB_ID': sabdab_id,
        'redundancy_group': redundancy_group, 'is_primary_representative': is_primary,
        'resolution': row.get('resolution'),
        'total_vh_delta_sasa_due_to_vl': None,
        'n_interface_residues_1': None, 'n_interface_residues_5': None, 'n_interface_residues_10': None,
        'vh_positions_observed': None, 'fraction_vh_positions_with_interface_5': None,
        'negative_delta_count': None,
        'sasa_qc_pass': False,
        'failure_or_review_reason': None,
    }

    if lchain is None:
        summary['failure_or_review_reason'] = 'NO_LCHAIN_NOT_A_PAIRED_STRUCTURE'
        return {'residue_rows': [], 'structure_summary': summary}

    path = PROJECT_ROOT / row['extracted_structure_path']
    try:
        s = MMCIFStructure(str(path))
        model_id = s.get_model()
    except (mmCIFParseError, FileNotFoundError) as e:
        summary['failure_or_review_reason'] = f'PARSE_ERROR: {e}'
        return {'residue_rows': [], 'structure_summary': summary}

    vh_atoms = extract_variable_domain_atoms(s, hchain, model_id)
    vl_atoms = extract_variable_domain_atoms(s, lchain, model_id)

    if not vh_atoms:
        summary['failure_or_review_reason'] = 'NO_VH_ATOMS_EXTRACTED'
        return {'residue_rows': [], 'structure_summary': summary}
    if not vl_atoms:
        summary['failure_or_review_reason'] = 'NO_VL_ATOMS_EXTRACTED'
        return {'residue_rows': [], 'structure_summary': summary}

    per_res_alone, total_alone = calc_sasa(vh_atoms)
    per_res_with_vl, total_with_vl = calc_sasa(vh_atoms + vl_atoms)

    # Build per-residue rows: one per unique (auth_seq_id, ins_code) in vh_atoms
    residue_keys = sorted(set((a['auth_seq_id'], a['ins_code']) for a in vh_atoms))
    residue_meta = {}
    for a in vh_atoms:
        key = (a['auth_seq_id'], a['ins_code'])
        if key not in residue_meta:
            residue_meta[key] = a['residue_name']

    residue_rows = []
    total_delta = 0.0
    n_interface_1 = n_interface_5 = n_interface_10 = 0
    n_negative = 0

    for (auth_seq_id, ins_code) in residue_keys:
        resnum_str = f"{auth_seq_id}{ins_code}"
        key = (hchain, resnum_str)
        sasa_alone = per_res_alone.get(key)
        sasa_with_vl = per_res_with_vl.get(key)

        if sasa_alone is None or sasa_with_vl is None:
            continue

        delta = sasa_alone - sasa_with_vl
        total_delta += delta
        if delta < 0:
            n_negative += 1

        interface_1 = delta >= 1.0
        interface_5 = delta >= 5.0
        interface_10 = delta >= 10.0
        n_interface_1 += int(interface_1)
        n_interface_5 += int(interface_5)
        n_interface_10 += int(interface_10)

        residue_rows.append({
            'PDB': pdb, 'INSTANCE': instance, 'SABDAB_ID': sabdab_id,
            'redundancy_group': redundancy_group, 'is_primary_representative': is_primary,
            'Hchain': hchain, 'Lchain': lchain,
            'imgt_position': auth_seq_id, 'imgt_number': auth_seq_id, 'insertion_code': ins_code,
            'region': get_region(auth_seq_id),
            'residue_name': residue_meta[(auth_seq_id, ins_code)],
            'sasa_vh_alone': sasa_alone,
            'sasa_vh_with_vl': sasa_with_vl,
            'delta_sasa_vl_raw': delta,
            'interface_dSASA_1': interface_1,
            'interface_dSASA_5': interface_5,
            'interface_dSASA_10': interface_10,
            'atom_completeness_flag': True,
            'nonstandard_residue_flag': residue_meta[(auth_seq_id, ins_code)] not in
                {'ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS',
                 'MET','PHE','PRO','SER','THR','TRP','TYR','VAL'},
            'sasa_warning': None,
        })

    n_observed = len(residue_rows)
    summary['total_vh_delta_sasa_due_to_vl'] = total_delta
    summary['n_interface_residues_1'] = n_interface_1
    summary['n_interface_residues_5'] = n_interface_5
    summary['n_interface_residues_10'] = n_interface_10
    summary['vh_positions_observed'] = n_observed
    summary['fraction_vh_positions_with_interface_5'] = n_interface_5 / n_observed if n_observed else None
    summary['negative_delta_count'] = n_negative
    summary['sasa_qc_pass'] = n_observed > 0

    return {'residue_rows': residue_rows, 'structure_summary': summary}


def run_full_cohort(master: pd.DataFrame) -> tuple:
    all_residue_rows = []
    all_summaries = []
    n = len(master)
    systematic_errors = 0

    for idx, row in master.iterrows():
        try:
            result = process_one_structure(row)
        except Exception as e:
            logger.error(f"UNEXPECTED error on {row.get('INSTANCE','?')}: {e}")
            logger.error(traceback.format_exc())
            systematic_errors += 1
            result = {
                'residue_rows': [],
                'structure_summary': {
                    'PDB': row.get('PDB'), 'INSTANCE': row.get('INSTANCE'),
                    'SABDAB_ID': row.get('SABDAB_ID'),
                    'failure_or_review_reason': f'UNEXPECTED_EXCEPTION: {e}',
                    'sasa_qc_pass': False,
                },
            }
            if systematic_errors >= 20:
                logger.error("STOPPING: too many unexpected exceptions — systematic error suspected.")
                raise

        all_residue_rows.extend(result['residue_rows'])
        all_summaries.append(result['structure_summary'])

        if (idx + 1) % 100 == 0 or (idx + 1) == n:
            pd.DataFrame(all_residue_rows).to_parquet(
                CHECKPOINT_DIR / f"residues_{idx+1:05d}.parquet", index=False)
            pd.DataFrame(all_summaries).to_parquet(
                CHECKPOINT_DIR / f"summaries_{idx+1:05d}.parquet", index=False)
            logger.info(f"Checkpoint: {idx+1}/{n} structures processed")

    return pd.DataFrame(all_residue_rows), pd.DataFrame(all_summaries)


def main():
    logger.info("=" * 80)
    logger.info("STAGE 3 — EMPIRICAL VH-VL INTERFACE DERIVATION BY ΔSASA")
    logger.info("=" * 80)

    master = pd.read_parquet(PROJECT_ROOT / "data" / "processed" / "stage3_all_instances_with_redundancy.parquet")
    logger.info(f"Processing {len(master)} paired structures (all-instance sensitivity set)...")

    residue_df, summary_df = run_full_cohort(master)

    logger.info(f"Total residue rows: {len(residue_df)}")
    logger.info(f"Total structure summaries: {len(summary_df)}")
    logger.info(f"QC pass: {summary_df['sasa_qc_pass'].sum()}/{len(summary_df)}")

    return residue_df, summary_df


if __name__ == "__main__":
    main()
