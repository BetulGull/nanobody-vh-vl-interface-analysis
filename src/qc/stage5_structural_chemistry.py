"""
Stage 5 — Native 3D Solvent Exposure: Combined Residue Table

Reuses Stage-4's already-validated mapping-status long table (sequence
presence/absence/ambiguity, amino acid, KD, interface weights) as the join
backbone -- no reimplementation of that logic. Adds structural SASA on top:
  - VHH_NATIVE: newly computed (native VHH variable domain, antigen/other
    chains excluded, CDR3 intact) using the SAME validated Stage-3 engine.
  - VH_WITHOUT_VL_GEOMETRIC / VH_WITH_VL_PAIRED: reused directly from frozen
    Stage-3 per-residue SASA (sasa_vh_alone / sasa_vh_with_vl).

SASA only exists where coordinates exist: PRESENT_OBSERVED residues may get
a real SASA value; PRESENT_UNRESOLVED / ABSENT_CONFIRMED / AMBIGUOUS_UNRESOLVED
residues get sasa = NaN (never 0 -- 0 means "resolved with zero exposure",
not "no coordinates").
"""

import sys
from pathlib import Path
import pandas as pd
import numpy as np
import logging

sys.path.insert(0, str(Path(__file__).parent.parent))

from qc.mmcif_qc import MMCIFStructure
from qc.stage3_variable_domain import extract_variable_domain_atoms
from qc.stage3_sasa import calc_sasa

PROJECT_ROOT = Path(__file__).parent.parent.parent
PROCESSED_DIR = PROJECT_ROOT / "data" / "processed"
LOG_DIR = PROJECT_ROOT / "logs"
LOG_DIR.mkdir(exist_ok=True)

logging.basicConfig(
    level=logging.INFO, format='%(asctime)s - %(message)s',
    handlers=[logging.FileHandler(LOG_DIR / "stage5_surface_chemistry.log"), logging.StreamHandler()],
)
logger = logging.getLogger(__name__)


def compute_vhh_native_sasa(vhh_master: pd.DataFrame) -> dict:
    """Returns {INSTANCE: {(auth_seq_id, ins_code): sasa}} for every VHH representative."""
    sasa_lookup = {}
    n_errors = 0
    for _, row in vhh_master.iterrows():
        instance = row['INSTANCE']
        hchain = row['Hchain']
        path = PROJECT_ROOT / row['extracted_structure_path']
        try:
            s = MMCIFStructure(str(path))
            model_id = s.get_model()
            vhh_atoms = extract_variable_domain_atoms(s, hchain, model_id)
            per_res, _ = calc_sasa(vhh_atoms)
            # per_res keyed by (chain, resnum_str); convert to (auth_seq_id, ins_code)
            structure_lookup = {}
            for a in vhh_atoms:
                key_res = (hchain, f"{a['auth_seq_id']}{a['ins_code']}")
                if key_res in per_res:
                    structure_lookup[(a['auth_seq_id'], a['ins_code'])] = per_res[key_res]
            sasa_lookup[instance] = structure_lookup
        except Exception as e:
            n_errors += 1
            logger.error(f"Error computing VHH native SASA for {instance}: {e}")
            sasa_lookup[instance] = {}

    logger.info(f"VHH native SASA computed for {len(vhh_master)} representatives, {n_errors} errors")
    return sasa_lookup


def build_combined_table(long_df: pd.DataFrame, vhh_sasa_lookup: dict,
                          stage3_residues: pd.DataFrame, vhh_meta: pd.DataFrame,
                          vh_meta: pd.DataFrame) -> pd.DataFrame:
    vhh_meta_lookup = vhh_meta.set_index('INSTANCE')[['PDB', 'Hchain']].to_dict('index')
    vh_meta_lookup = vh_meta.set_index('INSTANCE')[['PDB', 'Hchain']].to_dict('index')

    # Stage-3 frozen VH residue lookup: (INSTANCE, auth_seq_id, ins_code) -> (alone, with_vl)
    stage3_lookup = {}
    for _, r in stage3_residues.iterrows():
        key = (r['INSTANCE'], r['imgt_position'], r['insertion_code'])
        stage3_lookup[key] = (r['sasa_vh_alone'], r['sasa_vh_with_vl'])

    rows = []
    for _, r in long_df.iterrows():
        group = r['group']
        instance = r['representative_INSTANCE']
        pos, ins = r['imgt_position'], r['insertion_code']
        coord_observed = r['coordinate_observed_if_available']

        if group == 'VHH':
            meta = vhh_meta_lookup.get(instance, {})
            sasa_native = vhh_sasa_lookup.get(instance, {}).get((pos, ins)) if coord_observed else None
            sasa_warning = None if coord_observed else f"no_coordinates_state={r['sequence_state']}"

            rows.append({
                'group': 'VHH', 'structural_state': 'VHH_NATIVE',
                'sequence_group_id': r['sequence_group_id'], 'SABDAB_ID': r['SABDAB_ID'],
                'INSTANCE': instance, 'PDB': meta.get('PDB'), 'Hchain': meta.get('Hchain'),
                'imgt_position': pos, 'imgt_number': pos, 'insertion_code': ins, 'region': r['region'],
                'mapping_status': r['sequence_state'], 'amino_acid': r['amino_acid'],
                'kyte_doolittle': r['kyte_doolittle'],
                'primary_binary_interface': r['primary_binary_interface'],
                'interface_frequency_5': r['interface_frequency_5'],
                'mean_delta_sasa': r['mean_delta_sasa'], 'median_delta_sasa': r['median_delta_sasa'],
                'coordinate_observed': coord_observed,
                'sasa': sasa_native,
                'positive_KD': r['residue_positive_KD'],
                'noncanonical_residue_flag': (r['amino_acid'] is not None) and (r['amino_acid'] not in set('ARNDCQEGHILKMFPSTWYV')),
                'sasa_warning': sasa_warning,
            })
        else:  # VH
            meta = vh_meta_lookup.get(instance, {})
            key = (instance, pos, ins)
            alone, with_vl = stage3_lookup.get(key, (None, None))

            for structural_state, sasa_val in [
                ('VH_WITHOUT_VL_GEOMETRIC', alone if coord_observed else None),
                ('VH_WITH_VL_PAIRED', with_vl if coord_observed else None),
            ]:
                sasa_warning = None if coord_observed else f"no_coordinates_state={r['sequence_state']}"
                rows.append({
                    'group': 'VH', 'structural_state': structural_state,
                    'sequence_group_id': r['sequence_group_id'], 'SABDAB_ID': r['SABDAB_ID'],
                    'INSTANCE': instance, 'PDB': meta.get('PDB'), 'Hchain': meta.get('Hchain'),
                    'imgt_position': pos, 'imgt_number': pos, 'insertion_code': ins, 'region': r['region'],
                    'mapping_status': r['sequence_state'], 'amino_acid': r['amino_acid'],
                    'kyte_doolittle': r['kyte_doolittle'],
                    'primary_binary_interface': r['primary_binary_interface'],
                    'interface_frequency_5': r['interface_frequency_5'],
                    'mean_delta_sasa': r['mean_delta_sasa'], 'median_delta_sasa': r['median_delta_sasa'],
                    'coordinate_observed': coord_observed,
                    'sasa': sasa_val,
                    'positive_KD': r['residue_positive_KD'],
                    'noncanonical_residue_flag': (r['amino_acid'] is not None) and (r['amino_acid'] not in set('ARNDCQEGHILKMFPSTWYV')),
                    'sasa_warning': sasa_warning,
                })

    return pd.DataFrame(rows)


def main():
    long_df = pd.read_parquet(PROCESSED_DIR / "stage4_interface_projection_long.parquet")
    stage3_residues = pd.read_parquet(PROCESSED_DIR / "stage3_vh_vl_delta_sasa_residues.parquet")
    vhh_meta = pd.read_parquet(PROCESSED_DIR / "stage4_vhh_sequence_representatives.parquet")
    vh_meta = pd.read_parquet(PROCESSED_DIR / "stage3_paired_representatives.parquet")

    logger.info(f"Computing native VHH SASA for {vhh_meta['INSTANCE'].nunique()} representatives...")
    vhh_sasa_lookup = compute_vhh_native_sasa(vhh_meta)

    logger.info("Building combined Stage-5 residue table...")
    combined = build_combined_table(long_df, vhh_sasa_lookup, stage3_residues, vhh_meta, vh_meta)

    out_path = PROCESSED_DIR / "stage5_interface_surface_chemistry_residues.parquet"
    combined.to_parquet(out_path, index=False)
    logger.info(f"Wrote {out_path}: {len(combined)} rows")

    logger.info("\nRow counts by structural_state:")
    logger.info(combined['structural_state'].value_counts().to_string())

    logger.info("\nSASA availability by structural_state (non-null count):")
    logger.info(combined.groupby('structural_state')['sasa'].apply(lambda x: x.notna().sum()).to_string())

    return combined


if __name__ == "__main__":
    main()
