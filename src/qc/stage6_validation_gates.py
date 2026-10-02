"""
Stage 6 — Validation Gates: CDR3-Removal Correctness + Translated-CDR3 Control
(spec §14-15)
"""

import sys
from pathlib import Path
import copy
import pandas as pd
import logging

sys.path.insert(0, str(Path(__file__).parent.parent))

from qc.mmcif_qc import MMCIFStructure
from qc.stage3_variable_domain import extract_variable_domain_atoms
from qc.stage3_sasa import calc_sasa
from qc.stage6_cdr3_removal import extract_no_cdr3_atoms, CDR3_START, CDR3_END, get_shieldable_interface_positions

PROJECT_ROOT = Path(__file__).parent.parent.parent
AUDITS_DIR = PROJECT_ROOT / "results" / "audits"

logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger(__name__)

TRANSLATE_DISTANCE = 100000.0
TRANSLATED_TOLERANCE = 0.5  # Å²


def validate_removal(instance: str, master: pd.DataFrame, category: str) -> dict:
    row = master[master['INSTANCE'] == instance].iloc[0]
    hchain = row['Hchain']
    path = PROJECT_ROOT / row['extracted_structure_path']

    s = MMCIFStructure(str(path))
    model_id = s.get_model()
    native_atoms = extract_variable_domain_atoms(s, hchain, model_id)
    nocdr3_atoms = extract_no_cdr3_atoms(native_atoms)

    n_cdr3_in_native = sum(1 for a in native_atoms if CDR3_START <= a['auth_seq_id'] <= CDR3_END)
    n_cdr3_in_nocdr3 = sum(1 for a in nocdr3_atoms if CDR3_START <= a['auth_seq_id'] <= CDR3_END)
    n_noncdr3_native = len(native_atoms) - n_cdr3_in_native
    n_noncdr3_removed_state = len(nocdr3_atoms)

    # Coordinate identity of retained (non-CDR3) atoms
    native_noncdr3 = [a for a in native_atoms if not (CDR3_START <= a['auth_seq_id'] <= CDR3_END)]
    coords_identical = (
        len(native_noncdr3) == len(nocdr3_atoms) and
        all((a['x'], a['y'], a['z']) == (b['x'], b['y'], b['z']) for a, b in zip(native_noncdr3, nocdr3_atoms))
    )

    per_res_native, total_native = calc_sasa(native_atoms)
    per_res_nocdr3, total_nocdr3 = calc_sasa(nocdr3_atoms)
    per_res_native_2, total_native_2 = calc_sasa(native_atoms)
    per_res_nocdr3_2, total_nocdr3_2 = calc_sasa(nocdr3_atoms)

    native_reproducible = abs(total_native - total_native_2) < 1e-9
    nocdr3_reproducible = abs(total_nocdr3 - total_nocdr3_2) < 1e-9

    return {
        'category': category, 'INSTANCE': instance, 'Hchain': hchain,
        'n_cdr3_atoms_in_native': n_cdr3_in_native,
        'n_cdr3_atoms_remaining_in_nocdr3_state': n_cdr3_in_nocdr3,
        'n_noncdr3_atoms_native': n_noncdr3_native,
        'n_atoms_in_nocdr3_state': n_noncdr3_removed_state,
        'no_atoms_outside_cdr3_removed': n_noncdr3_native == n_noncdr3_removed_state,
        'zero_cdr3_atoms_in_nocdr3_state': n_cdr3_in_nocdr3 == 0,
        'retained_coords_identical': coords_identical,
        'native_sasa_reproducible': native_reproducible,
        'nocdr3_sasa_reproducible': nocdr3_reproducible,
        'total_native_sasa': total_native,
        'total_nocdr3_sasa': total_nocdr3,
    }


def validate_translated_control(instance: str, master: pd.DataFrame, shieldable_positions: set) -> dict:
    row = master[master['INSTANCE'] == instance].iloc[0]
    hchain = row['Hchain']
    path = PROJECT_ROOT / row['extracted_structure_path']

    s = MMCIFStructure(str(path))
    model_id = s.get_model()
    native_atoms = extract_variable_domain_atoms(s, hchain, model_id)
    nocdr3_atoms = extract_no_cdr3_atoms(native_atoms)

    cdr3_atoms = [a for a in native_atoms if CDR3_START <= a['auth_seq_id'] <= CDR3_END]
    noncdr3_atoms = [a for a in native_atoms if not (CDR3_START <= a['auth_seq_id'] <= CDR3_END)]

    cdr3_translated = copy.deepcopy(cdr3_atoms)
    for a in cdr3_translated:
        a['x'] += TRANSLATE_DISTANCE

    combined_translated = noncdr3_atoms + cdr3_translated

    per_res_nocdr3, _ = calc_sasa(nocdr3_atoms)
    per_res_translated, _ = calc_sasa(combined_translated)

    diffs = []
    for pos, ins in shieldable_positions:
        key = (hchain, f"{pos}{ins}")
        v1 = per_res_nocdr3.get(key)
        v2 = per_res_translated.get(key)
        if v1 is not None and v2 is not None:
            diffs.append(abs(v1 - v2))

    max_diff = max(diffs) if diffs else None
    mean_diff = sum(diffs) / len(diffs) if diffs else None

    return {
        'INSTANCE': instance, 'Hchain': hchain,
        'n_shieldable_positions_compared': len(diffs),
        'max_abs_diff_nocdr3_vs_translated': max_diff,
        'mean_abs_diff_nocdr3_vs_translated': mean_diff,
        'control_pass': (max_diff is not None) and (max_diff < TRANSLATED_TOLERANCE),
    }


def main():
    vhh_master = pd.read_parquet(PROJECT_ROOT / "data/processed/stage4_vhh_sequence_representatives.parquet")
    primary_map = pd.read_parquet(PROJECT_ROOT / "data/processed/stage3_empirical_vh_vl_interface_primary.parquet")
    _, _, shieldable = get_shieldable_interface_positions(primary_map)
    shieldable_positions = set(zip(shieldable['imgt_position'], shieldable['insertion_code']))

    cases = [
        ('pdb_00003k3q-A', 'long_CDR3_insertion_rich'),
        ('pdb_00009imx-C', 'low_coverage_short_CDR3'),
        ('pdb_00006zbp-FFF', 'unusual_multichar_chain'),
        ('pdb_00006ul6-B1', 'unusual_multichar_chain_2'),
        ('pdb_00008qf4-E', 'high_resolution'),
        ('pdb_00001zmy-A', 'low_resolution'),
        ('pdb_00007m1h-G', 'low_coverage_2'),
        ('pdb_00003zkx-C', 'incomplete_coordinate_case'),
        ('pdb_00006oca-D', 'long_CDR3_2'),
        ('pdb_00009gv3-E', 'long_CDR3_3'),
    ]

    removal_results = [validate_removal(inst, vhh_master, cat) for inst, cat in cases]
    removal_df = pd.DataFrame(removal_results)
    removal_out = AUDITS_DIR / "stage6_cdr3_removal_validation.csv"
    removal_df.to_csv(removal_out, index=False)
    logger.info(f"Wrote {removal_out}")
    logger.info(removal_df[['category', 'INSTANCE', 'zero_cdr3_atoms_in_nocdr3_state',
                              'no_atoms_outside_cdr3_removed', 'retained_coords_identical',
                              'native_sasa_reproducible', 'nocdr3_sasa_reproducible']].to_string(index=False))

    removal_gate = (
        removal_df['zero_cdr3_atoms_in_nocdr3_state'].all() and
        removal_df['no_atoms_outside_cdr3_removed'].all() and
        removal_df['retained_coords_identical'].all() and
        removal_df['native_sasa_reproducible'].all() and
        removal_df['nocdr3_sasa_reproducible'].all()
    )
    logger.info(f"\nCDR3 REMOVAL VALIDATION GATE: {'PASS' if removal_gate else 'FAIL'}")

    translated_results = [validate_translated_control(inst, vhh_master, shieldable_positions) for inst, _ in cases[:6]]
    translated_df = pd.DataFrame(translated_results)
    translated_out = AUDITS_DIR / "stage6_translated_cdr3_control.csv"
    translated_df.to_csv(translated_out, index=False)
    logger.info(f"\nWrote {translated_out}")
    logger.info(translated_df.to_string(index=False))

    translated_gate = translated_df['control_pass'].all()
    logger.info(f"\nTRANSLATED-CDR3 CONTROL GATE: {'PASS' if translated_gate else 'FAIL'}")

    return removal_gate and translated_gate, removal_df, translated_df


if __name__ == "__main__":
    ok, _, _ = main()
    sys.exit(0 if ok else 1)
