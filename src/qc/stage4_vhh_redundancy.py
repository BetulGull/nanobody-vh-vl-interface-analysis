"""
Stage 4 — VHH Redundancy Verification and Deterministic Representative Selection

Mirrors src/qc/stage3_redundancy.py but for the single-chain strict natural
VHH cohort (no VL/pairing dimension).
"""

import sys
from pathlib import Path
import pandas as pd
import logging

sys.path.insert(0, str(Path(__file__).parent.parent))

from qc.mmcif_qc import MMCIFStructure

PROJECT_ROOT = Path(__file__).parent.parent.parent
AUDITS_DIR = PROJECT_ROOT / "results" / "audits"
PROCESSED_DIR = PROJECT_ROOT / "data" / "processed"

logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger(__name__)

THREE_TO_ONE = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q',
    'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
    'MET': 'M', 'MSE': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T',
    'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
}


def build_variable_domain_sequence(structure: MMCIFStructure, chain_id: str, model_id: str) -> dict:
    scheme = structure.build_residue_scheme(chain_id, model_id)
    seq = {}
    for r in scheme:
        if not r['observed'] or r['auth_seq_id'] is None:
            continue
        if not (1 <= r['auth_seq_id'] <= 128):
            continue
        pos_str = f"{r['auth_seq_id']}{r['ins_code']}"
        seq[pos_str] = THREE_TO_ONE.get(r['mon_id'], 'X')
    return seq


def sequences_consistent(seq_a: dict, seq_b: dict) -> bool:
    shared = set(seq_a) & set(seq_b)
    if not shared:
        return True
    return all(seq_a[pos] == seq_b[pos] for pos in shared)


def main():
    master = pd.read_parquet(PROCESSED_DIR / "stage2_vhh_strict_primary.parquet")
    logger.info(f"Building VHH variable-domain sequences for {len(master)} structures...")

    seq_cache = {}
    coverage_cache = {}
    errors = []

    for idx, row in master.iterrows():
        instance = row['INSTANCE']
        hchain = row['Hchain']
        path = PROJECT_ROOT / row['extracted_structure_path']
        try:
            s = MMCIFStructure(str(path))
            model_id = s.get_model()
            vh_seq = build_variable_domain_sequence(s, hchain, model_id)
            seq_cache[instance] = vh_seq
            coverage_cache[instance] = len(vh_seq)
        except Exception as e:
            errors.append({'INSTANCE': instance, 'error': str(e)})
            seq_cache[instance] = {}
            coverage_cache[instance] = 0

    logger.info(f"Sequence extraction complete. Errors: {len(errors)}")

    # --- Verify SABDAB_ID sequence consistency ---
    consistency_rows = []
    redundancy_group_override = {}

    for sabdab_id, group in master.groupby('SABDAB_ID'):
        instances = group['INSTANCE'].tolist()
        if len(instances) < 2:
            continue
        seqs = [seq_cache[i] for i in instances]
        consistent = all(sequences_consistent(seqs[0], seqs[j]) for j in range(1, len(seqs)))
        consistency_rows.append({
            'SABDAB_ID': sabdab_id, 'n_instances': len(instances),
            'sequence_consistent': consistent,
        })
        if not consistent:
            for i in instances:
                seq = seq_cache[i]
                seq_str = ';'.join(f"{k}:{v}" for k, v in sorted(seq.items()))
                redundancy_group_override[i] = f"SEQ::{hash(seq_str) & 0xffffffff}"

    consistency_df = pd.DataFrame(consistency_rows)
    consistency_out = AUDITS_DIR / "stage4_vhh_redundancy_validation.csv"
    consistency_df.to_csv(consistency_out, index=False)
    n_inconsistent = (~consistency_df['sequence_consistent']).sum() if len(consistency_df) else 0
    logger.info(f"Wrote {consistency_out}: {len(consistency_df)} groups checked, {n_inconsistent} inconsistent")

    master = master.copy()
    master['redundancy_group'] = master.apply(
        lambda r: redundancy_group_override.get(r['INSTANCE'], r['SABDAB_ID']), axis=1
    )
    master['vh_coverage'] = master['INSTANCE'].map(coverage_cache)

    def select_representative(group: pd.DataFrame) -> str:
        g = group.copy()
        g['resolution_sort'] = pd.to_numeric(g.get('resolution', pd.Series([None] * len(g))), errors='coerce')
        g['r_free_sort'] = pd.to_numeric(g.get('r_free', pd.Series([None] * len(g))), errors='coerce')
        g['r_factor_sort'] = pd.to_numeric(g.get('r_factor', pd.Series([None] * len(g))), errors='coerce')
        g_sorted = g.sort_values(
            by=['vh_coverage', 'resolution_sort', 'r_free_sort', 'r_factor_sort', 'INSTANCE'],
            ascending=[False, True, True, True, True], na_position='last',
        )
        return g_sorted.iloc[0]['INSTANCE']

    representative_map = {}
    for group_key, group in master.groupby('redundancy_group'):
        representative_map[group_key] = select_representative(group)

    master['is_primary_representative'] = master.apply(
        lambda r: representative_map[r['redundancy_group']] == r['INSTANCE'], axis=1
    )

    n_groups = master['redundancy_group'].nunique()
    n_reps = master['is_primary_representative'].sum()
    logger.info(f"VHH redundancy groups: {n_groups}; representatives: {n_reps}")
    assert n_groups == n_reps

    rep_df = master[master['is_primary_representative']].copy()
    rep_df.to_parquet(PROCESSED_DIR / "stage4_vhh_sequence_representatives.parquet", index=False)
    rep_df.to_csv(PROCESSED_DIR / "stage4_vhh_sequence_representatives.csv", index=False)
    logger.info(f"Wrote stage4_vhh_sequence_representatives: {len(rep_df)} rows")

    master.to_parquet(PROCESSED_DIR / "stage4_vhh_all_instances_with_redundancy.parquet", index=False)

    selection_audit = master[[
        'INSTANCE', 'PDB', 'SABDAB_ID', 'redundancy_group', 'Hchain',
        'vh_coverage', 'is_primary_representative',
    ]].copy()
    selection_audit.to_csv(AUDITS_DIR / "stage4_vhh_representative_selection.csv", index=False)
    logger.info(f"Wrote stage4_vhh_representative_selection.csv")

    return master, rep_df


if __name__ == "__main__":
    main()
