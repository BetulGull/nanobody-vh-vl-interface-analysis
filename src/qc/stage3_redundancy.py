"""
Stage 3 — Redundancy Verification and Deterministic Representative Selection

Verifies SABDAB_ID sequence-consistency for the paired cohort using actual
IMGT-numbered VH/VL sequences (not assumed), then selects one representative
structure per validated redundancy group using a fixed, non-ΔSASA-based
hierarchy.
"""

import sys
from pathlib import Path
from typing import Dict, Optional
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


def normalize_lchain(raw):
    if pd.isna(raw) or str(raw).strip() == 'NA':
        return None
    return raw


def build_variable_domain_sequence(structure: MMCIFStructure, chain_id: str, model_id: str) -> Dict[str, str]:
    """Returns dict {imgt_position_string ('105', '112A'): one_letter_residue}."""
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


def sequences_consistent(seq_a: Dict[str, str], seq_b: Dict[str, str]) -> bool:
    """Consistent iff every SHARED IMGT position has the same residue identity."""
    shared = set(seq_a) & set(seq_b)
    if not shared:
        return True  # nothing to contradict; conservative
    return all(seq_a[pos] == seq_b[pos] for pos in shared)


def coverage_score(seq: Dict[str, str]) -> int:
    return len(seq)


def main():
    master = pd.read_parquet(PROCESSED_DIR / "stage2_paired_primary.parquet")
    logger.info(f"Building VH/VL variable-domain sequences for {len(master)} paired structures...")

    seq_cache = {}
    coverage_cache = {}
    errors = []

    for idx, row in master.iterrows():
        instance = row['INSTANCE']
        hchain = row['Hchain']
        lchain = normalize_lchain(row['Lchain'])
        path = PROJECT_ROOT / row['extracted_structure_path']

        try:
            s = MMCIFStructure(str(path))
            model_id = s.get_model()
            vh_seq = build_variable_domain_sequence(s, hchain, model_id)
            vl_seq = build_variable_domain_sequence(s, lchain, model_id) if lchain else {}
            seq_cache[instance] = (vh_seq, vl_seq)
            coverage_cache[instance] = (coverage_score(vh_seq), coverage_score(vl_seq))
        except Exception as e:
            errors.append({'INSTANCE': instance, 'error': str(e)})
            seq_cache[instance] = ({}, {})
            coverage_cache[instance] = (0, 0)

    logger.info(f"Sequence extraction complete. Errors: {len(errors)}")

    # --- Verify SABDAB_ID sequence consistency ---
    consistency_rows = []
    redundancy_group_override = {}  # INSTANCE -> override group key, only set when inconsistent

    for sabdab_id, group in master.groupby('SABDAB_ID'):
        instances = group['INSTANCE'].tolist()
        if len(instances) < 2:
            continue

        vh_seqs = [seq_cache[i][0] for i in instances]
        vl_seqs = [seq_cache[i][1] for i in instances]

        vh_consistent = all(sequences_consistent(vh_seqs[0], vh_seqs[j]) for j in range(1, len(vh_seqs)))
        vl_consistent = all(sequences_consistent(vl_seqs[0], vl_seqs[j]) for j in range(1, len(vl_seqs)))

        consistency_rows.append({
            'SABDAB_ID': sabdab_id,
            'n_instances': len(instances),
            'vh_sequence_consistent': vh_consistent,
            'vl_sequence_consistent': vl_consistent,
            'group_valid_for_sabdab_id_key': vh_consistent and vl_consistent,
        })

        if not (vh_consistent and vl_consistent):
            # Fall back to exact-sequence-pair grouping for this SABDAB_ID's members only
            for i in instances:
                vh_seq, vl_seq = seq_cache[i]
                vh_str = ';'.join(f"{k}:{v}" for k, v in sorted(vh_seq.items()))
                vl_str = ';'.join(f"{k}:{v}" for k, v in sorted(vl_seq.items()))
                redundancy_group_override[i] = f"SEQPAIR::{hash((vh_str, vl_str)) & 0xffffffff}"

    consistency_df = pd.DataFrame(consistency_rows)
    consistency_out = AUDITS_DIR / "stage3_sabdab_id_sequence_consistency.csv"
    consistency_df.to_csv(consistency_out, index=False)
    logger.info(f"Wrote {consistency_out}")
    n_inconsistent = (~consistency_df['group_valid_for_sabdab_id_key']).sum() if len(consistency_df) else 0
    logger.info(f"SABDAB_ID groups checked: {len(consistency_df)}; inconsistent: {n_inconsistent}")

    # --- Assign final redundancy_group key ---
    master = master.copy()
    master['redundancy_group'] = master.apply(
        lambda r: redundancy_group_override.get(r['INSTANCE'], r['SABDAB_ID']), axis=1
    )
    master['vh_coverage'] = master['INSTANCE'].map(lambda i: coverage_cache[i][0])
    master['vl_coverage'] = master['INSTANCE'].map(lambda i: coverage_cache[i][1])

    # --- Deterministic representative selection ---
    def select_representative(group: pd.DataFrame) -> str:
        g = group.copy()
        g['resolution_sort'] = pd.to_numeric(g.get('resolution', pd.Series([None]*len(g))), errors='coerce')
        g['r_free_sort'] = pd.to_numeric(g.get('r_free', pd.Series([None]*len(g))), errors='coerce')
        g['r_factor_sort'] = pd.to_numeric(g.get('r_factor', pd.Series([None]*len(g))), errors='coerce')

        g_sorted = g.sort_values(
            by=['vh_coverage', 'vl_coverage', 'resolution_sort', 'r_free_sort', 'r_factor_sort', 'INSTANCE'],
            ascending=[False, False, True, True, True, True],
            na_position='last',
        )
        return g_sorted.iloc[0]['INSTANCE']

    representative_map = {}
    for group_key, group in master.groupby('redundancy_group'):
        rep = select_representative(group)
        representative_map[group_key] = rep

    master['is_primary_representative'] = master.apply(
        lambda r: representative_map[r['redundancy_group']] == r['INSTANCE'], axis=1
    )

    n_groups = master['redundancy_group'].nunique()
    n_representatives = master['is_primary_representative'].sum()
    logger.info(f"Redundancy groups: {n_groups}; representatives selected: {n_representatives}")
    assert n_groups == n_representatives, "One representative per group required"

    # --- Outputs ---
    rep_df = master[master['is_primary_representative']].copy()
    rep_df.to_parquet(PROCESSED_DIR / "stage3_paired_representatives.parquet", index=False)
    rep_df.to_csv(PROCESSED_DIR / "stage3_paired_representatives.csv", index=False)
    logger.info(f"Wrote stage3_paired_representatives (.parquet/.csv): {len(rep_df)} rows")

    # Full annotated table (all instances, redundancy_group + is_primary_representative
    # columns attached) -- input for the full ΔSASA orchestrator run.
    master.to_parquet(PROCESSED_DIR / "stage3_all_instances_with_redundancy.parquet", index=False)
    logger.info(f"Wrote stage3_all_instances_with_redundancy.parquet: {len(master)} rows")

    selection_audit = master[[
        'INSTANCE', 'PDB', 'SABDAB_ID', 'redundancy_group', 'Hchain', 'Lchain',
        'vh_coverage', 'vl_coverage', 'is_primary_representative',
    ]].copy()
    selection_audit_out = AUDITS_DIR / "stage3_representative_selection.csv"
    selection_audit.to_csv(selection_audit_out, index=False)
    logger.info(f"Wrote {selection_audit_out}")

    return master, rep_df, errors


if __name__ == "__main__":
    main()
