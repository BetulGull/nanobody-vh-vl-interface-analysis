"""
Stage 7 — Build interaction long-table and summary tables from raw Arpeggio
JSON output + BSA results, for the full 533-VHH cohort.
"""
import sys
import json
from pathlib import Path
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from stage7_bsa import CDR3_START, CDR3_END, NONADJACENT_PACKING_TARGET, FULL_SHIELDABLE_TARGET

ROOT = Path("/Users/betulgul/Desktop/staj")
JSON_DIR = ROOT / "data/interim/stage7/cohort_json"
TABLES_DIR = ROOT / "results/tables"
TABLES_DIR.mkdir(parents=True, exist_ok=True)

# The atom-atom possible contact vocabulary (source-confirmed, stage7_arpeggio_schema.md)
ATOM_ATOM_LABELS = {'clash', 'covalent', 'vdw_clash', 'vdw', 'proximal', 'hbond', 'weak_hbond',
                     'xbond', 'ionic', 'metal_complex', 'aromatic', 'hydrophobic', 'carbonyl',
                     'polar', 'weak_polar'}


def load_instances():
    df = pd.read_parquet(ROOT / "data/processed/stage4_vhh_sequence_representatives.parquet")
    return df[['INSTANCE', 'PDB', 'Hchain']]


def classify_record(contact_labels):
    is_covalent = 'covalent' in contact_labels
    is_proximal_only = contact_labels == ['proximal']
    is_specific = any(l not in ('proximal', 'covalent') for l in contact_labels)
    return is_covalent, is_proximal_only, is_specific


def build_long_table():
    instances = load_instances()
    long_rows = []
    ring_type_seen = set()

    for _, inst_row in instances.iterrows():
        instance = inst_row['INSTANCE']
        json_path = JSON_DIR / f"{instance}.json"
        if not json_path.exists():
            continue
        with open(json_path) as f:
            data = json.load(f)

        for rec in data:
            if rec['type'] != 'atom-atom':
                ring_type_seen.add(rec['type'])
                continue  # ring-based interactions handled separately if ever observed

            bgn, end = rec['bgn'], rec['end']
            bgn_pos, end_pos = bgn['auth_seq_id'], end['auth_seq_id']

            bgn_is_cdr3 = CDR3_START <= bgn_pos <= CDR3_END
            end_is_cdr3 = CDR3_START <= end_pos <= CDR3_END
            bgn_is_full_target = bgn_pos in FULL_SHIELDABLE_TARGET
            end_is_full_target = end_pos in FULL_SHIELDABLE_TARGET

            # keep only CDR3<->target crossing records (either direction), excludes CDR3-CDR3 and target-target
            if bgn_is_cdr3 and end_is_full_target:
                cdr3_atom, target_atom = bgn, end
                cdr3_pos, target_pos = bgn_pos, end_pos
            elif end_is_cdr3 and bgn_is_full_target:
                cdr3_atom, target_atom = end, bgn
                cdr3_pos, target_pos = end_pos, bgn_pos
            else:
                continue

            contact = rec['contact']
            is_covalent, is_proximal_only, is_specific = classify_record(contact)
            is_seq_adjacent = abs(cdr3_pos - target_pos) <= 1
            in_nonadj_target = target_pos in NONADJACENT_PACKING_TARGET
            include_primary_noncovalent = is_specific and not is_covalent and in_nonadj_target

            long_rows.append({
                'INSTANCE': instance, 'PDB': inst_row['PDB'],
                'cdr3_imgt_position': cdr3_pos, 'cdr3_ins_code': cdr3_atom['pdbx_PDB_ins_code'].strip(),
                'cdr3_atom_id': cdr3_atom['auth_atom_id'], 'cdr3_resname': cdr3_atom['label_comp_id'],
                'target_imgt_position': target_pos, 'target_ins_code': target_atom['pdbx_PDB_ins_code'].strip(),
                'target_atom_id': target_atom['auth_atom_id'], 'target_resname': target_atom['label_comp_id'],
                'distance': rec['distance'],
                'contact_labels': ';'.join(contact),
                'is_covalent': is_covalent,
                'is_proximal_only': is_proximal_only,
                'is_specific_interaction': is_specific,
                'is_sequential_adjacent': is_seq_adjacent,
                'in_nonadjacent_target_set': in_nonadj_target,
                'include_primary_noncovalent': include_primary_noncovalent,
            })

    long_df = pd.DataFrame(long_rows)
    long_df.to_csv(TABLES_DIR / "stage7_interaction_records_long.csv", index=False)
    print(f"Long table: {len(long_df)} records across {long_df['INSTANCE'].nunique()} instances")
    print(f"Non-atom-atom record types observed (should be empty per smoke test): {ring_type_seen}")
    return long_df


if __name__ == '__main__':
    build_long_table()
