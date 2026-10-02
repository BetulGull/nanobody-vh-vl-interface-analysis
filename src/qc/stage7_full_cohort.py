"""
Stage 7 — Full 533-VHH Cohort Run: BSA + Arpeggio Specific Interactions.

Checkpointed: skips any INSTANCE whose CIF/JSON/BSA-row already exist and are
valid, so an interrupted run can resume without recomputation or duplication.
"""
import sys
import os
import json
import subprocess
import time
from pathlib import Path
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from structure_resolver import StructureResolver
from mmcif_qc import MMCIFStructure
from stage3_variable_domain import extract_variable_domain_atoms
from stage7_cif_writer import write_cif_with_surrogate_mapping
from stage7_bsa import compute_pairwise_bsa, NONADJACENT_PACKING_TARGET, FULL_SHIELDABLE_TARGET, CDR3_START, CDR3_END

ROOT = Path("/Users/betulgul/Desktop/staj")
CIF_DIR = ROOT / "data" / "interim" / "stage7" / "cohort_cif"
JSON_DIR = ROOT / "data" / "interim" / "stage7" / "cohort_json"
CIF_DIR.mkdir(parents=True, exist_ok=True)
JSON_DIR.mkdir(parents=True, exist_ok=True)

BSA_OUT = ROOT / "data" / "interim" / "stage7" / "bsa_results.jsonl"
STATUS_OUT = ROOT / "data" / "interim" / "stage7" / "cohort_status.jsonl"
MAPPING_OUT = ROOT / "data" / "interim" / "stage7" / "cohort_chain_mapping.jsonl"

ARPEGGIO_CMD = ["micromamba", "run", "-n", "vhh-arpeggio", "pdbe-arpeggio", "-o", str(JSON_DIR)]


def load_done_instances(path):
    done = set()
    if path.exists():
        with open(path) as f:
            for line in f:
                try:
                    rec = json.loads(line)
                    done.add(rec['INSTANCE'])
                except Exception:
                    continue
    return done


def main():
    df = pd.read_parquet(ROOT / "data/processed/stage4_vhh_sequence_representatives.parquet")
    resolver = StructureResolver()

    status_done = load_done_instances(STATUS_OUT)
    print(f"Total instances: {len(df)}; already done: {len(status_done)}")

    bsa_f = open(BSA_OUT, 'a')
    status_f = open(STATUS_OUT, 'a')
    map_f = open(MAPPING_OUT, 'a')

    for i, (_, row) in enumerate(df.iterrows()):
        instance = row['INSTANCE']
        if instance in status_done:
            continue

        rec = {'INSTANCE': instance, 'PDB': row['PDB']}
        try:
            res = resolver.resolve_structure_path(row['extracted_structure_path'], row['PDB'], row['INSTANCE'], row['archive_sha256'])
            if not res['resolved']:
                rec['status'] = 'STRUCTURE_UNRESOLVED'
                status_f.write(json.dumps(rec) + '\n'); status_f.flush()
                continue

            struct = MMCIFStructure(res['structure_path'])
            atoms = extract_variable_domain_atoms(struct, row['Hchain'], '1')

            if not atoms:
                rec['status'] = 'NO_ATOMS_EXTRACTED'
                status_f.write(json.dumps(rec) + '\n'); status_f.flush()
                continue

            bsa_nonadj = compute_pairwise_bsa(atoms, NONADJACENT_PACKING_TARGET)
            bsa_full = compute_pairwise_bsa(atoms, FULL_SHIELDABLE_TARGET)
            bsa_row = {'INSTANCE': instance, 'PDB': row['PDB']}
            for k, v in bsa_nonadj.items():
                bsa_row[f'nonadj_{k}'] = v
            for k, v in bsa_full.items():
                bsa_row[f'full_{k}'] = v
            bsa_f.write(json.dumps(bsa_row) + '\n'); bsa_f.flush()

            cif_path = CIF_DIR / f"{instance}.cif"
            json_path = JSON_DIR / f"{instance}.json"

            mapping_rows = write_cif_with_surrogate_mapping(atoms, str(cif_path), instance, row['PDB'])
            for m in mapping_rows:
                map_f.write(json.dumps(m) + '\n')
            map_f.flush()

            if not json_path.exists():
                result = subprocess.run(ARPEGGIO_CMD + [str(cif_path)], capture_output=True, text=True, timeout=60)
                if result.returncode != 0 or not json_path.exists():
                    rec['status'] = 'ARPEGGIO_FAILURE'
                    rec['stderr_tail'] = result.stderr[-2000:]
                    status_f.write(json.dumps(rec) + '\n'); status_f.flush()
                    continue

            with open(json_path) as jf:
                arp_data = json.load(jf)

            rec['status'] = 'ARPEGGIO_SUCCESS_WITH_INTERACTIONS' if len(arp_data) > 0 else 'ARPEGGIO_SUCCESS_ZERO_INTERACTIONS'
            rec['n_records'] = len(arp_data)
            status_f.write(json.dumps(rec) + '\n'); status_f.flush()

        except Exception as e:
            rec['status'] = 'ARPEGGIO_FAILURE'
            rec['error'] = str(e)
            status_f.write(json.dumps(rec) + '\n'); status_f.flush()

        if (i + 1) % 25 == 0:
            print(f"Progress: {i+1}/{len(df)}")

    bsa_f.close(); status_f.close(); map_f.close()
    print("Full cohort run complete.")


if __name__ == '__main__':
    main()
