import sys
import json
from pathlib import Path
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from structure_resolver import StructureResolver
from mmcif_qc import MMCIFStructure
from stage3_variable_domain import extract_variable_domain_atoms
from stage9_antigen_contacts import get_antigen_atoms_as_list
from stage9_arpeggio_pipeline import run_arpeggio_antibody_antigen, extract_inter_records

ROOT = Path("/Users/betulgul/Desktop/staj")
WORK_DIR = ROOT / "data/interim/stage9/arpeggio_cohort"
OUT_STATUS = ROOT / "data/interim/stage9/arpeggio_status.jsonl"
OUT_INTER = ROOT / "data/interim/stage9/arpeggio_inter_long.jsonl"
OUT_STATUS.parent.mkdir(parents=True, exist_ok=True)


def load_done(path):
    done = set()
    if path.exists():
        with open(path) as f:
            for line in f:
                try:
                    done.add(json.loads(line)['INSTANCE'])
                except Exception:
                    continue
    return done


def normalize_lchain(raw):
    if pd.isna(raw) or str(raw).strip() == 'NA':
        return None
    return raw


def process(row, group, resolver):
    instance = row['INSTANCE']
    rec = {'INSTANCE': instance, 'PDB': row['PDB'], 'group': group}
    res = resolver.resolve_structure_path(row['extracted_structure_path'], row['PDB'], instance, row['archive_sha256'])
    if not res['resolved']:
        rec['status'] = 'STRUCTURE_NOT_EVALUABLE'
        return rec, []
    try:
        struct = MMCIFStructure(res['structure_path'])
        ab_atoms = extract_variable_domain_atoms(struct, row['Hchain'], '1')
        if group == 'VH':
            lchain = normalize_lchain(row['Lchain'])
            if lchain:
                ab_atoms += extract_variable_domain_atoms(struct, lchain, '1')
        antigen_chains = [c for c in str(row['antigen_chains']).split(';') if c]
        if not antigen_chains:
            rec['status'] = 'NO_VALID_PROTEIN_ANTIGEN'
            return rec, []
        ag_atoms = get_antigen_atoms_as_list(struct, antigen_chains, '1')
        if not ag_atoms:
            rec['status'] = 'ANTIGEN_MAPPING_FAILURE'
            return rec, []

        result = run_arpeggio_antibody_antigen(ab_atoms, ag_atoms, str(WORK_DIR / instance), instance)
        if result['status'] != 'ARPEGGIO_SUCCESS':
            rec['status'] = 'ARPEGGIO_FAILURE'
            rec['error'] = result.get('stderr', '')[-500:]
            return rec, []

        inter = extract_inter_records(result['data'], result['mapping'], result['ab_surrogate'])
        for r in inter:
            r['INSTANCE'] = instance
            r['PDB'] = row['PDB']
            r['group'] = group
        rec['status'] = 'ARPEGGIO_SUCCESS'
        rec['n_records'] = result['n_records']
        rec['n_inter_records'] = len(inter)
        return rec, inter
    except Exception as e:
        rec['status'] = 'ARPEGGIO_FAILURE'
        rec['error'] = str(e)
        return rec, []


def main():
    vhh_df = pd.read_parquet(ROOT / "data/processed/stage4_vhh_sequence_representatives.parquet")
    vh_df = pd.read_parquet(ROOT / "data/processed/stage3_paired_representatives.parquet")
    vh_df = vh_df[vh_df['is_primary_representative'] == True]

    resolver = StructureResolver()
    done = load_done(OUT_STATUS)
    print(f"VHH N={len(vhh_df)}, VH N={len(vh_df)}, already done={len(done)}")

    status_f = open(OUT_STATUS, 'a')
    inter_f = open(OUT_INTER, 'a')

    total = len(vhh_df) + len(vh_df)
    i = 0
    for group, df in [('VHH', vhh_df), ('VH', vh_df)]:
        for _, row in df.iterrows():
            i += 1
            if row['INSTANCE'] in done:
                continue
            rec, inter = process(row, group, resolver)
            status_f.write(json.dumps(rec) + '\n')
            status_f.flush()
            for r in inter:
                inter_f.write(json.dumps(r) + '\n')
            inter_f.flush()
            if i % 100 == 0:
                print(f"Progress: {i}/{total}")

    status_f.close()
    inter_f.close()
    print("Done.")


if __name__ == '__main__':
    main()
