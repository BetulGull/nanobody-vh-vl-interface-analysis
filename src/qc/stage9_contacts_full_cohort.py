"""
Stage 9 — Full-cohort antibody-antigen contact computation.
Reuses the frozen Stage-2 Gordon-verified antigen-engaged cohorts (533 VHH,
1429 paired VH+VL) directly -- no further cohort reduction, per the
Stage-2 reuse audit (stage9_antigen_engagement_method_definition.md).
"""
import sys
import json
from pathlib import Path
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from structure_resolver import StructureResolver
from mmcif_qc import MMCIFStructure
from stage3_variable_domain import extract_variable_domain_atoms
from stage9_antigen_contacts import get_antigen_atoms_as_list, compute_residue_contacts

ROOT = Path("/Users/betulgul/Desktop/staj")
OUT_LONG = ROOT / "data/interim/stage9/contacts_long.jsonl"
OUT_STATUS = ROOT / "data/interim/stage9/cohort_status.jsonl"
OUT_LONG.parent.mkdir(parents=True, exist_ok=True)


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


def process_vhh(row, resolver):
    instance = row['INSTANCE']
    rec = {'INSTANCE': instance, 'PDB': row['PDB'], 'group': 'VHH'}
    res = resolver.resolve_structure_path(row['extracted_structure_path'], row['PDB'], instance, row['archive_sha256'])
    if not res['resolved']:
        rec['status'] = 'STRUCTURE_NOT_EVALUABLE'
        return rec, []
    try:
        struct = MMCIFStructure(res['structure_path'])
        ab_atoms = extract_variable_domain_atoms(struct, row['Hchain'], '1')
        antigen_chains = [c for c in str(row['antigen_chains']).split(';') if c]
        if not antigen_chains:
            rec['status'] = 'NO_VALID_PROTEIN_ANTIGEN'
            return rec, []
        ag_atoms = get_antigen_atoms_as_list(struct, antigen_chains, '1')
        if not ag_atoms:
            rec['status'] = 'ANTIGEN_MAPPING_FAILURE'
            return rec, []
        contact_rows = compute_residue_contacts(ab_atoms, ag_atoms, 'VHH')
        for c in contact_rows:
            c['INSTANCE'] = instance
            c['PDB'] = row['PDB']
            c['group'] = 'VHH'
        rec['status'] = 'ANTIGEN_ENGAGEMENT_SUCCESS' if any(c['contact_at_4p5A'] for c in contact_rows) else 'ANTIGEN_NOT_CONTACTING'
        rec['n_antigen_chains'] = len(antigen_chains)
        return rec, contact_rows
    except Exception as e:
        rec['status'] = 'ANTIBODY_MAPPING_FAILURE'
        rec['error'] = str(e)
        return rec, []


def process_vh(row, resolver):
    instance = row['INSTANCE']
    rec = {'INSTANCE': instance, 'PDB': row['PDB'], 'group': 'VH'}
    res = resolver.resolve_structure_path(row['extracted_structure_path'], row['PDB'], instance, row['archive_sha256'])
    if not res['resolved']:
        rec['status'] = 'STRUCTURE_NOT_EVALUABLE'
        return rec, []
    try:
        struct = MMCIFStructure(res['structure_path'])
        vh_atoms = extract_variable_domain_atoms(struct, row['Hchain'], '1')
        lchain = normalize_lchain(row['Lchain'])
        vl_atoms = extract_variable_domain_atoms(struct, lchain, '1') if lchain else []
        antigen_chains = [c for c in str(row['antigen_chains']).split(';') if c]
        if not antigen_chains:
            rec['status'] = 'NO_VALID_PROTEIN_ANTIGEN'
            return rec, []
        ag_atoms = get_antigen_atoms_as_list(struct, antigen_chains, '1')
        if not ag_atoms:
            rec['status'] = 'ANTIGEN_MAPPING_FAILURE'
            return rec, []
        contact_rows = compute_residue_contacts(vh_atoms, ag_atoms, 'VH')
        if vl_atoms:
            contact_rows += compute_residue_contacts(vl_atoms, ag_atoms, 'VL')
        for c in contact_rows:
            c['INSTANCE'] = instance
            c['PDB'] = row['PDB']
            c['group'] = 'VH'
        rec['status'] = 'ANTIGEN_ENGAGEMENT_SUCCESS' if any(c['contact_at_4p5A'] for c in contact_rows) else 'ANTIGEN_NOT_CONTACTING'
        rec['n_antigen_chains'] = len(antigen_chains)
        return rec, contact_rows
    except Exception as e:
        rec['status'] = 'ANTIBODY_MAPPING_FAILURE'
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
    long_f = open(OUT_LONG, 'a')

    total = len(vhh_df) + len(vh_df)
    i = 0
    for group, df, fn in [('VHH', vhh_df, process_vhh), ('VH', vh_df, process_vh)]:
        for _, row in df.iterrows():
            i += 1
            if row['INSTANCE'] in done:
                continue
            rec, contact_rows = fn(row, resolver)
            status_f.write(json.dumps(rec) + '\n')
            status_f.flush()
            for c in contact_rows:
                long_f.write(json.dumps(c) + '\n')
            long_f.flush()
            if i % 100 == 0:
                print(f"Progress: {i}/{total}")

    status_f.close()
    long_f.close()
    print("Done.")


if __name__ == '__main__':
    main()
