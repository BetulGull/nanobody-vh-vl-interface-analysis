import sys
import json
from pathlib import Path
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from structure_resolver import StructureResolver
from mmcif_qc import MMCIFStructure
from stage3_variable_domain import extract_variable_domain_atoms
from stage8_patch_scoring import write_minimal_pdb, run_psa, parse_psa_output, build_residue_registry, build_vicinity
from stage8_patch_topology import compute_topology, HYDROPHOBIC_AA, THRESHOLDS, PRIMARY_THRESHOLD

ROOT = Path("/Users/betulgul/Desktop/staj")
OUT = ROOT / "data/interim/stage8/topology_results.jsonl"
OUT.parent.mkdir(parents=True, exist_ok=True)

FULL_EMPIRICAL_INTERFACE = {40, 42, 44, 49, 50, 52, 55, 66, 67, 68, 69, 103,
                             107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 118, 119, 120}


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


def process_instance(instance, pdb, hchain, extracted_structure_path, archive_sha256, group, resolver):
    rec = {'INSTANCE': instance, 'PDB': pdb, 'group': group}
    res = resolver.resolve_structure_path(extracted_structure_path, pdb, instance, archive_sha256)
    if not res['resolved']:
        rec['status'] = 'STRUCTURE_UNRESOLVED'
        return rec
    try:
        struct = MMCIFStructure(res['structure_path'])
        atoms = extract_variable_domain_atoms(struct, hchain, '1')
        if not atoms:
            rec['status'] = 'NO_ATOMS_EXTRACTED'
            return rec

        import tempfile
        with tempfile.TemporaryDirectory() as td:
            pdb_path = f"{td}/t.pdb"
            write_minimal_pdb(atoms, pdb_path)
            stdout = run_psa(pdb_path)
            psa_data = parse_psa_output(stdout)
        residues = build_residue_registry(atoms, psa_data)

        exposed_hydrophobic = {k for k, r in residues.items() if r.get('is_surface_exposed') and r['residue_name'] in HYDROPHOBIC_AA}
        core_keys = {k for k in residues if k[0] in FULL_EMPIRICAL_INTERFACE}
        vicinity_keys = build_vicinity(residues, core_keys)
        local_exposed_hydrophobic = exposed_hydrophobic & vicinity_keys

        for thr in THRESHOLDS:
            g = compute_topology(residues, exposed_hydrophobic, thr)
            for k, v in g.items():
                rec[f'global_thr{thr}_{k}'] = v
            l = compute_topology(residues, local_exposed_hydrophobic, thr)
            for k, v in l.items():
                rec[f'local_thr{thr}_{k}'] = v

        rec['status'] = 'OK'
        rec['n_exposed_hydrophobic_global'] = len(exposed_hydrophobic)
        rec['n_exposed_hydrophobic_local'] = len(local_exposed_hydrophobic)
        return rec
    except Exception as e:
        rec['status'] = 'TOPOLOGY_FAILURE'
        rec['error'] = str(e)
        return rec


def main():
    vhh_df = pd.read_parquet(ROOT / "data/processed/stage4_vhh_sequence_representatives.parquet")
    vh_df = pd.read_parquet(ROOT / "data/processed/stage3_paired_representatives.parquet")
    vh_df = vh_df[vh_df['is_primary_representative'] == True]

    resolver = StructureResolver()
    done = load_done(OUT)
    print(f"VHH N={len(vhh_df)}, VH N={len(vh_df)}, already done={len(done)}")

    out_f = open(OUT, 'a')
    total = len(vhh_df) + len(vh_df)
    i = 0
    for group, df in [('VHH', vhh_df), ('VH', vh_df)]:
        for _, row in df.iterrows():
            i += 1
            if row['INSTANCE'] in done:
                continue
            rec = process_instance(row['INSTANCE'], row['PDB'], row['Hchain'],
                                    row['extracted_structure_path'], row['archive_sha256'], group, resolver)
            out_f.write(json.dumps(rec) + '\n')
            out_f.flush()
            if i % 200 == 0:
                print(f"Progress: {i}/{total}")

    out_f.close()
    print("Done.")


if __name__ == '__main__':
    main()
