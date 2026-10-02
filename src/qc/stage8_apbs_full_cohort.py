"""
Stage 8 — Full-cohort APBS run: VHH (533) + conventional VH (1429), primary
FULL COHORT decision per stage8_apbs_full_cohort_decision.md.

Samples potential at solvent-exposed former-interface atoms (core + vicinity,
reusing the exposure classification already computed in the patch branch)
and at a non-interface framework control set, plus a global-domain exposed
sample. Failures tracked explicitly, never converted to zero potential.
"""
import sys
import json
from pathlib import Path
import pandas as pd
import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from structure_resolver import StructureResolver
from mmcif_qc import MMCIFStructure
from stage3_variable_domain import extract_variable_domain_atoms
from stage8_apbs_pipeline import run_apbs_pipeline
from stage8_apbs_sampling import sample_potentials
from stage8_patch_scoring import write_minimal_pdb, run_psa, parse_psa_output, build_residue_registry, EXPOSURE_THRESHOLD_PCT

ROOT = Path("/Users/betulgul/Desktop/staj")
WORK_DIR = ROOT / "data/interim/stage8/apbs_cohort"
OUT = ROOT / "data/interim/stage8/apbs_results.jsonl"
OUT.parent.mkdir(parents=True, exist_ok=True)

FULL_EMPIRICAL_INTERFACE = {40, 42, 44, 49, 50, 52, 55, 66, 67, 68, 69, 103,
                             107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 118, 119, 120}
CDR_POSITIONS = set(range(27, 39)) | set(range(56, 66)) | set(range(105, 118))
VHH_RANGE = set(range(1, 129))


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


def describe(vals):
    vals = np.array([v for v in vals if v is not None], dtype=float)
    n_missing = sum(1 for v in vals if v is None) if False else None
    if len(vals) == 0:
        return {'n': 0, 'mean': None, 'median': None, 'sd': None, 'q25': None, 'q75': None,
                'fraction_positive': None, 'fraction_negative': None,
                'positive_mean': None, 'negative_mean': None, 'absolute_mean': None}
    pos = vals[vals > 0]
    neg = vals[vals < 0]
    return {
        'n': len(vals), 'mean': float(np.mean(vals)), 'median': float(np.median(vals)), 'sd': float(np.std(vals)),
        'q25': float(np.percentile(vals, 25)), 'q75': float(np.percentile(vals, 75)),
        'fraction_positive': float(len(pos) / len(vals)), 'fraction_negative': float(len(neg) / len(vals)),
        'positive_mean': float(np.mean(pos)) if len(pos) else None,
        'negative_mean': float(np.mean(neg)) if len(neg) else None,
        'absolute_mean': float(np.mean(np.abs(vals))),
    }


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

        apbs_result = run_apbs_pipeline(atoms, str(WORK_DIR / instance), instance)
        rec['pdb2pqr_runtime_seconds'] = apbs_result.get('pdb2pqr_runtime_seconds')
        rec['apbs_runtime_seconds'] = apbs_result.get('apbs_runtime_seconds')
        if apbs_result['status'] != 'PDB2PQR_SUCCESS_APBS_SUCCESS':
            rec['status'] = apbs_result['status']
            rec['error'] = apbs_result.get('stderr', '')[-500:]
            return rec
        rec['molecule_contained_in_grid'] = apbs_result['molecule_contained_in_grid']
        grid = apbs_result['grid']

        import tempfile
        with tempfile.TemporaryDirectory() as td:
            pdb_path = f"{td}/exp.pdb"
            write_minimal_pdb(atoms, pdb_path)
            stdout = run_psa(pdb_path)
            psa_data = parse_psa_output(stdout)
        residues = build_residue_registry(atoms, psa_data)

        exposed_keys = {k for k, r in residues.items() if r.get('is_surface_exposed')}
        interface_keys = {k for k in exposed_keys if k[0] in FULL_EMPIRICAL_INTERFACE}
        noninterface_fw_keys = {k for k in exposed_keys if k[0] not in FULL_EMPIRICAL_INTERFACE and k[0] not in CDR_POSITIONS}
        global_keys = exposed_keys

        def coords_for(keys):
            coords = []
            for k in keys:
                for a in residues[k]['atoms']:
                    coords.append((a['x'], a['y'], a['z']))
            return coords

        if apbs_result['status'] == 'PDB2PQR_SUCCESS_APBS_SUCCESS':
            pot_interface = sample_potentials(grid, coords_for(interface_keys))
            pot_fw = sample_potentials(grid, coords_for(noninterface_fw_keys))
            pot_global = sample_potentials(grid, coords_for(global_keys))

            desc_interface = describe(pot_interface)
            desc_fw = describe(pot_fw)
            desc_global = describe(pot_global)
            for k, v in desc_interface.items():
                rec[f'interface_potential_{k}'] = v
            for k, v in desc_fw.items():
                rec[f'framework_potential_{k}'] = v
            for k, v in desc_global.items():
                rec[f'global_potential_{k}'] = v

        rec['status'] = 'PDB2PQR_SUCCESS_APBS_SUCCESS'
        return rec
    except Exception as e:
        rec['status'] = 'APBS_FAILURE'
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
            if i % 50 == 0:
                print(f"Progress: {i}/{total}")

    out_f.close()
    print("Done.")


if __name__ == '__main__':
    main()
