"""
Stage 8 — Full-cohort PSH/PPC/PNC patch scoring: VHH (533, VHH_NATIVE) and
conventional VH (1429, VH_WITHOUT_VL_GEOMETRIC -- Hchain-only extraction,
same as Stage 5's counterfactual definition, VL never loaded).

Computes BOTH:
  - ADAPTED former-interface descriptors (PRIMARY): core = Stage-3 28-position
    empirical interface, vicinity = core + 4A-exposed-neighbors.
  - STANDARD TNP-context descriptors (SECONDARY): core = TNP's own internal
    CDR definitions (IMGT H1=25-40,H2=54-67,H3=103-119), for comparability
    with the published TNP context. Never substituted for the primary.
"""
import sys
import json
from pathlib import Path
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from structure_resolver import StructureResolver
from mmcif_qc import MMCIFStructure
from stage3_variable_domain import extract_variable_domain_atoms
from stage8_patch_scoring import compute_patch_for_target

ROOT = Path("/Users/betulgul/Desktop/staj")
OUT = ROOT / "data/interim/stage8/patch_results.jsonl"
OUT.parent.mkdir(parents=True, exist_ok=True)
COMPONENTS_OUT = ROOT / "data/interim/stage8/patch_components.jsonl"

FULL_EMPIRICAL_INTERFACE = {40, 42, 44, 49, 50, 52, 55, 66, 67, 68, 69, 103,
                             107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 118, 119, 120}
TNP_STANDARD_CORE = set(range(25, 41)) | set(range(54, 68)) | set(range(103, 120))
CDR_POSITIONS = set(range(27, 39)) | set(range(56, 66)) | set(range(105, 118))
VHH_DOMAIN_RANGE = set(range(1, 129))


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
        return rec, None
    try:
        struct = MMCIFStructure(res['structure_path'])
        atoms = extract_variable_domain_atoms(struct, hchain, '1')
        if not atoms:
            rec['status'] = 'NO_ATOMS_EXTRACTED'
            return rec, None

        adapted = compute_patch_for_target(atoms, FULL_EMPIRICAL_INTERFACE)
        standard = compute_patch_for_target(atoms, TNP_STANDARD_CORE)

        from stage8_patch_scoring import compute_patch_scores
        exposed_keys = {k for k, r in adapted['residues'].items() if r.get('is_surface_exposed')}
        noninterface_fw_keys = {k for k in exposed_keys if k[0] not in FULL_EMPIRICAL_INTERFACE and k[0] not in CDR_POSITIONS}
        global_keys = exposed_keys
        framework_scores = compute_patch_scores(adapted['residues'], noninterface_fw_keys)
        global_scores = compute_patch_scores(adapted['residues'], global_keys)

        n_missing = sum(1 for r in adapted['residues'].values() if r['exposure_status'] != 'OK')
        rec['status'] = 'OK'
        rec['n_residues'] = len(adapted['residues'])
        rec['n_exposed'] = sum(1 for r in adapted['residues'].values() if r.get('is_surface_exposed'))
        rec['n_psa_missing'] = n_missing
        rec['n_core'] = len(adapted['core_keys'])
        rec['n_vicinity'] = len(adapted['vicinity_keys'])
        rec['former_interface_core_PSH'] = adapted['core_scores']['PSH']
        rec['former_interface_core_PPC'] = adapted['core_scores']['PPC']
        rec['former_interface_core_PNC'] = adapted['core_scores']['PNC']
        rec['former_interface_vicinity_PSH'] = adapted['vicinity_scores']['PSH']
        rec['former_interface_vicinity_PPC'] = adapted['vicinity_scores']['PPC']
        rec['former_interface_vicinity_PNC'] = adapted['vicinity_scores']['PNC']
        rec['standard_tnp_core_PSH'] = standard['core_scores']['PSH']
        rec['standard_tnp_core_PPC'] = standard['core_scores']['PPC']
        rec['standard_tnp_core_PNC'] = standard['core_scores']['PNC']
        rec['standard_tnp_vicinity_PSH'] = standard['vicinity_scores']['PSH']
        rec['standard_tnp_vicinity_PPC'] = standard['vicinity_scores']['PPC']
        rec['standard_tnp_vicinity_PNC'] = standard['vicinity_scores']['PNC']
        rec['noninterface_framework_PSH'] = framework_scores['PSH']
        rec['noninterface_framework_PPC'] = framework_scores['PPC']
        rec['noninterface_framework_PNC'] = framework_scores['PNC']
        rec['n_noninterface_framework_exposed'] = len(noninterface_fw_keys)
        rec['global_PSH_like'] = global_scores['PSH']
        rec['global_positive_patch'] = global_scores['PPC']
        rec['global_negative_patch'] = global_scores['PNC']
        rec['n_global_exposed'] = len(global_keys)

        component_rows = []
        for key, r in adapted['residues'].items():
            component_rows.append({
                'INSTANCE': instance, 'PDB': pdb, 'group': group,
                'imgt_position': key[0], 'ins_code': key[1], 'amino_acid': r['residue_name'],
                'primary_empirical_interface': key[0] in FULL_EMPIRICAL_INTERFACE,
                'interface_vicinity_member': key in adapted['vicinity_keys'],
                'sidechain_asa_abs': r.get('sidechain_asa_abs'), 'relative_sidechain_asa': r.get('relative_sidechain_asa'),
                'is_surface_exposed': r.get('is_surface_exposed'), 'exposure_status': r.get('exposure_status'),
                'kyte_doolittle_norm': r.get('kyte_doolittle_norm'), 'charge_class': r.get('charge_class'),
            })
        return rec, component_rows
    except Exception as e:
        rec['status'] = 'PATCH_FAILURE'
        rec['error'] = str(e)
        return rec, None


def main():
    vhh_df = pd.read_parquet(ROOT / "data/processed/stage4_vhh_sequence_representatives.parquet")
    vh_df = pd.read_parquet(ROOT / "data/processed/stage3_paired_representatives.parquet")
    vh_df = vh_df[vh_df['is_primary_representative'] == True]

    print(f"VHH N={len(vhh_df)}, VH N={len(vh_df)}")
    resolver = StructureResolver()
    done = load_done(OUT)
    print(f"Already done: {len(done)}")

    out_f = open(OUT, 'a')
    comp_f = open(COMPONENTS_OUT, 'a')

    total = len(vhh_df) + len(vh_df)
    i = 0
    for group, df in [('VHH', vhh_df), ('VH', vh_df)]:
        for _, row in df.iterrows():
            i += 1
            if row['INSTANCE'] in done:
                continue
            rec, components = process_instance(row['INSTANCE'], row['PDB'], row['Hchain'],
                                                 row['extracted_structure_path'], row['archive_sha256'], group, resolver)
            out_f.write(json.dumps(rec) + '\n')
            out_f.flush()
            if components:
                for c in components:
                    comp_f.write(json.dumps(c) + '\n')
                comp_f.flush()
            if i % 100 == 0:
                print(f"Progress: {i}/{total}")

    out_f.close(); comp_f.close()
    print("Done.")


if __name__ == '__main__':
    main()
