"""
Stage 12 — write aligned PDB files for the two Raybould-style exemplars,
with per-atom B-factor encoding the frozen Stage-8 normalized hydropathy
value ONLY for residues that are (a) within the frozen 28-position
FULL_EMPIRICAL_INTERFACE and (b) surface-exposed per the frozen Stage-8
classification. All other atoms get a sentinel B-factor (-1) so PyMOL can
render them white. No new hydropathy/exposure values are computed here --
all values are read directly from the already-frozen Stage-8 per-residue
components table (see stage12_raybould_style_surface_residue_mapping.csv).
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from structure_resolver import StructureResolver
from mmcif_qc import MMCIFStructure
from stage3_variable_domain import extract_variable_domain_atoms
from stage12_raybould_alignment import get_ca_by_imgt, kabsch

ROOT = Path("/Users/betulgul/Desktop/staj")
OUT_DIR = ROOT / "data/interim/stage12_raybould_structures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

FULL_EMPIRICAL_INTERFACE = {40, 42, 44, 49, 50, 52, 55, 66, 67, 68, 69, 103,
                             107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 118, 119, 120}
SENTINEL = -1.0


def write_pdb(atoms_by_chain, mapping_csv_path, out_path, chain_label_map):
    """atoms_by_chain: list of (chain_letter_for_pdb, atoms_list, hydropathy_lookup or None)."""
    lines = []
    serial = 1
    for chain_letter, atoms, hydro_lookup in atoms_by_chain:
        for a in atoms:
            resname = a['residue_name']
            resseq = a['auth_seq_id']
            icode = a['ins_code'] if a['ins_code'] else ''
            key = (resseq, icode)
            bfactor = SENTINEL
            if hydro_lookup is not None and key in hydro_lookup:
                bfactor = hydro_lookup[key]
            line = (
                f"ATOM  {serial:5d} {a['atom_name']:<4s} {resname:<3s} {chain_letter:1s}"
                f"{resseq:4d}{icode:1s}   "
                f"{a['x']:8.3f}{a['y']:8.3f}{a['z']:8.3f}{a['occupancy']:6.2f}{bfactor:6.2f}          "
                f"{a['element']:>2s}"
            )
            lines.append(line)
            serial += 1
    lines.append("END")
    with open(out_path, 'w') as f:
        f.write('\n'.join(lines) + '\n')


def build_hydro_lookup(residue_df, instance):
    """From the frozen residue mapping table: only interface-exposed residues get a real value."""
    sub = residue_df[residue_df['INSTANCE'] == instance]
    lookup = {}
    for _, r in sub.iterrows():
        if r['is_primary_empirical_interface'] and r['surface_exposed']:
            icode = '' if pd.isna(r.get('ins_code')) else r.get('ins_code', '')
            # residue mapping table doesn't carry ins_code column name exactly; reconstruct from IMGT_position + look at raw comp table
    return lookup


def main():
    residue_map = pd.read_csv(ROOT / "results/structures/stage12_raybould_style_surface_residue_mapping.csv")
    # need ins_code too -- reload from the frozen patch components table directly (already has ins_code)
    comp = pd.read_parquet(ROOT / "data/processed/stage8_patch_components.parquet")
    comp = comp[comp['INSTANCE'].isin(['pdb_00008f6v-A', 'pdb_00006wo3-H-L'])]

    def hydro_lookup_for(instance):
        sub = comp[comp['INSTANCE'] == instance]
        lookup = {}
        for _, r in sub.iterrows():
            if bool(r['primary_empirical_interface']) and bool(r['is_surface_exposed']):
                icode = r['ins_code'] if pd.notna(r['ins_code']) and r['ins_code'] else ''
                lookup[(int(r['imgt_position']), icode)] = float(r['kyte_doolittle_norm'])
        return lookup

    vhh_df = pd.read_parquet(ROOT / "data/processed/stage4_vhh_sequence_representatives.parquet")
    vh_df = pd.read_parquet(ROOT / "data/processed/stage3_paired_representatives.parquet")
    vh_df = vh_df[vh_df['is_primary_representative'] == True]

    vhh_row = vhh_df[vhh_df.INSTANCE == 'pdb_00008f6v-A'].iloc[0]
    vh_row = vh_df[vh_df.INSTANCE == 'pdb_00006wo3-H-L'].iloc[0]

    resolver = StructureResolver()
    res1 = resolver.resolve_structure_path(vhh_row['extracted_structure_path'], vhh_row['PDB'], vhh_row['INSTANCE'], vhh_row['archive_sha256'])
    res2 = resolver.resolve_structure_path(vh_row['extracted_structure_path'], vh_row['PDB'], vh_row['INSTANCE'], vh_row['archive_sha256'])
    struct1 = MMCIFStructure(res1['structure_path'])
    struct2 = MMCIFStructure(res2['structure_path'])

    vhh_atoms = extract_variable_domain_atoms(struct1, vhh_row['Hchain'], '1')
    vh_atoms = extract_variable_domain_atoms(struct2, vh_row['Hchain'], '1')
    vl_atoms = extract_variable_domain_atoms(struct2, vh_row['Lchain'], '1')

    # alignment: recompute R,t exactly as stage12_raybould_alignment.py (framework CA only)
    ca1 = get_ca_by_imgt(vhh_atoms)
    ca2 = get_ca_by_imgt(vh_atoms)
    common = sorted(set(ca1.keys()) & set(ca2.keys()))
    P = np.array([ca2[p] for p in common])  # VH (to be moved)
    Q = np.array([ca1[p] for p in common])  # VHH (reference)
    R, t, rmsd = kabsch(P, Q)
    print(f"Alignment RMSD: {rmsd:.3f} A over {len(common)} positions")

    def transform_atoms(atoms):
        out = []
        for a in atoms:
            xyz = np.array([a['x'], a['y'], a['z']])
            new_xyz = R @ xyz + t
            a2 = dict(a)
            a2['x'], a2['y'], a2['z'] = new_xyz
            out.append(a2)
        return out

    vh_atoms_aligned = transform_atoms(vh_atoms)
    vl_atoms_aligned = transform_atoms(vl_atoms)

    vhh_hydro = hydro_lookup_for('pdb_00008f6v-A')
    vh_hydro = hydro_lookup_for('pdb_00006wo3-H-L')

    write_pdb([('A', vhh_atoms, vhh_hydro)], None, OUT_DIR / "vhh_exemplar.pdb", {})
    write_pdb([('H', vh_atoms_aligned, vh_hydro), ('L', vl_atoms_aligned, None)], None, OUT_DIR / "vh_exemplar_with_vl_aligned.pdb", {})
    write_pdb([('H', vh_atoms_aligned, vh_hydro)], None, OUT_DIR / "vh_exemplar_only_aligned.pdb", {})

    print("Wrote PDB files to", OUT_DIR)
    print("VHH interface-exposed hydro residues:", len(vhh_hydro))
    print("VH interface-exposed hydro residues:", len(vh_hydro))


if __name__ == '__main__':
    main()
