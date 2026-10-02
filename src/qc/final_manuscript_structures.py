"""
Final manuscript figures -- shared structure-extraction/writing helpers.
Reuses the exact validated Stage-3/12 machinery (IMGT-numbered variable-domain
atom extraction, common-framework Kabsch alignment) -- no new numbering logic.
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import gemmi

sys.path.insert(0, str(Path(__file__).parent))
from structure_resolver import StructureResolver
from mmcif_qc import MMCIFStructure
from stage3_variable_domain import extract_variable_domain_atoms
from stage12_raybould_alignment import get_ca_by_imgt, kabsch, FR_POSITIONS

ROOT = Path("/Users/betulgul/Desktop/staj")
SENTINEL = -1.0


def _row_lookup():
    vhh_df = pd.read_parquet(ROOT / "data/processed/stage4_vhh_sequence_representatives.parquet")
    vh_df = pd.read_parquet(ROOT / "data/processed/stage3_paired_representatives.parquet")
    vh_df = vh_df[vh_df['is_primary_representative'] == True]
    rows = {}
    for _, r in vhh_df.iterrows():
        rows[r['INSTANCE']] = r
    for _, r in vh_df.iterrows():
        rows[r['INSTANCE']] = r
    return rows


ROWS = None


def get_row(instance):
    global ROWS
    if ROWS is None:
        ROWS = _row_lookup()
    return ROWS[instance]


def resolve(instance):
    row = get_row(instance)
    resolver = StructureResolver()
    res = resolver.resolve_structure_path(row['extracted_structure_path'], row['PDB'], row['INSTANCE'], row['archive_sha256'])
    return MMCIFStructure(res['structure_path']), row


def get_variable_domain_atoms(instance, chain_which='H'):
    struct, row = resolve(instance)
    chain_id = row['Hchain'] if chain_which == 'H' else row['Lchain']
    return extract_variable_domain_atoms(struct, chain_id, '1')


def get_full_chain_atoms(instance, chain_id):
    """Full polymer atoms (e.g. antigen), no IMGT restriction, no renumbering."""
    row = get_row(instance)
    resolver = StructureResolver()
    res_info = resolver.resolve_structure_path(row['extracted_structure_path'], row['PDB'], row['INSTANCE'], row['archive_sha256'])
    gstruct = gemmi.read_structure(res_info['structure_path'])
    gstruct.setup_entities()
    model = gstruct[0]
    scheme_atoms = []
    chain = None
    for c in model:
        if c.name == chain_id:
            chain = c
            break
    if chain is None:
        raise KeyError(f"chain {chain_id} not found in {instance}")
    for res in chain:
        info = gemmi.find_tabulated_residue(res.name)
        if info is None or not info.is_amino_acid():
            continue
        for atom in res:
            if atom.element.is_hydrogen:
                continue
            scheme_atoms.append({
                'auth_seq_id': res.seqid.num, 'ins_code': res.seqid.icode.strip(),
                'residue_name': res.name, 'atom_name': atom.name, 'element': atom.element.name,
                'x': atom.pos.x, 'y': atom.pos.y, 'z': atom.pos.z,
                'occupancy': atom.occ if atom.occ else 1.0,
            })
    return scheme_atoms


def build_bfactor_lookup(instance, comp_df, source_col='kyte_doolittle_norm', filter_interface_exposed=True):
    sub = comp_df[comp_df['INSTANCE'] == instance]
    lookup = {}
    for _, r in sub.iterrows():
        if filter_interface_exposed and not (bool(r['primary_empirical_interface']) and bool(r['is_surface_exposed'])):
            continue
        icode = r['ins_code'] if pd.notna(r['ins_code']) and r['ins_code'] else ''
        val = r[source_col]
        if pd.notna(val):
            lookup[(int(r['imgt_position']), icode)] = float(val)
    return lookup


def write_pdb(atoms_by_chain, out_path):
    """atoms_by_chain: list of (chain_letter, atoms_list, bfactor_lookup_or_None)."""
    lines = []
    serial = 1
    for chain_letter, atoms, bfac_lookup in atoms_by_chain:
        for a in atoms:
            resname = a['residue_name']
            resseq = a['auth_seq_id']
            icode = a['ins_code'] if a['ins_code'] else ''
            key = (resseq, icode)
            bfactor = SENTINEL
            if bfac_lookup is not None and key in bfac_lookup:
                bfactor = bfac_lookup[key]
            elif bfac_lookup is None:
                bfactor = 0.0
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


def align_to_reference(mobile_atoms, ref_common_ca):
    """Kabsch-align mobile_atoms (variable-domain, IMGT-numbered) onto a reference
    set of IMGT-framework CA coordinates. Returns (R, t, rmsd, n_common, transform_fn)."""
    mobile_ca = get_ca_by_imgt(mobile_atoms)
    common = sorted(set(mobile_ca.keys()) & set(ref_common_ca.keys()))
    P = np.array([mobile_ca[p] for p in common])
    Q = np.array([ref_common_ca[p] for p in common])
    R, t, rmsd = kabsch(P, Q)

    def transform(atoms):
        out = []
        for a in atoms:
            xyz = np.array([a['x'], a['y'], a['z']])
            new_xyz = R @ xyz + t
            a2 = dict(a)
            a2['x'], a2['y'], a2['z'] = new_xyz
            out.append(a2)
        return out
    return R, t, rmsd, len(common), transform
