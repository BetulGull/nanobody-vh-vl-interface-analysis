"""
Stage 12 — Raybould-style figure: common-framework CA alignment between the
selected VHH and VH exemplars. Uses IMGT framework positions only (FR1-FR4),
never CDR3/antigen/VL/the empirical interface itself as the alignment target.
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from structure_resolver import StructureResolver
from mmcif_qc import MMCIFStructure
from stage3_variable_domain import extract_variable_domain_atoms

FR_POSITIONS = set(range(1, 27)) | set(range(39, 56)) | set(range(66, 105)) | set(range(118, 129))


def get_ca_by_imgt(atoms):
    d = {}
    for a in atoms:
        if a['atom_name'] == 'CA' and a['auth_seq_id'] in FR_POSITIONS and not a['ins_code']:
            d[a['auth_seq_id']] = np.array([a['x'], a['y'], a['z']])
    return d


def kabsch(P, Q):
    """Superpose P onto Q. Returns R, t, rmsd such that Q ~= P @ R.T + t"""
    Pc = P - P.mean(axis=0)
    Qc = Q - Q.mean(axis=0)
    H = Pc.T @ Qc
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1, 1, d])
    R = Vt.T @ D @ U.T
    P_rot = (R @ Pc.T).T
    rmsd = np.sqrt(np.mean(np.sum((P_rot - Qc) ** 2, axis=1)))
    t = Q.mean(axis=0) - R @ P.mean(axis=0)
    return R, t, rmsd


def main():
    ROOT = Path("/Users/betulgul/Desktop/staj")
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

    atoms1 = extract_variable_domain_atoms(struct1, vhh_row['Hchain'], '1')
    atoms2 = extract_variable_domain_atoms(struct2, vh_row['Hchain'], '1')

    ca1 = get_ca_by_imgt(atoms1)
    ca2 = get_ca_by_imgt(atoms2)
    common = sorted(set(ca1.keys()) & set(ca2.keys()))

    P = np.array([ca1[p] for p in common])
    Q = np.array([ca2[p] for p in common])
    R, t, rmsd = kabsch(P, Q)

    audit = pd.DataFrame([{
        'VHH_PDB': vhh_row['PDB'], 'VHH_INSTANCE': vhh_row['INSTANCE'], 'VHH_chain': vhh_row['Hchain'],
        'VH_PDB': vh_row['PDB'], 'VH_INSTANCE': vh_row['INSTANCE'], 'VH_chain': vh_row['Hchain'], 'VL_chain': vh_row['Lchain'],
        'imgt_positions_used': ';'.join(str(p) for p in common), 'n_aligned_positions': len(common),
        'rmsd_angstrom': rmsd, 'alignment_target': 'common IMGT framework (FR1-FR4) C-alpha atoms',
    }])
    audit.to_csv(ROOT / "results/structures/stage12_raybould_style_alignment_audit.csv", index=False)
    print(f"N aligned FR positions: {len(common)}, RMSD: {rmsd:.3f} A")
    print(common)

    np.save(ROOT / "data/interim/stage12_alignment_R.npy", R)
    np.save(ROOT / "data/interim/stage12_alignment_t.npy", t)
    return R, t, rmsd, common


if __name__ == '__main__':
    main()
