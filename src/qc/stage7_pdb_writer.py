"""
Stage 7 — Minimal PDB Writer with Reversible Identity Mapping

Arpeggio requires PDB (or mmCIF via gemmi) input. PDB format has strict
single-character chain ID and limited residue-number field width, which is
exactly the failure class already found in Stage 3 (FreeSASA silently
truncating multi-character chain labels). This writer assigns a
deterministic, collision-safe single-character surrogate chain ID and
preserves insertion codes in the standard PDB icode column, while returning
an explicit reversible mapping table.
"""

import string
from typing import List, Dict, Tuple


def write_pdb_with_surrogate_mapping(atoms: List[Dict], out_path: str, instance: str, pdb: str) -> List[Dict]:
    """
    atoms: list of atom dicts as produced by extract_variable_domain_atoms
           (auth_asym_id, label_seq_id, auth_seq_id, ins_code, residue_name,
           atom_name, element, x, y, z, occupancy, altloc)
    Returns: mapping_rows, one per distinct (auth_asym_id) -> surrogate chain,
             plus per-residue identity preserved via auth_seq_id/ins_code
             (both fit PDB's native residue-number/icode columns unchanged,
             so only the CHAIN id needs surrogate treatment here).
    """
    surrogate_map: Dict[str, str] = {}
    available = iter(string.ascii_uppercase)
    for a in atoms:
        if a['auth_asym_id'] not in surrogate_map:
            surrogate_map[a['auth_asym_id']] = next(available)

    lines = []
    atom_serial = 1
    for a in atoms:
        surrogate_chain = surrogate_map[a['auth_asym_id']]
        resname = a['residue_name']
        resseq = a['auth_seq_id']
        icode = a['ins_code'] if a['ins_code'] else ''
        atom_name = a['atom_name']
        # PDB ATOM record column layout (1-indexed, fixed width)
        line = (
            f"ATOM  {atom_serial:5d} {atom_name:<4s} {resname:<3s} {surrogate_chain:1s}"
            f"{resseq:4d}{icode:1s}   "
            f"{a['x']:8.3f}{a['y']:8.3f}{a['z']:8.3f}{a['occupancy']:6.2f}{0.0:6.2f}          "
            f"{a['element']:>2s}"
        )
        lines.append(line)
        atom_serial += 1
    lines.append("END")

    with open(out_path, 'w') as f:
        f.write('\n'.join(lines) + '\n')

    mapping_rows = [
        {
            'INSTANCE': instance, 'PDB': pdb,
            'original_auth_asym_id': orig, 'temporary_chain_id': surrogate,
        }
        for orig, surrogate in surrogate_map.items()
    ]
    return mapping_rows
