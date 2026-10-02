"""
Stage 7 — Minimal mmCIF Writer for Arpeggio Input

Arpeggio 1.4.4's internal reader (arpeggio.core.protein_reader) calls
gemmi.cif.read() directly on the input filename regardless of nominal file
extension, i.e. it actually requires mmCIF-format input (confirmed
empirically via a real PDB-format smoke test failure, not assumed).
Using mmCIF input conveniently sidesteps PDB's single-character chain-ID
and fixed-width limitations entirely -- but a residue-identity mapping
table is still generated for the audit trail and for full round-trip
verification, per instruction.
"""

import string
from typing import List, Dict


def write_cif_with_surrogate_mapping(atoms: List[Dict], out_path: str, instance: str, pdb: str) -> List[Dict]:
    surrogate_map: Dict[str, str] = {}
    available = iter(string.ascii_uppercase)
    for a in atoms:
        if a['auth_asym_id'] not in surrogate_map:
            surrogate_map[a['auth_asym_id']] = next(available)

    unique_residues = sorted(set(a['residue_name'] for a in atoms))
    chem_comp_lines = [
        "loop_",
        "_chem_comp.id",
        "_chem_comp.type",
        "_chem_comp.name",
    ]
    for res in unique_residues:
        # All residues here are canonical/noncanonical amino acids from the
        # VHH variable domain (never water/ion/heterogen -- those are already
        # excluded by extract_variable_domain_atoms). "L-PEPTIDE LINKING" is
        # the standard mmCIF _chem_comp.type for amino acids, mapped by
        # Arpeggio's config.ComponentType.from_chem_comp_type to type 1
        # (Polypeptide) -- confirmed by direct inspection of
        # arpeggio/core/config.py, not assumed.
        chem_comp_lines.append(f"{res} 'L-peptide linking' {res}")
    chem_comp_lines.append("#")

    lines = [
        "data_stage7input",
        "#",
    ] + chem_comp_lines + [
        "loop_",
        "_atom_site.group_PDB",
        "_atom_site.id",
        "_atom_site.type_symbol",
        "_atom_site.label_atom_id",
        "_atom_site.label_alt_id",
        "_atom_site.label_comp_id",
        "_atom_site.label_asym_id",
        "_atom_site.label_entity_id",
        "_atom_site.label_seq_id",
        "_atom_site.pdbx_PDB_ins_code",
        "_atom_site.Cartn_x",
        "_atom_site.Cartn_y",
        "_atom_site.Cartn_z",
        "_atom_site.occupancy",
        "_atom_site.B_iso_or_equiv",
        "_atom_site.pdbx_formal_charge",
        "_atom_site.auth_seq_id",
        "_atom_site.auth_asym_id",
        "_atom_site.pdbx_PDB_model_num",
    ]

    atom_serial = 1
    for a in atoms:
        surrogate_chain = surrogate_map[a['auth_asym_id']]
        icode = a['ins_code'] if a['ins_code'] else '?'
        lines.append(
            f"ATOM {atom_serial} {a['element']} {a['atom_name']} . {a['residue_name']} "
            f"{surrogate_chain} 1 {atom_serial} {icode} "
            f"{a['x']:.3f} {a['y']:.3f} {a['z']:.3f} {a['occupancy']:.2f} 0.0 ? "
            f"{a['auth_seq_id']} {surrogate_chain} 1"
        )
        atom_serial += 1

    with open(out_path, 'w') as f:
        f.write('\n'.join(lines) + '\n')

    mapping_rows = [
        {'INSTANCE': instance, 'PDB': pdb, 'original_auth_asym_id': orig, 'temporary_chain_id': surrogate}
        for orig, surrogate in surrogate_map.items()
    ]
    return mapping_rows
