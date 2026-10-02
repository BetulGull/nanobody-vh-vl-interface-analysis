"""
Stage 9 — Antibody-Antigen Arpeggio Secondary Interaction Layer.

New inter-partner use case vs Stage 7 (intradomain only): combined CIF with
both antibody (surrogate chain, first) and antigen (surrogate chain(s),
after) atoms, run with `-s /<ab_surrogate>//` to select the antibody as
"ligand" so interacting_entities distinguishes INTER (antibody-antigen,
the interaction of interest) from INTRA_SELECTION (within antibody) and
INTRA_NON_SELECTION (within antigen).
"""
import sys
import subprocess
import json
from pathlib import Path
from typing import Dict, List

sys.path.insert(0, str(Path(__file__).parent))
from stage7_cif_writer import write_cif_with_surrogate_mapping
from stage9_antigen_contacts import get_antigen_atoms_as_list, get_region

ARPEGGIO_CMD = ["micromamba", "run", "-n", "vhh-arpeggio", "pdbe-arpeggio"]


def run_arpeggio_antibody_antigen(ab_atoms: List[Dict], ag_atoms_raw: List[Dict], work_dir: str, instance: str, timeout: int = 60) -> Dict:
    work = Path(work_dir)
    work.mkdir(parents=True, exist_ok=True)
    cif_path = work / f"{instance}.cif"
    json_path = work / f"{instance}.json"

    ag_atoms = []
    for a in ag_atoms_raw:
        ag_atoms.append({
            'auth_asym_id': a['antigen_chain'], 'label_seq_id': 0, 'auth_seq_id': a['auth_seq_id'],
            'ins_code': a['ins_code'], 'residue_name': a['residue_name'], 'atom_name': a['atom_name'],
            'element': a['atom_name'][0], 'x': a['x'], 'y': a['y'], 'z': a['z'], 'occupancy': 1.0, 'altloc': '',
        })
    combined = ab_atoms + ag_atoms
    mapping = write_cif_with_surrogate_mapping(combined, str(cif_path), instance, instance)
    ab_surrogate = mapping[0]['temporary_chain_id']  # antibody atoms always written first

    if not json_path.exists():
        result = subprocess.run(ARPEGGIO_CMD + ["-o", ".", "-s", f"/{ab_surrogate}//", "--", cif_path.name],
                                 capture_output=True, text=True, timeout=timeout, cwd=str(work))
        if result.returncode != 0 or not json_path.exists():
            return {'status': 'ARPEGGIO_FAILURE', 'stderr': result.stderr[-1000:], 'mapping': mapping}

    with open(json_path) as f:
        data = json.load(f)

    return {'status': 'ARPEGGIO_SUCCESS', 'n_records': len(data), 'data': data, 'mapping': mapping, 'ab_surrogate': ab_surrogate}


def extract_inter_records(arpeggio_data: List[Dict], mapping: List[Dict], ab_surrogate: str) -> List[Dict]:
    """Filter to INTER records (antibody<->antigen), translate surrogate chain back to real IDs."""
    reverse_map = {m['temporary_chain_id']: m['original_auth_asym_id'] for m in mapping}
    rows = []
    for r in arpeggio_data:
        if r.get('interacting_entities') != 'INTER' or r.get('type') != 'atom-atom':
            continue
        bgn, end = r['bgn'], r['end']
        if bgn['auth_asym_id'] == ab_surrogate:
            ab_atom, ag_atom = bgn, end
        elif end['auth_asym_id'] == ab_surrogate:
            ab_atom, ag_atom = end, bgn
        else:
            continue
        contact = r['contact']
        rows.append({
            'antibody_imgt_position': ab_atom['auth_seq_id'], 'antibody_ins_code': ab_atom['pdbx_PDB_ins_code'].strip(),
            'antibody_region': get_region(ab_atom['auth_seq_id']),
            'antigen_chain': reverse_map.get(ag_atom['auth_asym_id'], ag_atom['auth_asym_id']),
            'antigen_residue_id': ag_atom['auth_seq_id'], 'antigen_ins_code': ag_atom['pdbx_PDB_ins_code'].strip(),
            'distance': r['distance'], 'contact_labels': ';'.join(contact),
            'is_proximal_only': contact == ['proximal'],
            'is_specific_interaction': any(l not in ('proximal', 'covalent') for l in contact),
            'is_covalent': 'covalent' in contact,
        })
    return rows
