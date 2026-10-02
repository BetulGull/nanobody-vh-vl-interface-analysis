"""
Stage 9 — Antibody-Antigen Residue-Level Contact Computation

Primary geometric contact definition (source-verified, stage2_structure_qc.py /
contact_search.py, reused): minimum heavy-atom distance <= threshold. Primary
threshold 4.5A (Gordon-verified); sensitivity thresholds 4.0A and 5.0A also
computed for every pair.

Uses the frozen Stage-2 AntigenQuality (protein-polymer-only, sibling-antibody-
and multi-instance-masked) for antigen chain identification, and the same
IMGT region boundaries used throughout Stages 3-8.
"""
import sys
from pathlib import Path
from typing import List, Dict, Set, Tuple
import numpy as np
from scipy.spatial import cKDTree

sys.path.insert(0, str(Path(__file__).parent))
from stage3_variable_domain import extract_variable_domain_atoms
from antigen_qc import AntigenQuality

IMGT_REGIONS = {
    'FR1': (1, 26), 'CDR1': (27, 38), 'FR2': (39, 55), 'CDR2': (56, 65),
    'FR3': (66, 104), 'CDR3': (105, 117), 'FR4': (118, 128),
}


def get_region(imgt_number: int) -> str:
    for region, (start, end) in IMGT_REGIONS.items():
        if start <= imgt_number <= end:
            return region
    return 'UNKNOWN'


def get_antigen_atoms_as_list(structure, antigen_chains: List[str], model_id: str) -> List[Dict]:
    aq = AntigenQuality()
    atom_map = aq.get_antigen_atoms(structure, antigen_chains, model_id)

    # AntigenQuality.get_antigen_atoms filters chains at the entity level
    # (protein-polymer only, via identify_antigen_chains) but does NOT filter
    # individual atom records within a legitimate protein-antigen chain --
    # so water/heteroatoms sharing the same author chain ID (a real, common
    # mmCIF pattern: crystallographic waters are often assigned to a nearby
    # polymer chain's auth_asym_id) were being silently included as "antigen"
    # atoms. Discovered via a real Arpeggio smoke test showing HOH residues
    # in the antigen selection (see stage9_arpeggio_validation notes) and
    # confirmed directly against raw _atom_site rows: HOH entries in an
    # antigen chain carry group_PDB=='HETATM' while genuine polymer residues
    # carry group_PDB=='ATOM'. Filtered here at the atom-record level, never
    # touching the frozen Stage-2 AntigenQuality entity-level chain filter.
    resname_lookup = {}
    valid_atom_keys = set()
    for row in structure.data.get('_atom_site', []):
        chain = row.get('_atom_site.auth_asym_id')
        if chain not in antigen_chains or row.get('_atom_site.pdbx_PDB_model_num') != model_id:
            continue
        resnum = int(row.get('_atom_site.auth_seq_id', '0'))
        icode = row.get('_atom_site.pdbx_PDB_ins_code', '').strip()
        icode = '' if icode in ('?', '.', '') else icode
        resname_lookup[(chain, resnum, icode)] = row.get('_atom_site.label_comp_id', '').strip()
        if row.get('_atom_site.group_PDB') == 'ATOM':
            valid_atom_keys.add((chain, resnum, icode))

    # altloc selection: same highest-occupancy/alphabetical-tiebreak policy as
    # extract_variable_domain_atoms (MMCIFStructure.select_altloc), applied
    # here per (chain, residue, atom_name) group -- required to avoid duplicate
    # atoms (e.g. downstream Arpeggio input construction rejects duplicates).
    by_key: Dict = {}
    for chain, chain_atoms in atom_map.items():
        for a in chain_atoms:
            res_key = (chain, a['residue_number'], a['insertion_code'])
            if res_key not in valid_atom_keys:
                continue  # HETATM (water, ion, crystallization additive) -- not a real antigen residue
            k = (chain, a['residue_number'], a['insertion_code'], a['name'])
            by_key.setdefault(k, []).append(a)

    atoms = []
    for (chain, resnum, icode, atom_name), conformers in by_key.items():
        conformers.sort(key=lambda x: (-x['occupancy'], x['altloc']))
        a = conformers[0]
        key = (chain, resnum, icode)
        atoms.append({
            'antigen_chain': chain, 'auth_seq_id': resnum, 'ins_code': icode,
            'atom_name': atom_name, 'x': a['x'], 'y': a['y'], 'z': a['z'],
            'residue_name': resname_lookup.get(key, ''),
        })
    return atoms


def compute_residue_contacts(antibody_atoms: List[Dict], antigen_atoms: List[Dict],
                              domain_label: str, thresholds=(4.0, 4.5, 5.0)) -> List[Dict]:
    """
    For every antibody residue x antigen residue pair, compute minimum heavy-atom
    distance. Returns one row per PAIR THAT IS WITHIN THE LARGEST THRESHOLD
    (others are simply absent, not zeroed) -- with contact_at_4A/4p5A/5A booleans.
    """
    if not antibody_atoms or not antigen_atoms:
        return []

    max_thr = max(thresholds)
    ab_coords = np.array([[a['x'], a['y'], a['z']] for a in antibody_atoms])
    ag_coords = np.array([[a['x'], a['y'], a['z']] for a in antigen_atoms])

    ab_tree = cKDTree(ab_coords)
    ag_tree = cKDTree(ag_coords)
    pairs = ab_tree.query_ball_tree(ag_tree, r=max_thr)

    pair_min_dist: Dict[Tuple, float] = {}
    for i, ag_idx_list in enumerate(pairs):
        if not ag_idx_list:
            continue
        ab_atom = antibody_atoms[i]
        ab_key = (ab_atom['auth_asym_id'], ab_atom['auth_seq_id'], ab_atom['ins_code'])
        for j in ag_idx_list:
            ag_atom = antigen_atoms[j]
            ag_key = (ag_atom['antigen_chain'], ag_atom['auth_seq_id'], ag_atom['ins_code'])
            d = float(np.linalg.norm(ab_coords[i] - ag_coords[j]))
            key = (ab_key, ag_key)
            if key not in pair_min_dist or d < pair_min_dist[key]:
                pair_min_dist[key] = d

    rows = []
    ab_lookup = {}
    for a in antibody_atoms:
        k = (a['auth_asym_id'], a['auth_seq_id'], a['ins_code'])
        if k not in ab_lookup:
            ab_lookup[k] = a
    ag_lookup = {}
    for a in antigen_atoms:
        k = (a['antigen_chain'], a['auth_seq_id'], a['ins_code'])
        if k not in ag_lookup:
            ag_lookup[k] = a

    for (ab_key, ag_key), dist in pair_min_dist.items():
        ab_a = ab_lookup[ab_key]
        ag_a = ag_lookup[ag_key]
        rows.append({
            'antibody_domain': domain_label,
            'antibody_imgt_position': ab_key[1], 'antibody_ins_code': ab_key[2],
            'antibody_region': get_region(ab_key[1]), 'antibody_residue': ab_a['residue_name'],
            'antigen_chain': ag_key[0], 'antigen_residue_id': ag_key[1], 'antigen_ins_code': ag_key[2],
            'antigen_residue_name': ag_a.get('residue_name', ''),
            'minimum_heavy_atom_distance': dist,
            'contact_at_4A': dist <= 4.0,
            'contact_at_4p5A': dist <= 4.5,
            'contact_at_5A': dist <= 5.0,
        })
    return rows
