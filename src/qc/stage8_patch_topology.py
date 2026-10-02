"""
Stage 8 — PROJECT-SPECIFIC SURFACE PATCH TOPOLOGY (explicitly not TAP/TNP).

Builds a spatial residue graph over solvent-exposed hydrophobic residues
(same exposure rule/threshold as the source-verified PSH machinery), using
the TAP/TNP source-confirmed 4.0A all-heavy-atom edge criterion as the
PRIMARY threshold (reused without contradiction, since the source code
already defines this exact spatial radius for its own neighborhood
construction) plus two sensitivity thresholds (3.5A, 4.5A) predeclared
BEFORE inspecting any group result.
"""
import sys
from pathlib import Path
from typing import List, Dict, Set, Tuple
sys.path.insert(0, str(Path(__file__).parent))
from stage8_patch_scoring import KD_RAW, _min_heavy_atom_distance

HYDROPHOBIC_AA = {'ALA', 'VAL', 'LEU', 'ILE', 'PHE', 'MET', 'TRP', 'CYS'}  # KD > 0, source-consistent
THRESHOLDS = [3.5, 4.0, 4.5]
PRIMARY_THRESHOLD = 4.0


def _union_find_components(nodes: List, edges: List[Tuple]) -> List[Set]:
    parent = {n: n for n in nodes}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    for a, b in edges:
        union(a, b)

    comps = {}
    for n in nodes:
        r = find(n)
        comps.setdefault(r, set()).add(n)
    return list(comps.values())


def compute_topology(residues: Dict, key_set: Set, threshold: float, sasa_lookup: Dict = None) -> Dict:
    """key_set: residues to consider (already restricted to exposed + hydrophobic by caller)."""
    nodes = list(key_set)
    edges = []
    keys_list = list(key_set)
    for i, k1 in enumerate(keys_list):
        for k2 in keys_list[i + 1:]:
            d = _min_heavy_atom_distance(residues[k1]['atoms'], residues[k2]['atoms'])
            if d < threshold:
                edges.append((k1, k2))

    components = _union_find_components(nodes, edges)
    if not components:
        return {
            'number_of_hydrophobic_surface_components': 0, 'largest_patch_residue_count': 0,
            'largest_patch_exposed_sasa': 0.0, 'largest_patch_KD_area_sum': 0.0,
            'total_hydrophobic_patch_sasa': 0.0, 'fraction_hydrophobic_surface_in_largest_patch': None,
            'mean_component_size': 0.0,
        }

    total_sasa = sum(residues[k].get('sidechain_asa_abs') or 0.0 for k in key_set)
    comp_sizes = [len(c) for c in components]
    largest = max(components, key=len)
    largest_sasa = sum(residues[k].get('sidechain_asa_abs') or 0.0 for k in largest)
    largest_kd_sum = sum((residues[k].get('sidechain_asa_abs') or 0.0) * KD_RAW.get(residues[k]['residue_name'], 0.0) for k in largest)

    return {
        'number_of_hydrophobic_surface_components': len(components),
        'largest_patch_residue_count': len(largest),
        'largest_patch_exposed_sasa': largest_sasa,
        'largest_patch_KD_area_sum': largest_kd_sum,
        'total_hydrophobic_patch_sasa': total_sasa,
        'fraction_hydrophobic_surface_in_largest_patch': (largest_sasa / total_sasa) if total_sasa > 0 else None,
        'mean_component_size': sum(comp_sizes) / len(comp_sizes),
    }
