"""
Stage 8 — APBS OpenDX potential grid reader + exposed-atom potential sampling.

Sampling scheme (explicit, per Section 28): trilinear interpolation of the
APBS potential grid at the (x,y,z) coordinates of solvent-exposed
former-interface side-chain atoms. This is atom-coordinate sampling, NOT a
full molecular-surface mesh integral -- outputs are named
"...potential sampled at exposed interface atoms", never "surface potential".
"""
import numpy as np
from typing import Tuple, List, Dict


def read_dx(path: str) -> Dict:
    with open(path) as f:
        lines = f.readlines()

    counts = None
    origin = None
    deltas = []
    data_start = None
    n_items = None
    for i, line in enumerate(lines):
        if line.startswith('object 1 class gridpositions'):
            counts = tuple(int(x) for x in line.split()[-3:])
        elif line.startswith('origin'):
            origin = tuple(float(x) for x in line.split()[1:4])
        elif line.startswith('delta'):
            deltas.append(tuple(float(x) for x in line.split()[1:4]))
        elif line.startswith('object 3 class array'):
            n_items = int(line.split('items')[1].split()[0])
            data_start = i + 1
            break

    values = []
    for line in lines[data_start:]:
        if line.startswith('attribute'):
            break
        values.extend(float(x) for x in line.split())
    values = np.array(values[:n_items], dtype=float).reshape(counts)

    delta = np.array([deltas[0][0], deltas[1][1], deltas[2][2]])
    origin = np.array(origin)

    return {'values': values, 'origin': origin, 'delta': delta, 'counts': counts}


def trilinear_sample(grid: Dict, xyz: Tuple[float, float, float]) -> float:
    origin, delta, values, counts = grid['origin'], grid['delta'], grid['values'], grid['counts']
    rel = (np.array(xyz) - origin) / delta
    i0 = np.floor(rel).astype(int)
    frac = rel - i0

    if np.any(i0 < 0) or np.any(i0 >= np.array(counts) - 1):
        return None  # outside grid -- explicit None, never silently 0

    v000 = values[i0[0], i0[1], i0[2]]
    v100 = values[i0[0] + 1, i0[1], i0[2]]
    v010 = values[i0[0], i0[1] + 1, i0[2]]
    v001 = values[i0[0], i0[1], i0[2] + 1]
    v110 = values[i0[0] + 1, i0[1] + 1, i0[2]]
    v101 = values[i0[0] + 1, i0[1], i0[2] + 1]
    v011 = values[i0[0], i0[1] + 1, i0[2] + 1]
    v111 = values[i0[0] + 1, i0[1] + 1, i0[2] + 1]

    fx, fy, fz = frac
    v00 = v000 * (1 - fx) + v100 * fx
    v01 = v001 * (1 - fx) + v101 * fx
    v10 = v010 * (1 - fx) + v110 * fx
    v11 = v011 * (1 - fx) + v111 * fx
    v0 = v00 * (1 - fy) + v10 * fy
    v1 = v01 * (1 - fy) + v11 * fy
    return float(v0 * (1 - fz) + v1 * fz)


def sample_potentials(grid: Dict, coords: List[Tuple[float, float, float]]) -> List[float]:
    return [trilinear_sample(grid, c) for c in coords]
