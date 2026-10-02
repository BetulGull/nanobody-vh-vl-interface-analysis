"""
Stage 8 — PDB2PQR + APBS pipeline: single-structure end-to-end run, potential
sampling at solvent-exposed former-interface atoms.
"""
import subprocess
import time
import shutil
from pathlib import Path
from typing import Dict, List, Tuple
import numpy as np

import sys
sys.path.insert(0, str(Path(__file__).parent))
from stage7_pdb_writer import write_pdb_with_surrogate_mapping
from stage8_apbs_sampling import read_dx, sample_potentials

PH = 7.4
PDB2PQR_CMD = ["micromamba", "run", "-n", "vhh-elec", "pdb2pqr30", "--ff", "PARSE",
               "--with-ph", str(PH), "--titration-state-method", "propka", "--keep-chain"]
APBS_CMD = ["micromamba", "run", "-n", "vhh-elec", "apbs"]


def run_apbs_pipeline(atoms: List[Dict], work_dir: str, instance: str, timeout: int = 300) -> Dict:
    work = Path(work_dir)
    work.mkdir(parents=True, exist_ok=True)
    pdb_path = work / f"{instance}.pdb"
    pqr_path = work / f"{instance}.pqr"
    in_path = work / f"{instance}.in"

    mapping = write_pdb_with_surrogate_mapping(atoms, str(pdb_path), instance, instance)

    t0 = time.time()
    result = subprocess.run(PDB2PQR_CMD + ["--apbs-input", in_path.name, pdb_path.name, pqr_path.name],
                             capture_output=True, text=True, timeout=timeout, cwd=str(work))
    t1 = time.time()
    pdb2pqr_runtime = t1 - t0
    if result.returncode != 0 or not pqr_path.exists():
        return {'status': 'PDB2PQR_FAILURE', 'stderr': result.stderr[-2000:], 'pdb2pqr_runtime_seconds': pdb2pqr_runtime}

    n_pqr_atoms = sum(1 for l in open(pqr_path) if l.startswith('ATOM'))

    t2 = time.time()
    apbs_result = subprocess.run(APBS_CMD + [in_path.name], capture_output=True, text=True, timeout=timeout, cwd=str(work))
    t3 = time.time()
    apbs_runtime = t3 - t2

    dx_path = work / f"{instance}.pqr.dx"
    if apbs_result.returncode != 0 or not dx_path.exists():
        return {'status': 'APBS_FAILURE', 'stderr': apbs_result.stderr[-2000:],
                'pdb2pqr_runtime_seconds': pdb2pqr_runtime, 'apbs_runtime_seconds': apbs_runtime,
                'n_pqr_atoms': n_pqr_atoms}

    grid = read_dx(str(dx_path))

    mol_coords = []
    with open(pqr_path) as f:
        for line in f:
            if line.startswith('ATOM'):
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                mol_coords.append((x, y, z))
    mol_coords = np.array(mol_coords)
    grid_min = grid['origin']
    grid_max = grid['origin'] + (np.array(grid['counts']) - 1) * grid['delta']
    contained = bool(np.all(mol_coords >= grid_min) and np.all(mol_coords <= grid_max))

    return {
        'status': 'PDB2PQR_SUCCESS_APBS_SUCCESS',
        'pdb2pqr_runtime_seconds': pdb2pqr_runtime, 'apbs_runtime_seconds': apbs_runtime,
        'n_pqr_atoms': n_pqr_atoms, 'grid_counts': grid['counts'], 'grid_delta': tuple(grid['delta']),
        'grid_origin': tuple(grid['origin']), 'molecule_contained_in_grid': contained,
        'mapping': mapping, 'pqr_path': str(pqr_path), 'dx_path': str(dx_path), 'grid': grid,
    }
