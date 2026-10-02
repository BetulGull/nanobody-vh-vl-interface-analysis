"""
Stage 8 — TAP-style PSH/PPC/PNC Patch Scoring (source-faithful)

Formula and thresholds recovered directly from external/TNP's
Hydrophobicity_and_Charge_Assigner.py, which the TNP paper explicitly states
reproduces methods from TAP (Raybould et al. 2019 PNAS) -- see
results/audits/stage8_tap_source_gate.md for full derivation and evidence.

Relative side-chain ASA is obtained by running the actual author-provided
`psa` binary (external/TNP/bin/psa_mac) directly, since its internal
reference max-ASA table is compiled into the binary and not independently
recoverable as text -- using the real tool is more source-faithful than
approximating its denominator with a different published table.
"""
import subprocess
import string
from pathlib import Path
from typing import List, Dict, Set, Tuple

PSA_BIN = "/Users/betulgul/Desktop/staj/external/TNP/bin/psa_mac"

# Kyte & Doolittle (1982) raw values, confirmed byte-identical to
# external/TNP/lib/python/theraprofnano/Hydrophobicity_and_Charge_Profiler/dat/Hydrophobics.txt column 0
KD_RAW = {
    'ILE': 4.5, 'VAL': 4.2, 'LEU': 3.8, 'PHE': 2.8, 'CYS': 2.5, 'MET': 1.9, 'ALA': 1.8,
    'GLY': -0.4, 'THR': -0.7, 'SER': -0.8, 'TRP': -0.9, 'TYR': -1.3, 'PRO': -1.6,
    'HIS': -3.2, 'GLU': -3.5, 'GLN': -3.5, 'ASP': -3.5, 'ASN': -3.5, 'LYS': -3.9, 'ARG': -4.5,
}
_kd_min, _kd_max = min(KD_RAW.values()), max(KD_RAW.values())
KD_NORM = {aa: (v - _kd_min) / (_kd_max - _kd_min) + 1.0 for aa, v in KD_RAW.items()}

# ASP/GLU=-1, ARG/LYS=+1, HIS=+0.1 (pH 7.4 Henderson-Hasselbalch approximation, source-confirmed)
CHARGE_CLASS = {'ASP': -1.0, 'GLU': -1.0, 'ARG': 1.0, 'LYS': 1.0, 'HIS': 0.1}

EXPOSURE_THRESHOLD_PCT = 7.5  # source-confirmed: side-chain relative ASA >= 7.5%
VICINITY_DISTANCE_A = 4.0     # source-confirmed: exposed-neighbor inclusion radius
PATCH_PAIR_CUTOFF_A = 7.5     # source-confirmed: pairwise patch-scoring distance cutoff
SALTBRIDGE_CUTOFF_A = 3.2     # source-confirmed: donor-acceptor salt-bridge neutralization cutoff

DONOR_ATOMS = {'LYS': ['NZ'], 'ARG': ['NH1', 'NH2', 'NE']}
ACCEPTOR_ATOMS = {'ASP': ['OD1', 'OD2'], 'GLU': ['OE1', 'OE2']}


def write_minimal_pdb(atoms: List[Dict], out_path: str) -> Dict[str, str]:
    """Single fixed chain 'H'; reversible surrogate mapping still recorded for
    the original (possibly multi-character) auth_asym_id, per project audit policy."""
    surrogate_map: Dict[str, str] = {}
    available = iter(string.ascii_uppercase)
    for a in atoms:
        if a['auth_asym_id'] not in surrogate_map:
            surrogate_map[a['auth_asym_id']] = next(available)

    lines = []
    serial = 1
    for a in atoms:
        icode = a['ins_code'] if a['ins_code'] else ''
        line = (
            f"ATOM  {serial:5d} {a['atom_name']:<4s} {a['residue_name']:<3s} H"
            f"{a['auth_seq_id']:4d}{icode:1s}   "
            f"{a['x']:8.3f}{a['y']:8.3f}{a['z']:8.3f}{a['occupancy']:6.2f}{0.0:6.2f}          "
            f"{a['element']:>2s}"
        )
        lines.append(line)
        serial += 1
    lines.append("END")
    with open(out_path, 'w') as f:
        f.write('\n'.join(lines) + '\n')
    return surrogate_map


def run_psa(pdb_path: str, timeout: int = 30) -> str:
    result = subprocess.run([PSA_BIN, "-t", str(pdb_path)], capture_output=True, text=True, timeout=timeout)
    return result.stdout


def _split_resnum(field: str) -> Tuple[int, str]:
    field = field.strip()
    if field and field[-1].isalpha():
        return int(field[:-1]), field[-1]
    return int(field), ''


def parse_psa_output(stdout: str) -> Dict[Tuple[int, str], Dict]:
    """Returns {(auth_seq_id, ins_code): {'restype', 'sidechain_asa_abs', 'sidechain_rel_pct'}}"""
    result = {}
    for line in stdout.splitlines():
        if not line.startswith('ACCESS'):
            continue
        resnum_field = line[6:12]
        restype = line[12:18].strip()
        try:
            resnum, icode = _split_resnum(resnum_field)
            side_abs = float(line[55:61])
            side_pct = float(line[61:67])
        except (ValueError, IndexError):
            continue
        result[(resnum, icode)] = {'restype': restype, 'sidechain_asa_abs': side_abs, 'sidechain_rel_pct': side_pct}
    return result


def build_residue_registry(atoms: List[Dict], psa_data: Dict[Tuple[int, str], Dict]) -> Dict[Tuple[int, str], Dict]:
    """One entry per real residue present in `atoms`, joined with psa exposure data.
    Missing PSA lookup (should not occur for atoms actually written) is tracked explicitly, never zeroed."""
    residues = {}
    for a in atoms:
        key = (a['auth_seq_id'], a['ins_code'] if a['ins_code'] else '')
        if key not in residues:
            residues[key] = {
                'imgt_position': a['auth_seq_id'], 'ins_code': a['ins_code'] if a['ins_code'] else '',
                'residue_name': a['residue_name'], 'atoms': [],
            }
        residues[key]['atoms'].append(a)

    for key, r in residues.items():
        psa_entry = psa_data.get(key)
        if psa_entry is None:
            r['sidechain_asa_abs'] = None
            r['relative_sidechain_asa'] = None
            r['is_surface_exposed'] = None
            r['exposure_status'] = 'PSA_LOOKUP_MISSING'
            continue
        r['sidechain_asa_abs'] = psa_entry['sidechain_asa_abs']
        r['relative_sidechain_asa'] = psa_entry['sidechain_rel_pct']
        r['is_surface_exposed'] = psa_entry['sidechain_rel_pct'] >= EXPOSURE_THRESHOLD_PCT
        r['exposure_status'] = 'OK'

    for key, r in residues.items():
        aa = r['residue_name']
        r['kyte_doolittle_norm'] = KD_NORM.get(aa)
        r['charge_class'] = CHARGE_CLASS.get(aa, 0.0)
    return residues


def _min_heavy_atom_distance(atoms1: List[Dict], atoms2: List[Dict]) -> float:
    mindist = float('inf')
    for a1 in atoms1:
        for a2 in atoms2:
            d = ((a1['x'] - a2['x']) ** 2 + (a1['y'] - a2['y']) ** 2 + (a1['z'] - a2['z']) ** 2) ** 0.5
            if d < mindist:
                mindist = d
    return mindist


def apply_saltbridge_neutralization(residues: Dict[Tuple[int, str], Dict]) -> None:
    keys = list(residues.keys())
    for k1 in keys:
        r1 = residues[k1]
        if r1['residue_name'] not in DONOR_ATOMS:
            continue
        donor_atoms = [a for a in r1['atoms'] if a['atom_name'] in DONOR_ATOMS[r1['residue_name']]]
        if not donor_atoms:
            continue
        for k2 in keys:
            if k1 == k2:
                continue
            r2 = residues[k2]
            if r2['residue_name'] not in ACCEPTOR_ATOMS:
                continue
            acceptor_atoms = [a for a in r2['atoms'] if a['atom_name'] in ACCEPTOR_ATOMS[r2['residue_name']]]
            if not acceptor_atoms:
                continue
            d = _min_heavy_atom_distance(donor_atoms, acceptor_atoms)
            if d < SALTBRIDGE_CUTOFF_A:
                r1['charge_class'] = 0.0
                r1['kyte_doolittle_norm'] = KD_NORM['GLY']
                r2['charge_class'] = 0.0
                r2['kyte_doolittle_norm'] = KD_NORM['GLY']


def build_vicinity(residues: Dict[Tuple[int, str], Dict], core_keys: Set[Tuple[int, str]]) -> Set[Tuple[int, str]]:
    """core + any other exposed residue within VICINITY_DISTANCE_A (min all-heavy-atom distance) of a core residue."""
    exposed_keys = {k for k, r in residues.items() if r.get('is_surface_exposed')}
    core_exposed = core_keys & exposed_keys
    vicinity = set(core_exposed)
    for k2 in exposed_keys:
        if k2 in vicinity:
            continue
        for k1 in core_exposed:
            d = _min_heavy_atom_distance(residues[k1]['atoms'], residues[k2]['atoms'])
            if d < VICINITY_DISTANCE_A:
                vicinity.add(k2)
                break
    return vicinity


def compute_patch_scores(residues: Dict[Tuple[int, str], Dict], vicinity: Set[Tuple[int, str]]) -> Dict[str, float]:
    """PSH/PPC/PNC = sum over ordered pairs (r1,r2) in vicinity, r1!=r2, d<=7.5A, of score1*score2/d^2."""
    psh = ppc = pnc = 0.0
    vicinity_list = list(vicinity)
    for i, k1 in enumerate(vicinity_list):
        for k2 in vicinity_list:
            if k1 == k2:
                continue
            r1, r2 = residues[k1], residues[k2]
            d = _min_heavy_atom_distance(r1['atoms'], r2['atoms'])
            if d > PATCH_PAIR_CUTOFF_A or d == 0:
                continue
            if r1['kyte_doolittle_norm'] is not None and r2['kyte_doolittle_norm'] is not None:
                psh += (r1['kyte_doolittle_norm'] * r2['kyte_doolittle_norm']) / (d * d)
            c1, c2 = r1['charge_class'], r2['charge_class']
            if c1 > 0 and c2 > 0:
                ppc += (c1 * c2) / (d * d)
            elif c1 < 0 and c2 < 0:
                pnc += (c1 * c2) / (d * d)
    return {'PSH': psh, 'PPC': ppc, 'PNC': pnc}


def compute_patch_for_target(atoms: List[Dict], core_positions: Set[int]) -> Dict:
    """Full pipeline: write PDB, run psa, parse, build registry, salt-bridge neutralize,
    vicinity, patch scores -- for a given 'core' IMGT-position set (any insertion code at
    that integer position counts as core)."""
    import tempfile
    with tempfile.TemporaryDirectory() as td:
        pdb_path = f"{td}/input.pdb"
        write_minimal_pdb(atoms, pdb_path)
        stdout = run_psa(pdb_path)
        psa_data = parse_psa_output(stdout)

    residues = build_residue_registry(atoms, psa_data)
    apply_saltbridge_neutralization(residues)

    core_keys = {k for k in residues if k[0] in core_positions}
    vicinity_keys = build_vicinity(residues, core_keys)

    core_scores = compute_patch_scores(residues, core_keys & {k for k, r in residues.items() if r.get('is_surface_exposed')})
    vicinity_scores = compute_patch_scores(residues, vicinity_keys)

    return {
        'residues': residues,
        'core_keys': core_keys,
        'vicinity_keys': vicinity_keys,
        'core_scores': core_scores,
        'vicinity_scores': vicinity_scores,
    }
