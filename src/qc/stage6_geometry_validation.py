"""
Stage 6 — Gordon/Zhang Geometry Validation (spec §23, §25)

Validates rigid-body translation and rotation invariance on deterministic
real structures before full-cohort geometry calculation.
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd
import copy
import logging

sys.path.insert(0, str(Path(__file__).parent.parent))

from qc.mmcif_qc import MMCIFStructure
from qc.stage6_geometry import compute_cdr3_length, compute_gordon_rho_compactness, compute_zhang_fr2_cdr3_dis

PROJECT_ROOT = Path(__file__).parent.parent.parent
AUDITS_DIR = PROJECT_ROOT / "results" / "audits"

logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger(__name__)


def rotation_matrix(axis, theta):
    axis = axis / np.linalg.norm(axis)
    a = np.cos(theta / 2.0)
    b, c, d = -axis * np.sin(theta / 2.0)
    return np.array([
        [a*a+b*b-c*c-d*d, 2*(b*c+a*d), 2*(b*d-a*c)],
        [2*(b*c-a*d), a*a+c*c-b*b-d*d, 2*(c*d+a*b)],
        [2*(b*d+a*c), 2*(c*d-a*b), a*a+d*d-b*b-c*c],
    ])


def make_transformed_structure(structure: MMCIFStructure, translate: np.ndarray, rot: np.ndarray) -> MMCIFStructure:
    """Return a deep-copied structure with all _atom_site coordinates rigidly transformed."""
    transformed = copy.deepcopy(structure)
    for row in transformed.data.get('_atom_site', []):
        x = float(row['_atom_site.Cartn_x'])
        y = float(row['_atom_site.Cartn_y'])
        z = float(row['_atom_site.Cartn_z'])
        new_xyz = rot @ np.array([x, y, z]) + translate
        row['_atom_site.Cartn_x'] = str(new_xyz[0])
        row['_atom_site.Cartn_y'] = str(new_xyz[1])
        row['_atom_site.Cartn_z'] = str(new_xyz[2])
    return transformed


def validate_case(instance: str, master: pd.DataFrame, category: str) -> dict:
    row = master[master['INSTANCE'] == instance].iloc[0]
    hchain = row['Hchain']
    path = PROJECT_ROOT / row['extracted_structure_path']

    s_orig = MMCIFStructure(str(path))
    model_id = s_orig.get_model()

    length_info = compute_cdr3_length(s_orig, hchain, model_id)
    gordon_orig = compute_gordon_rho_compactness(s_orig, hchain, model_id, length_info['cdr3_sequence_length'])
    zhang_orig = compute_zhang_fr2_cdr3_dis(s_orig, hchain, model_id)

    # Rigid transform: translate by (1000, -500, 250), rotate 37 degrees about an arbitrary axis
    translate = np.array([1000.0, -500.0, 250.0])
    rot = rotation_matrix(np.array([0.3, 0.7, 0.2]), np.radians(37))
    s_transformed = make_transformed_structure(s_orig, translate, rot)

    gordon_transformed = compute_gordon_rho_compactness(s_transformed, hchain, model_id, length_info['cdr3_sequence_length'])
    zhang_transformed = compute_zhang_fr2_cdr3_dis(s_transformed, hchain, model_id)

    rho_diff = (
        abs(gordon_orig['cdr3_rho'] - gordon_transformed['cdr3_rho'])
        if gordon_orig['cdr3_rho'] is not None and gordon_transformed['cdr3_rho'] is not None else None
    )
    compactness_diff = (
        abs(gordon_orig['cdr3_compactness'] - gordon_transformed['cdr3_compactness'])
        if gordon_orig['cdr3_compactness'] is not None and gordon_transformed['cdr3_compactness'] is not None else None
    )
    fr2_diff = (
        abs(zhang_orig['FR2_CDR3_Dis'] - zhang_transformed['FR2_CDR3_Dis'])
        if zhang_orig['FR2_CDR3_Dis'] is not None and zhang_transformed['FR2_CDR3_Dis'] is not None else None
    )

    return {
        'category': category, 'INSTANCE': instance,
        'cdr3_sequence_length': length_info['cdr3_sequence_length'],
        'gordon_valid_original': gordon_orig['gordon_geometry_valid'],
        'rho_original': gordon_orig['cdr3_rho'], 'rho_transformed': gordon_transformed['cdr3_rho'],
        'rho_invariance_diff': rho_diff,
        'compactness_original': gordon_orig['cdr3_compactness'], 'compactness_transformed': gordon_transformed['cdr3_compactness'],
        'compactness_invariance_diff': compactness_diff,
        'zhang_valid_original': zhang_orig['zhang_metric_valid'],
        'FR2_CDR3_Dis_original': zhang_orig['FR2_CDR3_Dis'], 'FR2_CDR3_Dis_transformed': zhang_transformed['FR2_CDR3_Dis'],
        'FR2_CDR3_Dis_invariance_diff': fr2_diff,
    }


def main():
    vhh_master = pd.read_parquet(PROJECT_ROOT / "data/processed/stage4_vhh_sequence_representatives.parquet")

    cases = [
        ('pdb_00009imx-C', 'short_CDR3'),
        ('pdb_00003k3q-A', 'long_CDR3_insertion_rich'),
        ('pdb_00006oca-D', 'long_CDR3_2'),
        ('pdb_00007m1h-G', 'medium_CDR3_low_coverage'),
        ('pdb_00008qf4-E', 'high_resolution'),
        ('pdb_00001zmy-A', 'low_resolution'),
        ('pdb_00003zkx-C', 'low_coverage'),
        ('pdb_00009gv3-E', 'long_CDR3_3'),
        ('pdb_00006zbp-FFF', 'unusual_chain_id'),
        ('pdb_00006ul6-B1', 'unusual_chain_id_2'),
    ]

    results = [validate_case(inst, vhh_master, cat) for inst, cat in cases]
    df = pd.DataFrame(results)
    out_path = AUDITS_DIR / "stage6_gordon_geometry_validation.csv"
    df.to_csv(out_path, index=False)
    logger.info(f"Wrote {out_path}")
    logger.info(df[['category', 'INSTANCE', 'gordon_valid_original', 'rho_invariance_diff',
                     'compactness_invariance_diff', 'zhang_valid_original', 'FR2_CDR3_Dis_invariance_diff']].to_string(index=False))

    gate_pass = (
        df.loc[df['gordon_valid_original'], 'rho_invariance_diff'].max() < 1e-6 and
        df.loc[df['gordon_valid_original'], 'compactness_invariance_diff'].max() < 1e-6 and
        (df.loc[df['zhang_valid_original'], 'FR2_CDR3_Dis_invariance_diff'].max() < 1e-6
         if df['zhang_valid_original'].any() else True)
    )
    logger.info(f"\nGEOMETRY RIGID-BODY INVARIANCE GATE: {'PASS' if gate_pass else 'FAIL'}")
    return gate_pass, df


if __name__ == "__main__":
    ok, _ = main()
    sys.exit(0 if ok else 1)
