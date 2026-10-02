"""
Stage 5 — Reproduce Frozen Stage-3 VH-Alone SASA (spec §9, §38)

Independently recomputes VH-alone SASA for a deterministic sample of
conventional VH representatives using Stage-5 code paths (which are the SAME
Stage-3 extraction/SASA modules, reused not reimplemented), and compares
against the frozen Stage-3 per-residue values. If this disagrees materially,
STOP -- do not silently replace Stage-3 values.
"""

import sys
from pathlib import Path
import pandas as pd
import logging

sys.path.insert(0, str(Path(__file__).parent.parent))

from qc.mmcif_qc import MMCIFStructure
from qc.stage3_variable_domain import extract_variable_domain_atoms
from qc.stage3_sasa import calc_sasa

PROJECT_ROOT = Path(__file__).parent.parent.parent
AUDITS_DIR = PROJECT_ROOT / "results" / "audits"

logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger(__name__)

TOLERANCE = 1e-6  # Å² -- same underlying deterministic calculation


def main():
    vh_master = pd.read_parquet(PROJECT_ROOT / "data/processed/stage3_paired_representatives.parquet")
    frozen_residues = pd.read_parquet(PROJECT_ROOT / "data/processed/stage3_vh_vl_delta_sasa_residues.parquet")

    sample = vh_master.head(12)
    rows = []

    for _, row in sample.iterrows():
        instance = row['INSTANCE']
        hchain = row['Hchain']
        path = PROJECT_ROOT / row['extracted_structure_path']

        s = MMCIFStructure(str(path))
        model_id = s.get_model()
        vh_atoms = extract_variable_domain_atoms(s, hchain, model_id)
        per_res, total_recomputed = calc_sasa(vh_atoms)

        frozen_rows = frozen_residues[frozen_residues['INSTANCE'] == instance]
        frozen_total = frozen_rows['sasa_vh_alone'].sum()

        # Per-residue comparison
        n_matched = 0
        n_mismatched = 0
        max_diff = 0.0
        for _, frow in frozen_rows.iterrows():
            key = (hchain, f"{frow['imgt_position']}{frow['insertion_code']}")
            recomputed_val = per_res.get(key)
            if recomputed_val is None:
                continue
            diff = abs(recomputed_val - frow['sasa_vh_alone'])
            max_diff = max(max_diff, diff)
            if diff < TOLERANCE:
                n_matched += 1
            else:
                n_mismatched += 1

        rows.append({
            'INSTANCE': instance,
            'frozen_total_vh_alone_sasa': frozen_total,
            'recomputed_total_vh_alone_sasa': total_recomputed,
            'total_diff': abs(frozen_total - total_recomputed),
            'n_residues_matched': n_matched,
            'n_residues_mismatched': n_mismatched,
            'max_per_residue_diff': max_diff,
        })

    df = pd.DataFrame(rows)
    out_path = AUDITS_DIR / "stage5_stage3_vh_sasa_reproduction.csv"
    df.to_csv(out_path, index=False)
    logger.info(f"Wrote {out_path}")
    logger.info(df.to_string(index=False))

    gate_pass = (df['total_diff'] < TOLERANCE).all() and (df['n_residues_mismatched'] == 0).all()
    logger.info(f"\nSTAGE-3 VH SASA REPRODUCTION GATE: {'PASS' if gate_pass else 'FAIL'}")
    return gate_pass, df


if __name__ == "__main__":
    ok, _ = main()
    sys.exit(0 if ok else 1)
