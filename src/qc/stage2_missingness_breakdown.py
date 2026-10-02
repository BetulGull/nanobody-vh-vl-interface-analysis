"""
Stage-2 Final Semantic Audit — Granular CDR-Missingness Breakdown

Replaces the vague "163 flags changed" figure with actual biological
missingness counts: residues, loops, chains, INSTANCE rows, and a per-CDR
(H1/H2/H3/L1/L2/L3) breakdown.
"""

import sys
from pathlib import Path
import pandas as pd
import logging

sys.path.insert(0, str(Path(__file__).parent.parent))

from qc.mmcif_qc import MMCIFStructure
from qc.cdr_qc import CDRQuality

PROJECT_ROOT = Path(__file__).parent.parent.parent
AUDITS_DIR = PROJECT_ROOT / "results" / "audits"

logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger(__name__)


def normalize_lchain(raw):
    if pd.isna(raw) or str(raw).strip() == 'NA':
        return None
    return raw


def main():
    master = pd.read_parquet(PROJECT_ROOT / "data" / "processed" / "stage2_structure_qc_master.parquet")
    cdrq = CDRQuality()

    residue_rows = []
    n_errors = 0

    for _, row in master.iterrows():
        if not row['qc_structure_resolved'] or not row['qc_hchain_exists']:
            continue

        hchain = row['Hchain']
        lchain = normalize_lchain(row['Lchain'])
        path = PROJECT_ROOT / row['extracted_structure_path']

        try:
            s = MMCIFStructure(str(path))
            model_id = s.get_model()
        except Exception as e:
            n_errors += 1
            continue

        h_missing = cdrq.get_missing_cdr_positions(s, hchain, 'H', model_id)
        for cdr_name, missing_list in h_missing.items():
            for m in missing_list:
                residue_rows.append({
                    'INSTANCE': row['INSTANCE'], 'PDB': row['PDB'], 'chain': hchain,
                    'chain_role': 'H', 'CDR': cdr_name,
                    'label_seq_id': m['label_seq_id'], 'mon_id': m['mon_id'], 'reason': m['reason'],
                })

        if lchain is not None and row['qc_lchain_exists']:
            l_missing = cdrq.get_missing_cdr_positions(s, lchain, 'L', model_id)
            for cdr_name, missing_list in l_missing.items():
                for m in missing_list:
                    residue_rows.append({
                        'INSTANCE': row['INSTANCE'], 'PDB': row['PDB'], 'chain': lchain,
                        'chain_role': 'L', 'CDR': cdr_name,
                        'label_seq_id': m['label_seq_id'], 'mon_id': m['mon_id'], 'reason': m['reason'],
                    })

    residue_df = pd.DataFrame(residue_rows)
    out_path = AUDITS_DIR / "stage2_cdr_missingness_breakdown.csv"
    residue_df.to_csv(out_path, index=False)
    logger.info(f"Wrote {out_path} ({len(residue_df)} missing-residue rows)")
    logger.info(f"Parse errors during recomputation: {n_errors}")

    if len(residue_df) == 0:
        logger.info("No missing CDR residues found.")
        return residue_df

    total_missing_residues = len(residue_df)
    affected_loops = residue_df.drop_duplicates(subset=['INSTANCE', 'CDR'])
    n_affected_loops = len(affected_loops)
    affected_chains = residue_df.drop_duplicates(subset=['INSTANCE', 'chain'])
    n_affected_chains = len(affected_chains)
    n_affected_instances = residue_df['INSTANCE'].nunique()

    per_cdr_residue_counts = residue_df['CDR'].value_counts().to_dict()
    per_cdr_loop_counts = affected_loops['CDR'].value_counts().to_dict()

    summary = {
        'total_missing_real_CDR_residues': total_missing_residues,
        'affected_CDR_loops (INSTANCE,CDR pairs)': n_affected_loops,
        'affected_antibody_chains (INSTANCE,chain pairs)': n_affected_chains,
        'affected_INSTANCE_rows': n_affected_instances,
    }
    for cdr in ['H1', 'H2', 'H3', 'L1', 'L2', 'L3']:
        summary[f'residues_missing_{cdr}'] = per_cdr_residue_counts.get(cdr, 0)
        summary[f'loops_affected_{cdr}'] = per_cdr_loop_counts.get(cdr, 0)

    summary_df = pd.DataFrame([summary])
    summary_out = AUDITS_DIR / "stage2_cdr_missingness_summary.csv"
    summary_df.to_csv(summary_out, index=False)
    logger.info(f"Wrote {summary_out}")
    logger.info("\n" + summary_df.T.to_string(header=False))

    return residue_df, summary_df


if __name__ == "__main__":
    main()
