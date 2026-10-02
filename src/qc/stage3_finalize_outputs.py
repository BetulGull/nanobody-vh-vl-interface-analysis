"""
Stage 3 — Finalize Outputs from Checkpoints

Loads the final (highest-numbered) checkpoint files and writes the official
Stage-3 core outputs.
"""

import sys
from pathlib import Path
import pandas as pd
import logging

PROJECT_ROOT = Path(__file__).parent.parent.parent
CHECKPOINT_DIR = PROJECT_ROOT / "data" / "interim" / "stage3" / "checkpoints"
PROCESSED_DIR = PROJECT_ROOT / "data" / "processed"

logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger(__name__)


def load_final_checkpoint(prefix: str) -> pd.DataFrame:
    checkpoints = sorted(CHECKPOINT_DIR.glob(f"{prefix}_*.parquet"))
    if not checkpoints:
        raise FileNotFoundError(f"No checkpoints for {prefix}")
    final = checkpoints[-1]
    df = pd.read_parquet(final)
    logger.info(f"Loaded {final.name}: {len(df)} rows")
    return df


def main():
    residue_df = load_final_checkpoint("residues")
    summary_df = load_final_checkpoint("summaries")

    residue_out = PROCESSED_DIR / "stage3_vh_vl_delta_sasa_residues.parquet"
    residue_df.to_parquet(residue_out, index=False)
    logger.info(f"Wrote {residue_out}: {len(residue_df)} rows")

    summary_out_parquet = PROCESSED_DIR / "stage3_structure_delta_sasa_summary.parquet"
    summary_out_csv = PROCESSED_DIR / "stage3_structure_delta_sasa_summary.csv"
    summary_df.to_parquet(summary_out_parquet, index=False)
    summary_df.to_csv(summary_out_csv, index=False)
    logger.info(f"Wrote {summary_out_parquet} / .csv: {len(summary_df)} rows")

    n_pass = summary_df['sasa_qc_pass'].sum()
    logger.info(f"QC pass: {n_pass}/{len(summary_df)}")

    failures = summary_df[summary_df['sasa_qc_pass'] != True]
    failures_out = PROJECT_ROOT / "results" / "audits" / "stage3_sasa_failures.csv"
    failures.to_csv(failures_out, index=False)
    logger.info(f"Wrote {failures_out}: {len(failures)} failed/reviewed structures")

    if len(failures) > 0:
        logger.info("Failure reasons:")
        logger.info(failures['failure_or_review_reason'].value_counts().to_string())

    return residue_df, summary_df


if __name__ == "__main__":
    main()
