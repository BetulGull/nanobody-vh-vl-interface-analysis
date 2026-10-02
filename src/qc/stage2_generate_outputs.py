"""
Stage-2 Output Generation

Loads final checkpoints, merges with input metadata, builds the master QC
table, tier cohorts, failure-reason tables, antigen chain map, coverage
summary, and cohort funnel.
"""

import sys
from pathlib import Path
import pandas as pd
import numpy as np
import logging

sys.path.insert(0, str(Path(__file__).parent.parent))

PROJECT_ROOT = Path(__file__).parent.parent.parent
CHECKPOINT_DIR = PROJECT_ROOT / "data" / "interim" / "stage2" / "checkpoints"
PROCESSED_DIR = PROJECT_ROOT / "data" / "processed"
TABLES_DIR = PROJECT_ROOT / "results" / "tables"
AUDITS_DIR = PROJECT_ROOT / "results" / "audits"

for d in [PROCESSED_DIR, TABLES_DIR, AUDITS_DIR]:
    d.mkdir(parents=True, exist_ok=True)

logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger(__name__)


def load_final_checkpoint(prefix: str, expected_n: int) -> pd.DataFrame:
    """Load the highest-numbered checkpoint for a cohort (the full run)."""
    checkpoints = sorted(CHECKPOINT_DIR.glob(f"{prefix}_*.parquet"))
    if not checkpoints:
        raise FileNotFoundError(f"No checkpoints found for prefix {prefix}")
    final = checkpoints[-1]
    df = pd.read_parquet(final)
    logger.info(f"Loaded {final.name}: {len(df)} rows")
    if len(df) != expected_n:
        logger.warning(f"Expected {expected_n} rows, got {len(df)}")
    return df


def main():
    metadata = pd.read_parquet(PROJECT_ROOT / "data" / "metadata" / "master_stage1_resolved.parquet")

    vhh_input = metadata[metadata["vhh_candidate_structural_qc"] == True].copy().reset_index(drop=True)
    paired_input = metadata[metadata["primary_paired_candidate_metadata"] == True].copy().reset_index(drop=True)

    vhh_qc = load_final_checkpoint("vhh_checkpoint", len(vhh_input))
    paired_qc = load_final_checkpoint("paired_checkpoint", len(paired_input))

    # Merge QC results back onto full metadata rows (preserve all Stage-1 columns)
    vhh_merged = vhh_input.merge(vhh_qc, on="INSTANCE", how="left", suffixes=("", "_qc"))
    paired_merged = paired_input.merge(paired_qc, on="INSTANCE", how="left", suffixes=("", "_qc"))

    # Master table = union of both cohorts (INSTANCE should not overlap between VHH/paired pools)
    master = pd.concat([vhh_merged, paired_merged], ignore_index=True)
    dup_instances = master["INSTANCE"].value_counts()
    dup_instances = dup_instances[dup_instances > 1]
    if len(dup_instances) > 0:
        logger.warning(f"Found {len(dup_instances)} INSTANCE rows present in both VHH and paired pools (unexpected)")

    master_path_parquet = PROCESSED_DIR / "stage2_structure_qc_master.parquet"
    master_path_csv = PROCESSED_DIR / "stage2_structure_qc_master.csv"
    master.to_parquet(master_path_parquet, index=False)
    master.to_csv(master_path_csv, index=False)
    logger.info(f"Wrote master table: {len(master)} rows -> {master_path_parquet.name}")

    # --- VHH tier tables ---
    vhh_merged.to_parquet(PROCESSED_DIR / "stage2_vhh_qc_all.parquet", index=False)

    vhh_strict_primary = vhh_merged[
        (vhh_merged["primary_vhh_candidate_strict"] == True) &
        (vhh_merged["qc_stage2_core_structure_pass"] == True) &
        (vhh_merged["qc_gordon_cdr_antigen_contact"] == True)
    ].copy()
    vhh_strict_primary.to_parquet(PROCESSED_DIR / "stage2_vhh_strict_primary.parquet", index=False)
    logger.info(f"VHH strict-primary structural cohort: {len(vhh_strict_primary)} rows")

    vhh_secondary_structural = vhh_merged[
        (vhh_merged["primary_vhh_candidate_strict"] != True) &
        (vhh_merged["qc_stage2_core_structure_pass"] == True)
    ].copy()
    vhh_secondary_structural.to_parquet(PROCESSED_DIR / "stage2_vhh_secondary_structural.parquet", index=False)
    logger.info(f"VHH secondary/sensitivity structural cohort: {len(vhh_secondary_structural)} rows")

    # --- Paired tier tables ---
    paired_merged.to_parquet(PROCESSED_DIR / "stage2_paired_qc_all.parquet", index=False)

    paired_primary = paired_merged[
        (paired_merged["qc_stage2_core_structure_pass"] == True) &
        (paired_merged["qc_gordon_cdr_antigen_contact"] == True)
    ].copy()
    paired_primary.to_parquet(PROCESSED_DIR / "stage2_paired_primary.parquet", index=False)
    logger.info(f"Paired primary structural cohort: {len(paired_primary)} rows")

    # --- Antigen chain map ---
    antigen_rows = []
    for _, row in master.iterrows():
        if row.get("antigen_chains"):
            for achain in str(row["antigen_chains"]).split(";"):
                if achain:
                    antigen_rows.append({
                        "INSTANCE": row["INSTANCE"],
                        "PDB": row["PDB"],
                        "cohort": row.get("cohort", ""),
                        "antigen_chain": achain,
                        "n_antigen_chains_total": row.get("n_antigen_chains", 0),
                    })
    antigen_map_df = pd.DataFrame(antigen_rows)
    antigen_map_df.to_parquet(PROCESSED_DIR / "stage2_antigen_chain_map.parquet", index=False)
    antigen_map_df.to_csv(PROCESSED_DIR / "stage2_antigen_chain_map.csv", index=False)
    logger.info(f"Antigen chain map: {len(antigen_map_df)} chain-level rows")

    # --- Failure reasons table ---
    failure_rows = master[master["failure_reasons"].notna() & (master["failure_reasons"] != "")]
    failure_reasons_expanded = []
    for _, row in failure_rows.iterrows():
        for reason in str(row["failure_reasons"]).split(";"):
            if reason:
                failure_reasons_expanded.append({
                    "INSTANCE": row["INSTANCE"],
                    "PDB": row["PDB"],
                    "cohort": row.get("cohort", ""),
                    "failure_reason": reason,
                })
    failure_df = pd.DataFrame(failure_reasons_expanded)
    failure_df.to_csv(TABLES_DIR / "stage2_failure_reasons.csv", index=False)
    logger.info(f"Failure reasons: {len(failure_df)} reason-level rows from {len(failure_rows)} failed/flagged structures")

    # Failure reason frequency summary
    if len(failure_df) > 0:
        freq = failure_df["failure_reason"].value_counts()
        logger.info("Top failure reasons:")
        for reason, count in freq.head(10).items():
            logger.info(f"  {reason}: {count}")

    # --- Failed structures (core structure pass = False) ---
    failed_structures = master[master["qc_stage2_core_structure_pass"] != True].copy()
    failed_structures.to_csv(AUDITS_DIR / "stage2_failed_structures.csv", index=False)
    logger.info(f"Failed structures (core pass=False): {len(failed_structures)}")

    # --- Missing residues table (from CDR completeness flags — position-level not tracked
    #     per-position in master row; report which CDRs were incomplete per structure) ---
    missing_rows = []
    for _, row in master.iterrows():
        for cdr_flag, cdr_name in [
            ("qc_h_cdr1_complete", "H1"), ("qc_h_cdr2_complete", "H2"), ("qc_h_cdr3_complete", "H3"),
            ("qc_l_cdr1_complete", "L1"), ("qc_l_cdr2_complete", "L2"), ("qc_l_cdr3_complete", "L3"),
        ]:
            val = row.get(cdr_flag)
            if val is False:
                missing_rows.append({"INSTANCE": row["INSTANCE"], "PDB": row["PDB"], "cdr": cdr_name})
    missing_df = pd.DataFrame(missing_rows)
    missing_df.to_csv(TABLES_DIR / "stage2_missing_residues.csv", index=False)
    logger.info(f"Missing/incomplete CDR records: {len(missing_df)}")

    # --- Coverage summary ---
    coverage_rows = []
    for cohort_name, df in [("VHH", vhh_merged), ("PAIRED", paired_merged)]:
        n = len(df)
        coverage_rows.append({
            "cohort": cohort_name,
            "n_input": n,
            "n_structure_resolved": (df["qc_structure_resolved"] == True).sum(),
            "n_hchain_found": (df["qc_hchain_exists"] == True).sum(),
            "n_h_cdr_complete": (df["qc_h_all_cdr_complete"] == True).sum(),
            "n_has_antigen": (df["qc_has_candidate_protein_antigen"] == True).sum(),
            "n_gordon_contact": (df["qc_gordon_cdr_antigen_contact"] == True).sum(),
            "n_core_structure_pass": (df["qc_stage2_core_structure_pass"] == True).sum(),
            "n_gordon_complex_pass": (df["qc_stage2_gordon_complex_pass"] == True).sum(),
        })
    coverage_df = pd.DataFrame(coverage_rows)
    coverage_df.to_csv(TABLES_DIR / "stage2_structural_coverage_summary.csv", index=False)
    logger.info("\nCoverage summary:")
    logger.info(coverage_df.to_string(index=False))

    # --- Cohort funnel ---
    funnel_rows = [
        {"stage": "Stage-1 VHH structural-QC pool", "n": len(vhh_input)},
        {"stage": "Stage-2 VHH core structure pass", "n": (vhh_merged["qc_stage2_core_structure_pass"] == True).sum()},
        {"stage": "Stage-2 VHH Gordon complex pass", "n": (vhh_merged["qc_stage2_gordon_complex_pass"] == True).sum()},
        {"stage": "Stage-2 VHH strict-primary final", "n": len(vhh_strict_primary)},
        {"stage": "Stage-1 paired primary pool", "n": len(paired_input)},
        {"stage": "Stage-2 paired core structure pass", "n": (paired_merged["qc_stage2_core_structure_pass"] == True).sum()},
        {"stage": "Stage-2 paired Gordon complex pass", "n": (paired_merged["qc_stage2_gordon_complex_pass"] == True).sum()},
        {"stage": "Stage-2 paired primary final", "n": len(paired_primary)},
    ]
    funnel_df = pd.DataFrame(funnel_rows)
    funnel_df.to_csv(TABLES_DIR / "stage2_cohort_funnel.csv", index=False)
    logger.info("\nCohort funnel:")
    logger.info(funnel_df.to_string(index=False))

    return {
        "master": master,
        "vhh_merged": vhh_merged,
        "paired_merged": paired_merged,
        "vhh_strict_primary": vhh_strict_primary,
        "vhh_secondary_structural": vhh_secondary_structural,
        "paired_primary": paired_primary,
        "coverage_df": coverage_df,
        "funnel_df": funnel_df,
        "failure_df": failure_df,
    }


if __name__ == "__main__":
    main()
