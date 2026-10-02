"""
Stage-2 Closure Audit — Independent vs Sequential Counts, H/L CDR Contingency
Table, and Direct Boolean-Intersection Final Cohort Counts.
"""

import sys
from pathlib import Path
import pandas as pd
import logging

sys.path.insert(0, str(Path(__file__).parent.parent))

PROJECT_ROOT = Path(__file__).parent.parent.parent
TABLES_DIR = PROJECT_ROOT / "results" / "tables"

logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger(__name__)


def main():
    master = pd.read_parquet(PROJECT_ROOT / "data" / "processed" / "stage2_structure_qc_master.parquet")
    vhh = master[master["cohort"] == "VHH"]
    paired = master[master["cohort"] == "PAIRED"]

    # --- Section 3: independent vs sequential counts ---
    def independent_vs_sequential(df, label):
        n_total = len(df)
        n_core = (df["qc_stage2_core_structure_pass"] == True).sum()
        n_gordon_independent = (df["qc_gordon_cdr_antigen_contact"] == True).sum()
        n_both = ((df["qc_stage2_core_structure_pass"] == True) & (df["qc_gordon_cdr_antigen_contact"] == True)).sum()
        return {
            "cohort": label,
            "N_total": n_total,
            "N_core_structure_pass": n_core,
            "N_gordon_contact_pass_independent_of_core": n_gordon_independent,
            "N_both_core_and_gordon": n_both,
        }

    vhh_counts = independent_vs_sequential(vhh, "VHH")
    vhh_counts["N_strict_primary_metadata_and_core_and_gordon"] = (
        (vhh["primary_vhh_candidate_strict"] == True) &
        (vhh["qc_stage2_core_structure_pass"] == True) &
        (vhh["qc_gordon_cdr_antigen_contact"] == True)
    ).sum()

    paired_counts = independent_vs_sequential(paired, "PAIRED")
    paired_counts["N_final_paired_primary"] = (
        (paired["qc_stage2_core_structure_pass"] == True) &
        (paired["qc_gordon_cdr_antigen_contact"] == True)
    ).sum()

    reconciliation_df = pd.DataFrame([vhh_counts, paired_counts])
    out1 = TABLES_DIR / "stage2_independent_vs_sequential_counts.csv"
    reconciliation_df.to_csv(out1, index=False)
    logger.info(f"Wrote {out1}")
    logger.info(reconciliation_df.to_string(index=False))

    # --- Section 4: H/L CDR completeness contingency table (paired cohort) ---
    h_complete = paired["qc_h_all_cdr_complete"] == True
    l_complete = paired["qc_l_all_cdr_complete"] == True

    contingency = {
        "both_H_and_L_complete": (h_complete & l_complete).sum(),
        "H_incomplete_only": ((~h_complete) & l_complete).sum(),
        "L_incomplete_only": (h_complete & (~l_complete)).sum(),
        "both_incomplete": ((~h_complete) & (~l_complete)).sum(),
        "H_all_CDR_complete_total": h_complete.sum(),
        "H_any_CDR_incomplete_total": (~h_complete).sum(),
        "L_all_CDR_complete_total": l_complete.sum(),
        "L_any_CDR_incomplete_total": (~l_complete).sum(),
        "n_paired_total": len(paired),
    }
    contingency_df = pd.DataFrame([contingency])
    out2 = TABLES_DIR / "stage2_hl_cdr_contingency_table.csv"
    contingency_df.to_csv(out2, index=False)
    logger.info(f"\nWrote {out2}")
    logger.info(contingency_df.to_string(index=False))

    # Sanity check: the four contingency cells must sum to n_paired_total
    assert (contingency["both_H_and_L_complete"] + contingency["H_incomplete_only"] +
            contingency["L_incomplete_only"] + contingency["both_incomplete"]) == len(paired), \
        "Contingency cells do not sum to total paired count"

    # --- Section 5: direct boolean-intersection final cohort counts ---
    strict_vhh_mask = (
        (master["primary_vhh_candidate_strict"] == True) &
        (master["qc_stage2_core_structure_pass"] == True) &
        (master["qc_gordon_cdr_antigen_contact"] == True)
    )
    n_strict_vhh = strict_vhh_mask.sum()

    paired_final_mask = (
        (master["cohort"] == "PAIRED") &
        (master["qc_stage2_core_structure_pass"] == True) &
        (master["qc_gordon_cdr_antigen_contact"] == True)
    )
    n_paired_final = paired_final_mask.sum()

    vhh_secondary_mask = (
        (master["cohort"] == "VHH") &
        (master["primary_vhh_candidate_strict"] != True) &
        (master["qc_stage2_core_structure_pass"] == True)
    )
    n_vhh_secondary = vhh_secondary_mask.sum()

    final_counts = {
        "strict_VHH_final_N_direct_intersection": int(n_strict_vhh),
        "paired_primary_final_N_direct_intersection": int(n_paired_final),
        "vhh_secondary_sensitivity_final_N_direct_intersection": int(n_vhh_secondary),
    }
    logger.info(f"\nDirect boolean-intersection final counts: {final_counts}")

    # Cross-check against the parquet tier files on disk (must match exactly)
    vhh_strict_saved = pd.read_parquet(PROJECT_ROOT / "data/processed/stage2_vhh_strict_primary.parquet")
    paired_primary_saved = pd.read_parquet(PROJECT_ROOT / "data/processed/stage2_paired_primary.parquet")
    vhh_secondary_saved = pd.read_parquet(PROJECT_ROOT / "data/processed/stage2_vhh_secondary_structural.parquet")

    logger.info(f"Saved-file cross-check: strict_vhh={len(vhh_strict_saved)}, "
                f"paired_primary={len(paired_primary_saved)}, vhh_secondary={len(vhh_secondary_saved)}")

    assert len(vhh_strict_saved) == n_strict_vhh, "strict VHH count mismatch vs saved file"
    assert len(paired_primary_saved) == n_paired_final, "paired primary count mismatch vs saved file"
    assert len(vhh_secondary_saved) == n_vhh_secondary, "VHH secondary count mismatch vs saved file"
    logger.info("Cross-check PASS: direct intersection counts match saved tier files exactly")

    return reconciliation_df, contingency_df, final_counts


if __name__ == "__main__":
    main()
