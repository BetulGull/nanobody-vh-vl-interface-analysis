"""Stage 1 — generate summary tables and the QC audit markdown."""
import json
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path("/Users/betulgul/Desktop/staj")


def main():
    df = pd.read_parquet(PROJECT_ROOT / "data" / "metadata" / "master_stage1.parquet")

    tables_dir = PROJECT_ROOT / "results" / "tables"
    audits_dir = PROJECT_ROOT / "results" / "audits"

    # --- stage1_counts.csv ---
    counts = []
    for cohort in ["SD-H", "FAB", "FV"]:
        sub = df[df["source_cohort"] == cohort]
        counts.append({
            "source_cohort": cohort,
            "n_rows": len(sub),
            "unique_PDB": sub["PDB"].nunique(),
            "unique_SABDAB_ID": sub["SABDAB_ID"].nunique(),
            "n_xray": int(sub["qc_method_xray"].sum()),
            "n_resolution_le_3": int(sub["qc_resolution_le_3"].sum()),
            "n_xray_and_resolution_le_3": int((sub["qc_method_xray"] & sub["qc_resolution_le_3"]).sum()),
            "n_resolution_missing": int(sub["qc_resolution_numeric"].isna().sum()),
            "n_rfactor_missing": int((sub["r_factor"] == "NA").sum() + sub["r_factor"].isna().sum()),
            "n_rfree_missing": int((sub["r_free"] == "NA").sum() + sub["r_free"].isna().sum()),
        })
    counts_df = pd.DataFrame(counts)
    counts_df.to_csv(tables_dir / "stage1_counts.csv", index=False)

    # --- species counts (SD-H) ---
    sdh = df[df["source_cohort"] == "SD-H"]
    species_counts = sdh.groupby(["normalized_species", "species_evidence_field", "camelid_species_flag"]).size().reset_index(name="n")
    species_counts = species_counts.sort_values("n", ascending=False)
    species_counts.to_csv(tables_dir / "stage1_species_counts.csv", index=False)

    # --- SD-H classification counts ---
    sdh_class_counts = sdh["sdh_classification"].value_counts().reset_index()
    sdh_class_counts.columns = ["sdh_classification", "n"]
    sdh_class_counts.to_csv(tables_dir / "stage1_sdh_classification.csv", index=False)

    # --- construct counts (FAB/FV) ---
    paired = df[df["source_cohort"].isin(["FAB", "FV"])]
    construct_counts = paired.groupby(["source_cohort", "construct_bucket"]).size().reset_index(name="n")
    construct_counts.to_csv(tables_dir / "stage1_construct_counts.csv", index=False)

    # --- exclusion reasons ---
    reason_rows = []
    for cohort in ["SD-H", "FAB", "FV"]:
        sub = df[df["source_cohort"] == cohort]
        reason_counter = Counter()
        for reasons_str in sub["exclusion_reasons"]:
            for r in reasons_str.split(";"):
                reason_counter[r] += 1
        for reason, n in reason_counter.items():
            reason_rows.append({"source_cohort": cohort, "exclusion_reason": reason, "n_rows": n})
    exclusion_df = pd.DataFrame(reason_rows).sort_values(["source_cohort", "n_rows"], ascending=[True, False])
    exclusion_df.to_csv(tables_dir / "stage1_exclusion_reasons.csv", index=False)

    # --- redundancy inventory ---
    redundancy_rows = []
    for cohort in ["SD-H", "FAB", "FV", "ALL"]:
        sub = df if cohort == "ALL" else df[df["source_cohort"] == cohort]
        per_id = sub.groupby("SABDAB_ID").size()
        redundancy_rows.append({
            "group": cohort,
            "n_instance_rows": len(sub),
            "unique_PDB": sub["PDB"].nunique(),
            "unique_SABDAB_ID": sub["SABDAB_ID"].nunique(),
            "unique_HEAVY_ID": sub["HEAVY_ID"].nunique(),
            "instances_per_SABDAB_ID_median": float(per_id.median()),
            "instances_per_SABDAB_ID_mean": float(per_id.mean()),
            "instances_per_SABDAB_ID_max": int(per_id.max()),
            "instances_per_SABDAB_ID_p95": float(np.percentile(per_id, 95)),
        })
    redundancy_df = pd.DataFrame(redundancy_rows)
    redundancy_df.to_csv(tables_dir / "stage1_redundancy_inventory.csv", index=False)

    top20 = df.groupby("SABDAB_ID").size().sort_values(ascending=False).head(20)
    top20_df = top20.reset_index()
    top20_df.columns = ["SABDAB_ID", "n_instances"]
    top20_df.to_csv(tables_dir / "stage1_top20_repeated_sabdab_ids.csv", index=False)

    # --- known example sanity checks ---
    known_examples = {}
    for pdb_short, expected in [
        ("1bzq", "natural Camelus dromedarius VHH candidate"),
        ("2vyr", "human single-domain VH; must NOT enter primary natural VHH candidate cohort"),
        ("27uh", "synthetic/de-novo category"),
        ("1g6v", "should fail resolution <=3.0 A"),
    ]:
        matches = sdh[sdh["PDB"].str.contains(pdb_short, case=False, na=False)]
        if len(matches) == 0:
            known_examples[pdb_short] = {"present": False}
        else:
            r = matches.iloc[0]
            known_examples[pdb_short] = {
                "present": True,
                "n_rows": len(matches),
                "expected": expected,
                "actual_sdh_classification": r["sdh_classification"],
                "actual_normalized_species": r["normalized_species"],
                "actual_resolution": r["qc_resolution_numeric"],
                "actual_resolution_le_3": bool(r["qc_resolution_le_3"]),
                "actual_primary_vhh_candidate": bool(r["primary_vhh_candidate_metadata"]),
            }

    with open(audits_dir / "stage1_known_example_checks.json", "w") as f:
        json.dump(known_examples, f, indent=2, default=str)

    # --- headline numbers for the audit doc / final report ---
    headline = {
        "total_rows": len(df),
        "sdh_rows": len(sdh),
        "fab_rows": len(df[df["source_cohort"] == "FAB"]),
        "fv_rows": len(df[df["source_cohort"] == "FV"]),
        "sdh_classification_counts": sdh["sdh_classification"].value_counts().to_dict(),
        "primary_vhh_candidate_count": int(df["primary_vhh_candidate_metadata"].sum()),
        "fv_construct_counts": df[df["source_cohort"] == "FV"]["construct_bucket"].value_counts().to_dict(),
        "fab_construct_counts": df[df["source_cohort"] == "FAB"]["construct_bucket"].value_counts().to_dict(),
        "primary_paired_candidate_overall": int(df["primary_paired_candidate_metadata"].sum()),
        "primary_paired_candidate_fab": int(df[df["source_cohort"] == "FAB"]["primary_paired_candidate_metadata"].sum()),
        "primary_paired_candidate_fv": int(df[df["source_cohort"] == "FV"]["primary_paired_candidate_metadata"].sum()),
        "unique_pdb_overall": df["PDB"].nunique(),
        "unique_sabdab_id_overall": df["SABDAB_ID"].nunique(),
        "n_rfactor_missing_overall": int((df["r_factor"] == "NA").sum()),
        "n_rfree_missing_overall": int((df["r_free"] == "NA").sum()),
        "duplicate_archive_rows_with_flag": int(df["content_duplicate_flag"].sum()),
    }
    with open(audits_dir / "stage1_headline_numbers.json", "w") as f:
        json.dump(headline, f, indent=2, default=str)

    print(json.dumps(headline, indent=2, default=str))


if __name__ == "__main__":
    main()
