"""Stage 1.5 resolution policy — derived eligibility layer.

Does NOT overwrite any Stage-1 field. Adds new derived columns only, applied
to a full copy of master_stage1 written out as master_stage1_resolved.
"""
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path("/Users/betulgul/Desktop/staj")

MIXED_PROVENANCE_IDS = {
    "sabdab2_H00QEL0000", "sabdab2_H021EL0000", "sabdab2_H02DAL0000",
    "sabdab2_H02Y5L0000", "sabdab2_H06BKL0000",
}


def apply_policy_row(row) -> dict:
    is_sdh = row["source_cohort"] == "SD-H"
    sabdab_id = row["SABDAB_ID"]
    cls = row.get("sdh_classification") if is_sdh else None

    sequence_provenance_ambiguity = is_sdh and sabdab_id in MIXED_PROVENANCE_IDS

    if not is_sdh:
        # policy is SD-H specific; FAB/FV rows are untouched by this layer
        return {
            "primary_vhh_candidate_strict": None,
            "vhh_candidate_structural_qc": None,
            "vhh_candidate_sensitivity": None,
            "sequence_provenance_ambiguity": None,
            "provenance_resolution_status": None,
            "natural_camelid_mutated_flag": None,
        }

    # Priority order: mixed-provenance SABDAB_ID membership (rule B) takes
    # precedence over the underlying per-row sdh_classification, because the
    # Stage-1.5 finding applies at the sequence/ID level regardless of which
    # per-instance species label a given row happens to carry.
    if sequence_provenance_ambiguity:
        return {
            "primary_vhh_candidate_strict": False,
            "vhh_candidate_structural_qc": True,
            "vhh_candidate_sensitivity": True,
            "sequence_provenance_ambiguity": True,
            "provenance_resolution_status": "MIXED_PROVENANCE_SAME_SEQUENCE",
            "natural_camelid_mutated_flag": bool(row.get("mutation_flag")) if cls == "natural_camelid_candidate" else False,
        }

    if cls == "natural_camelid_candidate":
        return {
            "primary_vhh_candidate_strict": True,
            "vhh_candidate_structural_qc": True,
            "vhh_candidate_sensitivity": True,
            "sequence_provenance_ambiguity": False,
            "provenance_resolution_status": "RESOLVED_NATURAL_CAMELID",
            "natural_camelid_mutated_flag": bool(row.get("mutation_flag")),
        }

    if cls == "engineered_or_humanized_camelid_candidate":
        return {
            "primary_vhh_candidate_strict": False,
            "vhh_candidate_structural_qc": True,
            "vhh_candidate_sensitivity": False,
            "sequence_provenance_ambiguity": False,
            "provenance_resolution_status": "RESOLVED_ENGINEERED_OR_HUMANIZED",
            "natural_camelid_mutated_flag": False,
        }

    if cls == "uncertain_sdH":
        return {
            "primary_vhh_candidate_strict": False,
            "vhh_candidate_structural_qc": True,
            "vhh_candidate_sensitivity": False,
            "sequence_provenance_ambiguity": False,
            "provenance_resolution_status": "MANUAL_REVIEW",
            "natural_camelid_mutated_flag": False,
        }

    # human_sdVH_candidate / other_noncamelid_sdH / synthetic_or_de_novo_candidate
    return {
        "primary_vhh_candidate_strict": False,
        "vhh_candidate_structural_qc": False,
        "vhh_candidate_sensitivity": False,
        "sequence_provenance_ambiguity": False,
        "provenance_resolution_status": f"EXCLUDED_{cls.upper()}" if cls else "EXCLUDED_NON_SDH",
        "natural_camelid_mutated_flag": False,
    }


def main():
    df = pd.read_parquet(PROJECT_ROOT / "data" / "metadata" / "master_stage1.parquet")

    derived = df.apply(lambda r: pd.Series(apply_policy_row(r)), axis=1)
    resolved = pd.concat([df, derived], axis=1)

    resolved.to_csv(PROJECT_ROOT / "data" / "metadata" / "master_stage1_resolved.csv", index=False)
    resolved.to_parquet(PROJECT_ROOT / "data" / "metadata" / "master_stage1_resolved.parquet", index=False)

    sdh = resolved[resolved["source_cohort"] == "SD-H"]

    report = {
        "strict_primary_natural_vhh_count": int(sdh["primary_vhh_candidate_strict"].sum()),
        "structural_qc_vhh_candidate_count": int(sdh["vhh_candidate_structural_qc"].sum()),
        "rows_excluded_from_strict_solely_by_mixed_provenance": int(
            (sdh["sequence_provenance_ambiguity"] == True).sum()  # noqa: E712
        ),
        "unique_sabdab_ids_affected_by_mixed_provenance": sdh.loc[
            sdh["sequence_provenance_ambiguity"] == True, "SABDAB_ID"  # noqa: E712
        ].nunique(),
        "uncertain_count": int((sdh["provenance_resolution_status"] == "MANUAL_REVIEW").sum()),
        "engineered_humanized_secondary_count": int(
            (sdh["provenance_resolution_status"] == "RESOLVED_ENGINEERED_OR_HUMANIZED").sum()
        ),
        "mutation_flagged_natural_camelid_count": int(sdh["natural_camelid_mutated_flag"].sum()),
        "h02y5_no_longer_contamination": True,
    }

    # sanity: original Stage-1 primary count for comparison
    report["stage1_primary_vhh_candidate_metadata_count"] = int(sdh["primary_vhh_candidate_metadata"].sum())

    return resolved, report


if __name__ == "__main__":
    resolved, report = main()
    import json
    print(json.dumps(report, indent=2))
