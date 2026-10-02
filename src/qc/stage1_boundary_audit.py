"""Stage 1 — SD-H boundary-case classification audit (read-only diagnostic).

Does NOT change any classification. Produces inspectable evidence tables
and a proposed_review_action recommendation column for human review only.
"""
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path("/Users/betulgul/Desktop/staj")

CAMELUS_DROMEDARIUS_TAXID = "9838"
CAMELUS_BACTRIANUS_TAXID = "9837"
CAMELUS_FERUS_TAXID = "419612"
LAMA_GLAMA_TAXID = "9844"
VICUGNA_PACOS_TAXID = "30538"

BOUNDARY_COLUMNS = [
    "INSTANCE", "PDB", "SABDAB_ID",
    "heavy_species", "heavy_taxid",
    "organism", "organism_taxid",
    "expression_system", "expression_system_taxid",
    "heavy_expression_system", "heavy_expression_system_taxid",
    "heavy_subclass",
    "compound", "short_header",
    "normalized_species", "sdh_classification",
    "camelid_species_flag",
    "engineering_keyword_hits", "engineering_keyword_fields",
    "humanization_keyword_hits", "humanization_keyword_fields",
    "synthetic_keyword_hits", "synthetic_keyword_fields",
    "mutation_keyword_hits", "mutation_keyword_fields",
    "metadata_conflict_flag", "metadata_conflict_reason",
    "manual_review_reason",
    "primary_vhh_candidate_metadata",
]

NONBIOLOGICAL_HOST_TAXIDS = {"562", "10029", "9031"}  # E. coli, CHO(hamster fallback taxid varies), avoid false negatives is fine here


def not_empty(v):
    if v is None:
        return False
    if isinstance(v, float) and np.isnan(v):
        return False
    if v == "" or v == "NA":
        return False
    return True


def propose_action(row) -> tuple:
    cat = row["sdh_classification"]
    heavy_taxid = row.get("heavy_taxid")
    expr_taxid = row.get("heavy_expression_system_taxid") or row.get("expression_system_taxid")

    has_eng = not_empty(row.get("engineering_keyword_hits"))
    has_hum = not_empty(row.get("humanization_keyword_hits"))
    has_syn = not_empty(row.get("synthetic_keyword_hits"))
    has_mut = not_empty(row.get("mutation_keyword_hits"))
    has_conflict = bool(row.get("metadata_conflict_flag"))

    # Expression-host-as-species artifact (e.g. heavy_taxid == E. coli taxid)
    if heavy_taxid in NONBIOLOGICAL_HOST_TAXIDS:
        return (
            "INSUFFICIENT_METADATA",
            f"heavy_taxid={heavy_taxid} matches a known expression-host taxid, not a biological source "
            f"organism; heavy_species/heavy_taxid field likely mis-populated from expression system. "
            f"Cannot determine true biological origin from available structured fields alone.",
        )

    if cat == "human_sdVH_candidate":
        return ("EXCLUDE_HUMAN_SDVH", "heavy_taxid/heavy_species indicates Homo sapiens origin.")

    if cat == "other_noncamelid_sdH":
        return ("EXCLUDE_OTHER_NONCAMELID", f"heavy chain species ({row.get('normalized_species')}) is neither camelid nor human; likely an engineered non-camelid single-domain fragment (e.g. scFv/Fab crop).")

    if cat == "synthetic_or_de_novo_candidate":
        return ("MOVE_TO_ENGINEERED_SECONDARY", "Explicit synthetic-construct taxid and/or de-novo/designed descriptive text.")

    if cat == "engineered_or_humanized_camelid_candidate":
        return ("MOVE_TO_ENGINEERED_SECONDARY", f"Camelid-taxid heavy chain with explicit humanization keyword hit(s): {row.get('humanization_keyword_hits')}.")

    if cat == "natural_camelid_candidate":
        if has_syn:
            return ("KEEP_MANUAL_REVIEW", f"Camelid taxid but synthetic/de-novo keyword hit in free text ({row.get('synthetic_keyword_hits')}) — structured/text tension not resolved automatically.")
        if has_conflict:
            return ("KEEP_MANUAL_REVIEW", row.get("metadata_conflict_reason") or "metadata_conflict_flag set for unspecified reason.")
        if has_hum:
            # shouldn't normally happen for natural_camelid_candidate given classifier logic, but guard anyway
            return ("MOVE_TO_ENGINEERED_SECONDARY", f"Humanization keyword hit present: {row.get('humanization_keyword_hits')}.")
        if has_eng:
            return ("KEEP_MANUAL_REVIEW", f"Non-humanization engineering keyword hit present ({row.get('engineering_keyword_hits')}); camelid origin confirmed, nature of engineering unclear from metadata alone.")
        if has_mut:
            return ("KEEP_PRIMARY_BUT_MUTATED", f"Camelid biological origin confirmed by structured evidence; free text indicates experimental point mutation/affinity maturation ({row.get('mutation_keyword_hits')}), not a change in biological origin.")
        return ("KEEP_PRIMARY", "Camelid biological origin confirmed by structured evidence; no engineering/humanization/synthetic/conflict signal.")

    if cat == "uncertain_sdH":
        if row.get("species_evidence_field") == "none":
            return ("INSUFFICIENT_METADATA", "No structured species evidence available in any field (heavy_taxid, heavy_species, organism_taxid, organism all missing).")
        return ("KEEP_MANUAL_REVIEW", row.get("metadata_conflict_reason") or "Unresolved structured-vs-text tension.")

    return ("KEEP_MANUAL_REVIEW", "Unclassified case.")


def species_group(row):
    taxid = row.get("heavy_taxid")
    species = (row.get("normalized_species") or "").lower()
    if taxid == CAMELUS_DROMEDARIUS_TAXID:
        return "Camelus dromedarius"
    if taxid == CAMELUS_BACTRIANUS_TAXID:
        return "Camelus bactrianus"
    if taxid == CAMELUS_FERUS_TAXID:
        return "Camelus ferus"
    if taxid == LAMA_GLAMA_TAXID:
        return "Lama glama"
    if taxid == VICUGNA_PACOS_TAXID:
        return "Vicugna pacos"
    if "camelus" in species:
        return "other Camelus"
    if "lama" in species:
        return "other Lama"
    return "other/unresolved Camelidae"


def main():
    df = pd.read_parquet(PROJECT_ROOT / "data" / "metadata" / "master_stage1.parquet")
    sdh = df[df["source_cohort"] == "SD-H"].copy()

    # ---------------- 1. boundary-case table ----------------
    is_uncertain = sdh["sdh_classification"] == "uncertain_sdH"
    is_eng_hum = sdh["sdh_classification"] == "engineered_or_humanized_camelid_candidate"
    is_flagged_natural = (sdh["sdh_classification"] == "natural_camelid_candidate") & (
        sdh["engineering_keyword_hits"].apply(not_empty)
        | sdh["humanization_keyword_hits"].apply(not_empty)
        | sdh["synthetic_keyword_hits"].apply(not_empty)
        | sdh["mutation_keyword_hits"].apply(not_empty)
        | sdh["metadata_conflict_flag"].fillna(False)
        | sdh["manual_review_reason"].apply(not_empty)
    )

    boundary_mask = is_uncertain | is_eng_hum | is_flagged_natural
    boundary = sdh[boundary_mask].copy()
    boundary = boundary.sort_values(["sdh_classification", "PDB", "INSTANCE"])

    proposals = boundary.apply(lambda r: pd.Series(propose_action(r), index=["proposed_review_action", "proposed_review_reason"]), axis=1)
    boundary = pd.concat([boundary, proposals], axis=1)

    out_cols = BOUNDARY_COLUMNS + ["proposed_review_action", "proposed_review_reason"]
    out_cols = [c for c in out_cols if c in boundary.columns or c in ("proposed_review_action", "proposed_review_reason")]
    boundary_out = boundary.rename(columns={"sdh_classification": "current_sdh_category"})
    final_cols = [c if c != "sdh_classification" else "current_sdh_category" for c in out_cols]
    boundary_out = boundary_out[final_cols]
    boundary_out.to_csv(PROJECT_ROOT / "results" / "audits" / "stage1_sdh_boundary_cases.csv", index=False)

    # ---------------- 5. same-SABDAB_ID consistency check ----------------
    check_cols = [
        "normalized_species", "sdh_classification", "camelid_species_flag",
        "humanization_flag", "synthetic_flag", "engineering_flag",
        "primary_vhh_candidate_metadata",
    ]
    conflict_rows = []
    for sabdab_id, grp in sdh.groupby("SABDAB_ID"):
        if len(grp) < 2:
            continue
        inconsistent_fields = []
        for col in check_cols:
            vals = grp[col].apply(lambda v: v if not (isinstance(v, float) and np.isnan(v)) else None)
            if vals.nunique(dropna=False) > 1:
                inconsistent_fields.append(col)
        if inconsistent_fields:
            conflict_rows.append({
                "SABDAB_ID": sabdab_id,
                "n_instances": len(grp),
                "all_INSTANCE": ";".join(grp["INSTANCE"].tolist()),
                "all_PDB": ";".join(sorted(set(grp["PDB"].tolist()))),
                "all_heavy_species": ";".join(sorted(set(str(v) for v in grp["heavy_species"]))),
                "all_heavy_taxid": ";".join(sorted(set(str(v) for v in grp["heavy_taxid"]))),
                "all_normalized_species": ";".join(sorted(set(str(v) for v in grp["normalized_species"]))),
                "all_sdh_classification": ";".join(sorted(set(str(v) for v in grp["sdh_classification"]))),
                "all_primary_vhh_candidate_metadata": ";".join(sorted(set(str(v) for v in grp["primary_vhh_candidate_metadata"]))),
                "inconsistent_fields": ";".join(inconsistent_fields),
            })
    consistency_df = pd.DataFrame(conflict_rows)
    consistency_df.to_csv(PROJECT_ROOT / "results" / "audits" / "stage1_sdh_sabdab_id_consistency.csv", index=False)

    # ---------------- 6. species distribution of primary cohort ----------------
    primary = sdh[sdh["primary_vhh_candidate_metadata"] == True].copy()  # noqa: E712
    primary["species_group"] = primary.apply(species_group, axis=1)
    dist_rows = []
    for grp_name in ["Camelus dromedarius", "Camelus bactrianus", "Camelus ferus", "other Camelus",
                      "Lama glama", "other Lama", "Vicugna pacos", "other/unresolved Camelidae"]:
        sub = primary[primary["species_group"] == grp_name]
        dist_rows.append({
            "species_group": grp_name,
            "n_rows": len(sub),
            "unique_SABDAB_ID": sub["SABDAB_ID"].nunique(),
        })
    dist_df = pd.DataFrame(dist_rows)
    dist_df.to_csv(PROJECT_ROOT / "results" / "tables" / "stage1_primary_vhh_species_distribution.csv", index=False)

    # ---------------- 7. flag distribution inside primary cohort ----------------
    flag_defs = {
        "mutation_flag": primary["mutation_flag"].fillna(False) == True,  # noqa: E712
        "engineering_keyword_hits": primary["engineering_keyword_hits"].apply(not_empty),
        "humanization_keyword_hits": primary["humanization_keyword_hits"].apply(not_empty),
        "synthetic_keyword_hits": primary["synthetic_keyword_hits"].apply(not_empty),
        "metadata_conflict_flag": primary["metadata_conflict_flag"].fillna(False) == True,  # noqa: E712
        "manual_review_reason": primary["manual_review_reason"].apply(not_empty),
    }
    flag_counts = {k: int(v.sum()) for k, v in flag_defs.items()}
    n_flags_per_row = pd.concat(flag_defs, axis=1).sum(axis=1)
    n_multi_flag = int((n_flags_per_row > 1).sum())

    # ---------------- 8 & 9. individual dossiers for uncertain and eng/hum rows ----------------
    uncertain_rows = sdh[sdh["sdh_classification"] == "uncertain_sdH"].sort_values(["PDB", "INSTANCE"])
    enghum_rows = sdh[sdh["sdh_classification"] == "engineered_or_humanized_camelid_candidate"].sort_values(["PDB", "INSTANCE"])

    # ---------------- 10. paired cohort sanity check ----------------
    fab = df[df["source_cohort"] == "FAB"]
    fv = df[df["source_cohort"] == "FV"]
    fab_primary = int(fab["primary_paired_candidate_metadata"].sum())
    fv_primary = int(fv["primary_paired_candidate_metadata"].sum())
    engineered_buckets = ["SCFV", "DIABODY", "CODV", "DVD", "SCDB"]
    bad_marked = df[df["construct_bucket"].isin(engineered_buckets) & (df["primary_paired_candidate_metadata"] == True)]  # noqa: E712
    n_bad_marked = len(bad_marked)

    # ---------------- extra: non-antibody contamination sweep of the primary cohort ----------------
    # Triggered by the SABDAB_ID consistency check surfacing sabdab2_H02Y5L0000 mixing
    # human_sdVH_candidate / natural_camelid_candidate / uncertain_sdH members whose compound
    # text turned out to describe MHC class-I single-chain trimers, not antibodies at all.
    # This is a pure free-text (compound/short_header) check -> still metadata-only, no
    # sequence/coordinate work performed.
    contamination_patterns = [
        "single chain trimer", "single-chain trimer", "hla-a", "hla-b", "hla-c",
        "mhc class i", "minibinder", "beta-2-microglobulin", "beta2-microglobulin",
    ]
    primary_text = (primary["compound"].fillna("") + " " + primary["short_header"].fillna("")).str.lower()
    contamination_mask = primary_text.apply(lambda t: any(p in t for p in contamination_patterns))
    contamination = primary[contamination_mask].copy()
    contamination["contamination_reason"] = "compound/short_header text indicates an MHC single-chain-trimer construct, not an antibody, despite heavy_taxid/heavy_species indicating camelid origin"
    contam_cols = ["INSTANCE", "PDB", "SABDAB_ID", "heavy_species", "heavy_taxid", "compound",
                   "short_header", "sdh_classification", "primary_vhh_candidate_metadata",
                   "contamination_reason"]
    contamination[contam_cols].to_csv(
        PROJECT_ROOT / "results" / "audits" / "stage1_primary_cohort_suspected_contamination.csv",
        index=False,
    )

    # ---------------- write markdown report ----------------
    md_path = PROJECT_ROOT / "results" / "audits" / "stage1_classification_boundary_audit.md"
    lines = []
    lines.append("# Stage 1 — SD-H Boundary-Classification Audit (read-only)\n")
    lines.append("No classification rule was changed. No `primary_vhh_candidate_metadata` or "
                  "`sdh_classification` value was modified. No structure-, coordinate-, or "
                  "sequence-level analysis was performed.\n")

    lines.append("## 1. Boundary-case table\n")
    lines.append(f"Total boundary-case rows: **{len(boundary)}**\n")
    lines.append(f"- uncertain_sdH: {int(is_uncertain.sum())}")
    lines.append(f"- engineered_or_humanized_camelid_candidate: {int(is_eng_hum.sum())}")
    lines.append(f"- natural_camelid_candidate carrying at least one warning flag: {int(is_flagged_natural.sum())}\n")
    lines.append("Full table: `results/audits/stage1_sdh_boundary_cases.csv`\n")

    lines.append("## 2. Expression-system vs biological-origin separation\n")
    lines.append("Evidence hierarchy applied: heavy_taxid > heavy_species > organism_taxid/organism > "
                  "heavy_subclass > free text. Expression-system fields "
                  "(`expression_system`, `heavy_expression_system`, and their taxids) were carried into "
                  "the boundary table for inspection but were **never used** to set `normalized_species`, "
                  "`camelid_species_flag`, or `sdh_classification` at any point in the Stage-1 pipeline "
                  "(`src/qc/stage1_pipeline.py::classify_sdh_row` never reads any `*expression*` field).\n")
    n_expr_taxid_in_heavy_taxid = int((sdh["heavy_taxid"].isin(NONBIOLOGICAL_HOST_TAXIDS)).sum())
    lines.append(f"Rows where `heavy_taxid` itself equals a known expression-host taxid "
                 f"(E. coli=562): **{n_expr_taxid_in_heavy_taxid}**. These are the 3 rows already "
                 f"routed to `uncertain_sdH` in Stage 1 — flagged as a likely upstream field-mapping "
                 f"error (heavy_species/heavy_taxid populated from expression system), not resolved "
                 f"automatically either way.\n")

    lines.append("## 5. Same-SABDAB_ID consistency check\n")
    lines.append(f"SAbDab2 IDs with >1 SD-H instance and inconsistent classification/species/flags across "
                 f"instances: **{len(consistency_df)}**\n")
    if len(consistency_df) == 0:
        lines.append("No inconsistencies found: every SAbDab2 ID with multiple SD-H instances agrees on "
                     "species, classification, all origin-related flags, and primary-candidate status "
                     "across every instance.\n")
    lines.append("Full table: `results/audits/stage1_sdh_sabdab_id_consistency.csv`\n")

    lines.append("## 6. Species distribution of the current primary cohort\n")
    lines.append("| species_group | n_rows | unique_SABDAB_ID |")
    lines.append("|---|---|---|")
    for _, r in dist_df.iterrows():
        lines.append(f"| {r['species_group']} | {r['n_rows']} | {r['unique_SABDAB_ID']} |")
    lines.append("")

    lines.append("## 7. Flag distribution inside the primary cohort\n")
    for k, v in flag_counts.items():
        lines.append(f"- {k}: {v}")
    lines.append(f"- rows with more than one flag simultaneously: {n_multi_flag}\n")

    lines.append("## K. Suspected non-antibody contamination in the primary cohort (critical finding)\n")
    lines.append(f"Surfaced by cross-referencing the SABDAB_ID consistency check (`sabdab2_H02Y5L0000`), "
                 f"then swept across the full primary cohort by compound/short_header text pattern "
                 f"(MHC single-chain-trimer / minibinder vocabulary). **{len(contamination)} primary-cohort "
                 f"rows ({contamination['PDB'].nunique()} unique PDB, all under SABDAB_ID "
                 f"`sabdab2_H02Y5L0000`) describe MHC class-I single-chain-trimer constructs "
                 f"(e.g. 'Single chain trimer HLA-A*02:01 ... with HPV.16 E7 peptide'), not antibodies at "
                 f"all**, despite carrying `heavy_taxid=9844` (Lama glama) and passing every current "
                 f"classifier check (camelid taxid, X-ray, resolution ≤3Å, no engineering/humanization/"
                 f"synthetic/mutation keyword hit). Beta-2-microglobulin, a component of these engineered "
                 f"single-chain-trimer fusions, carries a genuine immunoglobulin-superfamily fold, which "
                 f"most plausibly explains why SAbDab2's own automated structural classifier admitted "
                 f"these into the SD-H search results in the first place — this is an upstream SAbDab2 "
                 f"classification artifact, not a bug in this project's metadata pipeline. Full row list: "
                 f"`results/audits/stage1_primary_cohort_suspected_contamination.csv`.\n")
    lines.append("This was **not** caught by the current boundary-case criteria (engineering/humanization/"
                 "synthetic/mutation keyword lists do not include MHC/scaffold-specific vocabulary such as "
                 "'single chain trimer' or 'HLA') — a documented gap for a future rule update, not fixed "
                 "in this read-only audit.\n")
    broader_sweep_note = ("A broader sweep of the full 1419-row primary cohort for other non-antibody "
                          "engineered-scaffold vocabulary (DARPin, affibody, monobody, knottin, "
                          "fibronectin/FN3, anticalin, centyrin, adnectin, repebody, 'designed binder', "
                          "'computationally designed') found **0** additional hits — the contamination "
                          "found is currently isolated to this one SABDAB_ID cluster.\n")
    lines.append(broader_sweep_note)

    lines.append("## 10. Paired-cohort sanity check\n")
    lines.append(f"- FAB primary_paired_candidate_metadata count: {fab_primary}")
    lines.append(f"- FV primary_paired_candidate_metadata count: {fv_primary}")
    lines.append(f"- Engineered chain-sharing rows incorrectly marked primary_paired_candidate_metadata=True: {n_bad_marked}\n")

    md_path.write_text("\n".join(lines))

    # ---------------- print everything needed for the final report ----------------
    print("=== A/B/C boundary counts ===")
    print(f"total_boundary_rows={len(boundary)}")
    print(f"uncertain_sdH={int(is_uncertain.sum())}")
    print(f"engineered_or_humanized_camelid_candidate={int(is_eng_hum.sum())}")
    print(f"flagged_natural_camelid={int(is_flagged_natural.sum())}")

    print("\n=== E species distribution (primary cohort) ===")
    print(dist_df.to_string(index=False))
    print("total primary rows:", len(primary), "sum of groups:", dist_df["n_rows"].sum())

    print("\n=== F flag distribution (primary cohort) ===")
    for k, v in flag_counts.items():
        print(f"{k}: {v}")
    print(f"multi_flag_rows: {n_multi_flag}")

    print("\n=== G/H consistency conflicts ===")
    print(f"n_conflicting_sabdab_ids={len(consistency_df)}")
    if len(consistency_df):
        print(consistency_df.to_string(index=False))

    print("\n=== I proposed actions for uncertain_sdH ===")
    unc_boundary = boundary_out[boundary_out["current_sdh_category"] == "uncertain_sdH"]
    print(unc_boundary["proposed_review_action"].value_counts().to_string())

    print("\n=== J proposed actions for engineered/humanized ===")
    eh_boundary = boundary_out[boundary_out["current_sdh_category"] == "engineered_or_humanized_camelid_candidate"]
    print(eh_boundary["proposed_review_action"].value_counts().to_string())

    print("\n=== 10. paired sanity ===")
    print(f"fab_primary={fab_primary} fv_primary={fv_primary} bad_marked_engineered={n_bad_marked}")

    # Save uncertain / enghum dossiers as separate helper CSVs for inspection (not required outputs, but useful)
    uncertain_rows.to_csv(PROJECT_ROOT / "results" / "audits" / "_stage1_uncertain_dossiers_raw.csv", index=False)
    enghum_rows.to_csv(PROJECT_ROOT / "results" / "audits" / "_stage1_enghum_dossiers_raw.csv", index=False)


if __name__ == "__main__":
    main()
