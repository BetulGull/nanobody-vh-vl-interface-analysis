"""Stage 1 — metadata ingestion and cohort classification.

Metadata-level processing only. No SASA/geometry/interaction/electrostatics
calculation happens here.
"""
import csv
import re
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path("/Users/betulgul/Desktop/staj")

SUMMARY_CSVS = {
    "SD-H": "sabdab_search_summary_20260901T161115.csv",
    "FAB": "sabdab_search_summary_20260901T161202.csv",
    "FV": "sabdab_search_summary_20260901T161254.csv",
}

CANONICAL_FIELDS = [
    "INSTANCE", "PDB", "SABDAB_ID", "HEAVY_ID", "LIGHT_ID", "Hchain", "Lchain",
    "model", "antigen_chain", "antigen_type", "antigen_name", "organism",
    "organism_taxid", "heavy_species", "light_species", "heavy_taxid",
    "light_taxid", "method", "resolution", "r_free", "r_factor", "type",
    "chainsharing_construct", "chainsharing_construct_partners",
    "non_chainsharing_construct", "non_chainsharing_construct_partners",
    "heavy_subclass", "light_subclass", "light_ctype", "compound",
    "short_header", "authors", "date",
]

# ---------------------------------------------------------------------------
# Species / taxid reference tables (structured evidence, priority 1-2)
# ---------------------------------------------------------------------------
CAMELID_TAXIDS = {
    "9838": "Camelus dromedarius",
    "9837": "Camelus bactrianus",
    "419612": "Camelus ferus",
    "9844": "Lama glama",
    "30538": "Vicugna pacos",
    "9835": "Camelidae (family-level)",
    "9839": "Lama (genus-level)",
    "1579311": "Camelidae mixed library (pooled camelid donors)",
}
HUMAN_TAXID = "9606"
SYNTHETIC_TAXID = "32630"

CAMELID_TEXT_PATTERNS = [
    "camelus dromedarius", "camelus bactrianus", "camelus ferus", "camelus spp",
    "lama glama", "lama spp", "vicugna pacos", "camelidae",
]

ENGINEERING_KEYWORDS = ["engineered", "framework graft", "graft", "minimizer"]
HUMANIZATION_KEYWORDS = ["humanized", "humanised"]
SYNTHETIC_KEYWORDS = ["synthetic", "de novo", "de-novo", "designed"]
MUTATION_KEYWORDS = ["mutant", "mutation", "variant", "affinity-matured", "affinity matured"]
GERMLINE_REVERSION_KEYWORDS = ["germline reversion", "germline-reverted", "germlined"]

CHAINSHARING_MAP = {
    "SCFV": "construct_scfv",
    "DIABODY": "construct_diabody",
    "CODV": "construct_codv",
    "DVD": "construct_dvd",
    "SCDB": "construct_scdb",
}


def na_or(val):
    return None if val in (None, "", "NA") else val


def load_csv(path: Path):
    with open(path, newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def load_manifest_and_matches():
    matches_path = PROJECT_ROOT / "data" / "metadata" / "archive_structure_to_summary_matches.tsv"
    with open(matches_path, newline="") as f:
        matches = list(csv.DictReader(f, delimiter="\t"))
    return matches


def load_archive_sha256_map():
    crosswalk_path = PROJECT_ROOT / "results" / "audits" / "archive_summary_crosswalk.tsv"
    mapping = {}
    with open(crosswalk_path, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row["archive"] != "NOT_PRESENT":
                mapping[row["archive"]] = row["archive_sha256"]
    return mapping


DUPLICATE_ARCHIVE = "sabdab_search_structures_20260901T151141.tgz"


def build_master_raw():
    matches = load_manifest_and_matches()
    sha_map = load_archive_sha256_map()
    # key: (PDB, Hchain, Lchain-or-NA) -> list of match rows (may include duplicate archive)
    match_index = defaultdict(list)
    for m in matches:
        key = (m["parsed_PDB"], m["parsed_Hchain"], m["parsed_Lchain"])
        match_index[key].append(m)

    master_rows = []
    for cohort, csv_name in SUMMARY_CSVS.items():
        rows = load_csv(PROJECT_ROOT / csv_name)
        for r in rows:
            key = (r["PDB"], r["Hchain"], r["Lchain"] if r["Lchain"] != "NA" else "NA")
            candidates = match_index.get(key, [])
            archives = sorted(set(c["archive"] for c in candidates))
            is_dup_only = DUPLICATE_ARCHIVE in archives
            non_dup_candidates = [c for c in candidates if c["archive"] != DUPLICATE_ARCHIVE]
            primary = non_dup_candidates[0] if non_dup_candidates else (candidates[0] if candidates else None)

            row = dict(r)
            row["source_cohort"] = cohort
            row["source_summary_csv"] = csv_name
            row["source_archive"] = ";".join(archives) if archives else ""
            row["source_archive_member"] = primary["original_member_path"] if primary else ""
            row["extracted_structure_path"] = primary["extracted_path"] if primary else ""
            row["archive_sha256"] = ";".join(sha_map.get(a, "") for a in archives) if archives else ""
            row["content_duplicate_flag"] = DUPLICATE_ARCHIVE in archives
            row["n_matched_archives"] = len(archives)
            master_rows.append(row)

    return pd.DataFrame(master_rows)


# ---------------------------------------------------------------------------
# QC flags
# ---------------------------------------------------------------------------
def normalize_method(method_raw):
    if not method_raw or method_raw == "NA":
        return None
    m = method_raw.strip().upper()
    if "XRAY" in m.replace("-", "").replace(" ", "") or "X-RAY" in m:
        return "XRAY"
    return m


def coerce_resolution(res_raw):
    if not res_raw or res_raw == "NA":
        return None
    try:
        return float(res_raw)
    except ValueError:
        return None


def add_qc_flags(df: pd.DataFrame) -> pd.DataFrame:
    df["method_normalized"] = df["method"].apply(normalize_method)
    df["qc_method_xray"] = df["method_normalized"] == "XRAY"
    df["qc_resolution_numeric"] = df["resolution"].apply(coerce_resolution)
    df["qc_resolution_le_3"] = df["qc_resolution_numeric"].apply(
        lambda v: (v is not None) and (v <= 3.0)
    )
    df["qc_contains_protein_annotation"] = df["antigen_type"].apply(
        lambda v: isinstance(v, str) and "PROTEIN" in v.upper()
    )
    df["qc_expected_source_type"] = df["type"] == df["source_cohort"]
    df["qc_has_heavy_metadata"] = df["Hchain"].apply(lambda v: na_or(v) is not None)
    df["qc_has_light_metadata"] = df["Lchain"].apply(lambda v: na_or(v) is not None)
    df["qc_basic_metadata_complete"] = (
        df["qc_method_xray"]
        & df["qc_resolution_numeric"].notna()
        & df["qc_has_heavy_metadata"]
    )
    return df


# ---------------------------------------------------------------------------
# SD-H classification
# ---------------------------------------------------------------------------
def keyword_hits(text_fields: dict, keywords):
    hits = []
    fields_hit = []
    for field, text in text_fields.items():
        if not text or text == "NA":
            continue
        low = text.lower()
        for kw in keywords:
            if kw in low:
                hits.append(kw)
                fields_hit.append(field)
    return hits, fields_hit


def classify_sdh_row(row) -> dict:
    heavy_taxid = na_or(row.get("heavy_taxid"))
    heavy_species = na_or(row.get("heavy_species"))
    organism_taxid = na_or(row.get("organism_taxid"))
    organism = na_or(row.get("organism"))
    heavy_subclass = na_or(row.get("heavy_subclass"))
    compound = row.get("compound") or ""
    short_header = row.get("short_header") or ""

    text_fields = {
        "compound": compound,
        "short_header": short_header,
        "heavy_subclass": heavy_subclass or "",
    }

    eng_hits, eng_fields = keyword_hits(text_fields, ENGINEERING_KEYWORDS)
    hum_hits, hum_fields = keyword_hits(text_fields, HUMANIZATION_KEYWORDS)
    syn_hits, syn_fields = keyword_hits(text_fields, SYNTHETIC_KEYWORDS)
    mut_hits, mut_fields = keyword_hits(text_fields, MUTATION_KEYWORDS)
    germ_hits, germ_fields = keyword_hits(text_fields, GERMLINE_REVERSION_KEYWORDS)

    # --- structured species evidence, priority order ---
    species_evidence_field = None
    species_evidence_value = None
    normalized_species = None
    camelid_species_flag = False
    is_human = False
    is_synthetic_taxid = False

    if heavy_taxid is not None:
        species_evidence_field = "heavy_taxid"
        species_evidence_value = heavy_taxid
        if heavy_taxid in CAMELID_TAXIDS:
            camelid_species_flag = True
            normalized_species = CAMELID_TAXIDS[heavy_taxid]
        elif heavy_taxid == HUMAN_TAXID:
            is_human = True
            normalized_species = "Homo sapiens"
        elif heavy_taxid == SYNTHETIC_TAXID:
            is_synthetic_taxid = True
            normalized_species = "synthetic construct"
        else:
            normalized_species = heavy_species or f"taxid:{heavy_taxid}"
    elif heavy_species is not None:
        species_evidence_field = "heavy_species"
        species_evidence_value = heavy_species
        low = heavy_species.lower()
        if any(p in low for p in CAMELID_TEXT_PATTERNS):
            camelid_species_flag = True
            normalized_species = heavy_species
        elif low == "homo sapiens":
            is_human = True
            normalized_species = "Homo sapiens"
        elif low == "synthetic construct":
            is_synthetic_taxid = True
            normalized_species = "synthetic construct"
        else:
            normalized_species = heavy_species
    elif organism_taxid is not None or organism is not None:
        species_evidence_field = "organism_taxid" if organism_taxid else "organism"
        species_evidence_value = organism_taxid or organism
        val = organism_taxid or organism
        if organism_taxid in CAMELID_TAXIDS or (organism and any(p in organism.lower() for p in CAMELID_TEXT_PATTERNS)):
            camelid_species_flag = True
        normalized_species = organism or f"taxid:{organism_taxid}"
    else:
        species_evidence_field = "none"
        species_evidence_value = None
        normalized_species = None

    # --- metadata conflict detection ---
    metadata_conflict_flag = False
    conflict_reasons = []

    # (a) structured evidence present but implausible as a real antibody source
    #     organism (E. coli, other clear non-mammalian expression hosts) while
    #     descriptive text explicitly uses camelid single-domain terminology.
    NONBIOLOGICAL_HOST_TAXIDS = {"562"}  # E. coli
    text_all = f"{compound} {short_header}".lower()
    text_says_nanobody_vhh = ("nanobody" in text_all) or ("vhh" in text_all) or (re.search(r"\bnb\d", text_all) is not None)

    if heavy_taxid in NONBIOLOGICAL_HOST_TAXIDS:
        metadata_conflict_flag = True
        conflict_reasons.append(
            f"heavy_taxid={heavy_taxid} indicates an expression host (E. coli), not a biological antibody-source "
            f"organism; likely heavy_species/heavy_taxid populated from expression_system by upstream error"
            + (" (descriptive text uses nanobody/VHH terminology)" if text_says_nanobody_vhh else "")
        )

    if species_evidence_field == "none":
        metadata_conflict_flag = True
        conflict_reasons.append("no structured species evidence available (heavy_taxid, heavy_species, organism_taxid, organism all missing)")

    if is_human and text_says_nanobody_vhh:
        # Not a conflict per the explicit rule (structured Homo sapiens wins), but recorded for transparency.
        conflict_reasons.append(
            "INFO: heavy_taxid=9606 (Homo sapiens) but descriptive text uses nanobody/VHH terminology; "
            "structured species evidence takes priority per project rule, classified as human_sdVH_candidate"
        )

    # --- primary classification bucket ---
    if metadata_conflict_flag and species_evidence_field in ("none",) :
        sdh_class = "uncertain_sdH"
    elif heavy_taxid in NONBIOLOGICAL_HOST_TAXIDS:
        sdh_class = "uncertain_sdH"
    elif is_human:
        sdh_class = "human_sdVH_candidate"
    elif is_synthetic_taxid:
        sdh_class = "synthetic_or_de_novo_candidate"
    elif camelid_species_flag:
        if hum_hits:
            sdh_class = "engineered_or_humanized_camelid_candidate"
        elif syn_hits:
            # camelid-origin text also flagged synthetic/designed keywords -> conflict between
            # structured camelid provenance and a stronger synthetic/de-novo textual claim
            metadata_conflict_flag = True
            conflict_reasons.append(
                "structured species evidence indicates camelid origin, but descriptive text contains a "
                "synthetic/de-novo keyword; sent to manual review rather than resolved automatically"
            )
            sdh_class = "uncertain_sdH"
        else:
            sdh_class = "natural_camelid_candidate"
    else:
        sdh_class = "other_noncamelid_sdH"

    manual_review_reason = "; ".join(conflict_reasons) if conflict_reasons else None

    return {
        "species_evidence_field": species_evidence_field,
        "species_evidence_value": species_evidence_value,
        "normalized_species": normalized_species,
        "camelid_species_flag": camelid_species_flag,
        "engineering_keyword_hits": ";".join(eng_hits) if eng_hits else None,
        "engineering_keyword_fields": ";".join(sorted(set(eng_fields))) if eng_fields else None,
        "humanization_keyword_hits": ";".join(hum_hits) if hum_hits else None,
        "humanization_keyword_fields": ";".join(sorted(set(hum_fields))) if hum_fields else None,
        "synthetic_keyword_hits": ";".join(syn_hits) if syn_hits else None,
        "synthetic_keyword_fields": ";".join(sorted(set(syn_fields))) if syn_fields else None,
        "mutation_keyword_hits": ";".join(mut_hits) if mut_hits else None,
        "mutation_keyword_fields": ";".join(sorted(set(mut_fields))) if mut_fields else None,
        "germline_reversion_keyword_hits": ";".join(germ_hits) if germ_hits else None,
        "engineering_flag": bool(eng_hits),
        "humanization_flag": bool(hum_hits),
        "synthetic_flag": bool(syn_hits) or is_synthetic_taxid,
        "mutation_flag": bool(mut_hits),
        "metadata_conflict_flag": metadata_conflict_flag,
        "metadata_conflict_reason": manual_review_reason,
        "manual_review_reason": manual_review_reason if (metadata_conflict_flag or sdh_class == "uncertain_sdH") else None,
        "sdh_classification": sdh_class,
    }


def apply_sdh_classification(df: pd.DataFrame) -> pd.DataFrame:
    is_sdh = df["source_cohort"] == "SD-H"
    classifications = df[is_sdh].apply(lambda r: pd.Series(classify_sdh_row(r)), axis=1)
    for col in classifications.columns:
        if col not in df.columns:
            df[col] = None
        df.loc[is_sdh, col] = classifications[col]

    df["primary_vhh_candidate_metadata"] = False
    mask = (
        is_sdh
        & (df["sdh_classification"] == "natural_camelid_candidate")
        & (df["qc_method_xray"])
        & (df["qc_resolution_le_3"])
    )
    df.loc[mask, "primary_vhh_candidate_metadata"] = True
    return df


# ---------------------------------------------------------------------------
# Paired (FAB/FV) construct classification
# ---------------------------------------------------------------------------
def classify_construct(row) -> dict:
    cs = na_or(row.get("chainsharing_construct"))
    flags = {f"construct_{k}": False for k in
             ["conventional", "scfv", "diabody", "codv", "dvd", "scdb", "other_engineered", "unknown"]}

    if cs is None:
        flags["construct_conventional"] = True
        bucket = "conventional"
    elif cs.upper() in CHAINSHARING_MAP:
        flags[CHAINSHARING_MAP[cs.upper()]] = True
        bucket = cs.upper()
    else:
        flags["construct_other_engineered"] = True
        bucket = f"other:{cs}"

    return {**flags, "construct_bucket": bucket}


def apply_paired_classification(df: pd.DataFrame) -> pd.DataFrame:
    is_paired = df["source_cohort"].isin(["FAB", "FV"])
    df["paired_metadata_candidate"] = (
        is_paired
        & df["qc_has_heavy_metadata"]
        & df["qc_has_light_metadata"]
    )

    construct_cols = df[is_paired].apply(lambda r: pd.Series(classify_construct(r)), axis=1)
    for col in construct_cols.columns:
        if col not in df.columns:
            df[col] = None
        df.loc[is_paired, col] = construct_cols[col]

    is_engineered_chainsharing = df["construct_bucket"].isin(["SCFV", "DIABODY", "CODV", "DVD", "SCDB"]) | df["construct_bucket"].astype(str).str.startswith("other:")

    df["primary_paired_candidate_metadata"] = (
        is_paired
        & df["qc_has_heavy_metadata"]
        & df["qc_has_light_metadata"]
        & df["qc_method_xray"]
        & df["qc_resolution_le_3"]
        & (~is_engineered_chainsharing.fillna(False))
    )
    return df


# ---------------------------------------------------------------------------
# Exclusion reasons (applies across cohorts, additive, non-destructive)
# ---------------------------------------------------------------------------
def compute_exclusion_reasons(row) -> str:
    reasons = []
    if not row.get("qc_method_xray"):
        reasons.append("FAIL_NOT_XRAY")
    if row.get("qc_resolution_numeric") is None:
        reasons.append("FAIL_RESOLUTION_MISSING")
    elif not row.get("qc_resolution_le_3"):
        reasons.append("FAIL_RESOLUTION_GT_3")
    if not row.get("qc_has_heavy_metadata"):
        reasons.append("FAIL_MISSING_HEAVY")

    cohort = row.get("source_cohort")
    if cohort in ("FAB", "FV"):
        if not row.get("qc_has_light_metadata"):
            reasons.append("FAIL_MISSING_LIGHT")
        bucket = row.get("construct_bucket")
        if bucket not in (None, "conventional") :
            reasons.append("FAIL_ENGINEERED_CHAINSHARING")
    if cohort == "SD-H":
        cls = row.get("sdh_classification")
        if cls == "other_noncamelid_sdH":
            reasons.append("FAIL_NONCAMELID_SDH")
        elif cls == "synthetic_or_de_novo_candidate":
            reasons.append("FAIL_SYNTHETIC_SDH")
        elif cls == "engineered_or_humanized_camelid_candidate":
            reasons.append("FAIL_HUMANIZED_SDH")
        elif cls == "uncertain_sdH":
            reasons.append("REVIEW_UNCERTAIN_SDH")
        elif cls == "human_sdVH_candidate":
            reasons.append("FAIL_NONCAMELID_SDH")

    return ";".join(reasons) if reasons else "NONE"


def main():
    print("Building master_raw...")
    df = build_master_raw()
    print(f"master_raw rows: {len(df)}")

    df.to_csv(PROJECT_ROOT / "data" / "metadata" / "master_raw.csv", index=False)
    df.to_parquet(PROJECT_ROOT / "data" / "metadata" / "master_raw.parquet", index=False)

    print("Adding QC flags...")
    df = add_qc_flags(df)

    print("Classifying SD-H rows...")
    df = apply_sdh_classification(df)

    print("Classifying paired FAB/FV rows...")
    df = apply_paired_classification(df)

    print("Computing exclusion reasons...")
    df["exclusion_reasons"] = df.apply(compute_exclusion_reasons, axis=1)

    df.to_csv(PROJECT_ROOT / "data" / "metadata" / "master_stage1.csv", index=False)
    df.to_parquet(PROJECT_ROOT / "data" / "metadata" / "master_stage1.parquet", index=False)
    print(f"Wrote master_stage1 with {len(df)} rows, {len(df.columns)} columns")

    return df


if __name__ == "__main__":
    main()
