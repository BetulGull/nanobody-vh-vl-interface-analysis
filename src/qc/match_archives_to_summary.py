"""Match every downloaded (extracted) structure to official SAbDab2 summary
CSV rows, using original archive-member identity (via the extraction
manifest) rather than any renamed/collision-disambiguated filename.

Archive -> summary is the primary direction: every downloaded structure
should ideally find a match. Summary rows with no downloaded structure are
informational only (summary_only_rows), since the summaries are known
supersets.
"""
import csv
import re
from collections import defaultdict
from pathlib import Path

PROJECT_ROOT = Path("/Users/betulgul/Desktop/staj")

SUMMARY_CSVS = [
    "sabdab_search_summary_20260901T161115.csv",
    "sabdab_search_summary_20260901T161202.csv",
    "sabdab_search_summary_20260901T161254.csv",
]

ARCHIVE_LABELS = {
    "sabdab_search_structures_20260901T143831.tgz": "archive_143831",
    "sabdab_search_structures_20260901T151240.tgz": "archive_151240",
    "sabdab_search_structures_20260901T150654.tgz": "archive_150654",
    "sabdab_search_structures_20260901T151141.tgz": "archive_151141_dup",
    "sabdab_search_structures_20260901T145541.tgz": "archive_145541",
    "sabdab_search_structures_20260901T150035.tgz": "archive_150035",
}

# member basename pattern: pdb_XXXXX_<chain1>_<chain2-or-+>_abag.cif
MEMBER_RE = re.compile(r"^(pdb_[0-9a-z]+)_(.+?)_(.+?)_abag\.cif$")


def parse_member_name(name: str):
    """Return (pdb, chain1, chain2_or_None) from an original archive member basename."""
    base = Path(name).name
    m = MEMBER_RE.match(base)
    if not m:
        return None
    pdb, chain1, chain2 = m.groups()
    if chain2 == "+":
        chain2 = None
    return pdb, chain1, chain2


def load_summary_csv(path: Path):
    with open(path, newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        rows = list(reader)
    return rows


def build_index(rows, csv_name):
    """Key: (PDB, Hchain, Lchain-or-None) -> list of (row, csv_name)."""
    index = defaultdict(list)
    for r in rows:
        pdb = r["PDB"]
        hchain = r["Hchain"]
        lchain = r["Lchain"] if r["Lchain"] not in ("NA", "") else None
        key = (pdb, hchain, lchain)
        index[key].append((r, csv_name))
    return index


def main():
    manifest_path = PROJECT_ROOT / "data" / "raw_extracted" / "extraction_manifest.tsv"
    with open(manifest_path, newline="") as f:
        manifest_rows = list(csv.DictReader(f, delimiter="\t"))

    summaries = {}
    combined_index = defaultdict(list)
    for csv_name in SUMMARY_CSVS:
        path = PROJECT_ROOT / csv_name
        rows = load_summary_csv(path)
        summaries[csv_name] = rows
        idx = build_index(rows, csv_name)
        for k, v in idx.items():
            combined_index[k].extend(v)

    # dedupe manifest to one row per unique original (archive, member_path):
    # the duplicate archive (151141) shares content with 150654, but is still
    # a separate archive-level provenance record, so keep it separate here.
    matched_rows = []
    unmatched_rows = []

    per_archive_stats = defaultdict(lambda: {"total": 0, "matched": 0, "matched_types": defaultdict(int)})

    for row in manifest_rows:
        archive = row["original_archive"]
        member_path = row["original_member_path"]
        parsed = parse_member_name(member_path)

        per_archive_stats[archive]["total"] += 1

        if parsed is None:
            unmatched_rows.append({
                "archive": archive,
                "original_member_path": member_path,
                "extracted_path": row["extracted_path"],
                "parsed_PDB": "",
                "parsed_Hchain": "",
                "parsed_Lchain": "",
                "possible_reason": "filename_did_not_match_expected_pattern",
            })
            continue

        pdb, chain1, chain2 = parsed
        key = (pdb, chain1, chain2)
        candidates = combined_index.get(key, [])

        if not candidates:
            # try case-insensitive chain match as a diagnostic fallback (not accepted as a real match)
            reason = "no_exact_PDB_Hchain_Lchain_match_in_any_summary"
            fallback_hint = None
            for (kpdb, khc, klc), vals in combined_index.items():
                if kpdb == pdb and khc.lower() == (chain1 or "").lower() and (klc or "").lower() == (chain2 or "").lower():
                    fallback_hint = "case_insensitive_chain_match_exists"
                    break
            if fallback_hint:
                reason = "case_sensitive_chain_ID_mismatch"
            unmatched_rows.append({
                "archive": archive,
                "original_member_path": member_path,
                "extracted_path": row["extracted_path"],
                "parsed_PDB": pdb,
                "parsed_Hchain": chain1,
                "parsed_Lchain": chain2 or "NA",
                "possible_reason": reason,
            })
            continue

        per_archive_stats[archive]["matched"] += 1
        for (csv_row, csv_name) in candidates:
            per_archive_stats[archive]["matched_types"][csv_row["type"]] += 1

        matched_rows.append({
            "archive": archive,
            "original_member_path": member_path,
            "extracted_path": row["extracted_path"],
            "parsed_PDB": pdb,
            "parsed_Hchain": chain1,
            "parsed_Lchain": chain2 or "NA",
            "matched_csv": candidates[0][1],
            "matched_type": candidates[0][0]["type"],
            "matched_INSTANCE": candidates[0][0]["INSTANCE"],
            "matched_SABDAB_ID": candidates[0][0]["SABDAB_ID"],
        })

    # ---- print summary stats ----
    print("=== Per-archive match stats (against ALL summary CSVs combined) ===")
    for archive, stats in per_archive_stats.items():
        total = stats["total"]
        matched = stats["matched"]
        frac = matched / total if total else 0
        print(f"{ARCHIVE_LABELS.get(archive, archive)}: total={total} matched={matched} "
              f"unmatched={total - matched} match_fraction={frac:.4f} "
              f"matched_types={dict(stats['matched_types'])}")

    # ---- summary-only rows per CSV (informational) ----
    print("\n=== Summary-only rows per CSV (rows with no matched archive structure) ===")
    matched_keys = set()
    for row in matched_rows:
        matched_keys.add((row["parsed_PDB"], row["parsed_Hchain"], row["parsed_Lchain"] if row["parsed_Lchain"] != "NA" else None))

    for csv_name, rows in summaries.items():
        n_only = 0
        for r in rows:
            key = (r["PDB"], r["Hchain"], r["Lchain"] if r["Lchain"] not in ("NA", "") else None)
            if key not in matched_keys:
                n_only += 1
        print(f"{csv_name}: total_rows={len(rows)} summary_only_rows={n_only}")

    # ---- write outputs ----
    out_unmatched = PROJECT_ROOT / "results" / "audits" / "unmatched_downloaded_structures.tsv"
    with open(out_unmatched, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "archive", "original_member_path", "extracted_path",
            "parsed_PDB", "parsed_Hchain", "parsed_Lchain", "possible_reason",
        ], delimiter="\t")
        writer.writeheader()
        writer.writerows(unmatched_rows)
    print(f"\nWrote {len(unmatched_rows)} unmatched rows to {out_unmatched}")

    out_matched = PROJECT_ROOT / "data" / "metadata" / "archive_structure_to_summary_matches.tsv"
    with open(out_matched, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "archive", "original_member_path", "extracted_path",
            "parsed_PDB", "parsed_Hchain", "parsed_Lchain",
            "matched_csv", "matched_type", "matched_INSTANCE", "matched_SABDAB_ID",
        ], delimiter="\t")
        writer.writeheader()
        writer.writerows(matched_rows)
    print(f"Wrote {len(matched_rows)} matched rows to {out_matched}")


if __name__ == "__main__":
    main()
