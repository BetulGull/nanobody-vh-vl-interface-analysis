"""Build the authoritative full extraction manifest.

For every member of every source .tgz archive, records exactly which
extracted file on disk holds its content — including members that landed
under a __casecollisionN filename because of APFS case-insensitivity.

Output columns:
  original_archive, archive_sha256, original_member_path, original_filename,
  extracted_path, content_sha256, case_collision_flag
"""
import csv
import hashlib
import os
import tarfile
from pathlib import Path

PROJECT_ROOT = Path("/Users/betulgul/Desktop/staj")

# (archive filename, dest dir relative to PROJECT_ROOT)
ARCHIVE_DEST_MAP = [
    ("sabdab_search_structures_20260901T143831.tgz", "data/raw_extracted/sdh"),
    ("sabdab_search_structures_20260901T151240.tgz", "data/raw_extracted/fv"),
    ("sabdab_search_structures_20260901T150654.tgz", "data/raw_extracted/fab"),
    ("sabdab_search_structures_20260901T151141.tgz", "data/raw_extracted/fab"),  # duplicate of 150654
    ("sabdab_search_structures_20260901T145541.tgz", "data/raw_extracted/unknown/archive_145541"),
    ("sabdab_search_structures_20260901T150035.tgz", "data/raw_extracted/unknown/archive_150035"),
]


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    return sha256_bytes(path.read_bytes())


def build_disk_index(dest_root: Path):
    """Map content_sha256 -> list of relative paths currently on disk under dest_root."""
    index = {}
    for p in dest_root.rglob("*.cif"):
        rel = p.relative_to(dest_root)
        h = sha256_file(p)
        index.setdefault(h, []).append(str(rel))
    return index


def main():
    out_path = PROJECT_ROOT / "data" / "raw_extracted" / "extraction_manifest.tsv"
    rows = []

    disk_index_cache = {}

    for archive_name, dest_rel in ARCHIVE_DEST_MAP:
        archive_path = PROJECT_ROOT / archive_name
        dest_root = PROJECT_ROOT / dest_rel
        archive_sha256 = sha256_file(archive_path)

        if dest_rel not in disk_index_cache:
            disk_index_cache[dest_rel] = build_disk_index(dest_root)
        disk_index = disk_index_cache[dest_rel]

        with tarfile.open(archive_path, "r:gz") as tf:
            for member in tf.getmembers():
                if not member.isfile():
                    continue
                data = tf.extractfile(member).read()
                content_hash = sha256_bytes(data)
                original_filename = Path(member.name).name

                direct_dest = dest_root / member.name
                candidates = disk_index.get(content_hash, [])

                if str(Path(member.name)) in candidates:
                    extracted_path = str(Path(member.name))
                    case_collision_flag = False
                elif candidates:
                    # content exists on disk, but not at the exact member.name path ->
                    # either a plain collision-remap, or (for the duplicate archive)
                    # content shared with the primary archive's extraction.
                    same_dir_candidates = [
                        c for c in candidates
                        if str(Path(c).parent) == str(Path(member.name).parent)
                    ]
                    chosen = same_dir_candidates[0] if same_dir_candidates else candidates[0]
                    extracted_path = chosen
                    case_collision_flag = (chosen != str(Path(member.name)))
                else:
                    extracted_path = "MISSING_ON_DISK"
                    case_collision_flag = False

                rows.append({
                    "original_archive": archive_name,
                    "archive_sha256": archive_sha256,
                    "original_member_path": member.name,
                    "original_filename": original_filename,
                    "extracted_path": f"{dest_rel}/{extracted_path}" if extracted_path != "MISSING_ON_DISK" else "MISSING_ON_DISK",
                    "content_sha256": content_hash,
                    "case_collision_flag": case_collision_flag,
                })

    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "original_archive", "archive_sha256", "original_member_path",
            "original_filename", "extracted_path", "content_sha256",
            "case_collision_flag",
        ], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    n_missing = sum(1 for r in rows if r["extracted_path"] == "MISSING_ON_DISK")
    n_collision = sum(1 for r in rows if r["case_collision_flag"])
    print(f"Wrote {len(rows)} rows to {out_path}")
    print(f"  case_collision_flag=True: {n_collision}")
    print(f"  MISSING_ON_DISK: {n_missing}")


if __name__ == "__main__":
    main()
