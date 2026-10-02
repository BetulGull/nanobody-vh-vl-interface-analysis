"""Safe, idempotent extraction of raw SAbDab2 .tgz archives.

Never overwrites an existing extracted file that differs from the archive
version; logs any conflict instead. Rejects archive members with unsafe
paths (absolute paths, path traversal).
"""
import hashlib
import os
import tarfile
from pathlib import Path


def _is_safe_member(member: tarfile.TarInfo, dest_root: Path) -> bool:
    if member.name.startswith("/"):
        return False
    dest_path = (dest_root / member.name).resolve()
    try:
        dest_path.relative_to(dest_root.resolve())
    except ValueError:
        return False
    return True


def _sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _actual_case_on_disk(path: Path):
    """Return the real on-disk filename matching `path.name` case-insensitively,
    or None if no such entry exists. Needed because APFS/HFS+ default to
    case-insensitive (but case-preserving) filesystems, under which archive
    members differing only by chain-letter case (e.g. `_B_` vs `_b_`) collide.
    """
    parent = path.parent
    if not parent.exists():
        return None
    target_lower = path.name.lower()
    for entry in os.listdir(parent):
        if entry.lower() == target_lower:
            return entry
    return None


def safe_extract(archive_path: str, dest_dir: str, log_path: str,
                  case_collision_manifest: str = None) -> dict:
    archive_path = Path(archive_path)
    dest_root = Path(dest_dir)
    dest_root.mkdir(parents=True, exist_ok=True)

    stats = {
        "extracted": 0,
        "skipped_existing_identical": 0,
        "conflicts_logged": 0,
        "rejected_unsafe_paths": 0,
        "case_collisions_remapped": 0,
    }

    manifest_f = None
    if case_collision_manifest:
        new_manifest = not Path(case_collision_manifest).exists()
        manifest_f = open(case_collision_manifest, "a")
        if new_manifest:
            manifest_f.write("archive\ttar_member_name\textracted_relative_path\n")

    with open(log_path, "a") as logf, tarfile.open(archive_path, "r:gz") as tf:
        for member in tf.getmembers():
            if not member.isfile():
                continue
            if not _is_safe_member(member, dest_root):
                stats["rejected_unsafe_paths"] += 1
                logf.write(f"REJECTED_UNSAFE_PATH\t{archive_path.name}\t{member.name}\n")
                continue

            dest_path = dest_root / member.name
            src_f = tf.extractfile(member)
            data = src_f.read()
            new_hash = _sha256_bytes(data)

            if dest_path.exists():
                existing_hash = _sha256_bytes(dest_path.read_bytes())
                if existing_hash == new_hash:
                    stats["skipped_existing_identical"] += 1
                    continue

                actual_name = _actual_case_on_disk(dest_path)
                if actual_name is not None and actual_name != dest_path.name:
                    # Case-insensitive filesystem collision (e.g. chain "B" vs "b"):
                    # this is genuinely different archive content, not a true
                    # duplicate path. Extract under a disambiguated filename and
                    # record the remap so the original tar member name is recoverable.
                    stem, suffix = dest_path.stem, dest_path.suffix
                    counter = 2
                    remapped_path = dest_path.parent / f"{stem}__casecollision{counter}{suffix}"
                    while remapped_path.exists() or _actual_case_on_disk(remapped_path) is not None:
                        counter += 1
                        remapped_path = dest_path.parent / f"{stem}__casecollision{counter}{suffix}"
                    remapped_path.write_bytes(data)
                    stats["case_collisions_remapped"] += 1
                    if manifest_f:
                        rel_path = remapped_path.relative_to(dest_root)
                        manifest_f.write(f"{archive_path.name}\t{member.name}\t{rel_path}\n")
                    continue

                stats["conflicts_logged"] += 1
                logf.write(
                    f"CONFLICT_NOT_OVERWRITTEN\t{archive_path.name}\t{member.name}\t"
                    f"existing_sha256={existing_hash}\tarchive_sha256={new_hash}\n"
                )
                continue

            dest_path.parent.mkdir(parents=True, exist_ok=True)
            dest_path.write_bytes(data)
            stats["extracted"] += 1
            if manifest_f:
                rel_path = dest_path.relative_to(dest_root)
                manifest_f.write(f"{archive_path.name}\t{member.name}\t{rel_path}\n")

    if manifest_f:
        manifest_f.close()

    return stats


if __name__ == "__main__":
    import argparse
    import json

    parser = argparse.ArgumentParser()
    parser.add_argument("archive")
    parser.add_argument("dest")
    parser.add_argument("--log", default="logs/extraction_conflicts.log")
    parser.add_argument("--case_manifest", default=None)
    args = parser.parse_args()

    result = safe_extract(args.archive, args.dest, args.log, args.case_manifest)
    print(json.dumps(result, indent=2))
