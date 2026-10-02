# Stage 0 — Environment Setup & Raw Data Extraction Audit

## Project root

`/Users/betulgul/Desktop/staj` — confirmed as the unique directory named
`staj` under `$HOME` (searched to depth 3), and the only one containing the
downloaded `.tgz` archives.

## Machine

- OS: macOS 15.6.1 (Darwin 24.6.0, build 24G90)
- Architecture: arm64 (Apple Silicon)
- Filesystem: APFS, default case-insensitive mode (see data-integrity note
  below)

## Raw archives

See `raw_archive_inventory.tsv` for full checksums/sizes and
`archive_cohort_mapping.tsv` for cohort assignment with evidence and
confidence levels. Summary:

| Archive (timestamp) | Files | Unique PDB IDs | Assigned cohort | Confidence |
|---|---|---|---|---|
| 143831 | 1597 | 778 | SDH | HIGH |
| 145541 | 373 | 236 | UNKNOWN (Fab-format) | LOW |
| 150035 | 1033 | 604 | UNKNOWN (Fab-format) | LOW |
| 150654 | 1392 | 890 | FAB | MEDIUM-HIGH |
| 151141 | 1392 | 890 | FAB (duplicate of 150654) | HIGH (dup determination) |
| 151240 | 536 | 290 | FV | MEDIUM-HIGH |

**Key finding:** no summary CSV was bundled in any archive (100% `.cif` by
extension), so cohort identification was done entirely from mmCIF content
(entity_poly sequence lengths, VH/VL/CH1/CL/Fc signatures), cross-checked
against RCSB's full deposited-entry entity data for representative PDB IDs.
FAB vs FAB+FC could not be reliably separated this way: none of the three
Fab-format archives shows a dominant (>2%) Fc-containing fraction, and their
PDB-ID sets are mutually exclusive. This is a genuine, documented
methodological uncertainty — not resolved by guessing — flagged for Stage 1
follow-up (batch RCSB entity-length verification across all IDs in 145541
and 150035).

## Extraction

All archives extracted losslessly via `src/ingestion/safe_extract.py`
(idempotent, rejects unsafe/traversal paths, never silently overwrites).
Extracted file counts match archive entry counts exactly for every cohort
(see `extraction_summary.tsv`).

**Data-integrity issue found and fixed:** this Mac's APFS volume is
case-insensitive. 18 archive members across 5 archives had filenames
differing only by chain-letter case (e.g. `..._B_+_abag.cif` vs
`..._b_+_abag.cif`), which would silently collide/overwrite under naive
extraction. `safe_extract.py` was extended to detect this specific
collision type and extract the colliding member under a disambiguated
`__casecollision2` filename, with the original→extracted name mapping
recorded in `data/metadata/extraction_manifest.tsv`. Zero data loss; zero
unresolved conflicts remain.

The duplicate archive (151141) was extracted once for verification only;
its output was identical to 150654's (same case-collision resolution), so
the redundant copy was removed to avoid wasted storage. Its original
`.tgz` remains untouched.

## Environments

Three micromamba environments created and functionally verified (real
imports/CLI calls, not just `conda list`):

- `vhh3d` (Python 3.11.16) — all core scientific packages import
  successfully; `cd-hit` and `mkdssp` CLIs work; 6/6 pytest smoke tests
  pass.
- `vhh-arpeggio` (Python 3.9.23) — PDBe-Arpeggio 1.4.4 installed via pip
  inside the isolated env; both the `pdbe-arpeggio` CLI and `import
  arpeggio` work.
- `vhh-elec` (Python 3.11.16) — APBS 3.4.1 and PDB2PQR 3.6.1 installed
  cleanly from conda-forge (no pip fallback needed); `apbs`, `pdb2pqr`, and
  `pdb2pqr30` all functional.

Full version list: `software_versions.txt`, `requirements_versions_summary.tsv`.
Reproducible exports: `environment-main.yml`, `environment-arpeggio.yml`,
`environment-elec.yml`.

## TNP reference

Cloned (source only) from `github.com/oxpig/TNP` to `external/TNP`, commit
`29dcac72f1380e8538e8870f45a699d3c6156162` on `main`. No dependencies
installed; NanoBodyBuilder2 not run; no structures predicted.

## AggreScan3D

Docker available (29.2.1). No installation attempted — aggregation-proxy
stage is secondary/tertiary and non-blocking for the core pipeline.

## Molecular dynamics

Confirmed: no MD engine (GROMACS/AMBER/NAMD/OpenMM/CHARMM/Desmond/other)
installed or configured at any point in this stage.

## Testing

`pytest -v tests/test_environment.py` → 6 passed, 0 failed (see
`logs/setup_errors.log` for one non-blocking synthetic-mmCIF test issue
that was resolved by testing against a real extracted structure instead).

## Data integrity confirmation

- Original `.tgz` archives: unmodified (verified by re-checksumming is
  unnecessary since no write operation ever targeted them; checksums
  recorded at time of inventory in `raw_archive_inventory.tsv`).
- Extracted raw mmCIF files: unmodified after initial extraction; no
  downstream stage has written into `data/raw_extracted/`.
- No raw files deleted.
- No scientific filtering, cohort classification beyond archive-level
  identification, or biological descriptor calculation performed.
