# VHH structural mechanism project

## Primary biological question

How can VHH (single-domain camelid) antibody domains remain soluble,
structurally viable, and functional in the complete absence of a VL domain,
given that in conventional antibodies the VL domain physically covers part
of the VH surface?

## Short mechanistic hypothesis

Loss of VL is not compensated by a single sequence motif. Compensation, if
present, may involve some combination of: VHH-specific sequence adaptation
at the empirically defined former VH–VL interface; how that sequence
chemistry translates into actual folded-surface exposure, hydrophobicity,
patch topology, and electrostatics; CDR3 geometric behavior relative to
that surface (orientation, compactness, and shielding are measured, not
assumed); and, where shielding occurs, its physical basis in intramolecular
packing (hydrogen bonds, hydrophobic, aromatic, cation-π, ionic, van der
Waals contacts). No compensatory, complementary, or independent relationship
between interface chemistry and CDR3 shielding is assumed in advance — all
are tested.

## Expected datasets

Four SAbDab2 export categories were expected: SD-H, FV, FAB, FAB+FC.
Six `.tgz` archives were actually downloaded (one exact duplicate, and two
whose FAB vs FAB+FC identity could not be confidently determined from crop
content alone). See `results/audits/archive_cohort_mapping.tsv` for the full
evidence-based mapping and confidence levels.

## Folder structure

```
data/
  raw_extracted/   # read-only extraction of the original .tgz archives (see README_RAW_DATA.md)
  metadata/        # extraction manifests, chain-mapping tables
  interim/         # intermediate derived data (not final)
  processed/       # analysis-ready derived data
src/
  ingestion/ qc/ interfaces/ geometry/ physicochemistry/
  interactions/ antigen/ aggregation/ statistics/ utils/
tests/             # unit and smoke tests
configs/           # config.yaml — analysis parameters, IMGT windows, thresholds
results/
  tables/ figures/ structures/ audits/
logs/              # setup and pipeline logs
external/          # third-party reference source (e.g. TNP); not tracked in this repository, see below
references/
notebooks/
```

## Environments

Three isolated micromamba environments (created via `environment-*.yml`):

- `vhh3d` — main scientific Python environment (Python 3.11): numpy, pandas,
  scipy, scikit-learn, statsmodels, matplotlib, biopython, gemmi, pyarrow,
  freesasa, cd-hit, mkdssp, etc.
- `vhh-arpeggio` — isolated environment (Python 3.9) for PDBe-Arpeggio
  interaction analysis.
- `vhh-elec` — isolated environment (Python 3.11) for APBS / PDB2PQR
  continuum electrostatics.

See `software_versions.txt` and `requirements_versions_summary.tsv` for
exact verified versions.

## Data integrity

All original `.tgz` archives and everything under `data/raw_extracted/` are
immutable and read-only. No pipeline stage may overwrite them. See
`data/raw_extracted/README_RAW_DATA.md` for details, including a documented
case-insensitive-filesystem extraction issue and its resolution.

## Molecular dynamics

Molecular dynamics is explicitly out of scope for the current pipeline. No
MD software has been installed. It may be considered only as future work.

## Current pipeline stage

**Stage 0 (environment setup and raw-data extraction) only.** No scientific
filtering, cohort classification, SASA/geometry/interaction/electrostatics
calculation, or statistical analysis has been performed. No biological or
structural result exists yet in this repository.

## External dependency: TNP

Stage 8 (`src/qc/stage8_patch_scoring.py`) uses the `psa` binary and parameter files from the Therapeutic Nanobody Profiler (TNP). Its source is not redistributed here. To reproduce Stage 8, clone it into `external/`:

```bash
git clone https://github.com/oxpig/TNP.git external/TNP
```

The scoring script expects the binary at `external/TNP/bin/psa_mac`; adjust `PSA_BIN` in `src/qc/stage8_patch_scoring.py` for your platform and install location. TNP is distributed under its own license.
