# Stage 1 — Metadata Ingestion & Cohort Classification QC

## Scope change carried in from this stage's kickoff

FAB+FC has been deliberately removed from the primary project scope
(project-design decision, not a data gap). Primary raw cohorts are now
**SD-H, FAB, FV**. FAB+FC remains an optional future-expansion dataset;
its absence is not treated as a failure at any stage from here on.
`configs/config.yaml` updated (`cohorts.primary_raw_cohorts` /
`cohorts.optional_future_expansion`).

## Input data

Three official SAbDab2 summary CSVs (SD-H 1597 rows, FAB 2798 rows, FV 536
rows) combined into one master table, joined to archive/extraction
provenance via `data/metadata/archive_structure_to_summary_matches.tsv`
(itself built from the full `data/raw_extracted/extraction_manifest.tsv`).
The content-duplicate FAB archive (`...151141.tgz`) is represented only in
provenance columns (`source_archive`, `content_duplicate_flag`) — it
contributes **zero** additional rows to the master table.

## Outputs produced

- `data/metadata/master_raw.csv` / `.parquet` — combined CSVs + provenance
  columns, before any QC/classification columns.
- `data/metadata/master_stage1.csv` / `.parquet` — full table with QC
  flags, SD-H classification, construct classification, exclusion
  reasons. **4931 rows, 96 columns.** No row physically removed at any
  point.
- `results/tables/stage1_counts.csv`, `stage1_species_counts.csv`,
  `stage1_sdh_classification.csv`, `stage1_construct_counts.csv`,
  `stage1_exclusion_reasons.csv`, `stage1_redundancy_inventory.csv`,
  `stage1_top20_repeated_sabdab_ids.csv`.
- `results/audits/stage1_known_example_checks.json`,
  `stage1_headline_numbers.json`.

## Basic structural-quality revalidation

All three cohorts were re-checked from scratch (method string normalized,
resolution coerced to numeric) rather than trusting the original SAbDab2
web-search filters. Result: **100% of all 4931 rows are X-ray with
resolution ≤ 3.0 Å** in every cohort (SD-H 1597/1597, FAB 2798/2798, FV
536/536). The originally stated search criteria held up under
revalidation; no violations found. R-factor/R-free missingness recorded
separately (5 / 24 missing overall) and not used as a filter at this
stage, per instruction.

## SD-H classification — evidence hierarchy applied

Priority order used per row: `heavy_taxid` → `heavy_species` →
`organism_taxid`/`organism` → `heavy_subclass` → free text
(`compound`/`short_header`). Structured evidence (taxid/species fields)
was checked for internal self-consistency first (`heavy_taxid` vs
`organism_taxid`: 0 disagreements found across all 1597 rows), then used
to set `normalized_species` / `camelid_species_flag` before any keyword
search of free text ran.

Six-way classification (1597 SD-H rows):

| Bucket | n |
|---|---|
| natural_camelid_candidate | 1419 |
| human_sdVH_candidate | 91 |
| synthetic_or_de_novo_candidate | 61 |
| uncertain_sdH (manual review) | 21 |
| engineered_or_humanized_camelid_candidate | 3 |
| other_noncamelid_sdH | 2 |

`primary_vhh_candidate_metadata = True` for all 1419
`natural_camelid_candidate` rows (all of them already pass X-ray +
resolution ≤ 3 Å, so no rows were dropped by the structural-quality
filter at this step).

### Notable edge cases surfaced (not silently resolved)

- **3 rows with `heavy_taxid=562` (E. coli)** — implausible as a
  biological antibody-source organism; `compound` text uses
  "Nanobody"/`Nb`-prefixed naming in 1 of the 3, consistent with
  `heavy_species`/`heavy_taxid` having been populated from the expression
  system rather than the source organism upstream. Classified
  `uncertain_sdH`, `metadata_conflict_flag=True`, full reasoning preserved
  in `manual_review_reason`. Not silently classified either way.
- **2 rows with fully missing structured species evidence** (`heavy_taxid`,
  `heavy_species`, `organism_taxid`, `organism` all `NA`) — classified
  `uncertain_sdH` rather than inferred from text alone.
- **6 rows: camelid-taxid heavy chain (`Lama glama`, taxid 9844) but
  `short_header=DE NOVO PROTEIN` / `compound` uses "artificial"/"de novo"
  language** (e.g. `pdb_00005bop`, "artificial nanobody octarellinV.1") —
  genuine structured-vs-text tension between camelid germline framework
  and an explicitly artificial/designed scaffold. Flagged
  `metadata_conflict_flag=True`, routed to `uncertain_sdH` rather than
  assumed natural.
- **1 row, `mus musculus` heavy chain** (`pdb_00006zqk`, compound
  "HER2-binding scFv-Fab fusion 841") — a heavy-chain-only crop of an
  engineered mouse scFv-Fab fusion that satisfied SAbDab2's structural
  SD-H criteria, not a biological single-domain antibody. Correctly
  routed to `other_noncamelid_sdH`, exactly the scenario the project brief
  warned SD-H can contain.
- **3 rows explicitly humanized camelid VHHs** (`pdb_00003eba`,
  `pdb_00006rnk`, `pdb_00006z10`) — camelid-taxid heavy chain +
  `humanized` keyword hit in `compound`. This is not treated as a
  conflict (camelid origin + humanization modification are expected to
  co-occur); routed to `engineered_or_humanized_camelid_candidate`.
- **`mutant`/`variant`/`affinity-matured` keyword hits never by themselves
  moved a row out of `natural_camelid_candidate`** — verified directly
  (e.g. camelid-taxid rows with a `mutant` compound-text hit remained
  `natural_camelid_candidate`, `metadata_conflict_flag=False`), per
  explicit instruction and covered by
  `tests/test_metadata.py::test_mutant_alone_does_not_exclude_natural_camelid`.
- **`camelidae mixed library` (taxid 1579311, 18 rows)** — a pooled/naive
  camelid phage-display library annotation (not a single wild-type
  individual, but still genuinely camelid-derived; `heavy_subclass` in
  several rows independently confirms `IGHV3 (Vicpac)` germline usage).
  Classified `natural_camelid_candidate` with no conflict, since all
  available evidence (taxid, species text, subclass) agrees on camelid
  origin — the "mixed library" qualifier describes selection provenance,
  not a synthetic/humanized/non-camelid identity.

## Known-example sanity checks

| PDB | Expected | Present | Actual classification | Match |
|---|---|---|---|---|
| 1BZQ | natural *C. dromedarius* VHH | Yes (4 rows) | `natural_camelid_candidate`, species=Camelus dromedarius, resolution 2.8 Å | ✅ |
| 2VYR | human sdVH, must not enter primary VHH cohort | Yes (8 rows) | `human_sdVH_candidate`, `primary_vhh_candidate_metadata=False` | ✅ |
| 27UH | synthetic/de-novo | Yes (4 rows) | `synthetic_or_de_novo_candidate`, species=synthetic construct | ✅ |
| 1G6V | should fail resolution ≤3 Å | **Not present** in this dataset | — | reported as absent, not forced |

All three present examples classify exactly as expected; no forced or
contradicted classification was needed. Full detail in
`results/audits/stage1_known_example_checks.json`.

## Paired (FAB/FV) construct classification

Construct architecture recomputed directly from the official
`chainsharing_construct` field (not hard-coded from the earlier Stage-0
observation):

**FV (536 rows):** conventional 294, SCFV 223, DIABODY 17, DVD 1, CODV 1
— exactly reproducing the Stage-0 estimate, now derived programmatically.

**FAB (2798 rows):** conventional 2798, 0 engineered chain-sharing
constructs.

`paired_metadata_candidate` requires only Hchain+Lchain present and
`type in {FAB,FV}` — it does **not** assume FAB/FV rows are automatically
usable paired antibodies (per instruction). `primary_paired_candidate_metadata`
additionally requires X-ray, resolution ≤3 Å, and no chain-sharing
construct: **3092 total** (FAB 2798, FV 294). Coordinate-level VH/VL
completeness is explicitly deferred to Stage 2 — nothing here claims a
row is a "usable paired VH–VL" structure yet.

## Protein antigen metadata flag

`metadata_contains_protein_antigen` set wherever `antigen_type` contains
the substring `PROTEIN` (including multi-value fields like
`PROTEIN|HAPTEN`). This is metadata-only; actual protein antigen chain
identification from mmCIF polymer/entity records remains a Stage 2 task,
as instructed — no chain-level zipping of `antigen_chain`/`antigen_type`
was performed.

## Pseudoreplication inventory (characterization only, no reduction)

| Group | Instance rows | Unique PDB | Unique SABDAB_ID | Unique HEAVY_ID | Instances/ID median | mean | max | p95 |
|---|---|---|---|---|---|---|---|---|
| SD-H | 1597 | 778 | 654 | 654 | 2.0 | 2.44 | 66 | 6.0 |
| FAB | 2798 | 1730 | 1360 | 1333 | 2.0 | 2.06 | 80 | 4.0 |
| FV | 536 | 290 | 245 | 228 | 2.0 | 2.19 | 16 | 6.0 |
| ALL | 4931 | 2734 | 2230 | 2185 | 2.0 | 2.21 | 80 | 5.0 |

Top 20 most-repeated SABDAB_IDs in `stage1_top20_repeated_sabdab_ids.csv`
(max: `sabdab2_H00DFL00BV` with 80 instance rows — a FAB entry, almost
certainly a large multi-copy asymmetric unit or a heavily-crystallized
antibody). No redundancy reduction or representative selection has been
performed — this is inventory only, as instructed.

## Duplicate archive policy verification

FAB row count in the master table is exactly **2798**, matching the
official FAB summary CSV row count precisely (not 4190, not any other
inflated number). `content_duplicate_flag=True` on exactly 1392 rows —
the FAB rows that trace back to structures present in *both* `150654.tgz`
and its duplicate `151141.tgz` — but these 1392 rows are represented once
each in the master table, with `source_archive` listing both archives.
Confirmed by `tests/test_metadata.py::test_duplicate_archive_does_not_duplicate_master_rows`.

## Design notes preserved for later mechanistic stages

- Later "VHH vs VH" comparisons will use **VH coordinates extracted from
  real paired VH–VL structures** (FAB/FV cohort), not an experimentally
  isolated VH. When VL is computationally omitted for SASA purposes, this
  is a **geometric exposure counterfactual only** — never to be described
  as an experimental isolated-VH state.
- Future empirical former-interface projection onto VHH will compute
  `interface_coverage_fraction = observed empirical-interface positions /
  expected empirical-interface positions`. No coverage cutoff has been
  chosen yet (`configs/config.yaml:
  interface_projection.coverage_threshold_primary` remains `null`).
- The future interface map will retain **both** a binary consensus
  interface **and** continuous per-position weights (interface frequency,
  mean ΔSASA, median ΔSASA) — never collapsed to binary-only. A weighted
  interface-sensitivity analysis is planned for later
  (`configs/config.yaml: former_interface.calculate_weighted_interface_sensitivity: true`).

## Testing

`pytest -v` (full suite, environment + metadata tests): **18 passed, 0
failed.**

## Data integrity confirmation

- No raw `.tgz` archive or extracted mmCIF file modified.
- No metadata row physically deleted at any point — every SD-H/FAB/FV row
  from the official summaries appears exactly once in `master_stage1`,
  carrying explicit `exclusion_reasons` rather than being dropped.
- No SASA, ΔSASA, CDR3 geometry, PSH/PPC/PNC, APBS, Arpeggio interaction,
  antigen-interface, aggregation, or final statistical calculation was
  performed at this stage.

## Status

**STAGE 1 STATUS: PASS**
