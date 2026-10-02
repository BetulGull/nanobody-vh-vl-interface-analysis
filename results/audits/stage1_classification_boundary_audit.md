# Stage 1 — SD-H Boundary-Classification Audit (read-only)

No classification rule was changed. No `primary_vhh_candidate_metadata` or `sdh_classification` value was modified. No structure-, coordinate-, or sequence-level analysis was performed.

## 1. Boundary-case table

Total boundary-case rows: **135**

- uncertain_sdH: 21
- engineered_or_humanized_camelid_candidate: 3
- natural_camelid_candidate carrying at least one warning flag: 111

Full table: `results/audits/stage1_sdh_boundary_cases.csv`

## 2. Expression-system vs biological-origin separation

Evidence hierarchy applied: heavy_taxid > heavy_species > organism_taxid/organism > heavy_subclass > free text. Expression-system fields (`expression_system`, `heavy_expression_system`, and their taxids) were carried into the boundary table for inspection but were **never used** to set `normalized_species`, `camelid_species_flag`, or `sdh_classification` at any point in the Stage-1 pipeline (`src/qc/stage1_pipeline.py::classify_sdh_row` never reads any `*expression*` field).

Rows where `heavy_taxid` itself equals a known expression-host taxid (E. coli=562): **3**. These are the 3 rows already routed to `uncertain_sdH` in Stage 1 — flagged as a likely upstream field-mapping error (heavy_species/heavy_taxid populated from expression system), not resolved automatically either way.

## 5. Same-SABDAB_ID consistency check

SAbDab2 IDs with >1 SD-H instance and inconsistent classification/species/flags across instances: **14**

Full table: `results/audits/stage1_sdh_sabdab_id_consistency.csv`

## 6. Species distribution of the current primary cohort

| species_group | n_rows | unique_SABDAB_ID |
|---|---|---|
| Camelus dromedarius | 92 | 41 |
| Camelus bactrianus | 25 | 16 |
| Camelus ferus | 0 | 0 |
| other Camelus | 0 | 0 |
| Lama glama | 774 | 273 |
| other Lama | 3 | 1 |
| Vicugna pacos | 474 | 225 |
| other/unresolved Camelidae | 51 | 24 |

## 7. Flag distribution inside the primary cohort

- mutation_flag: 88
- engineering_keyword_hits: 23
- humanization_keyword_hits: 0
- synthetic_keyword_hits: 0
- metadata_conflict_flag: 0
- manual_review_reason: 0
- rows with more than one flag simultaneously: 0

## K. Suspected non-antibody contamination in the primary cohort (critical finding)

Surfaced by cross-referencing the SABDAB_ID consistency check (`sabdab2_H02Y5L0000`), then swept across the full primary cohort by compound/short_header text pattern (MHC single-chain-trimer / minibinder vocabulary). **58 primary-cohort rows (8 unique PDB, all under SABDAB_ID `sabdab2_H02Y5L0000`) describe MHC class-I single-chain-trimer constructs (e.g. 'Single chain trimer HLA-A*02:01 ... with HPV.16 E7 peptide'), not antibodies at all**, despite carrying `heavy_taxid=9844` (Lama glama) and passing every current classifier check (camelid taxid, X-ray, resolution ≤3Å, no engineering/humanization/synthetic/mutation keyword hit). Beta-2-microglobulin, a component of these engineered single-chain-trimer fusions, carries a genuine immunoglobulin-superfamily fold, which most plausibly explains why SAbDab2's own automated structural classifier admitted these into the SD-H search results in the first place — this is an upstream SAbDab2 classification artifact, not a bug in this project's metadata pipeline. Full row list: `results/audits/stage1_primary_cohort_suspected_contamination.csv`.

This was **not** caught by the current boundary-case criteria (engineering/humanization/synthetic/mutation keyword lists do not include MHC/scaffold-specific vocabulary such as 'single chain trimer' or 'HLA') — a documented gap for a future rule update, not fixed in this read-only audit.

A broader sweep of the full 1419-row primary cohort for other non-antibody engineered-scaffold vocabulary (DARPin, affibody, monobody, knottin, fibronectin/FN3, anticalin, centyrin, adnectin, repebody, 'designed binder', 'computationally designed') found **0** additional hits — the contamination found is currently isolated to this one SABDAB_ID cluster.

## 10. Paired-cohort sanity check

- FAB primary_paired_candidate_metadata count: 2798
- FV primary_paired_candidate_metadata count: 294
- Engineered chain-sharing rows incorrectly marked primary_paired_candidate_metadata=True: 0
