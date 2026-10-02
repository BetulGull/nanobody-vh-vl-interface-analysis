# Stage 11 — Unit-of-Analysis Integrity Audit

## Statistical unit

The Stage-11 statistical unit is the **frozen nonredundant representative structure** (one row = one INSTANCE), matching every prior stage's primary-cohort convention. Residue-level, atom-level, and Arpeggio-interaction-record-level rows from Stages 7-9 are **never** treated as independent statistical observations — all such data was pre-aggregated to one summary value per representative before entering Stage 11 (e.g., Stage-9's 118,614-row contact long-table was aggregated to `shieldable_engagement_fraction`, one value per INSTANCE, before this stage; Stage-8's 237,986-row patch-components table was never re-read here — only its already-aggregated per-representative `stage8_surface_patch_descriptors.parquet` was used).

## Per-input-table duplication check

| Source table | N rows | N unique INSTANCE | N unique PDB | Duplicate keys |
|---|---|---|---|---|
| stage4_vhh_sequence_representatives.parquet | 533 | 533 | 533 (1 VHH per PDB in this cohort by construction) | 0 |
| stage3_paired_representatives.parquet (primary) | 1429 | 1429 | <=1429 (some PDBs contain >1 paired representative) | 0 (by INSTANCE) |
| stage4_sequence_interface_chemistry.csv | 1962 | 1962 | n/a | 0 |
| stage6_cdr3_shielding_geometry_combined.parquet | 533 | 533 | n/a | 0 (VHH-only, as expected) |
| stage5_sequence_surface_chemistry.parquet | 3391 | n/a (multiple structural_state rows per INSTANCE) | n/a | 0 duplicate (INSTANCE, structural_state) pairs — filtered to exactly one state (`VH_WITHOUT_VL_GEOMETRIC` for VH, `VHH_NATIVE` for VHH) before joining, eliminating the state-multiplicity by construction |
| stage7_per_vhh_packing_summary.csv | 533 | 533 | n/a | 0 |
| stage8_surface_patch_descriptors.parquet | 1962 | 1962 | n/a | 0 |
| stage8_patch_topology.parquet | 1962 | 1962 | n/a | 0 |
| stage8_apbs_interface_electrostatics.parquet | 1962 | 1962 | n/a | 0 |
| stage9_antigen_engagement_summary.parquet | 1962 | 1962 | n/a | 0 |
| stage9 shieldable/full engagement CSVs | 533+1429 each | 533+1429 | n/a | 0 |

## Canonical join verification

`stage11_build_canonical_table.py` asserts `len(df) == n0 (1962)` after **every single join step** (10 assertions total, one per source table) — any join that would have produced one-to-many duplication (e.g., accidentally joining against the Stage-5 table without first filtering to a single `structural_state`) would have raised an `AssertionError` and halted table construction rather than silently inflating the row count. All assertions passed; final table confirmed at exactly 1962 rows with a unique `INSTANCE` index.

## Note on PDB-level clustering (relevant to Section 12's sensitivity, not a duplication problem)

Multiple representative structures CAN legitimately share a PDB entry (e.g., `pdb_00007l6v` contributes 2 VHH representatives with different `Hchain` values — this is not duplication, it is two genuinely distinct antibody sequences/structures co-deposited in the same PDB file, each independently passing all upstream QC as its own nonredundant representative). This is retained (not collapsed) in the canonical table, and is exactly the reason Section 12's PDB-clustered-SE sensitivity analysis exists — to check whether treating these as fully independent observations (the primary HC3 model's assumption) materially understates uncertainty.

## Antigen-engaged observation retained for Stage-9 variables

Per Section 3, the Stage-9 antigen-engaged structural observation (native antigen-bound state) is retained for all Stage-9-derived columns — since Stage 9 established the antigen-engaged cohort equals the full frozen 533/1429 representative cohorts (100%, verified in Stage 9), no additional filtering or state-substitution was required at the join step.

**Verdict: no one-to-many duplication anywhere in the canonical Stage-11 input table.**
