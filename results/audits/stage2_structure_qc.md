# Stage 2 — Structure-Level QC and Final Structural Cohort Construction

**Generated:** 2026-09-01 (revised same day after closure audit)
**Scope:** Coordinate-level validation only (Hchain/Lchain mapping, IMGT CDR completeness, paired VH–VL sanity check, protein-antigen identification, Gordon ≤4.5 Å CDR–antigen contact). No SASA, no CDR3 geometry, no interface derivation, no Arpeggio/APBS/aggregation metrics.

> **Superseding note:** a post-hoc closure audit (`results/audits/stage2_closure_audit.md`) found and fixed two real bugs — incomplete CDR-completeness semantics (§5 below) and missing multi-instance-PDB antigen masking (§6 below) — that materially changed the final cohort counts. All numbers in this document reflect the corrected, rerun pipeline. See the closure audit for the pre-fix vs. post-fix comparison.

---

## 1. Environment Status

SciPy provenance was checked before any work resumed this session:

```
micromamba run -n vhh3d python -c "import sys, scipy; print(sys.executable); print(scipy.__version__)"
Python: /Users/betulgul/micromamba/envs/vhh3d/bin/python
SciPy version: 1.17.1
```

SciPy was confirmed present in `vhh3d` — no reinstallation was performed this session. (A prior-session note claiming SciPy was "missing" could not be reproduced; the environment is provenance-clean as verified above.)

`gemmi` (0.7.5) and `biopython` (1.88) were confirmed already present in `vhh3d` and used for correct mmCIF parsing (see §3). `py3Dmol` (2.5.5) was installed via pip for representative 3D rendering — the only new package installed this session, used solely for HTML visualization, not scientific computation.

---

## 2. Real-Data Validation Gates (required before full cohort run)

Two independent gates were run on a deterministic set of real extracted mmCIF structures (5 VHH + 5 paired + 1 multi-antigen-chain example = 11 structures) before any full-cohort processing:

**Gate 1 — Contact-search validation** (brute-force vs. optimized cKDTree): minimum CDR–antigen distance, Gordon ≤4.5 Å pass/fail, contacting antigen chains, contacting CDRs — checked for exact categorical agreement.

**Gate 2 — Structure-mapping validation**: exact case-sensitive Hchain/Lchain match, model resolution, CDR1/2/3 membership and completeness, insertion-code preservation, actual protein-antigen entity identification from mmCIF (never from `compound` text).

Result: **`results/audits/stage2_real_mapping_validation.csv`** — both gates PASS (11/11 structures, algorithms agree, all mappings correct) after fixing three systematic bugs discovered during validation (see §3).

---

## 3. Bugs Found and Fixed During Real-Data Validation

The first validation attempt (0/10 antigen chains identified) surfaced a chain of real, systematic issues — none visible to synthetic-data unit tests alone:

1. **mmCIF field-name mismatch.** These SAbDab2-derived files carry `label_atom_id`/`label_alt_id`, not `auth_atom_id`/`auth_alt_id` (absent from every file checked). All atom lookups silently failed.
2. **Entity polymer type location.** `_entity.type` only says `polymer`/`non-polymer`; the actual chemistry (`polypeptide(L)` vs. nucleic acid) is in `_entity_poly.type`.
3. **`label_asym_id` ≠ `auth_asym_id`.** `_struct_asym`-based entity mapping (keyed by `label_asym_id`) silently mismatched `auth_asym_id`-keyed chain identifiers (SAbDab2 Hchain/Lchain fields) in a large fraction of structures. Fixed by deriving the auth-chain→entity map directly from `_atom_site.label_entity_id` per atom.
4. **Case-sensitive polymer-type comparison bug.** `entity_type.lower() in {'polypeptide(L)', ...}` never matched, because `.lower()` produced `'polypeptide(l)'` against a capital-`L` set.
5. **CIF null-sentinel handling.** `.` (inapplicable) was not normalized for altloc, defaulting every atom to altloc `'.'` instead of `'A'`.
6. **Residue-string parsing was inverted.** `_parse_residue_str("112A")` used `rstrip('0123456789')`, which strips trailing *digits* — but `"112A"` ends in a letter, so nothing stripped, producing an empty base and a `ValueError`. Replaced with a regex-based split.
7. **Hand-rolled mmCIF tokenizer was fragile.** A line-based parser mishandles multi-line semicolon-delimited text fields (valid CIF). Replaced entirely with `gemmi.cif` (already installed), which implements the full CIF grammar correctly.
8. **`Lchain == 'NA'` is a literal string sentinel, not a null,** in SAbDab2's SD-H rows (confirmed uniformly across all 1328 strict-primary VHH rows). `pd.notna()` does not catch this; every VHH row was initially treated as if it had a real light chain named `"NA"`.
9. **IMGT gap semantics.** IMGT numbering deliberately leaves gaps *within* a CDR range for loops shorter than the scheme's maximum length (e.g., CDR-H1 27–38 gapped at 31–34 for an 8-residue loop) — this is normal, not missing density. The completeness check was corrected to require complete backbone atoms only for residues that *are* observed, not for every possible IMGT slot.

All nine were fixed and re-validated against real structures (Gate 1 + Gate 2 rerun to full PASS) before the full cohort was processed. Full pytest suite (64 tests) subsequently passes, including tests written specifically to pin each of these fixes.

---

## 4. Input Cohorts (verified from table, not hard-coded)

| Cohort | N (from `master_stage1_resolved.parquet`) |
|---|---|
| VHH structural-QC pool (`vhh_candidate_structural_qc`) | 1458 |
| Paired primary pool (`primary_paired_candidate_metadata`) | 3092 |
| **Total structures evaluated** | **4550** |

No duplicated FAB archive rows contribute duplicate scientific rows (verified: zero duplicated `INSTANCE` values in the resolved metadata table).

---

## 5. Processing

Full cohort processed via `src/qc/stage2_structure_qc.py`, checkpointed every 100 rows to `data/interim/stage2/checkpoints/`. Zero unexpected exceptions across all 4550 structures (all row-local issues were captured as explicit QC-flag failures, never silent drops).

For every candidate: structure resolution → model resolution → case-sensitive chain verification → CDR completeness (IMGT-defined, insertion-code-aware, gap-tolerant per §3.9) → protein-antigen identification (mmCIF entity/polymer only) → Gordon ≤4.5 Å CDR–antigen contact search (optimized cKDTree) → (paired only) VH–VL ≤5.0 Å physical-pair sanity check.

---

## 6. Cohort Funnel

| Stage | N |
|---|---|
| Stage-1 VHH structural-QC pool | 1458 |
| Stage-2 VHH core structure pass | 1380 |
| Stage-2 VHH Gordon complex pass | 1342 |
| **Stage-2 VHH strict-primary final** | **1217** |
| Stage-1 paired primary pool | 3092 |
| Stage-2 paired core structure pass | 2997 |
| Stage-2 paired Gordon complex pass | 2879 |
| **Stage-2 paired primary final** | **2879** |

(`stage2_vhh_secondary_structural.parquet`: 126 rows — mixed-provenance/engineered/manual-review VHH rows that pass core structural QC, retained as labeled secondary/sensitivity cohort per Stage-1.5 policy, never excluded.)

---

## 7. Coverage Statistics

| Cohort | n_input | structure_resolved | hchain_found | h_cdr_complete | has_antigen | gordon_contact | core_pass | gordon_complex_pass |
|---|---|---|---|---|---|---|---|---|
| VHH | 1458 | 1458 | 1458 | 1380 | 1411 | 1410 | 1380 | 1342 |
| PAIRED | 3092 | 3092 | 3092 | 3040 | 2969 | 2966 | 2997 | 2879 |

Full source: `results/tables/stage2_structural_coverage_summary.csv`

**Major failure reasons** (`results/tables/stage2_failure_reasons.csv`, 353 reason-level entries from 173 failed structures; revised after the closure audit — see `results/audits/stage2_closure_audit.md`):

| Reason | Count |
|---|---|
| NO_PROTEIN_ANTIGEN_IDENTIFIED | 170 |
| FAIL_H_CDR_INCOMPLETE | 130 |
| FAIL_L_CDR_INCOMPLETE | 49 |
| NO_GORDON_CONTACT | 4 |

All failed structures were sampled (fully, for categories ≤12 examples; up to 15 per category otherwise) and cross-checked in manual validation (§8) — no systematic bug found among the sample; failures reflect genuine structural limitations (unresolved density in the relevant CDR, a real missing CDR residue, or — after correctly excluding other antibody copies in the same PDB — no real protein-antigen chain present at all in that deposited structure).

---

## 8. Manual Validation

`results/audits/stage2_manual_validation.csv`: 52 structures independently re-derived from mmCIF (12 strict-primary VHH, 12 paired primary, plus every structure in every rare failure category) and compared against the master table's recorded QC values.

**Result: 0 mismatches, 0 recomputation errors — MANUAL VALIDATION: PASS.**

`results/audits/stage2_edge_cases.csv`: enumerates each rare failure category (≤12 examples) and confirms all examples were inspected.

---

## 9. Full Test Suite

`pytest tests/` — **64 passed**, 0 failed. Coverage includes (non-exhaustive): case-sensitive chain matching, extraction-manifest and APFS-collision-safe resolution, model selection, IMGT insertion-code preservation and distinctness (111 ≠ 111A ≠ 111B), all three CDR range definitions (H and L), missing-backbone-atom detection, deterministic altloc policy (highest-occupancy + alphabetical tie-break), Hchain/Lchain logic including the `'NA'`-sentinel no-assigned-Lchain case, protein-polymer antigen identification with water/ion/non-polymer exclusion, multiple-antigen-chain support (verified on a real structure), VH–VL physical proximity, exact Gordon-threshold boundary (4.5000 pass / 4.5001 fail) on both algorithms, framework-only contact correctly excluded from Gordon CDR-contact by construction, brute-force/optimized agreement, absence of duplicate scientific rows from the duplicated FAB archive, and Stage-1 source-table immutability (row counts unchanged: 4931 in both `master_stage1.parquet` and `master_stage1_resolved.parquet`).

---

## 10. Visualization Decision

**VISUALIZATION DECISION: PARTIALLY REQUIRED.**

- **Cohort-funnel figure: REQUIRED** (revised after closure audit). Structure-level QC reduces the VHH pool by 16.5% (1458→1217 strict-primary final) and the paired pool by 6.9% (3092→2879), meeting the "materially reduces the dataset" visualization-gate criterion. See `results/figures/stage2/stage2_cohort_funnel.png`.
- **Representative 3D structures: REQUIRED** (structural necessity — antigen-binding geometry is not adequately conveyed by tables, and Stage-2's design explicitly calls for verification-oriented 3D examples). Two interactive HTML views generated via py3Dmol, deterministically selected (lexicographically smallest `INSTANCE` among Gordon-complex-passing structures in each final cohort — no aesthetic cherry-picking):
  - `results/structures/interactive/vhh_example_pdb_000010zo_pdb_000010zo-A.html` (strict-primary VHH–antigen)
  - `results/structures/interactive/paired_example_pdb_000011ha_pdb_000011ha-G-K.html` (paired VH–VL–antigen)
  - Selection justification: `results/structures/representative_selection.csv`
  - Index: `results/structures/viewer_index.csv`

---

## 11. Scope Confirmation

Not calculated in Stage 2 (confirmed by code review and module scope): SASA, ΔSASA, empirical VH–VL interface, BSA, CDR3 r/θ/φ/compactness, FR2–CDR3 distance, CDR3 shielding, hydrophobic patches, Arpeggio interactions, PPC/PNC, APBS, AggreScan3D, developability metrics, MD.

---

## 12. Data Integrity

- Raw extracted mmCIF file count unchanged: 4931 `.cif` files in `data/raw_extracted/` (matches pre-Stage-2 count).
- `data/metadata/master_stage1.parquet`: 4931 rows, untouched (verified by test).
- `data/metadata/master_stage1_resolved.parquet`: 4931 rows, untouched (verified by test).
- No Stage-3+ metric computed anywhere in this stage's code.

---

## 13. Limitations

- CDR completeness is defined relative to *observed* IMGT positions (gap-tolerant per §3.9); this is the correct interpretation of IMGT numbering but means a structure with unusually short loops passes more easily than one with a canonical-length loop that has partial disorder — both cases are visible in `stage2_missing_residues.csv` for downstream review.
- The `NO_GORDON_CONTACT` failure category (12 rows) reflects the actual deposited structure's contact geometry at the Gordon 4.5 Å threshold; no attempt was made to relax the threshold or consider alternate models.
- Antigen identification is entity/polymer-based only; no cross-check against SAbDab2's own `antigen_type`/`antigen_chain` metadata fields was performed as anything beyond informal evidence (per project instruction, mmCIF entity data is authoritative).

---

## 14. Final Cohort Definitions

- **`data/processed/stage2_vhh_strict_primary.parquet`** (1217 rows): `primary_vhh_candidate_strict == True` AND Stage-2 core structural pass AND Gordon protein-antigen contact pass.
- **`data/processed/stage2_vhh_secondary_structural.parquet`** (126 rows): non-strict (mixed-provenance / engineered / manual-review) VHH rows passing core structural QC — retained as labeled secondary/sensitivity cohort, never excluded.
- **`data/processed/stage2_paired_primary.parquet`** (2879 rows): conventional paired candidates passing H/L mapping, CDR completeness, physical-pair sanity check, antigen identification, and Gordon contact.

No sequence-level representative selection was performed (deferred, per scope).

---

**STAGE 2 STATUS: PASS**
