# Stage 1.5 — Resolution Policy (derived eligibility layer)

This is a **new derived layer**, not a rewrite. `data/metadata/master_stage1.parquet`
(and `.csv`) are unchanged — verified: still 96 columns, 4931 rows,
`primary_vhh_candidate_metadata` sum still 1419. All policy output lives
in `data/metadata/master_stage1_resolved.{csv,parquet}` (102 columns —
the original 96 plus 6 new derived fields), with the same 4931 rows.

## 1. Contamination flag retracted

The Stage-1 "58 suspected contaminant rows" interpretation for
`sabdab2_H02Y5L0000` is retracted, per the Stage-1.5 chain-level
verification: the SabDab2 Hchain in every one of the 11 PDBs is a genuine
short VH/VHH-family domain, structurally and sequence-wise separate from
the MHC/HLA/beta-2-microglobulin antigen material. **This is no longer
treated as contamination anywhere downstream.** The original finding and
its correction remain preserved in
`results/audits/stage1_classification_boundary_audit.md` (the false
positive) and `results/audits/stage1_5_targeted_verification.md` (the
correction) as the historical record of how the false positive arose and
was resolved — not deleted or overwritten.

## 2–3. Derived eligibility fields and policy

New fields on `master_stage1_resolved`: `primary_vhh_candidate_strict`,
`vhh_candidate_structural_qc`, `vhh_candidate_sensitivity`,
`sequence_provenance_ambiguity`, `provenance_resolution_status`,
`natural_camelid_mutated_flag`.

**Policy priority order** (rule B is evaluated before the
classification-based rules A/C/D/E, since the Stage-1.5 finding applies at
the SABDAB_ID/sequence level regardless of which per-instance species
label an individual row happens to carry — a `human_sdVH_candidate`-labeled
row that shares a sequence-identical SABDAB_ID with confirmed
natural-camelid rows is governed by rule B, not by its own label under
rule E):

1. **Rule B** — `SABDAB_ID` in the 5 mixed-provenance set (108 rows, 5
   IDs, spanning `natural_camelid_candidate`, `human_sdVH_candidate`, and
   `uncertain_sdH` labels underneath): `primary_vhh_candidate_strict=False`,
   `vhh_candidate_structural_qc=True`, `vhh_candidate_sensitivity=True`,
   `provenance_resolution_status=MIXED_PROVENANCE_SAME_SEQUENCE`.
2. **Rule A** — remaining `natural_camelid_candidate` (1328 rows):
   `primary_vhh_candidate_strict=True`, structural_qc=True, sensitivity=True,
   `provenance_resolution_status=RESOLVED_NATURAL_CAMELID`.
3. **Rule C** — `engineered_or_humanized_camelid_candidate` (3 rows):
   `primary_vhh_candidate_strict=False`, structural_qc=True, sensitivity=False,
   `provenance_resolution_status=RESOLVED_ENGINEERED_OR_HUMANIZED`.
4. **Rule D** — `uncertain_sdH` (19 rows — 2 of the original 21 were
   absorbed into rule B above): `primary_vhh_candidate_strict=False`,
   structural_qc=True, sensitivity=False,
   `provenance_resolution_status=MANUAL_REVIEW`.
5. **Rule E** — `human_sdVH_candidate` / `other_noncamelid_sdH` /
   `synthetic_or_de_novo_candidate`: `primary_vhh_candidate_strict=False`,
   structural_qc=False, sensitivity=False,
   `provenance_resolution_status=EXCLUDED_<CATEGORY>`. Retained in the
   table for explicitly labeled secondary/control analyses, never in the
   primary cohort.

FAB/FV rows are outside the scope of this SD-H-specific policy; all six
new fields are `None` for them.

## 4. Mutated natural-camelid VHHs

`natural_camelid_mutated_flag` set wherever the row is a resolved natural
candidate (rule A) or a mixed-provenance natural-labeled instance (rule
B) and Stage-1's `mutation_flag` was already `True`. **88 rows.** Mutation
status is not used as an exclusion criterion anywhere in this policy — a
future sensitivity analysis can compare
`primary_vhh_candidate_strict & ~natural_camelid_mutated_flag` against
the full strict set.

## 5. SABDAB_ID semantics — recorded conclusion

Stage 1.5 found **zero amino-acid mismatches at any shared IMGT position**
across all five investigated conflicting IDs. This empirically supports
continuing to use `SABDAB_ID` as the principal sequence-level redundancy
grouping key for deduplication. However, `SABDAB_ID` groups **sequence
identity**, not **biological provenance** — the same ID can legitimately
span rows with different (and in the H06BKL0000 case, internally
self-contradictory) species/taxid annotation. Sequence identity and
provenance identity must be treated as separate axes from here forward.

## Headline counts

| Metric | Count |
|---|---|
| Strict primary natural-VHH candidates | **1328** |
| Structural-QC VHH candidates (strict + mixed-provenance + uncertain + engineered) | **1458** |
| Mixed-provenance total rows (all 5 SABDAB_IDs, any prior label) | **108** |
| Unique SABDAB_IDs affected by mixed provenance | **5** |
| Uncertain (manual review) | **19** |
| Engineered/humanized secondary | **3** |
| Mutation-flagged natural-camelid rows | **88** |
| H02Y5L0000 still considered contamination? | **No — retracted** |

### Correction: "108 rows excluded from strict" was imprecise

An earlier version of this report stated "Rows excluded from strict solely
by mixed provenance: 108." That conflated two different things: the total
size of the 5 mixed-provenance SABDAB_IDs (108 rows) is **not** the same
as how many of those rows were actually removed *from the strict primary
count*, since only a subset of the 108 carried `natural_camelid_candidate`
in the first place — the rest were already non-primary under Stage 1
(`human_sdVH_candidate` or `uncertain_sdH`) and so were never part of the
1419 to begin with. Corrected breakdown:

| Quantity | Value |
|---|---|
| `mixed_provenance_total_rows` | **108** |
| `mixed_provenance_previously_natural_rows` (were `natural_camelid_candidate`, i.e. part of the original 1419) | **91** |
| `mixed_provenance_previously_nonprimary_rows` (were `human_sdVH_candidate` or `uncertain_sdH`, never in the 1419) | **17** |
| `strict_primary_rows_removed_due_to_mixed_provenance` | **91** |

**1419 − 91 = 1328** strict primary natural-VHH candidates — confirmed.
The other 17 mixed-provenance rows are correctly captured under
`MIXED_PROVENANCE_SAME_SEQUENCE` too (rule B applies to all 108
regardless of prior label, per the documented priority order), but they
don't subtract from the 1419→1328 arithmetic since they were never in the
1419.

### Why the derived uncertain/manual-review count is 19, not 21

Verified directly against `master_stage1_resolved.parquet`:

- `sdh_classification == "uncertain_sdH"` (original Stage-1 field,
  unmodified): **21** rows.
- Of those 21, exactly **2** belong to one of the five mixed-provenance
  SABDAB_IDs (`sabdab2_H02Y5L0000`) and are therefore reassigned to the
  higher-priority `MIXED_PROVENANCE_SAME_SEQUENCE` derived status (rule B
  takes precedence over rule D per the documented priority order — see
  §2–3 above).
- Remaining rows with `provenance_resolution_status == "MANUAL_REVIEW"`:
  **19** = 21 − 2.

This does **not** alter `sdh_classification`, which still reads
`uncertain_sdH` for all 21 original rows — only the *derived*
`provenance_resolution_status` field differs for the 2 that overlap with
the mixed-provenance set.

### Decomposition of `vhh_candidate_structural_qc = 1458`

| `provenance_resolution_status` | Rows |
|---|---|
| RESOLVED_NATURAL_CAMELID | 1328 |
| MIXED_PROVENANCE_SAME_SEQUENCE | 108 |
| MANUAL_REVIEW | 19 |
| RESOLVED_ENGINEERED_OR_HUMANIZED | 3 |
| **Total** | **1458** |

1328 + 108 + 19 + 3 = 1458, confirmed exactly.

Sanity check: 1328 (strict) + 91 (mixed-provenance rows previously
labeled `natural_camelid_candidate`) = 1419 = the original Stage-1
`primary_vhh_candidate_metadata` count, confirming the derived layer is a
strict re-partition of the original primary set, not a re-derivation from
scratch.

## Confirmation

No additional sequence or structural analysis was performed in this
step — only application of the policy function to already-computed
Stage-1/Stage-1.5 fields.

## STAGE 1 / 1.5 FINAL STATUS: PASS
