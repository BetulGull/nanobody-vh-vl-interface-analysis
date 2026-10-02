# Stage 2 — Closure Audit

**Generated:** 2026-09-01
**Scope:** Verification-only pass over the completed Stage-2 pipeline before Stage 3. Two real implementation bugs were found and fixed; no Stage-3 metric was calculated.

---

## 1. CDR Completeness Semantics — Bug Found and Fixed

**Finding:** the original implementation could not distinguish "no real residue occupies this IMGT slot" (a legitimate IMGT numbering gap) from "a real residue occupies this position but has zero resolved coordinates" (genuine missing density), because it only ever examined positions that already appeared in `_atom_site` — a residue entirely absent from the coordinate block was invisible to it by construction.

**Fix:** `MMCIFStructure.build_residue_scheme()` (new) cross-references `_entity_poly_seq` (the true, gap-free 1..N polymer residue list, keyed by `label_seq_id`) against `_atom_site` to determine, for every REAL residue in the construct, whether it was observed and whether its backbone is complete. `CDRQuality.check_cdr_completeness()` now locates each CDR's real-residue span by anchoring on the nearest observed framework residues immediately flanking the CDR (by `label_seq_id` order) and requires every real residue in that span to have complete N/CA/C/O backbone — never requiring unused IMGT slots to be filled.

**Verification:** a full-cohort comparison of old vs. new logic (4550 structures) found 163 rows whose H or L completeness flag changed. Manual inspection of `pdb_00001kxt-B` H1 confirmed the new result is correct: two real residues (GLY, TYR) are genuinely absent from density there — the old logic silently missed this.

Targeted audit with all four required example categories: **`results/audits/stage2_cdr_completeness_semantics.csv`**

| Category | PDB / Instance | CDR | Result |
|---|---|---|---|
| A — legitimate IMGT gap, complete coordinates | pdb_00001bzq-K, H1 | 27-38 (31-34 unused slots) | `True` |
| B — genuine missing-coordinate residue | pdb_00001kxt-B, H1 | 27-38 (27,28 unobserved reals) | `False` |
| C — insertion-coded residue | pdb_00001bzq-K, H3 | 105-117 (112A present) | `True` |
| D — fully complete normal CDR | pdb_00001bzq-K, H2 | 56-65 | `True` |

**Answer to §1: yes, the implementation now correctly distinguishes the two cases; it did not before.**

---

## 2. Multi-Instance PDB Antigen Masking — Bug Found and Fixed

**Finding:** `identify_antigen_chains()` only excluded the CURRENT INSTANCE's own H/L chains. For any PDB deposited with multiple independent SAbDab2 antibody copies (e.g. two Fab copies in the asymmetric unit), the other copy's H/L chains — being non-self protein polymers — were misidentified as candidate protein antigen.

Scan of the (pre-fix) master table found **444 rows** with at least one antibody chain from another INSTANCE wrongly present in `antigen_chains`.

**Fix:** `build_pdb_antibody_chain_map()` (new, in `stage2_structure_qc.py`) builds PDB → {all H/L chains across every INSTANCE in both input pools} once; `process_one_structure()` now excludes this full set (own chains ∪ all-other-instances' chains) before calling `identify_antigen_chains()`. `stage2_manual_validation.py` was updated identically for consistency.

**Verification:** `results/audits/stage2_multiinstance_antigen_masking.csv` — re-checked all 3235 INSTANCEs across all 1244 multi-INSTANCE PDBs post-fix: **0 false antibody-as-antigen assignments** (down from 444 pre-fix).

---

## 3. Independent vs. Sequential Counts

`results/tables/stage2_independent_vs_sequential_counts.csv`

| Cohort | N total | N core pass | N Gordon pass (independent) | N both core+Gordon | N final |
|---|---|---|---|---|---|
| VHH | 1458 | 1380 | 1410 | 1342 | 1217 (strict-primary metadata + core + Gordon) |
| PAIRED | 3092 | 2997 | 2966 | 2879 | 2879 (= core+Gordon; no further metadata filter) |

Independent Gordon-pass and core-pass are computed from separate criteria on separate row subsets — a row can have Gordon contact without passing core structure (rare: contact geometry can be resolved even if, say, an L-CDR is separately incomplete), or vice versa. Sequential "core → Gordon" (`qc_stage2_gordon_complex_pass`) requires both simultaneously, which is why it is ≤ min(core, Gordon-independent) in both cohorts. VHH strict-primary final (1217) is smaller still because it additionally requires `primary_vhh_candidate_strict` (a Stage-1.5 metadata condition, not a Stage-2 structural one).

---

## 4. H/L CDR Completeness Contingency Table (Paired Cohort)

`results/tables/stage2_hl_cdr_contingency_table.csv`

| | L complete | L incomplete | Row total |
|---|---|---|---|
| **H complete** | 2997 | 46 | 3043 |
| **H incomplete** | 43 | 6 | 49 |
| **Column total** | 3040 | 52 | **3092** |

(Cells sum to 3092, the full paired cohort — verified programmatically.) H-all-CDR-complete = 3040/3092; L-all-CDR-complete = 3043/3092. Both-complete = 2997; H-only-incomplete = 46; L-only-incomplete = 43; both-incomplete = 6.

---

## 5. Recomputed Final Cohorts (Direct Boolean Intersection, Not Subtraction)

| Cohort | N (direct row-level intersection) | Cross-check vs. saved `.parquet` |
|---|---|---|
| Strict-primary VHH (`primary_vhh_candidate_strict` ∧ core-pass ∧ Gordon-pass) | **1217** | matches exactly |
| Paired primary (core-pass ∧ Gordon-pass) | **2879** | matches exactly |
| VHH secondary/sensitivity (non-strict ∧ core-pass) | **126** | matches exactly |

All three were computed directly from row-level boolean masks over the regenerated master table (`stage2_closure_reconciliation.py`) and asserted equal to the counts in the corresponding saved tier files — no subtraction-derived arithmetic used.

---

## 6. Regeneration Scope

Both bugs required correcting core QC logic (`cdr_qc.py`, `stage2_structure_qc.py`), so the full 4550-structure cohort was rerun (checkpointed, 0 unexpected exceptions) and all downstream aggregate outputs (master table, tier `.parquet`/`.csv` files, antigen chain map, failure reasons, coverage summary, cohort funnel, manual validation, representative 3D structures) were regenerated from that rerun. No row was manually patched.

**Cohort funnel changed** (previously reported → corrected):

| Stage | Previous | Corrected |
|---|---|---|
| VHH core structure pass | 1452 | 1380 |
| VHH Gordon complex pass | 1450 | 1342 |
| VHH strict-primary final | 1320 | **1217** |
| Paired core structure pass | 3082 | 2997 |
| Paired Gordon complex pass | 3072 | 2879 |
| Paired primary final | 3072 | **2879** |

The drop is expected and correct: many previously-"antigen-positive" structures had their only "antigen" chain removed once other-instance antibody chains were correctly excluded (`NO_PROTEIN_ANTIGEN_IDENTIFIED` failures rose from 0 to 170), and the CDR-completeness fix surfaced genuine missing-coordinate cases that were previously invisible (`FAIL_H_CDR_INCOMPLETE`/`FAIL_L_CDR_INCOMPLETE` rose accordingly).

---

## 7. Full Verification Suite Re-run

- Full 4550-structure cohort: 0 unexpected exceptions.
- Manual validation (independently re-derived from mmCIF, using the fixed pipeline logic): 73 structures sampled (12 strict-primary VHH, 12 paired primary, plus a representative sample of every failure category including the now-larger ones) — **0 mismatches, PASS**.
- Full pytest suite: **64/64 passed** (no test modifications were needed; all existing assertions remained valid against the corrected pipeline).
- Representative 3D structures regenerated from the corrected master table (selection rule unchanged: deterministic, lexicographically-smallest INSTANCE among Gordon-complex-passing rows).

---

## FINAL REPORT

**A. Is CDR completeness now scientifically correct (IMGT gaps vs. true missing coordinates)?**
Yes — was not correct before this audit; the fix (`build_residue_scheme` + anchor-based real-residue span lookup) is verified against `_entity_poly_seq` ground truth and confirmed correct on manually-inspected real examples.

**B. Number of real missing-coordinate CDR cases:** 163 chain-level completeness flags changed between old and new logic across the full cohort (each reflecting ≥1 genuinely missing real CDR residue); see `results/audits/stage2_missing_residues.csv`-equivalent detail available via `CDRQuality.get_missing_cdr_positions()` per-row.

**C. Number of multi-INSTANCE PDBs checked:** 1244 (3235 INSTANCEs).

**D. Number of false antibody-as-antigen assignments:** 444 found pre-fix → **0 post-fix**.

**E. Independent vs. sequential Gordon-contact counts:** VHH — 1410 independent / 1342 sequential (core∧Gordon) / 1217 final (+strict-primary metadata). Paired — 2966 independent / 2879 sequential (= final, no extra metadata filter). See §3 table.

**F. Paired H/L CDR completeness contingency table:** both complete 2997, H-only-incomplete 46, L-only-incomplete 43, both-incomplete 6 (sums to 3092). See §4.

**G. Recomputed boolean-intersection final Ns:** strict VHH = **1217**, paired primary = **2879**, VHH secondary/sensitivity = **126** — all cross-checked exactly against saved tier files.

**H. Did any Stage-2 rows or final cohorts change?** Yes — both bugs required a full cohort rerun; final cohort sizes changed materially (VHH strict-primary 1320→1217, paired primary 3072→2879) as documented in §6.

**I. pytest result:** 64/64 passed after the fix (full suite rerun, no test changes required).

---

**STAGE 2 CLOSURE AUDIT: CORRECTION REQUIRED (bugs found) → CORRECTIONS APPLIED AND VERIFIED → PASS**
