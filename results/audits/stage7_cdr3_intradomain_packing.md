# Stage 7 — CDR3–Former-Interface Intradomain Packing (Arpeggio)

**Cohort:** 533 native VHH representatives (frozen from `stage4_vhh_sequence_representatives.parquet`, verified count, never hard-coded).

---

## A. Real end-to-end Arpeggio smoke-test gate

See `results/audits/stage7_arpeggio_smoke_test.md` (all 16 checklist items PASS, including a real mmCIF-input-format fix and `_chem_comp` category fix discovered by direct source inspection, and explicit multi-character chain ID / insertion-code round-trip verification on `pdb_00006zbp-FFF`).

**GATE: PASS**

## B. Tool/environment versions

PDBe-Arpeggio 1.4.4, OpenBabel 3.1.0 (`vhh-arpeggio` environment), matching Stage-0 records, re-verified not trusted.

## C. Input construction

Same validated `extract_variable_domain_atoms` (Stages 3-6) restricted to the VHH `Hchain` only — no antigen, no light chain, no heteroatoms. Written via `src/qc/stage7_cif_writer.py` (minimal mmCIF + required `_chem_comp` category), with a per-instance reversible chain-identity mapping table (`data/interim/stage7/cohort_chain_mapping.jsonl`).

## D. Runtime benchmark

`results/audits/stage7_arpeggio_runtime_benchmark.csv` — 8 real structures, 2.34-3.09 s each (mean ~2.7 s). Used for planning only; did not alter the scientific cohort. Full 533-structure run completed in ~22 minutes.

## E. BSA convention

Implemented exactly as specified in `src/qc/stage7_bsa.py`: `SASA_A_alone`, `SASA_B_alone`, `SASA_AB_total` (native coordinates, validated FreeSASA engine reused unchanged from Stages 3/5/6), decomposed into `BSA_CDR3_side`, `BSA_target_side`, `BSA_pair_total`, with an exact consistency check (`BSA_CDR3_side + BSA_target_side == BSA_pair_total`) verified on every one of the 533 structures (`nonadj_bsa_consistency_ok` / `full_bsa_consistency_ok`, both 100% True). Always referred to as "pairwise mutual buried surface area", never "the BSA" alone, throughout code and output columns.

**10-structure BSA validation gate:** `results/audits/stage7_bsa_validation.csv` — PASS (10/10 real structures, all consistent, all positive).

## F. Covalent-adjacency target-set policy

`results/audits/stage7_covalent_adjacent_policy.md` — derived from actual sequence order (not assumption): IMGT 118 is the CDR3's immediate covalent C-terminal neighbor (excluded from the primary target set); IMGT 103 is one residue removed from CDR3 (104 intervenes) and is retained.

- **NONADJACENT_PACKING_TARGET** (primary, 14 positions): 40, 42, 44, 49, 50, 52, 55, 66, 67, 68, 69, 103, 119, 120
- **FULL_SHIELDABLE_TARGET** (sensitivity, 15 positions): adds 118

## G. Arpeggio real-structure validation gate

`results/audits/stage7_arpeggio_validation.csv` — 10/10 real structures PASS, all `ARPEGGIO_SUCCESS_WITH_INTERACTIONS`, 4218-4983 total interaction records each.

## H. Full 533-VHH cohort run

`data/interim/stage7/cohort_status.jsonl` — **533/533 `ARPEGGIO_SUCCESS_WITH_INTERACTIONS`**, zero failures, zero `TARGET_NOT_EVALUABLE`. No silent zeros anywhere: every instance's status is explicitly recorded.

## I. Arpeggio schema / interaction vocabulary

`results/audits/stage7_arpeggio_smoke_test.md` §12-13 — full real, source-confirmed atom-atom contact vocabulary (15 possible labels; 11 observed across the smoke-test sample). Ring/plane-based interactions (`plane-plane`/`atom-plane`, CARBONPI/CATIONPI/etc.) confirmed unconditionally computed by Arpeggio but **observed zero times across the entire 533-structure full-cohort parse** (`stage7_build_tables.py` output: `ring_type_seen = set()`) — a genuine empirical finding for this application context (whole native VHH domain, no antigen/ligand selection), not a tool limitation, and the parser handles these types generically should they ever appear.

## J. Output tables (all six required, in `results/tables/`)

1. `stage7_interaction_records_long.csv` — 89,725 atom-level CDR3↔target crossing records across all 533 instances (multi-label `contact_labels`, `is_covalent`/`is_proximal_only`/`is_specific_interaction`/`is_sequential_adjacent`/`include_primary_noncovalent` flags).
2. `stage7_per_vhh_packing_summary.csv` — 533 rows, BSA (both target sets) + specific-interaction atom-level and unique-residue-pair counts (parallel reporting, avoiding atom-count inflation — verified `n_unique_cdr3_target_respairs <= n_specific_interactions_atomlevel` for all rows) + normalized per-100Å²-BSA density.
3. `stage7_position_pair_interaction_map.csv` — 91 CDR3×target position pairs observed, with cohort-wide interaction frequency.
4. `stage7_target_position_summary.csv` / `stage7_cdr3_position_summary.csv` — per-position involvement summaries (14 target positions, 13 CDR3 positions ever involved in a specific interaction — position 116/117 never observed as a specific-interaction partner in this cohort, tracked not hidden).
5. `stage7_stage4to7_combined.csv` — 533×93 combined table joining Stage 4 (interface chemistry), Stage 6 (geometry + shielding), and Stage 7 (BSA + packing) columns, ready for a future Stage 11 formal-inference stage.

## K. Interaction classification discipline

- `is_proximal_only` records never simultaneously marked `is_specific_interaction` (verified by test, `contact == ['proximal']` exactly means no other label present).
- `is_covalent` records excluded from all "specific interaction" counts. Discovered (not assumed) that Arpeggio's `covalent` label captures **both** the trivial 117-118 backbone peptide bond **and** genuine Cys-Cys disulfide bonds at non-adjacent positions (observed covalent CDR3-target crossings at target positions 50 and 55, i.e., framework cysteines forming disulfides with CDR3 — a real, biologically plausible camelid-VHH feature, not a pipeline artifact; original test assumption of pure 118-adjacency was found incorrect and corrected).
- `include_primary_noncovalent` = specific AND not covalent AND target position in NONADJACENT_PACKING_TARGET (i.e., excludes both covalent bonds and the geometrically trivial 118 boundary).

## L. Descriptive (non-inferential) correlations

`results/tables/stage7_exploratory_correlations.csv` (Pearson + Spearman, explicitly exploratory, no formal-inference claims):

| x | y | n | Pearson r | Spearman r |
|---|---|---|---|---|
| BSA_pair_total | Stage-6 shielding_fraction | 533 | 0.967 | 0.949 |
| BSA_pair_total | Stage-6 delta_sasa_shield | 533 | 0.981 | 0.972 |
| n_unique_cdr3_target_respairs | shielding_fraction | 533 | 0.902 | 0.871 |
| cdr3_sequence_length | BSA_pair_total | 533 | 0.689 | 0.671 |
| cdr3_rho (Gordon reach) | BSA_pair_total | 533 | -0.191 | -0.174 |
| cdr3_compactness (Gordon) | BSA_pair_total | 533 | 0.846 | 0.836 |
| FR2_CDR3_Dis (Zhang) | BSA_pair_total | 533 | -0.824 | -0.815 |
| H_intrinsic_binary_mean_KD | BSA_pair_total | 533 | 0.183 | 0.156 |
| H_native_full_interface_sasa_weighted_mean_KD | BSA_pair_total | 533 | -0.138 | -0.144 |

Gordon θ/φ **not** computed against BSA (remain `UNRESOLVED_METHOD_DEPENDENCY`, per Stage 6). No p-value in this table is treated as confirmatory; formal inference is explicitly deferred to Stage 11.

**Interpretation boundary maintained:** these are geometric/structural correlations. None of the language above states or implies CDR3 packing "stabilizes," "causes," or "prevents" anything.

## M. High-coverage sensitivity analysis

`results/tables/stage7_high_coverage_sensitivity.csv` — 531/533 structures pass the ≥0.90 dual-coverage threshold (CDR3 coordinate coverage AND shieldable-interface coordinate coverage). BSA-vs-shielding correlation is essentially unchanged (r=0.969 high-coverage vs r=0.967 full cohort) — the primary finding is not an artifact of low-coverage structures. The primary 533-structure cohort was not altered by this check.

## N. Manual validation

`results/audits/stage7_manual_validation.csv` — 10 lowest-shielding, 10 near-median, 10 highest-shielding, plus 3 BSA-vs-shielding residual outliers (33 unique instances after de-duplication). All `nonadj_status == OK`; BSA and unique-pair counts increase monotonically with the shielding category as expected; no anomalies found.

## O. Tests

`tests/test_stage7_cdr3_packing.py` — 30 tests covering target-set definitions, BSA computation (synthetic + consistency), CIF-writer identity round-trip (reversible mapping, `_chem_comp` presence, insertion-code preservation, multi-chain distinctness), and full output-table structural/content checks. **30/30 PASS.** Full repository suite: **205/205 PASS** (up from the 175-test pre-Stage-7 baseline, zero regressions).

## P. Figures (visualization gate: justified by strong pre-computed correlations)

`results/figures/stage7/`:
1. `stage7_fig1_packing_vs_shielding.png` — BSA_pair_total vs shielding_fraction scatter (r=0.97).
2. `stage7_fig2_position_pair_heatmap.png` — CDR3×target position-pair interaction-frequency heatmap.
3. `stage7_fig3_composition_by_quartile.png` — interaction-type composition by predeclared shielding-fraction quartile (descriptive grouping only, no inferential test between quartiles).

## Q. 3D interactive viewer

`results/structures/interactive/stage7_cdr3_packing_pdb_00004grw_pdb_00004grw-E.html` — deterministic representative `pdb_00004grw-E` (simultaneously closest to cohort medians of both shielding_fraction and BSA_pair_total, not cherry-picked; selection recorded in `results/structures/stage7_representative_selection.csv`). Distinguishes: grey (VHH backbone), teal (CDR3), orange sticks (full NONADJACENT_PACKING_TARGET set regardless of detected interaction), red sticks (positions actually involved in a specific noncovalent interaction in this structure) — explicitly separating mere spatial proximity from typed noncovalent contact. `viewer_index.csv` updated.

## R. Terminological/conceptual discipline maintained throughout

- "Pairwise mutual buried surface area" (BSA) never conflated with Stage 6 "direct shielding" (a distinct, VH-vs-VH-without-CDR3 SASA counterfactual) or with Arpeggio "specific interactions" (atom-pair contact typing) — all three tracked as separate, distinctly-named columns throughout.
- "Proximal" contacts never counted as "specific interactions."
- No count is ever described as "energy."
- No stage-11 (formal inference), TAP-PSH/PPC/PNC, APBS, aggregation/developability, or MD metric was computed.
- All multi-state residue provenance (from Stage 4's PRESENT_OBSERVED/PRESENT_UNRESOLVED/ABSENT_CONFIRMED/AMBIGUOUS_UNRESOLVED model) implicitly respected — atoms only exist in this pipeline if `extract_variable_domain_atoms` found real observed coordinates; no 0 stand-in for missing data anywhere in the BSA or Arpeggio outputs (`TARGET_NOT_EVALUABLE` used explicitly where applicable, though it occurred 0/533 times in the real cohort).
- Case-sensitive, multi-character chain IDs and IMGT insertion codes preserved end-to-end with a full reversible audit trail (verified directly, not assumed, per §A).

---

## STAGE 7 STATUS: PASS
## VISUALIZATION DECISION: REQUIRED (fulfilled — 3 figures + 1 interactive 3D viewer, all justified by pre-computed strong correlations, none speculative)

STOP before Stage 8. No Stage-8+ metrics (TAP-PSH, PPC/PNC, APBS, aggregation/developability, MD) were computed in this stage.
