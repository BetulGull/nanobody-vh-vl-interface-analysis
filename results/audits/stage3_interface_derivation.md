# Stage 3 — Empirical VH-VL Interface Derivation by ΔSASA

**Generated:** 2026-09-01
**Scope:** Derive, purely geometrically, which VH surface positions are buried by VL across experimentally solved paired antibody structures. No VHH projection, no hydrophobicity, no CDR3 geometry, no shielding claims, no mechanistic interpretation.

---

## A. Stage-2 Paired Input Count

`data/processed/stage2_paired_primary.parquet`: **2879 rows** (verified from file, matches frozen Stage-2 authoritative count). All rows confirmed `construct_conventional == True`; zero SCFV/diabody/CODV/DVD/SCDB rows present.

## B. Number Successfully Processed by FreeSASA

**2879/2879** (100%) after fixing one implementation bug discovered mid-run (see below). Initial run: 2860/2879 passed, 19 failed.

**Bug found and fixed:** FreeSASA's `Structure.addAtom()` silently truncates any `chainLabel` to its first character (confirmed empirically — a real chain named `'HHH'` is echoed back by `result.residueAreas()` keyed under `'H'`). All 19 initial failures were PDB entries using extended/multi-character chain naming (e.g. `HHH`/`LLL`/`AAA` in `pdb_00006ss2`, `pdb_00007pa9`, etc.) where this silently broke the per-residue dictionary lookup. **Fix:** `stage3_sasa.py` now assigns each distinct real chain a guaranteed-unique single-character surrogate label per `calc_sasa()` call (at most 2 chains per call: VH alone, or VH+VL) and translates FreeSASA's output back to the true chain ID — never depending on FreeSASA's truncated echo. Full cohort was rerun after the fix: 2879/2879 pass.

## C. Number Failed/Reviewed and Reasons

**0** after the fix (`results/audits/stage3_sasa_failures.csv` is empty).

## D. Number of Primary Nonredundant Redundancy Groups

**1429** redundancy groups from 2879 paired structures. Redundancy key: `SABDAB_ID` (verified sequence-consistent for 706/718 multi-instance groups using actual IMGT-numbered VH/VL sequences); **12 SABDAB_ID groups were found sequence-inconsistent** and were NOT silently collapsed — those members were instead re-grouped by an exact VH+IMGT-sequence + VL+IMGT-sequence pairing key (`results/audits/stage3_sabdab_id_sequence_consistency.csv`), which is why 1429 > 1403 raw unique `SABDAB_ID` values.

## E. Representative Selection Rule

Deterministic, non-ΔSASA-based hierarchy applied per redundancy group: (1) highest VH variable-domain coordinate coverage, (2) highest VL coverage, (3) lowest resolution, (4) lowest R-free, (5) lowest R-factor, (6) lexical `INSTANCE` tie-break. `results/audits/stage3_representative_selection.csv`; outputs `data/processed/stage3_paired_representatives.{parquet,csv}` (1429 rows).

## F. Number of Exact VH IMGT Positions Observed

**214** distinct (imgt_position, insertion_code) combinations across the cohort (IMGT 1–128 base positions plus insertion-coded CDR3 variants).

## G. Primary Binary Empirical-Interface Position Count

**28** positions meet the primary definition (per-structure ΔSASA ≥ 5 Å², interface frequency ≥ 0.25, observation fraction ≥ 0.25) out of 214 observed VH positions.

## H. Primary Interface Position List

Full list: `results/tables/stage3_primary_interface_positions.csv` (28 rows, columns: imgt_position, region, n_observed, observation_fraction, interface_frequency_5, mean_delta_sasa, median_delta_sasa, vhh_tetrad_position).

## I. Tetrad Results (42/49/50/52)

| Position | n_observed | mean ΔSASA | median ΔSASA | freq≥1 | freq≥5 | freq≥10 | In primary binary interface |
|---|---|---|---|---|---|---|---|
| 42 | 1429 | 6.34 | 5.57 | 0.995 | 0.634 | 0.083 | True |
| 49 | 1429 | 29.58 | 27.17 | 0.996 | 0.965 | 0.937 | True |
| 50 | 1429 | 97.10 | 98.18 | 1.000 | 1.000 | 1.000 | True |
| 52 | 1429 | 78.69 | 80.61 | 0.999 | 0.999 | 0.999 | True |

Descriptive only, per spec §20 — annotation, not a predefined-interface assumption; the tetrad was not forced into the definition and was found to already satisfy it empirically at these positions.

## J. Region-Level Localization Summary

| Region | n positions | n primary consensus | mean freq≥5 | median freq≥5 | mean ΔSASA |
|---|---|---|---|---|---|
| FR1 | 25 | 0 | 0.001 | 0.000 | 0.01 |
| CDR1 | 17 | 0 | 0.005 | 0.000 | 0.09 |
| FR2 | 17 | 7 | 0.343 | 0.064 | 15.89 |
| CDR2 | 22 | 0 | 0.008 | 0.000 | 0.15 |
| FR3 | 58 | 5 | 0.071 | 0.000 | 1.79 |
| CDR3 | 63 | 13 | 0.158 | 0.000 | 7.29 |
| FR4 | 12 | 3 | 0.211 | 0.000 | 7.88 |

Descriptive localization only (per spec §21) — no mechanistic interpretation offered.

## K. Distribution of Total VH ΔSASA Due to VL (2879 structures)

mean 807.9 Å², median 801.6 Å², SD 89.8, min 567.2, 25th pct 745.0, 75th pct 861.2, max 1219.0 Å². `n_interface_residues_5` per structure: mean 19.4, median 19, range 12–30. This is geometric surface burial, never labeled as binding energy or interface energy (§22).

## L. Negative ΔSASA Audit

Of 348,144 total residue observations: **2 negative values**, minimum = −7.1×10⁻¹⁵ Å² (machine-epsilon-scale floating-point noise). Zero values below −0.1, −0.5, or −1.0 Å². No systematic negative burial detected — no investigation stop triggered. `results/audits/stage3_negative_delta_sasa.csv`.

## M. FreeSASA Validation Result

`results/audits/stage3_sasa_validation.csv` — **PASS** on 5 deterministic real structures (FAB with VH+CH1 same chain, paired FV, insertion-coded, missing-coordinate, cross-resolution): per-residue sums exactly match total SASA (Å² tolerance <1e-6), repeated calculations reproduce identically (tolerance <1e-9). CLI-vs-API comparison deliberately not attempted — documented rationale: a CLI round-trip would validate a PDB-text writer, not the SASA calculation itself; internal reproducibility + the geometric translated-VL control are more direct validations of the actual computation.

FreeSASA config used throughout: freesasa-python 2.2.1, algorithm=LeeRichards (explicit), probe radius=1.4 Å (explicit), classifier=default (NACCESS-derived; confirmed to natively recognize MSE — no MSE occurred within any antibody H/L chain across the full 2879-structure cohort after exhaustive search, only on excluded antigen chains).

## N. Translated-VL Geometric-Control Result

ΔSASA from a VL copy translated 100,000 Å away ≈ **0** for all 5 validation structures (largest deviation: 5.5×10⁻¹² Å², i.e. exactly floating-point noise). Confirms measured burial arises from genuine spatial occlusion by VL, not calculation bookkeeping artifacts. Same control re-verified in the automated test suite (`tests/test_stage3_interface.py::TestTranslatedVLControl`).

## O. Per-Residue Sum vs Total-SASA Consistency

Verified exactly (tolerance <1e-6 Å²) on all 5 validation structures: `sum_i[SASA_i(VH alone) − SASA_i(VH+VL)]` equals `total_VH_alone − total_VH_within_VH+VL_calc` in every case (`total_delta_consistent = True` for all rows in `stage3_sasa_validation.csv`).

## P. Threshold-Sensitivity Results

`results/tables/stage3_interface_threshold_sensitivity.csv` — full 3×3×3 grid (27 combinations). Consensus position count ranges from 14 (ΔSASA≥10, freq≥0.50, obs≥0.50 — strictest) to 39 (ΔSASA≥1, freq≥0.10, obs≥0.10 — most permissive). Jaccard similarity to the primary definition ranges 0.50–0.93 across the grid; all four tetrad positions are included at ΔSASA thresholds 1 and 5 across every frequency/observation combination, dropping to 3/4 at ΔSASA≥10 with the strictest frequency (position 42's freq≥10 = 0.083, below 0.10). The primary definition (5 Å², 0.25, 0.25) was fixed before generating this table, not chosen from it.

## Q. Primary vs. All-Instance Redundancy Sensitivity

`results/tables/stage3_redundancy_sensitivity.csv`: Pearson r = **0.998**, Spearman r = **0.988** (interface_frequency_5, primary vs. all-instance, n=214 positions compared); mean absolute frequency difference = 0.004, max = 0.167; **binary-interface Jaccard = 1.000** (identical 28-position sets from both maps). Repeated crystal structures do **not** materially distort the primary map.

## R. pytest Result

**89/89 passed** (full repository suite, including 25 new Stage-3-specific tests covering all 23 required areas from spec §36).

## S. Visualization Paths

- Empirical interface profile (MAIN_CANDIDATE): `results/figures/stage3/stage3_interface_profile.png` / `.svg` (4770×2066, 300 dpi)
- Threshold sensitivity (SUPPLEMENTARY_SENSITIVITY): `results/figures/stage3/stage3_threshold_sensitivity.png` / `.svg` (4471×1484, 300 dpi)
- Interactive 3D VH-VL interface viewer: `results/structures/interactive/stage3_empirical_vh_vl_interface_pdb_00009piz_pdb_00009piz-H-I.html` — representative structure selected deterministically as closest to cohort medians of `total_vh_delta_sasa_due_to_vl` (799.3 vs. median 801.6) and `n_interface_residues_5` (19 vs. median 19), resolution as tie-break. The viewer clearly distinguishes structure-specific burial values (reported in its header text) from the population-level `interface_frequency_5` coloring projected onto that one representative VH.
- Figure manifest: `results/figures/stage3/figure_manifest.csv`; README: `results/figures/stage3/README.md`

## T. Confirmations

- **Stage-2 data unchanged:** verified by test (`stage2_paired_primary.parquet` = 2879 rows, `stage2_vhh_strict_primary.parquet` = 1217 rows, both unchanged).
- **No VHH projection performed.**
- **No hydrophobicity analysis performed.**
- **No CDR3 geometry calculated.**
- **No shielding claims made.**
- **No Arpeggio mechanistic analysis performed.**
- **No APBS performed.**
- **No aggregation/developability analysis performed.**
- **No MD performed.**

---

## Scientific Limitations (recorded per spec §44)

1. ΔSASA is geometric surface burial, not binding free energy.
2. These are experimentally solved paired VH–VL conformations only.
3. The derived map describes how frequently VH positions are buried by VL across this structural dataset — a population-level empirical frequency, not a per-structure guarantee.
4. It does not prove the same positions would behave identically in an isolated VH (e.g. VHH).
5. It does not yet explain VHH solubility or stability.
6. It does not yet show that CDR3 shields the interface.
7. It does not yet measure hydrophobicity.
8. Antigen was deliberately excluded from the SASA context throughout, to isolate VH–VL burial specifically from any antigen-driven occlusion.

---

## STAGE 3 STATUS: PASS

**VISUALIZATION DECISION: REQUIRED** — satisfied (MAIN_CANDIDATE interface profile + SUPPLEMENTARY_SENSITIVITY threshold heatmap + interactive 3D viewer, all generated and validated).
