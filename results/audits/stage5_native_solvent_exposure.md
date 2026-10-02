# Stage 5 — Native 3D Solvent Exposure and Former-Interface Hydrophobic Burden

**Generated:** 2026-09-02
**Scope:** H_native only (structural solvent exposure). No CDR3 removal, no shielding metrics, no patches, no electrostatics.

---

## A. Primary VHH Representative N
**533** (`data/processed/stage4_vhh_sequence_representatives.parquet`, frozen, verified from file).

## B. Conventional VH Representative N
**1429** (`data/processed/stage3_paired_representatives.parquet`, frozen, reused directly).

## C. VHH SASA Validation Result
**PASS** (`results/audits/stage5_vhh_sasa_validation.csv`, 11 diverse real structures: lowest/highest coverage, longest/shortest CDR3-insertion depth, two unusual multi-character chain IDs (`FFF`, `B1`), resolution extremes 1.02–3.0 Å, multi-antigen-chain complex). All checks passed: antigen absent, constant domain excluded, CDR3 fully intact, coordinates identical across repeated extraction calls, SASA reproducible.

## D. Stage-3 VH SASA Reproduction Result
**PASS** (`results/audits/stage5_stage3_vh_sasa_reproduction.csv`, 12 structures). Total and per-residue differences all at floating-point noise level (~1e-12 Å²); 0 mismatches across all compared residues. Frozen Stage-3 values reused directly for production (not recomputed for the full cohort), consistent with instruction.

## E. Primary-Interface Coordinate Coverage Distributions
`results/tables/stage5_interface_structural_coverage.csv`. Zero sequences (any structural state) with `interface_coordinate_coverage_given_present < 0.5` — the empirical interface is well-resolved structurally across the cohort.

## F. High-Coverage Sensitivity Result
`results/tables/stage5_high_coverage_sensitivity.csv`. 532/533 VHH and 1429/1429 VH already meet the ≥0.90 coordinate-coverage bar. `interface_sasa_weighted_mean_KD` shifts by ≤0.001 for VHH (−1.8839→−1.8840) and is unchanged for VH (already 100% at/above threshold). **The primary result is not a coverage artifact.**

## G. VHH_NATIVE interface_total_sasa Descriptive Statistics
See `results/tables/stage5_three_state_surface_comparison.csv` for full N/mean/median/SD/Q25/Q75.

## H. VH_WITHOUT_VL_GEOMETRIC Equivalent
Same table, `structural_state == VH_WITHOUT_VL_GEOMETRIC`.

## I. VH_WITH_VL_PAIRED Equivalent
Same table, `structural_state == VH_WITH_VL_PAIRED`.

## J. VHH vs VH-Without-VL: interface_sasa_weighted_mean_KD
| Structural state | N | mean | median | SD |
|---|---|---|---|---|
| VHH_NATIVE | 533 | **−1.884** | −1.942 | 0.478 |
| VH_WITHOUT_VL_GEOMETRIC | 1429 | **−1.057** | −1.086 | 0.478 |

VHH's native solvent-exposed surface at the empirical former interface is substantially more hydrophilic (lower/more-negative KD) than the VL-removed conventional VH geometric reference.

## K. VHH vs VH-Without-VL: interface_positiveKD_sasa_fraction
VHH mean **0.139**, VH-without-VL mean **0.201** — a smaller fraction of VHH's exposed interface area is carried by positive-KD (hydrophobic) residues.

## L. Weighted-Interface Results (freq5-weighted, sensitivity)
Same direction and comparable magnitude to the binary descriptor (see `stage5_sequence_surface_chemistry.parquet` columns `weighted_interface_sasa_mean_KD_freq5`), confirming the result is not an artifact of the binary threshold choice.

## M. Intrinsic-vs-Structural Correlations
`results/tables/stage5_intrinsic_vs_structural_chemistry.csv`:

| Group | Structural state | N | Pearson r (binary) | Spearman r (binary) |
|---|---|---|---|---|
| VHH | VHH_NATIVE | 533 | **0.722** | 0.702 |
| VH | VH_WITHOUT_VL_GEOMETRIC | 1429 | **0.818** | 0.804 |

## N. Does Stage-4 Chemistry Separation Persist After Structural Exposure?
**It amplifies.** Intrinsic (Stage 4) group means: VHH −1.063, VH −0.696, separation = **0.367**. Structural (Stage 5) group means: VHH −1.884, VH −1.057, separation = **0.827** — more than double. Sequence-level chemistry correlates strongly with structural exposure in both groups (§M), and actual solvent exposure widens rather than narrows the VHH-vs-VH difference. No causal claim made.

## O. Position-Level Largest Structural Chemistry Differences
Top by |Δmean_sasa×KD| within the 28-position interface: position 50 (FR2, Δ=−766.7), 49 (FR2, Δ=−352.5), 115 (CDR3, Δ=−182.8), 52 (FR2, Δ=+149.0), 120 (FR4, Δ=−77.4) (`data/processed/stage5_position_surface_chemistry.csv`).

## P. Interface-vs-Non-Interface-Framework Control
`results/tables/stage5_interface_vs_noninterface_framework.csv` (descriptive only, no inferential test):

| Category | Structural state | SASA-weighted mean KD | positiveKD fraction |
|---|---|---|---|
| EMPIRICAL_INTERFACE | VHH_NATIVE | −1.890 | 0.138 |
| CONTROL_NONINTERFACE_FRAMEWORK | VHH_NATIVE | −1.536 | 0.159 |
| EMPIRICAL_INTERFACE | VH_WITHOUT_VL_GEOMETRIC | −1.064 | 0.201 |
| CONTROL_NONINTERFACE_FRAMEWORK | VH_WITHOUT_VL_GEOMETRIC | −1.429 | 0.161 |

Notable asymmetry, reported descriptively: VHH's empirical interface is *more* hydrophilic than its own framework control (−1.89 vs −1.54), while VH's empirical interface is *less* hydrophilic than its own framework control (−1.06 vs −1.43) — the two groups show opposite-direction interface-vs-framework relationships. No mechanistic interpretation offered.

## Q. Noncanonical-Residue Audit
`results/audits/stage5_noncanonical_residues.csv`: 28 observations within the primary interface (14 unique residue-positions × 2 VH structural states, all traced to the single VH structure `pdb_00008vxj-A-B` already identified in Stage 4 as having crystallographer-deposited `UNK` residues across its CDR3 core). No invented KD values; all excluded from KD-weighted denominators.

## R. Structural Manual-Validation Result
**PASS** (`results/audits/stage5_manual_validation.csv`, 10 VHH + 10 VH × 2 structural states = 30 rows). Zero instances of a zero-SASA value misattributed to an unresolved or absent residue.

## S. Full pytest Result
**144/144 passed** (full repository suite, including 25 new Stage-5-specific tests).

## T. Generated Figures
- MAIN_CANDIDATE: `results/figures/stage5/stage5_structural_chemistry_shift.png`/`.svg`
- MAIN_OR_SUPPLEMENTARY_MECHANISM: `results/figures/stage5/stage5_intrinsic_vs_structural.png`/`.svg`
- OPTIONAL_THREE_STATE: `results/figures/stage5/stage5_three_state_comparison.png`/`.svg` (generated — the ordering is non-obvious and materially aids interpretation)
- Manifest/README: `results/figures/stage5/figure_manifest.csv`, `results/figures/stage5/README.md`

## U. 3D Viewer Paths
- VHH_NATIVE: `results/structures/interactive/stage5_vhh_native_surface_pdb_00008ybp_pdb_00008ybp-B.html`
- VH_WITHOUT_VL_GEOMETRIC: `results/structures/interactive/stage5_vh_without_vl_surface_pdb_00006o9h_pdb_00006o9h-A-B.html`

Both selected deterministically (closest to cohort medians of `interface_sasa_weighted_mean_KD` and `interface_total_sasa`, ≥90% coordinate coverage, resolution tie-break — `results/structures/stage5_representative_selection.csv`). The VH viewer explicitly states VL is hidden only for the geometric reference, not experimentally removed.

## V. Limitations
1. SASA measures geometric solvent accessibility, not solubility.
2. SASA-weighted KD is a descriptive surface-chemistry metric, not free energy.
3. VHH structures are experimentally observed (often antigen-bound) conformations with antigen computationally removed for this analysis.
4. Conventional VH-without-VL is a same-coordinate geometric counterfactual, not an equilibrated isolated VH — never described otherwise.
5. Crystal packing may influence observed conformations/exposure.
6. Missing coordinates reduce structural coverage at some positions (never converted to zero SASA).
7. Stage 5 cannot and does not attribute any observed difference specifically to CDR3 — no CDR3-removal calculation was performed.

---

## Confirmations
- Stage-2 unchanged (1217 / 2879 verified). Stage-3 unchanged (214 positions / 28 binary interface verified). Stage-4 unchanged (1962 sequence-chemistry rows verified).
- No empirical-interface redefinition — the same frozen 28-position Stage-3 map used throughout.
- No CDR3 removal, no shielding metric, no ΔSASA_shield, no H_unshielded, no ΔH_CDR3 (verified by source-code grep in the test suite).
- No CDR3 geometry, no BSA, no Arpeggio, no TAP PSH, no PPC/PNC, no APBS, no aggregation/developability, no MD.

---

## STAGE 5 STATUS: PASS

**VISUALIZATION DECISION: REQUIRED** — satisfied (1 MAIN_CANDIDATE + 1 MAIN_OR_SUPPLEMENTARY_MECHANISM + 1 OPTIONAL_THREE_STATE figure + 2 interactive 3D viewers, all generated and validated).
