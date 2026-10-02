# Stage 6 — CDR3 Geometry, Direct Former-Interface Shielding, and H_unshielded Decomposition

**Generated:** 2026-09-02 (status corrected 2026-09-02, post-hoc literature-dependency re-check)
**Scope:** Same-coordinate no-CDR3 SASA counterfactual, direct shielding, Gordon ρ/compactness, Zhang FR2_CDR3_Dis. No atomistic interaction analysis (Stage 7).

> **STATUS CORRECTION:** an earlier version of this report incorrectly concluded "STAGE 6 STATUS: PASS" despite Gordon θ/φ being unresolved. Per the project's predeclared status policy (an unresolved published-method dependency must not be silently absorbed into an overall PASS), the corrected authoritative status is:
>
> **STAGE 6 CORE SHIELDING: PASS**
> **GORDON RHO/COMPACTNESS: PASS**
> **ZHANG FR2_CDR3_DIS: PASS**
> **GORDON THETA/PHI: UNRESOLVED METHOD DEPENDENCY**
> **OVERALL STAGE 6 STATUS: PARTIAL**
>
> A focused re-search for new Gordon source material (`results/audits/stage6_literature_dependency_audit.md`) found nothing new: Text S2 is confirmed hosted only at the publisher's external website, not present anywhere in this repository; the paper's code pointer (`github.com/oxpig`) names an organization, not a specific repository, and no such repository is cloned locally; `external/TNP` (already inspected) contains no θ/φ implementation. No substitute was fabricated. All frozen shielding/ρ/compactness/FR2_CDR3_Dis values were verified byte-identical before and after this re-check (§ "Frozen-value verification" below); 175/175 tests still pass; no Stage 2–5 file was touched.

---

## A. Primary VHH Representative N
**533** (`data/processed/stage4_vhh_sequence_representatives.parquet`, frozen, verified from file).

## B. Full Empirical Interface Position Count
**28** (Stage-3 primary consensus positions, frozen, unchanged).

## C. CDR3 Positions Within Full Empirical Interface
**13** (IMGT 107,108,109,110,111,111A,112,112A,112B,113,114,115,116 — all in CDR3/FR3 boundary region).

## D. Shieldable-Interface Position Count and List
**15** positions (28 − 13). All in FR2/FR3/FR4: IMGT 40,42,44,49,50,52,55 (FR2), 66,67,68,69,103 (FR3), 118,119,120 (FR4). Full list: derivable from `data/processed/stage6_position_shielding_summary.csv`.

## E. CDR3-Removal Validation Result
**PASS** (`results/audits/stage6_cdr3_removal_validation.csv`, 10 diverse structures). All checks passed: zero CDR3 atoms remain in the no-CDR3 state, zero non-CDR3 atoms removed, retained coordinates exactly identical, both SASA states reproducible.

## F. Translated-CDR3 Control Result
**PASS** (`results/audits/stage6_translated_cdr3_control.csv`, 6 structures). Max absolute difference between no-CDR3 removal and far-translation (100,000 Å) at shieldable positions: **0.0 Å²** exactly — confirms the shielding signal is pure spatial occlusion, not a bookkeeping artifact.

## G. Successfully Evaluated Shielding N
**533/533** (100%). Zero processing errors across the full cohort.

## H. Shielding Structural Coverage Distribution
`shieldable_interface_coordinate_coverage`: mean 0.999, median 1.0, min 0.80 (`results/tables/stage6_shielding_distribution.csv` context; full distribution in `stage6_vhh_cdr3_shielding.parquet`).

## I. CDR3 Coordinate Coverage Distribution
**1.0 for all 533 structures** — expected, since Stage-2 QC already required complete H-CDR3 backbone for strict-primary cohort membership.

## J. delta_sasa_shield Descriptive Statistics
N=533, mean **245.1 Å²**, median 271.7, SD 97.9, Q25 143.9, Q75 319.4, min 69.2, max 553.4 (`results/tables/stage6_shielding_distribution.csv`).

## K. shielding_fraction Descriptive Statistics
N=533, mean **0.215**, median 0.236, SD 0.086, Q25 0.123, Q75 0.282, min 0.062, max 0.483.

## L. Fraction of VHHs at Predeclared Shielding Thresholds
| Threshold | Fraction ≥ threshold |
|---|---|
| ≥1 Å² | 100% |
| ≥5 Å² | 100% |
| ≥10 Å² | 100% |
| ≥25 Å² | 100% |
| ≥50 Å² | 100% |
| ≥100 Å² | 93.4% |
| ≥200 Å² | 66.0% |
| ≥300 Å² | 37.0% |
| ≥500 Å² | 0.4% |

**Every single VHH in this cohort shows at least 50 Å² of direct CDR3 shielding of the shieldable former-interface**; no VHHs with negligible/zero shielding were found (none hidden — this is the complete distribution).

## M. Negative Shielding Audit
**Zero** negative `delta_sasa_shield` values at either residue level (0/7987 evaluable observations) or structure level (0/533) (`results/audits/stage6_negative_shielding.csv`). No investigation-stop triggered.

## N. H_native_shieldable vs H_unshielded Results
Native and no-CDR3 SASA-weighted mean KD, positive-KD fraction, and KD area sum are all reported per-structure on the matched shieldable target (`data/processed/stage6_vhh_cdr3_shielding.parquet`, columns `native_sasa_weighted_mean_KD` / `noCDR3_sasa_weighted_mean_KD` etc.) — never computed by subtracting Stage-5's full-interface H_native (different residue set) from this stage's non-CDR3-target H_unshielded.

## O. delta_positiveKD_sasa_CDR3
Reported per-structure (`stage6_vhh_cdr3_shielding.parquet`, column `delta_positiveKD_sasa_CDR3`) — the increase in positive-KD-carried exposed area if CDR3 were geometrically absent.

## P. delta_sasa_weighted_mean_KD_CDR3
Reported per-structure (column `delta_sasa_weighted_mean_KD_CDR3`).

## Q. Weighted-Interface Shielding Results
`weighted_delta_sasa_shield_freq5` (using Stage-3 `interface_frequency_5` as weight) reported per-structure; correlates with the binary `delta_sasa_shield` in the same direction as the geometry relationships (§Y), confirming the result is not an artifact of the binary threshold.

## R. Position-Level Strongest Shielding Sites
| IMGT position | Region | Mean ΔSASA_shield |
|---|---|---|
| 118 | FR4 | **87.8 Å²** |
| 42 | FR2 | **42.3 Å²** |
| 52 | FR2 | 24.0 |
| 66 | FR3 | 27.6 |
| 55 | FR2 | 21.1 |
| 40 | FR2 | 18.1 |
| 50 | FR2 | 15.0 |

Position 118 (immediately C-terminal to CDR3) and position 42 (VHH tetrad, annotation only — not privileged in the definition) show by far the largest mean occlusion. Full table: `data/processed/stage6_position_shielding_summary.csv`.

## S. Gordon Method-Definition Verification Result
ρ and compactness: **resolved and implemented** (main paper text + corroborating author code from the same OPIG lab, `external/TNP`). θ and φ: **genuinely unresolved** — the paper explicitly defers the exact R3 coordinate-frame construction to a supplementary "Text S2" not present in this repository, and no corroborating code was found for either. **Not fabricated; reported as `NaN` for every structure.** Full detail: `results/audits/stage6_gordon_method_definition.md`.

## T. Gordon Geometry Evaluable N
**533/533** (100%) for ρ and compactness. **0/533** for θ and φ (not implemented, by design, not by failure).

## U. rho/theta/phi/compactness Distributions
ρ: N=533, mean 12.90 Å, median 12.82, SD 1.90, range 5.44–18.30. Compactness: N=533, mean 1.243, median 1.327, SD 0.314, range 0.557–3.124. θ/φ: not computed (see §S).

## V. Zhang Method-Definition Verification Result
**Fully resolved, unambiguous** — exact computational protocol recovered directly from the paper's Materials and Methods section (quoted verbatim in `results/audits/stage6_zhang_method_definition.md`): minimum, over CDR3 residues excluding the first and last two in true sequence order, of the Euclidean distance between residue centroids (position 42 vs. each internal CDR3 residue).

## W. FR2_CDR3_Dis Evaluable N
**533/533** (100%).

## X. CDR3 Length vs Shielding Relationship
`cdr3_sequence_length` vs `delta_sasa_shield`: Pearson r=**0.671**, Spearman r=0.667 (N=533). vs `shielding_fraction`: r=0.660. Positive, substantial, but **not the strongest predictor** among the geometry descriptors (see §Y) — longer CDR3 associates with more shielding, but compactness and FR2_CDR3_Dis explain more variance.

## Y. Geometry vs Shielding Relationships
| Relationship | Pearson r | Spearman r |
|---|---|---|
| compactness vs delta_sasa_shield | **0.851** | 0.844 |
| compactness vs shielding_fraction | 0.836 | 0.840 |
| FR2_CDR3_Dis vs delta_sasa_shield | **−0.823** | −0.797 |
| FR2_CDR3_Dis vs shielding_fraction | −0.814 | −0.787 |
| ρ vs delta_sasa_shield | −0.227 | −0.196 |
| ρ vs shielding_fraction | −0.229 | −0.204 |

Compactness (positive) and FR2_CDR3_Dis (negative) are by far the strongest correlates of shielding — a more compact, more closely-packed CDR3 (shorter distance to FR2) shields substantially more. θ/φ relationships could not be evaluated (not implemented). All correlations reported as exploratory (`results/tables/stage6_geometry_shielding_correlations.csv`); p-values labeled exploratory/unadjusted, no formal inference performed.

## Z. FR2_CDR3_Dis vs Shielding Relationship
Covered in §Y above (r=−0.82 to −0.81) — closer CDR3-to-FR2 distance is strongly associated with greater shielding, consistent with (but not proof of) a "bent-down" CDR3 conformation occluding the former interface more.

## AA. H_intrinsic vs Shielding Relationship
`H_intrinsic_binary_mean_KD` vs `delta_sasa_shield`: r=**0.176** (weak positive, N=533). vs `shielding_fraction`: r=0.186. `H_native_shieldable` (this stage's own native chemistry on the matched target) vs shielding: r=**−0.34 to −0.39** (moderate negative — expected, since shielding directly reduces exposed area, and thus mechanically reduces exposed hydrophobicity at those same positions). No causal claim made; both are measured in the same native structure. Full detail: `results/tables/stage6_intrinsic_shielding_correlations.csv`.

## AB. High-Coverage Sensitivity Result
531/533 structures already meet ≥0.90 coverage on both shieldable-interface and CDR3 coordinates. `delta_sasa_shield` mean shifts by <0.5% (245.08→244.97); `shielding_fraction` mean shifts by <0.2% (0.2146→0.2142). **Primary result is not a coverage artifact.**

## AC. Manual-Validation Result
**PASS** (`results/audits/stage6_manual_validation.csv`, 10 low + 10 near-median + 10 high-shielding VHHs, deterministically selected by sorting on `shielding_fraction`). Zero geometry failures found in the full cohort (533/533 valid for both Gordon ρ/compactness and Zhang FR2_CDR3_Dis) — no "rare failures" existed to include as a category.

## AD. Full Repository pytest Result
**175/175 passed** (including 31 new Stage-6-specific tests).

## AE. Generated Figures
- MAIN_CANDIDATE: `results/figures/stage6/stage6_shielding_distribution.png`/`.svg`
- MAIN_CANDIDATE_OR_MECHANISM: `results/figures/stage6/stage6_geometry_vs_shielding.png`/`.svg`
- SUPPLEMENTARY_SPATIAL_MECHANISM: `results/figures/stage6/stage6_position_shielding.png`/`.svg`
- Manifest/README: `results/figures/stage6/figure_manifest.csv`, `results/figures/stage6/README.md`

## AF. 3D Viewer Path
`results/structures/interactive/stage6_cdr3_shielding_pdb_00007z1d_pdb_00007z1d-FFF.html` — deterministic representative (closest to cohort medians of `shielding_fraction` and `delta_sasa_shield`, ≥90% coverage, resolution tie-break; selection recorded in `results/structures/stage6_representative_selection.csv`). Tabbed viewer distinguishes native VHH (CDR3 visible, teal) from the same-coordinate no-CDR3 SASA counterfactual (CDR3 hidden), explicitly labeled as a geometric reference, not a physical mutant.

## AG. Limitations
1. The no-CDR3 state is a geometric SASA counterfactual, not a folded deletion mutant.
2. Structures are bound-state experimental conformations (antigen computationally removed for this analysis, as in Stage 5).
3. Crystal packing may influence observed CDR3 geometry.
4. Missing coordinates reduce geometry/shielding evaluability (though coverage was ≈100% in this cohort — see §H, §I).
5. Correlations (§X-AA) do not establish causality.
6. Stage 6 does not identify the atomistic interactions (H-bonds, hydrophobic contacts, aromatic/ionic contacts) that hold CDR3 near the former-interface surface — that is Stage 7.
7. Gordon θ and φ could not be reproduced (Text S2 unavailable) and are reported as absent, not approximated.

---

## Confirmations
- Stage-2 unchanged (1217 verified). Stage-3 unchanged (214 positions / 28 binary interface verified). Stage-4 unchanged (533 VHH representatives verified). Stage-5 unchanged (3391 sequence-surface-chemistry rows verified).
- Empirical VH-VL interface (Stage 3) never redefined — the 28-position map is the sole source of interface membership throughout.
- No physical no-CDR3 protein modeled; no coordinate relaxation, minimization, repacking, or MD anywhere.
- No atomistic interaction analysis (H-bonds, hydrophobic/aromatic/ionic contacts, van der Waals) performed — verified by source-code grep in the test suite (no Arpeggio import anywhere in Stage-6 code).
- No APBS, no aggregation/developability metrics, no CDR3-interface BSA as a packing metric.

---

## Frozen-value verification (post literature-dependency re-check)

Re-verified byte-identical against the values in this same report prior to the re-check:

| Descriptor | Value | Unchanged |
|---|---|---|
| `delta_sasa_shield` mean | 245.081 Å² | ✓ |
| `shielding_fraction` mean | 0.2146 | ✓ |
| `cdr3_rho` mean | 12.901 Å | ✓ |
| `cdr3_compactness` mean | 1.2433 | ✓ |
| `FR2_CDR3_Dis` mean | 7.951 Å | ✓ |
| N (VHH representatives) | 533 | ✓ |
| Stage 2 (`stage2_vhh_strict_primary`) | 1217 | ✓ unchanged |
| Stage 3 (primary binary interface) | 28 | ✓ unchanged |
| Stage 4 (`stage4_vhh_sequence_representatives`) | 533 | ✓ unchanged |
| Stage 5 (`stage5_sequence_surface_chemistry`) | 3391 rows | ✓ unchanged |
| Full repository pytest | 175/175 | ✓ unchanged (no code modified — documentation-only update) |

No scientific row was recomputed or patched. No Stage 2–5 source file was touched.

---

## STAGE 6 CORE SHIELDING: PASS
## GORDON RHO/COMPACTNESS: PASS
## ZHANG FR2_CDR3_DIS: PASS
## GORDON THETA/PHI: UNRESOLVED METHOD DEPENDENCY
## OVERALL STAGE 6 STATUS: PARTIAL

No new Gordon source material was found in `references/`, `external/`, `docs/`, or the project root beyond what was already inspected (`results/audits/stage6_literature_dependency_audit.md`). Text S2 — the paper's own stated location for the exact θ/φ definitions — is confirmed hosted only at the publisher's external website and is not present in this repository. No approximate substitute was implemented under the Gordon name. All other required Stage-6 scope — shielding, ρ, compactness, FR2_CDR3_Dis, all correlational analyses, validation gates (CDR3-removal, translated-CDR3 control, rigid-body invariance), figures, 3D viewer, manual validation, and the full test suite — is complete, validated, and passing.

**VISUALIZATION DECISION: REQUIRED** — satisfied (1 MAIN_CANDIDATE + 1 MAIN_CANDIDATE_OR_MECHANISM + 1 SUPPLEMENTARY_SPATIAL_MECHANISM figure + 1 interactive 3D native-vs-no-CDR3 viewer, all generated and validated).
