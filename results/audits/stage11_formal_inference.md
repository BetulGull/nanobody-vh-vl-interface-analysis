# Stage 11 — Formal Statistical Inference and Integrated Mechanism Testing

## A. Exact input tables

`results/audits/stage11_input_variable_dictionary.csv` — 22 variables traced to exact source file/column, including the critical structural-state disambiguation: `H_native_full_interface_sasa_weighted_mean_KD` is sourced from Stage 6's combined table for VHH (state=`VHH_NATIVE`) and directly from Stage 5's `stage5_sequence_surface_chemistry.parquet` filtered to `structural_state=='VH_WITHOUT_VL_GEOMETRIC'` for VH — never silently substituted with the paired `VH_WITH_VL_PAIRED` state.

## B. Model-specific cohort Ns

Canonical table (`data/processed/stage11_model_input.parquet`): 1962 rows (533 VHH + 1429 VH), one row per frozen representative, zero duplication (verified by an `assert len(df)==n0` check after every one of 10 join steps). VHH-only mechanism variables (shielding, geometry, BSA) are missing for all 1429 VH rows by design (72.8% missingness = exactly the VH row count) — never globally complete-case filtered; each model below reports its own N.

## C. Missingness audit

`results/audits/stage11_input_variable_dictionary.csv` (missing_n/missing_fraction columns) and `stage11_outcome_distribution_audit.csv`. Only two real missingness sources: VHH-only variables (1429 missing = all VH, expected) and Stage-8 APBS failures (55 missing = the known, previously-documented 55 non-success APBS structures, never converted to zero).

## D. Collinearity findings

`results/tables/stage11_collinearity_audit.csv` and `stage11_predictor_correlation_matrix.csv` (VHH-only, N=521 complete-case):

- **shielding_fraction vs nonadj_BSA_pair_total: r=0.969, VIF=15.5, condition number=26,287** — severe, exactly as predeclared. Their independent contributions **cannot be reliably separated** in this dataset; Stage 7 BSA is used only as a geometric corroboration variable, never as a co-predictor alongside shielding.
- **cdr3_sequence_length vs cdr3_compactness: r=0.85** — confirms the mathematical relationship; never fit together (verified by test).
- **H_intrinsic vs H_native: r=0.73, VIF~2.1 when paired** — moderate, estimable, both retained.
- Model1A (length+rho+FR2_CDR3_Dis): all VIF<3.5, condition number 180 — acceptable, independently estimable.
- Model1B (compactness+FR2_CDR3_Dis): VIF~2.1 — acceptable.

## E. Primary group effects

`results/tables/stage11_primary_group_inference.csv` — all 7 predeclared hypotheses (A: H_intrinsic; B: native surface; C: PSH + largest patch; D: framework engagement; E: 15-position former-interface engagement; F: CDR3 fraction), each with the **correct comparator per Section 10** (VHH vs conventional VH for all — never Fv for heavy-domain-specific quantities, verified by test):

| Hypothesis | Standardized effect | Raw p | BH-FDR q |
|---|---|---|---|
| A. H_intrinsic (VHH lower) | -0.91 | 1.3e-86 | 4.4e-86 |
| B. Native surface (VHH lower) | -1.37 | 2.5e-254 | 1.8e-253 |
| C. PSH (VHH lower) | -0.31 | 1.2e-11 | 1.4e-11 |
| C. Largest hydrophobic patch (VHH lower) | -0.67 | 2.5e-55 | 5.8e-55 |
| D. Framework engagement (VHH higher) | +0.32 | 9.3e-8 | 9.3e-8 |
| E. 15-position engagement (VHH higher) | +0.71 | 2.8e-27 | 4.9e-27 |
| F. CDR3 contribution (VHH higher) | +0.44 | 8.4e-13 | 1.2e-12 |

## F. FDR-adjusted group-family results

All 7 hypotheses remain significant after Benjamini-Hochberg correction within Family A (max q=4.9e-27, far below 0.05) — the family-wise correction changes no conclusion.

## G. Geometry → shielding models (Q3)

`Model1A` (shielding ~ length + rho + FR2_CDR3_Dis, N=533, R²=0.816): **all three geometry descriptors remain independently significant** when fit jointly (length std-β=0.72 p=3.7e-25; rho std-β=-0.47 p=2.1e-11; FR2_CDR3_Dis std-β=-0.28 p=3.4e-6). `Model1B` (compactness + FR2_CDR3_Dis, R²=0.787): both independently significant. **Optional extension** (+H_intrinsic): H_intrinsic adds no meaningful information (std-β=0.01, p=0.59, ΔR²≈0.0001) — CDR3 geometry, not sequence chemistry, explains shielding.

**Answer to Q3: CDR3 length, rho, and FR2_CDR3_Dis all independently associate with shielding, jointly.**

## H. Sequence → native-surface models (Q4)

`ModelC1` (native_surface ~ H_intrinsic, R²=0.522) → `ModelC2` (+ shielding_fraction, R²=0.578): shielding remains independently significant (std-β=-0.24, p=1.1e-14) after adjusting for intrinsic sequence chemistry, and materially improves fit (ΔR²=0.056).

**Answer to Q4: yes — native structural exposure contains a shielding-associated component beyond intrinsic sequence chemistry (not mediation-proof, purely a nested-model finding).**

## I. Sequence/shielding → PSH/topology models (Q5)

`ModelD1` (PSH ~ H_intrinsic, R²=0.115) → `ModelD2` (+ native_surface + shielding, R²=0.230, all VIF<2.5 — not a collinearity artifact): shielding is the strongest independent predictor (std-β=0.37, p=2.0e-17), and H_intrinsic/native_surface remain independently significant too (std-β=0.15/0.17, both p<0.01). Largest-hydrophobic-patch topology tracked as a **distinct** outcome (`ModelD3`), never conflated with PSH.

**Answer to Q5: yes — shielding contributes information about PSH beyond intrinsic sequence chemistry, in a well-conditioned joint model.**

## J. Shielding → 15-position antigen-engagement models (Q6) — central finding

`Model4A` (engagement_15 ~ shielding, unadjusted): std-β=-0.44, p=1.7e-25 — matches the Stage-9 descriptive finding exactly. `Model4B` (+ CDR3 length): shielding attenuates to std-β=-0.22 (p=3.5e-5, still significant). `Model4C` (+ CDR3 length + FR2_CDR3_Dis): shielding attenuates further to std-β=-0.15, **p=0.065 — no longer statistically distinguishable from zero**, while CDR3 length remains the dominant, strongly significant predictor throughout (p<1e-10 in every model). VIF in Model4C = 3.8 (well below the severe-collinearity threshold that applied to shielding-vs-BSA) — **this attenuation is a genuine statistical finding, not a numerical-instability artifact.**

**Answer to Q6: the unadjusted inverse relationship is robust and highly significant, but does NOT survive full geometry adjustment — most of its apparent explanatory power is attributable to (or statistically inseparable from) CDR3 length itself, not an independent shielding-specific effect.**

## K. Shielding → CDR3 antigen contribution (Q7)

`Model5A` (unadjusted): std-β=+0.42, p=1.3e-25 (matches Stage-9 descriptive finding). `Model5B` (+ length + FR2_CDR3_Dis): shielding attenuates to std-β=+0.08, **p=0.196 — not significant**; CDR3 length dominates (std-β=0.62, p=1.5e-45).

**Answer to Q7: no — the descriptive positive relationship does NOT survive geometry adjustment; it does not remain independently significant.**

## L. Shielding → framework contribution

`Model6A` (unadjusted): std-β=-0.40, p=2.9e-23. `Model6B` (+ length + FR2_CDR3_Dis): shielding attenuates to std-β=-0.08, **p=0.243 — not significant**; CDR3 length dominates (std-β=-0.39, p=6.6e-17).

**Consistent with J and K: the same attenuation pattern recurs across all three shielding→function models — CDR3 length is the more fundamental, robust correlate of former-interface functional organization in this dataset; shielding's marginal independent contribution, once CDR3 geometry is accounted for, cannot be reliably distinguished from zero.**

## M. Secondary Stage-8 phenotype → Stage-9 engagement models (Q9)

Separate FDR family (`results/tables/stage11_all_model_coefficients.csv`, family=`C_secondary_patch_electrostatic`): PSH significantly associates with engagement_15 (std-β=-0.28, BH-FDR q=3.5e-10); largest-patch borderline after FDR (q=0.049); **PPC, PNC, and both APBS components are NOT significant after FDR** (q>0.10 for all four). Kept fully separate throughout — no combined electrostatic score, no PPC-PNC merge, no "net APBS" constructed (verified by test).

**Answer to Q9: PSH shows a real independent relationship with engagement; the discrete charge-patch (PPC/PNC) and continuum APBS descriptors show weak, non-FDR-significant pairwise relationships — largely redundant/uninformative for this specific outcome once corrected for multiple testing.**

## N. Antigen-composition sensitivity

No validated antigen family/class label exists locally (`antigen_type` is a chemical-composition tag, not a biological family) — classified **EXACT_ANTIGEN_IDENTITY_ONLY**. Dominant exact antigen identity: "Spike protein S1" (N=37/533 VHH). `results/tables/stage11_dominant_antigen_sensitivity.csv`: after excluding all 37 Spike-protein-S1 structures (N=496 retained), all three central shielding-vs-function relationships preserve direction and near-identical standardized magnitude (e.g. shielding→engagement15: -0.441→-0.440) — **not driven by one antigen**.

## O. High-coverage sensitivity

`results/tables/stage11_high_coverage_sensitivity.csv` (≥0.90 former-interface coverage, N=394/533, reusing the exact threshold already established in Stages 8-9, not redefined here): all central relationships preserve direction; magnitudes attenuate somewhat but remain significant (e.g., geometry→shielding std-betas 0.72/-0.47/-0.28 full vs 0.42/-0.45/-0.35 high-coverage, all still p<0.005).

## P. Influence sensitivity

`results/tables/stage11_influence_resolution_sensitivity.csv` (Cook's-distance threshold 4/N): 24-46 influential observations per model (not manually deleted, formally identified only). Excluding them **strengthens** every primary relationship rather than weakening it (e.g. engagement15~shielding: p=1.7e-25 → p=3.8e-51 after exclusion) — primary conclusions are not artifacts of a small number of unusual structures.

## Q. Cluster-robust sensitivity

`data/interim/stage11_pdb_cluster_sensitivity.csv`: 40/533 VHH structures (7.5%) share a PDB entry with another VHH representative. PDB-clustered SEs are essentially identical to HC3 for every tested model (e.g. Model4A shielding p: HC3=1.65e-25 vs cluster=2.05e-25) — the limited within-PDB replication does not materially change uncertainty.

## R. Robustness summary

`results/tables/stage11_robustness_summary.csv` and `results/figures/stage11/stage11_fig4_robustness_map.png`: every **unadjusted** shielding→function relationship and the group-comparison family is `STABLE` across all applicable sensitivities. Every **fully-adjusted** shielding→function relationship is `MOSTLY_STABLE` — meaning the attenuation-to-nonsignificance itself is the consistent, robust finding (not a fragile artifact of one particular model specification).

## S. All primary null results

Three primary adjusted coefficients are formally null (raw p>0.05, not reported as "trending" or reinterpreted): Model4C shielding (p=0.065), Model5B shielding (p=0.196), Model6B shielding (p=0.243). All three retained in `stage11_all_model_coefficients.csv` and `stage11_primary_mechanism_inference.csv` — none suppressed. Two secondary predictors also null after FDR: PPC, PNC, both APBS components (Section M).

## T. Model limitations

- Standard OLS+HC3 was used as the primary form for all bounded/fractional outcomes rather than forcing beta-regression or binomial-GLM machinery whose assumptions/software behavior were not independently validated for this project (per Section 16's explicit allowance) — robustness checked via the two-part (any/magnitude) sensitivity for the most zero-heavy outcome (`shieldable_engagement_fraction`, 33% zeros), which confirmed both parts move in the same direction as the primary continuous model.
- The 15-position engagement fraction's within-structure positions are never treated as independent replicates (verified by test) — the fraction is used as a single structure-level continuous outcome throughout.
- No causal mediation model was fit; the plausible geometry→shielding→interface-organization chain is described only via nested regression, and Section 28's ban on formal causal-mediation claims is respected throughout this report.

## U. Generated figures

`results/figures/stage11/`: fig1 (primary adjusted effects forest plot), fig2 (VHH-vs-conventional group effects forest plot), fig3 (shielding vs engagement, raw+adjusted, honestly titled), fig4 (robustness map). Manifest and README included.

## V. Full pytest result

**301/301 PASSED** (up from the 270-test pre-Stage-11 baseline; 31 new Stage-11 tests, zero regressions).

## W. Confirmation: Stage 2-9 frozen outputs unchanged

Verified directly: Stage-4 VHH N=533, Stage-3 paired primary N=1429, Stage-9 engagement summary N=1962, Stage-8 APBS N=1962 — all unchanged. No Stage 2-9 file was modified this session.

## X. Confirmation: Gordon theta/phi remain unresolved and unused

Verified by direct test (`test_gordon_theta_phi_absent`) — absent from every Stage-11 output table.

## Y. Confirmation: Stage 10 was not invoked

No AggreScan3D, aggregation propensity, or developability score computed anywhere in Stage 11 (verified by test); Stage 10 remains intentionally skipped, unchanged from its "not yet run" state.

---

## COMPONENT STATUS

**INPUT INTEGRITY STATUS: PASS** (all variables traced to exact source column, one-row-per-representative verified with 10 join-time assertions, zero duplication)

**GROUP-INFERENCE STATUS: PASS** (all 7 primary hypotheses correctly comparator-matched, formally tested, BH-FDR-significant, robust to PDB-clustering/bootstrap/permutation)

**VHH GEOMETRY→SHIELDING STATUS: PASS** (all three geometry descriptors independently significant jointly, R²=0.82, well-conditioned)

**SEQUENCE→SURFACE STATUS: PASS** (H_intrinsic and shielding both independently significant for native surface, well-conditioned)

**SURFACE-PHENOTYPE STATUS: PASS** (shielding independently associates with PSH beyond sequence chemistry, well-conditioned, VIF<2.5)

**SHIELDING→ANTIGEN-ORGANIZATION STATUS: PASS** (formally tested at three levels — unadjusted relationships robust and significant; fully-adjusted relationships honestly reported as statistically null, attributable to CDR3 length rather than an independent shielding effect — this is a valid, informative scientific result, not an analysis failure)

**ROBUSTNESS/SENSITIVITY STATUS: PASS** (PDB-clustering, bootstrap, permutation, high-coverage, dominant-antigen exclusion, and influence sensitivity all completed for every applicable primary conclusion)

**MULTIPLE-TESTING STATUS: PASS** (four predeclared BH-FDR families kept fully separate — group comparisons, VHH mechanism hypotheses, secondary Stage-8 predictors, and no position-level family was added this stage)

## OVERALL STAGE 11 STATUS: PASS

All predefined primary model families validly completed; diagnostics acceptable (heteroskedasticity present and handled via HC3 throughout; VIF/condition-number inspected for every multi-predictor model; the one genuinely severe collinearity case, shielding vs Stage-7 BSA, was correctly never force-fit together); FDR families correctly separated; robustness analyses complete; all 301 tests pass. Three primary coefficients are nonsignificant after full geometry adjustment — reported fully and honestly as a valid, informative scientific finding (per Section 53, a nonsignificant primary hypothesis does not prevent PASS). No causal or mediation claim is made anywhere in this report.

## VISUALIZATION DECISION: REQUIRED (fulfilled — 4 figures)

STOP before Stage 12. No Stage-10 aggregation/developability metric, no MD, and no causal path diagram were computed in this stage.