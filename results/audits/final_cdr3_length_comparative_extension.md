# Final Figures — CDR3-Length × Group Comparative Extension

The one predeclared new scientific extension permitted by the final-figures task spec. No causal interpretation of β3 (`length_x_group`) is intended anywhere below — this is a descriptive comparison of length-slopes between groups.

## A. Model specification

For each outcome Y, one model form, fit once:

```
Y = β0 + β1·CDR3_length + β2·group_bin + β3·(CDR3_length × group_bin)
```

`group_bin` = 1 for VHH, 0 for VH. `CDR3_length` is the fresh Stage-4-derived per-region length (`final_regional_sequence_physicochemistry.parquet`), not the Stage-11 `cdr3_sequence_length` column (which is populated for VHH only — see `final_figure_data_lineage.csv`). Robust HC3 standard errors throughout (`statsmodels.OLS(cov_type='HC3')`, matching the Stage-11 `fit_ols_hc3` convention). No stepwise selection; no other model form tried. Source: [`src/qc/final_cdr3_length_group_interactions.py`](../../src/qc/final_cdr3_length_group_interactions.py).

`slope_VH = β1`, `slope_VHH = β1 + β3`. β3 (`beta_interaction`) is the coefficient of primary interest.

## B. Cohort Ns

VHH N=533, VH N=1429 for every outcome except the two APBS electrostatic-potential outcomes (`interface_potential_positive_mean`/`negative_mean`), where APBS coverage reduces N to 521 VHH / 1386 VH (unchanged from Stage 8/11 — no new missingness introduced here).

## C. Six predefined FDR families, corrected separately (never pooled)

| family | scope | n_outcomes | status |
|---|---|---|---|
| FAMILY1_regional_charge | sequence charge proxy (pH 7.4), by IMGT region | 10 | complete |
| FAMILY2_regional_hydropathy | mean Kyte-Doolittle, by IMGT region | 10 | complete |
| FAMILY3_surface_hydrophobic | H_native, PSH, largest hydrophobic patch | 3 | complete |
| FAMILY4_surface_electrostatics | PPC/PNC, APBS positive/negative mean | 4 | complete |
| FAMILY5_antigen_allocation | framework/CDR3/shieldable/full-interface engagement fractions | 4 | complete |
| FAMILY6_interaction_chemistry | Stage-9 Arpeggio interaction-category composition | — | **deferred** — variables not yet aggregated (needed for Figure 8; will be added as a genuinely separate 6th family, not pooled into FAMILY1-5, when built) |

Benjamini-Hochberg FDR applied within each family only (`stage11_stats.bh_fdr`), per Section 6/27 of the task spec.

## D. A real bug found and fixed before these results were finalized

`src/qc/final_regional_physicochemistry.py`'s `CHARGE_MAP` was initially written with three-letter residue codes (`'ASP'`, `'GLU'`, ...) while the source table's `amino_acid` column uses single-letter codes — every charge column computed as exactly 0.0 for all 1962 representatives. Caught via `.describe()` inspection (all-zero mean/std), root-caused, and fixed (`{'D': -1.0, 'E': -1.0, 'R': 1.0, 'K': 1.0, 'H': 0.1}`). Verified non-zero, sensible output before regenerating FAMILY1 and the downstream merged model-input table. FAMILY2-5 were computed from unaffected columns and are unchanged by the fix (cross-checked numerically before and after).

## E-I. Results by family

Full numeric results in [`results/tables/final_cdr3_length_group_interactions.csv`](../tables/final_cdr3_length_group_interactions.csv) (31 rows: family, outcome, N_VHH, N_VH, slope_VH, slope_VHH, beta_interaction, CI_low, CI_high, p_interaction, BH_FDR_q, interpretation_tier, sign_flip).

**FAMILY1 (regional charge)**: 9/10 outcomes significant after FDR (all q<0.03 except FR3_charge, NO_CLEAR_RELATIONSHIP). Both groups' regional charge becomes more negative with increasing CDR3 length, but VHH declines more steeply for CDR3 itself, the full CDR set, the full domain, and the 28-position former interface (VHH_SPECIFIC_OR_STRONGER); several individual framework/CDR2 regions instead show `sign_flip=True` (VH trends positive, VHH trends negative) — tiered VH_SPECIFIC_OR_STRONGER or SIMILAR_SLOPE by magnitude, per the predeclared 4-tier rule (Section G), with the sign difference reported as a transparency flag rather than a 5th invented category.

**FAMILY2 (regional hydropathy)**: 6/10 significant after FDR. CDR1 and all-CDR hydropathy interactions are VHH-specific with a sign flip (VH flat/slightly positive, VHH negative); the 28-position former-interface and FR2 hydropathy interactions are VHH_SPECIFIC_OR_STRONGER with matched sign (both positive, VHH steeper) — the interface region gets *more* hydropathic with CDR3 length specifically in VHH.

**FAMILY3 (surface hydrophobic)**: PSH and largest-patch interactions both significant (q<3e-6); both groups' PSH rises with CDR3 length (matched sign) but VH rises more steeply (slope_VH=2.56 vs slope_VHH=1.63, VH_SPECIFIC_OR_STRONGER) — consistent with, and not contradicting, the standing Stage-8/11 finding that shielding/PSH covary tightly with CDR3 length in both groups. `H_native_full_interface_sasa_weighted_mean_KD` shows NO_CLEAR_RELATIONSHIP (q=0.234).

**FAMILY4 (surface electrostatics)**: no interaction survives FDR (all q≥0.064) — NO_CLEAR_RELATIONSHIP throughout, consistent with the Stage-11 finding that electrostatic descriptors are the weakest predictors of interface behavior in this cohort. This null result is reported as-is, not suppressed.

**FAMILY5 (antigen allocation)**: all 4 outcomes VHH_SPECIFIC_OR_STRONGER and significant (q<3.5e-5, matched sign in 3/4). Framework-contact-fraction and shieldable-engagement-fraction both decline more steeply with CDR3 length in VHH than VH; CDR3-contact-fraction rises more steeply in VHH. This is the clearest, most consistent family: as VHH CDR3 lengthens, antigen engagement reallocates toward CDR3 and away from framework/shieldable-interface positions more strongly than in conventional VH — directionally coherent with the already-frozen Stage-9/11 finding that CDR3 length is the dominant driver of paratope composition.

## G. Interpretation-tier rule (predeclared, applied uniformly)

- `BH_FDR_q ≥ 0.05` → `NO_CLEAR_RELATIONSHIP`
- else `|slope_VHH| / |slope_VH| > 1.5` → `VHH_SPECIFIC_OR_STRONGER`
- else `|slope_VHH| / |slope_VH| < 0.67` → `VH_SPECIFIC_OR_STRONGER`
- else → `SIMILAR_SLOPE`

Applied on magnitude ratio only, exactly the four tiers predeclared in the task spec — no fifth category invented. A separate boolean `sign_flip` column (opposite-sign slopes with both magnitudes >1e-6) is reported alongside the tier for full transparency where relevant, since magnitude-ratio tiering alone can obscure a genuine direction reversal.

## H. Limitations

- FAMILY1/2 regional charge/hydropathy are sequence-level proxies (pH-7.4 net charge, Kyte-Doolittle), never mixed with H_native/PSH/patch-topology, which are structure/exposure-level measurements — kept in separate families per the task's explicit instruction.
- `FR2_CDR3_Dis`, the geometric covariate central to Stage-11's adjustment, is VHH-only in the frozen data (never computed for VH) and is therefore not usable as a covariate in this VHH-vs-VH interaction model; the model here is unadjusted for geometry by necessity, unlike the Stage-11 formal models it is not a replacement for.
- No causal claim is made about β3 anywhere in this document or downstream figures; it describes an association between CDR3 length, group, and the outcome, nothing more.
- FAMILY6 remains deferred pending Figure-8 Arpeggio composition aggregation; this document will be updated (not silently replaced) when it is added.

## I. Status

FAMILY1-5 complete, 31 outcomes tested, BH-FDR applied within-family. Full null results (FAMILY4, several FAMILY1/2 individual outcomes) are reported and retained, not suppressed, per the task's explicit null-result-visibility requirement.
