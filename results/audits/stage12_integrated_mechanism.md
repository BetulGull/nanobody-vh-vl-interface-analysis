# Stage 12 — Final Integrated Mechanism, Evidence Synthesis, and Paper-Ready Mechanistic Model

## A. Frozen stages used

Stages 2-9 and 11 (all PASS/FROZEN except Stage 6 = PARTIAL due to Gordon θ/φ, correctly carried forward as unresolved and unused throughout). Stage 10 intentionally skipped. Verified directly: Stage-4 VHH N=533, Stage-3 paired primary N=1429, Stage-8 N=1962, Stage-9 N=1962, Stage-11 model input N=1962 — all unchanged from their frozen values.

## B. Claim-evidence ledger status

`results/tables/stage12_claim_evidence_ledger.csv` — 12 candidate claims (A-L from the Stage-12 spec), each traced to an exact source stage/table/variable, each assigned an `evidence_type`, `stage11_adjusted_status`, `robustness_status`, `allowed_final_wording`, and `prohibited_overclaim`. None accepted automatically — every claim's numeric support (effect size, p/q-value, N) was pulled directly from the frozen Stage 3-11 output files, not asserted from memory.

## C. Strongest direct structural findings

- Claim A: 28-position empirical former interface (LEVEL 1, Stage 3), validated by translated-far-away control.
- Claim D: CDR3 shielding observed in all 533/533 VHH structures (LEVEL 1, Stage 6).
- Claim I: former-interface antigen accessibility retained in 66.4% of VHH structures at the 15-position level (LEVEL 1, Stage 9).
- Claim G: VHH electrostatic potential differs from VH in BOTH positive and negative components, never a uniform shift (LEVEL 1/SECONDARY, Stage 8).

## D. Strongest adjusted statistical findings

- Claim B/C: H_intrinsic → native surface, R²=0.52 alone, R²=0.58 with shielding added (LEVEL 2, Stage 11 ModelC1/C2) — the single strongest quantitative link in the whole mechanism.
- Claim L: sequence chemistry, native surface, AND shielding all three independently associate with PSH (LEVEL 2, Stage 11 ModelD2, all VIF<2.5).
- Claim K: CDR3 length independently associates with shielding (R²=0.82 jointly with rho/FR2_CDR3_Dis) and with all three antigen-organization outcomes even where shielding's own contribution disappears (LEVEL 2, Stage 11).
- All 7 Family-A group comparisons (H_intrinsic, native surface, PSH, largest patch, framework fraction, 15-position engagement, CDR3 fraction) are BH-FDR significant with q<1e-6, robust to PDB-clustering/bootstrap/permutation.

## E. Highly coupled relationships

- **Shielding ↔ Stage-7 BSA** (Claim E): r=0.969, VIF=15.5, condition number=26,287 — cannot be separated. Represented in the final mechanism as ONE module ("CDR3 geometric accommodation"), never as two independent arrows to the same downstream node (per Section 8's hard constraint, verified by test).
- **CDR3 length vs compactness vs rho**: compactness mathematically derived from the other two; never fit together in any Stage-11 model (verified by test, carried forward unchanged).

## F. Attenuated/null relationships

**This is the central corrective finding of the entire project and is preserved, not softened, in the final mechanism.**

Claim J: the Stage-9 descriptive relationships (shielding ↑ ↔ engagement/framework ↓, CDR3 fraction ↑) were robust and highly significant unadjusted (p<1e-22 in every case), but **after adjusting for CDR3 length (+FR2_CDR3_Dis)**:

| Model | Adjusted shielding std-β | p | BH-FDR q |
|---|---|---|---|
| 4C (engagement_15) | -0.146 | 0.065 | 0.077 |
| 5B (CDR3 fraction) | +0.082 | 0.196 | 0.212 |
| 6B (framework fraction) | -0.083 | 0.243 | 0.243 |

All three are formally nonsignificant, while CDR3 length remains independently significant (p<1e-10) in every one of these same models. **The final mechanism figure draws this relationship DASHED (descriptive only), never SOLID, and represents CDR3 geometry — not shielding specifically — as the variable independently associated with antigen-interface organization** (edge `CDR3_GEOMETRIC_ACCOMMODATION_GEOMETRY_AXIS → ANTIGEN_INTERFACE_ORGANIZATION`, figure_style=SOLID, vs the shielding-specific edge, figure_style=DASHED).

Also null: Stage-8 PPC, PNC, and both APBS components show no FDR-significant relationship with 15-position engagement (Stage 11 secondary family, q>0.10 for all four) — correctly omitted (`figure_style=OMIT`) from the final mechanism figure rather than force-included.

## G. Sequence-adaptation module (Module 2)

H_intrinsic differs VHH-vs-VH (std -0.91, q=4.4e-86) and independently predicts native surface chemistry (R²=0.52) and PSH (R²=0.115 alone). Fully supported, LEVEL 1+2.

## H. CDR3 geometric-accommodation module (Module 4)

Shielding observed in 100% of VHH; CDR3 length/rho/FR2_CDR3_Dis jointly explain 82% of shielding variance; shielding and Stage-7 BSA packing are represented as one coupled module (never two). This module also independently contributes to native surface (R² +0.056) and PSH (std-β=0.367, strongest of three predictors) — but its *specific* link to antigen-interface organization does not survive geometry adjustment (Section F).

## I. Surface-phenotype module (Module 5)

PSH and largest hydrophobic patch differ VHH-vs-VH (LEVEL 1+2, both q<1e-10); electrostatics differ in both components without a net directional shift (LEVEL 1, not formally Stage-11-tested as a primary hypothesis, reported as SECONDARY_CONTEXTUAL).

## J. Antigen-interface module (Module 6)

Framework fraction, CDR3 fraction, and 15-position engagement all formally differ VHH-vs-VH (LEVEL 1+2, all q<1e-6); framework contribution is not FR2-exclusive (FR1/FR3/FR4 all contribute, Stage 9); former interface remains partially, not universally, accessible (66.4% nonzero engagement).

## K. Final supported connections

`results/tables/stage12_mechanism_edges.csv` — 12 edges with `final_include=True`, spanning: VL absence → former interface; former interface → sequence adaptation AND CDR3 accommodation; sequence adaptation → native surface; CDR3 module → native surface (independent) and PSH (independent); CDR3 geometry axis → antigen organization (independent, via length); former interface → antigen organization (direct accessibility fact); PSH → antigen organization; antigen organization + PSH → formal group differences.

## L. Connections explicitly omitted

4 edges with `final_include=False`: PPC→engagement, PNC→engagement, APBS→engagement (all non-FDR-significant secondary predictors), and Gordon θ/φ→anything (unresolved, never used). None of these appear in the final mechanism figure.

## M. Unresolved Gordon theta/phi status

Confirmed absent from every Stage-12 output (verified by test `test_gordon_theta_phi_absent`) and explicitly listed as an `UNRESOLVED, not used` node in the edge table with `final_include=False`.

## N. Stage-10 skipped status

Confirmed: no AggreScan3D, aggregation-propensity, or developability-score column appears anywhere in any Stage-12 output (verified by test).

## O. Final mechanism figure

`results/figures/stage12/stage12_fig_final_mechanism.png` — every edge traces to `stage12_mechanism_edges.csv`; solid/dashed/double-line styles applied exactly per the evidence classification in Section D-F above; no decorative or unsupported connections.

## P. Position-level integration figure

`results/figures/stage12/stage12_fig_position_level_integration.png` — refines Stage-9's fig3 with tetrad positions explicitly marked, over the frozen 28-position (25 unique integer positions) former interface.

## Q. Formal-effect figure

`results/figures/stage12/stage12_fig_final_effect_evidence.png` — includes and color-codes (red) the three attenuated/null shielding coefficients (Model4C/5B/6B) alongside the significant ones (blue), never hiding the null results.

## R. Main-paper figure selection

`results/figures/stage12/main_figure_selection.csv` — 26 candidate figures from Stages 3-11 audited; 15 recommended for the main paper set (organized into 7 candidate multi-panel figure groups), 11 demoted to supplementary. Stage-9 fig4 (unadjusted shielding-vs-engagement) explicitly marked `retain=False` — superseded by the Stage-11/12 adjusted version, since presenting it alone without the adjustment context would be misleading.

## S. Supplementary figure selection

`results/figures/stage12/supplementary_figure_selection.csv` — 11 figures (threshold sensitivity, secondary electrostatics, robustness maps, compositional/topology detail) correctly demoted from main.

## T. Limitations

`results/tables/stage12_limitations.csv` — 12 limitations (bound-state bias, crystal packing, no experimental solubility/stability/affinity data, Stage-10 skip, Gordon θ/φ, cohort selection bias, antigen heterogeneity, shielding/BSA/length interdependence, non-causality, APBS atom-sampling limitation), each with `what_it_limits`/`what_it_does_not_invalidate`/`recommended_final_wording`.

## U. Final scientific answer

See Section 40 output below (verbatim final answer).

## V. Full pytest result

**331/331 PASSED** (up from the 301-test pre-Stage-12 baseline; 30 new Stage-12 tests, zero regressions across all 11 prior stages).

## W. Confirmation: Stage 2-11 scientific outputs unchanged

Verified directly (Section A) — no Stage 2-11 file was modified this session.

---

## FINAL SCIENTIFIC ANSWER

**Plain-language paragraph:** VHH single-domain antibodies lack the light chain (VL) that normally packs against a specific, identifiable surface on the heavy chain. This project found that VHH sequences carry systematically different amino-acid chemistry at exactly that surface, and that this chemical difference is realized in the actual folded, exposed structure — not just at the sequence level. The third CDR loop (CDR3) additionally folds back to physically cover part of that surface in essentially all VHH structures examined, and this covering is tightly linked to how the loop packs against the surface (the two could not be statistically told apart). Together, these changes leave VHH with a smaller, differently organized hydrophobic surface than a conventional antibody heavy chain would have. Rather than losing its ability to bind antigen, VHH spreads its antigen contacts differently — using more of its framework surface and more of its CDR3 loop than conventional antibodies do — and the surface that VL used to occupy remains involved in binding in most, though not all, of the structures studied. One earlier-stage finding that looked like "more loop-covering means less antigen contact" turned out, on closer statistical testing, to be better explained by the loop's basic length and geometry rather than by the covering itself — a reminder that the loop's overall shape, not just how much it covers, is the more fundamental structural factor.

**Technical paragraph:** VHH domains show a statistically robust, sequence-level chemical adaptation (H_intrinsic, standardized effect -0.91, q=4.4e-86) at the empirically-derived 28-position former VH-VL interface (Stage 3), which is realized in the native solvent-exposed structure (native-surface group difference q=1.8e-253) and independently predicted by H_intrinsic (R²=0.52), with CDR3 shielding contributing additional independent explanatory power (R²→0.58, p=1.1e-14). CDR3 shielding is observed in 100% of VHH structures and is jointly explained by CDR3 length, ρ, and FR2_CDR3_Dis (R²=0.82); shielding and Stage-7 CDR3-target packing (BSA) are too collinear to separate (r=0.97, VIF=15.5) and are represented as one coupled geometric-accommodation module, which also independently associates with former-interface PSH (std-β=0.37, p=2.0e-17). VHH shows significantly lower PSH and smaller hydrophobic patches, altered (non-directional) electrostatic potential magnitude, and significantly greater framework and CDR3 contribution to antigen contacts than conventional VH (all q<1e-6). The former interface retains antigen accessibility in 66.4% of VHH structures. Critically, the Stage-9 descriptive shielding-vs-antigen-organization relationships (robust unadjusted, p<1e-22) attenuate to statistical nonsignificance (p=0.065-0.243) after adjusting for CDR3 length, which remains independently significant throughout — indicating that CDR3 geometry broadly, rather than shielding specifically, is the variable independently linked to antigen-interface reorganization in this dataset.

**What is directly demonstrated:**
- A reproducible, empirically-derived 28-position former VH-VL interface, broader than the classical tetrad.
- VHH-vs-conventional-VH sequence chemistry differences at that interface, formally significant.
- Realization of that difference in the native folded/exposed structure.
- CDR3 shielding of part of the former interface in all 533 VHH structures examined.
- Severe statistical coupling between shielding and CDR3-target packing.
- Lower/differently-organized hydrophobic patch phenotype and altered (non-uniform) electrostatics in VHH.
- Greater framework and CDR3 antigen-contact participation in VHH than conventional VH.
- Partial (not universal) retained antigen accessibility of the former interface in VHH.
- Formal statistical adjustment showing CDR3 length, not shielding specifically, as the variable independently associated with antigen-interface organization.

**What remains inferential / unresolved:**
- Whether these structural/chemical adaptations causally produce higher experimental solubility, stability, or reduced aggregation — none directly measured.
- Whether antigen-contact reorganization confers any functional/affinity advantage — not measured.
- The exact physical mechanism by which shielding and packing jointly arise (their statistical inseparability does not resolve their physical relationship).
- Gordon θ/φ CDR3 orientation descriptors — unresolved method dependency, unused throughout.
- Generalizability beyond the SAbDab2-derived, antigen-engaged, nonredundant structural cohort used in this project.

---

## COMPONENT STATUS

**CLAIM-EVIDENCE INTEGRITY STATUS: PASS** (12/12 claims traced to exact frozen source, none accepted automatically)

**SEQUENCE-ADAPTATION SYNTHESIS STATUS: PASS**

**CDR3-GEOMETRY SYNTHESIS STATUS: PASS** (coupling correctly represented as one module, never two independent mechanisms)

**SURFACE-PHENOTYPE SYNTHESIS STATUS: PASS**

**ANTIGEN-ORGANIZATION SYNTHESIS STATUS: PASS** (including the critical attenuation correction)

**STATISTICAL-INTEGRATION STATUS: PASS** (Stage-11 results faithfully carried forward, null results preserved and visualized, no new inference performed)

**FINAL-MECHANISM FIGURE STATUS: PASS** (every edge traceable to `stage12_mechanism_edges.csv`; verified by test)

**REPRODUCIBILITY STATUS: PASS** (manifest complete, all results traceable to existing scripts/outputs; no expensive recomputation performed)

## OVERALL STAGE 12 STATUS: PASS

Every central claim maps to frozen evidence; unsupported causal arrows are absent; the critical Stage-11 attenuation result is faithfully represented (dashed edge, red-highlighted null coefficients, explicit prose); shielding/BSA coupling is represented as one module; sequence/surface/shielding distinctions are preserved throughout; comparator states (VHH-native vs VH-without-VL-geometric; VHH vs VH for heavy-domain quantities) are correct; limitations are explicit; the main-paper figure set is coherent and non-redundant; the reproducibility manifest is complete; full pytest passes (331/331); Stage 2-11 scientific outputs remain unchanged.

## VISUALIZATION DECISION: REQUIRED (fulfilled — 4 final integration figures)

STOP. Do not begin Stage 10 or any new computational stage.
