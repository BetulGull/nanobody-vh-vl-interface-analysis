# Stage 12 — Manuscript-Ready Results Section Outline

## 1. Empirical definition of the former VH-VL interface
**Primary figure:** FIG_S3_01. **Primary table:** `stage3_empirical_vh_vl_interface_primary.csv`.
**Central result:** 28 IMGT positions meet ΔSASA≥5Å² AND interface_frequency_5≥0.25 AND observation_fraction≥0.25; validated via translated-far-away geometric control.
**Allowed wording:** OBSERVED. **Caveat:** the interface is defined from paired VH-VL crystal structures and reflects that specific bound conformation.

## 2. VHH sequence adaptation at the former interface
**Primary figure:** FIG_S4_01 / stage12_fig_vhh_vs_conventional_summary (panel 1). **Primary table:** `stage11_primary_group_inference.csv` (row A).
**Central result:** H_intrinsic differs VHH vs VH, standardized effect -0.91, BH-FDR q=4.4e-86.
**Allowed wording:** OBSERVED + ASSOCIATED. **Caveat:** intrinsic (sequence-only) descriptor; does not by itself demonstrate structural realization.

## 3. Sequence adaptation is realized in the native solvent-exposed surface
**Primary figure:** FIG_S5_01/02. **Primary table:** `stage11_primary_group_inference.csv` (row B); `stage11_all_model_coefficients.csv` (ModelC1/C2).
**Central result:** native surface differs VHH vs VH (q=1.8e-253); H_intrinsic predicts native surface (R²=0.52); shielding adds independent information (R²→0.58).
**Allowed wording:** OBSERVED + ASSOCIATED. **Caveat:** R²=0.58 leaves substantial unexplained variance; not a full explanation.

## 4. CDR3 provides a coupled shielding/packing accommodation
**Primary figure:** FIG_S6_01/02, stage7_fig1. **Primary table:** `stage11_collinearity_audit.csv`, `stage11_all_model_coefficients.csv` (Model1A/1B).
**Central result:** shielding observed in all 533 VHH; CDR3 length/rho/FR2_CDR3_Dis independently associate with shielding jointly (R²=0.82); shielding and Stage-7 BSA are inseparable (VIF=15.5).
**Allowed wording:** OBSERVED + ASSOCIATED + CONSISTENT_WITH (coupled module). **Caveat:** shielding and packing must be presented as one module, never as independent parallel mechanisms.

## 5. VHH surface patch organization differs from conventional VH
**Primary figure:** stage8 fig1/fig4. **Primary table:** `stage11_primary_group_inference.csv` (rows C).
**Central result:** VHH shows lower PSH (q=1.4e-11) and smaller largest hydrophobic patch (q=5.8e-55); electrostatic potential differs in both positive and negative components (not a uniform shift).
**Allowed wording:** OBSERVED + ASSOCIATED (patch); OBSERVED only (electrostatics, not formally Stage-11-tested). **Caveat:** never say "VHH is more electropositive/less negative"; APBS is atom-sampled, not a molecular-surface integral.

## 6. Antigen-interface organization is redistributed in VHH
**Primary figure:** stage9 fig1, stage12_fig_position_level_integration. **Primary table:** `stage11_primary_group_inference.csv` (rows D, E, F).
**Central result:** VHH shows greater framework participation (q=9.3e-8), greater 15-position former-interface engagement (q=4.9e-27), greater CDR3 dominance (q=1.2e-12) than conventional VH.
**Allowed wording:** OBSERVED + ASSOCIATED. **Caveat:** framework participation ≠ higher affinity; not reducible to FR2 alone (FR1/FR3/FR4 all contribute).

## 7. Formal inference separates shared CDR3 geometry from independent shielding effects
**Primary figure:** stage11 fig1/fig3, stage12_fig_final_effect_evidence. **Primary table:** `stage11_primary_mechanism_inference.csv` (Model4/5/6 series).
**Central result:** unadjusted shielding-vs-function relationships are robust (p<1e-22); after adjusting for CDR3 length (+FR2_CDR3_Dis), shielding's independent contribution is not statistically distinguishable from zero (p=0.065-0.243) in all three tested outcomes, while CDR3 length remains dominant throughout.
**Allowed wording:** ASSOCIATED (unadjusted) / NOT_SUPPORTED (as independent adjusted effect). **Caveat:** this is the single most important corrective finding of the whole project — must never be simplified back to "shielding reduces engagement."

## 8. Integrated mechanistic model
**Primary figure:** stage12_fig_final_mechanism. **Primary table:** `stage12_mechanism_edges.csv`, `stage12_claim_evidence_ledger.csv`.
**Central result:** modular structure — sequence adaptation and CDR3 geometric accommodation (as one coupled module) both independently feed into surface patch phenotype, which together with the CDR3 geometry axis (not shielding per se) relates to antigen-interface organization; all differences are formally confirmed at the group level.
**Allowed wording:** CONSISTENT_WITH (integrated model), each component's own tier as listed in `stage12_claim_language_guide.csv`. **Caveat:** no causal claim; no experimental solubility/stability/affinity claim anywhere.
