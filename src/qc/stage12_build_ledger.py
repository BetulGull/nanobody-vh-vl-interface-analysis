"""Stage 12 — Claim-Evidence Ledger."""
import pandas as pd
from pathlib import Path

ROOT = Path("/Users/betulgul/Desktop/staj")

CLAIMS = [
dict(claim_id='A', claim_text='The conventional VH-VL interface is a reproducible empirical structural surface rather than only the classical FR2 tetrad.',
     mechanistic_level='MODULE1', source_stage=3, source_table='stage3_empirical_vh_vl_interface_primary.csv',
     source_variable_or_result='28-position primary_binary_interface set (delta-SASA-derived, translated-far-away control confirmed genuine occlusion)',
     evidence_type='DIRECT_STRUCTURAL', effect_direction='n/a', effect_size_or_summary='28 positions, spanning FR2/FR3/CDR3/FR4, not limited to IMGT 42/49/50/52',
     stage11_adjusted_status='not applicable (not a Stage-11 hypothesis)', robustness_status='validated on real-structure sample + translated-far-away geometric control (Stage 3)',
     allowed_final_wording='OBSERVED: the empirical former VH-VL interface spans 28 IMGT positions, of which only 4 are the classical tetrad',
     prohibited_overclaim='Do not say the interface IS the tetrad, or that the tetrad alone explains VL loss accommodation', main_or_secondary='main', notes='Frozen Stage-3 result, unchanged'),

dict(claim_id='B', claim_text='VHHs show distinct amino-acid chemistry at the former-interface relative to conventional VH.',
     mechanistic_level='MODULE2', source_stage=4, source_table='stage4_sequence_interface_chemistry.csv / stage11_primary_group_inference.csv',
     source_variable_or_result='H_intrinsic_binary_mean_KD, VHH vs VH', evidence_type='DIRECT_SEQUENCE', effect_direction='VHH lower (less hydrophobic, more charged-favoring) than VH',
     stage11_adjusted_status='Formally tested, Family A hypothesis A', robustness_status='BH-FDR q=4.4e-86, PDB-cluster p=2.9e-87, permutation p=0.0002 -- STABLE',
     allowed_final_wording='OBSERVED + statistically ASSOCIATED: VHH former-interface intrinsic sequence chemistry differs from conventional VH (standardized effect -0.91)',
     prohibited_overclaim='Do not claim this alone explains solubility', main_or_secondary='main', notes=''),

dict(claim_id='C', claim_text='The sequence-level chemical adaptation is reflected in the native solvent-exposed 3D surface.',
     mechanistic_level='MODULE3', source_stage=5, source_table='stage5_sequence_surface_chemistry.parquet / stage11_primary_group_inference.csv, stage11 ModelC1',
     source_variable_or_result='H_native_full_interface_sasa_weighted_mean_KD, VHH_NATIVE vs VH_WITHOUT_VL_GEOMETRIC; ModelC1 (native_surface ~ H_intrinsic)',
     evidence_type='ADJUSTED_ASSOCIATION', effect_direction='VHH lower native surface hydrophobicity than VH; H_intrinsic strongly predicts native surface (R2=0.52)',
     stage11_adjusted_status='Group: Family A hypothesis B (q=1.8e-253). Mechanism: ModelC1 R2=0.522, p=2.2e-128',
     robustness_status='STABLE (group diff); ModelC1 well-conditioned (single predictor)',
     allowed_final_wording='OBSERVED + ASSOCIATED: native exposed former-interface chemistry differs VHH-vs-VH and is strongly predicted by intrinsic sequence chemistry',
     prohibited_overclaim='Do not claim sequence chemistry fully explains native surface (R2=0.52, not 1.0; shielding adds independent information, see Claim L)', main_or_secondary='main', notes=''),

dict(claim_id='D', claim_text='CDR3 physically shields a substantial portion of the non-CDR3 former-interface in native VHH structures.',
     mechanistic_level='MODULE4', source_stage=6, source_table='stage6_cdr3_shielding_geometry_combined.parquet',
     source_variable_or_result='shielding_fraction (mean 0.215, N=533), delta_sasa_shield (mean 245.1 A2)',
     evidence_type='DIRECT_STRUCTURAL', effect_direction='positive shielding observed in all 533 VHH structures',
     stage11_adjusted_status='not a group-comparison hypothesis (VHH-only geometric counterfactual)', robustness_status='validated via same-coordinate no-CDR3 SASA counterfactual, translated-far-away control (Stage 6)',
     allowed_final_wording='OBSERVED: CDR3 removal (geometric SASA counterfactual) increases exposure at the 15-position shieldable former interface in all 533 VHH structures',
     prohibited_overclaim='Do not call the no-CDR3 state a physical mutant or claim CDR3 "replaces" VL', main_or_secondary='main', notes=''),

dict(claim_id='E', claim_text='CDR3 shielding is geometrically coupled to CDR3-former-interface packing.',
     mechanistic_level='MODULE4', source_stage=7, source_table='stage7_per_vhh_packing_summary.csv / stage11_collinearity_audit.csv',
     source_variable_or_result='nonadj_BSA_pair_total vs shielding_fraction correlation + VIF',
     evidence_type='HIGHLY_COUPLED_GEOMETRIC_EVIDENCE', effect_direction='positive correlation', effect_size_or_summary='r=0.969, VIF=15.5, condition number=26,287 (Stage 11)',
     stage11_adjusted_status='Formally quantified as inseparable (Stage 11 collinearity gate)', robustness_status='confirmed severe, not estimable as independent predictors',
     allowed_final_wording='OBSERVED + HIGHLY COUPLED: shielding and CDR3-target packing (BSA) are so strongly correlated that their independent contributions cannot be reliably separated in this dataset',
     prohibited_overclaim='Never present shielding and BSA as two independent parallel mechanisms or interpret a coefficient difference between them as evidence of distinct mechanisms', main_or_secondary='main', notes='Directly drives Section 8 constraint'),

dict(claim_id='F', claim_text='VHH former-interface hydrophobic surface organization differs from VH, including lower PSH / smaller hydrophobic patch organization.',
     mechanistic_level='MODULE5', source_stage=8, source_table='stage8_surface_patch_descriptors.parquet, stage8_patch_topology.parquet / stage11_primary_group_inference.csv',
     source_variable_or_result='former_interface_core_PSH, local_thr4.0_largest_patch_residue_count, VHH vs VH',
     evidence_type='ROBUST_ADJUSTED_ASSOCIATION' if False else 'DIRECT_STRUCTURAL', effect_direction='VHH lower PSH and smaller largest patch than VH',
     stage11_adjusted_status='Family A hypotheses C (PSH q=1.4e-11; largest patch q=5.8e-55)', robustness_status='STABLE, PDB-cluster and permutation confirmed',
     allowed_final_wording='OBSERVED + statistically ASSOCIATED: VHH shows lower former-interface PSH and smaller largest hydrophobic patch than conventional VH',
     prohibited_overclaim='Do not equate PSH/patch-size difference with lower aggregation propensity (not tested)', main_or_secondary='main', notes=''),

dict(claim_id='G', claim_text='VHH former-interface electrostatics differ from VH, but cannot be summarized as a simple positive-vs-negative shift.',
     mechanistic_level='MODULE5', source_stage=8, source_table='stage8_apbs_interface_electrostatics.parquet',
     source_variable_or_result='interface_potential_positive_mean, interface_potential_negative_mean',
     evidence_type='DIRECT_STRUCTURAL', effect_direction='VHH more positive in positive component (118.3 vs 115.4) AND more negative in negative component (-168.2 vs -157.3)',
     stage11_adjusted_status='not a Stage-11 formal hypothesis (secondary/contextual)', robustness_status='corrected in Stage-8 closure audit; re-verified in Stage 9/11 as separate components throughout',
     allowed_final_wording='OBSERVED: VHH shows stronger potential magnitude in BOTH the positive and negative components sampled at exposed former-interface atoms, not a uniform directional shift',
     prohibited_overclaim='Never write "VHH is more electropositive" or "less negative"; never construct positive-minus-negative net score; never call this a solvent-excluded molecular-surface integral', main_or_secondary='secondary', notes=''),

dict(claim_id='H', claim_text='VHH antigen-interface organization differs from conventional antibody organization, including greater framework participation.',
     mechanistic_level='MODULE6', source_stage=9, source_table='stage9_antigen_engagement_summary.parquet / stage11_primary_group_inference.csv',
     source_variable_or_result='framework_contact_fraction, cdr3_contact_fraction, VHH vs conventional VH',
     evidence_type='ROBUST_ADJUSTED_ASSOCIATION', effect_direction='VHH higher framework fraction (0.209 vs 0.163) and higher CDR3 fraction (0.470 vs 0.403) than VH',
     stage11_adjusted_status='Family A hypotheses D (q=9.3e-8) and F (q=1.2e-12)', robustness_status='STABLE across PDB-cluster/bootstrap/permutation',
     allowed_final_wording='OBSERVED + statistically ASSOCIATED: VHH shows significantly greater framework participation and CDR3 dominance in antigen contacts than conventional VH',
     prohibited_overclaim='Do not equate framework participation with higher affinity or "better" binding; do not reduce to FR2 alone (FR1/FR3/FR4 all contribute, Stage 9)', main_or_secondary='main', notes=''),

dict(claim_id='I', claim_text='The former VH-VL surface remains antigen-accessible in a subset of VHHs rather than becoming universally functionally occluded.',
     mechanistic_level='MODULE6', source_stage=9, source_table='stage9_former_interface_antigen_position_summary.csv, vhh_shieldable_interface_engagement.csv',
     source_variable_or_result='shieldable_engagement_fraction (mean 0.186, N=533; 33.6% structures show zero engagement, but 66.4% show nonzero)',
     evidence_type='DIRECT_STRUCTURAL', effect_direction='partial, heterogeneous engagement', effect_size_or_summary='mean 18.6% of 15-position interface contacts antigen; two-part sensitivity confirms both presence and magnitude vary',
     stage11_adjusted_status='descriptive distributional fact, not a formal hypothesis test', robustness_status='confirmed via manual validation (Stage 9), both high- and zero-engagement VHH cases genuinely exist',
     allowed_final_wording='OBSERVED: the former VH-VL surface participates in antigen engagement in a majority but not all VHH structures (66.4% nonzero at the 15-position level)',
     prohibited_overclaim='Do not say every former-interface residue remains functional, and do not say the interface is universally sealed', main_or_secondary='main', notes=''),

dict(claim_id='J', claim_text='The descriptive shielding-antigen-engagement relationship does NOT represent a robust independent shielding effect after CDR3 geometry adjustment.',
     mechanistic_level='MODULE7', source_stage=11, source_table='stage11_primary_mechanism_inference.csv',
     source_variable_or_result='Model4C/5B/6B shielding_fraction coefficients', effect_direction='attenuates toward null',
     effect_size_or_summary='Model4C std-beta=-0.146 p=0.065 (q=0.077); Model5B std-beta=0.082 p=0.196 (q=0.212); Model6B std-beta=-0.083 p=0.243 (q=0.243) -- all three formally nonsignificant after adjusting for CDR3 length (+FR2_CDR3_Dis)',
     evidence_type='HIGHLY_COUPLED_GEOMETRIC_EVIDENCE', stage11_adjusted_status='Formally tested and null in the fully-adjusted models (VIF<4, not a numerical-instability artifact)',
     robustness_status='attenuation itself confirmed MOSTLY_STABLE across sensitivities (Stage 11 Section R)',
     allowed_final_wording='ASSOCIATED (unadjusted) but NOT_SUPPORTED as an independent effect after adjustment: the raw shielding-vs-engagement relationship is robust and significant, but does not survive CDR3-geometry adjustment',
     prohibited_overclaim='Never draw shielding -> reduced antigen engagement as a standalone causal or strong independent arrow', main_or_secondary='main', notes='This is the single most important correction carried into the final mechanism figure'),

dict(claim_id='K', claim_text='CDR3 geometry is a shared organizing variable connecting shielding and antigen-interface redistribution.',
     mechanistic_level='MODULE7', source_stage=11, source_table='stage11_all_model_coefficients.csv',
     source_variable_or_result='cdr3_sequence_length as dominant predictor in Model4B/4C, Model5B, Model6B',
     effect_direction='CDR3 length remains significant (p<1e-10) in every adjusted model where shielding becomes null',
     evidence_type='ROBUST_ADJUSTED_ASSOCIATION', stage11_adjusted_status='Formally supported: CDR3 length independently significant in Model1A (shielding), Model4B/4C (engagement), Model5B (CDR3 fraction), Model6B (framework fraction)',
     robustness_status='STABLE (influence-exclusion, resolution-adjustment, high-coverage all preserve significance)',
     allowed_final_wording='ASSOCIATED: CDR3 length (and related geometry descriptors) independently associates with both shielding and antigen-interface organization measures, consistent with a shared underlying geometric axis',
     prohibited_overclaim='Do not say "CDR3 length causes the VHH phenotype" -- this is a statistical association across nested models, not a causal demonstration', main_or_secondary='main', notes=''),

dict(claim_id='L', claim_text='Sequence chemistry and shielding both retain independent association with native/surface patch phenotype where supported by Stage 11.',
     mechanistic_level='MODULE3_MODULE5', source_stage=11, source_table='stage11_all_model_coefficients.csv (ModelC2, ModelD2)',
     source_variable_or_result='ModelC2: H_intrinsic + shielding -> native_surface (both significant); ModelD2: H_intrinsic + native_surface + shielding -> PSH (all three significant)',
     effect_direction='all positive/negative as expected, all independently significant', effect_size_or_summary='ModelC2 R2=0.578 (vs 0.522 H_intrinsic alone); ModelD2 R2=0.230 (vs 0.115 H_intrinsic alone), VIF<2.5 throughout',
     evidence_type='ROBUST_ADJUSTED_ASSOCIATION', stage11_adjusted_status='Formally tested and supported, well-conditioned (not a collinearity artifact)',
     robustness_status='STABLE (VIF/condition-number acceptable)',
     allowed_final_wording='ASSOCIATED: shielding contributes independent information about native surface chemistry and PSH phenotype beyond intrinsic sequence chemistry alone',
     prohibited_overclaim='Do not imply shielding replaces sequence-level adaptation, or vice versa -- both are complementary, independently-supported forms of accommodation', main_or_secondary='main', notes=''),
]


def build():
    df = pd.DataFrame(CLAIMS)
    out = ROOT / "results/tables/stage12_claim_evidence_ledger.csv"
    df.to_csv(out, index=False)
    print(f"Wrote {out}, {len(df)} claims")
    return df


if __name__ == '__main__':
    build()
