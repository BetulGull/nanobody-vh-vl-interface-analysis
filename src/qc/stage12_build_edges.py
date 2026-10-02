"""Stage 12 — Mechanism Edge Table."""
import pandas as pd
from pathlib import Path

ROOT = Path("/Users/betulgul/Desktop/staj")

EDGES = [
dict(source_node='ABSENT_VL', target_node='EMPIRICAL_FORMER_VH_VL_SURFACE', relationship='creates',
     evidence_level='LEVEL 1', supporting_stage=3, supporting_result='28-position empirical interface, delta-SASA + translated-far-away control',
     adjusted_support='n/a', robustness='validated', collinearity_issue='none', causal_claim_allowed='structural fact only (VL absence -> no VL burial by definition)',
     figure_style='SOLID', final_include=True, notes='Module 1'),

dict(source_node='EMPIRICAL_FORMER_VH_VL_SURFACE', target_node='SEQUENCE_ADAPTATION', relationship='is the target of',
     evidence_level='LEVEL 1', supporting_stage=4, supporting_result='H_intrinsic VHH vs VH group difference, q=4.4e-86',
     adjusted_support='Stage 11 Family A', robustness='STABLE', collinearity_issue='none', causal_claim_allowed='observed difference, no causal claim',
     figure_style='SOLID', final_include=True, notes='Module 2'),

dict(source_node='EMPIRICAL_FORMER_VH_VL_SURFACE', target_node='CDR3_GEOMETRIC_ACCOMMODATION', relationship='is shielded/packed by',
     evidence_level='LEVEL 1', supporting_stage=6, supporting_result='shielding_fraction observed in all 533 VHH (Stage 6), coupled to Stage-7 BSA',
     adjusted_support='n/a (direct structural + Stage-11 collinearity confirms coupling)', robustness='validated (SASA counterfactual + translated-far-away control)',
     collinearity_issue='n/a for this edge', causal_claim_allowed='observed geometric fact', figure_style='SOLID', final_include=True, notes='Module 4'),

dict(source_node='SEQUENCE_ADAPTATION', target_node='NATIVE_SURFACE_CHEMISTRY', relationship='predicts',
     evidence_level='LEVEL 2', supporting_stage=11, supporting_result='ModelC1: native_surface ~ H_intrinsic, R2=0.522, p=2.2e-128',
     adjusted_support='yes (single predictor, well-conditioned)', robustness='STABLE', collinearity_issue='none (VIF=1.0)',
     causal_claim_allowed='"retained an independent association" only', figure_style='SOLID', final_include=True, notes='Module 2->3, strongest link in the whole mechanism'),

dict(source_node='CDR3_GEOMETRIC_ACCOMMODATION', target_node='NATIVE_SURFACE_CHEMISTRY', relationship='adds independent information to',
     evidence_level='LEVEL 2', supporting_stage=11, supporting_result='ModelC2: +shielding_fraction, incremental R2 +0.056, p=1.1e-14',
     adjusted_support='yes, adjusted for H_intrinsic', robustness='STABLE', collinearity_issue='none (VIF=1.04 in ModelC2)',
     causal_claim_allowed='"shielding retained independent association with native surface after adjustment for sequence chemistry"',
     figure_style='DASHED', final_include=True, notes='Complementary, not replacing sequence adaptation'),

dict(source_node='SEQUENCE_ADAPTATION', target_node='SURFACE_PATCH_ELECTROSTATIC_PHENOTYPE', relationship='predicts',
     evidence_level='LEVEL 2', supporting_stage=11, supporting_result='ModelD1/D2: PSH ~ H_intrinsic (+native_surface+shielding), all significant, VIF<2.5',
     adjusted_support='yes', robustness='STABLE', collinearity_issue='none', causal_claim_allowed='"independently associated with PSH"',
     figure_style='SOLID', final_include=True, notes='Module 2/3->5'),

dict(source_node='CDR3_GEOMETRIC_ACCOMMODATION', target_node='SURFACE_PATCH_ELECTROSTATIC_PHENOTYPE', relationship='adds independent information to',
     evidence_level='LEVEL 2', supporting_stage=11, supporting_result='ModelD2: shielding_fraction std-beta=0.367, p=2.0e-17, strongest of the three predictors',
     adjusted_support='yes, adjusted for H_intrinsic + native_surface', robustness='STABLE (VIF<2.5)', collinearity_issue='none in this specific model',
     causal_claim_allowed='"independently associated with PSH beyond intrinsic sequence chemistry"', figure_style='SOLID', final_include=True, notes='Module 4->5'),

dict(source_node='CDR3_GEOMETRIC_ACCOMMODATION', target_node='ANTIGEN_INTERFACE_ORGANIZATION', relationship='descriptively associated with (NOT independent after adjustment)',
     evidence_level='LEVEL 3/UNRESOLVED-INDEPENDENCE', supporting_stage=9, supporting_result='Stage-9 descriptive: shielding vs engagement_15/framework/CDR3-fraction, r=-0.29 to +0.42',
     adjusted_support='NO -- Stage 11 Model4C/5B/6B all nonsignificant (p=0.065-0.243) after CDR3-length+FR2_CDR3_Dis adjustment',
     robustness='attenuation confirmed MOSTLY_STABLE', collinearity_issue='shares variance with CDR3 length (not BSA-level severe, VIF<4, but attenuates fully)',
     causal_claim_allowed='NONE as independent link -- represent only as a descriptive/shared relationship, see Claim J',
     figure_style='DASHED', final_include=True, notes='CRITICAL: this is the corrected edge per Section 6 of the Stage-12 spec. Never SOLID.'),

dict(source_node='CDR3_GEOMETRIC_ACCOMMODATION_GEOMETRY_AXIS', target_node='ANTIGEN_INTERFACE_ORGANIZATION', relationship='independently associated with (via CDR3 length)',
     evidence_level='LEVEL 2', supporting_stage=11, supporting_result='CDR3 length significant in Model4B/4C (p<1e-10), Model5B (p=1.5e-45), Model6B (p=6.6e-17)',
     adjusted_support='yes', robustness='STABLE', collinearity_issue='moderate with shielding (VIF~2-4), not severe',
     causal_claim_allowed='"CDR3 geometry independently associates with antigen-interface organization measures"',
     figure_style='SOLID', final_include=True, notes='This is the edge that REPLACES the naive shielding->engagement arrow, per Claim K'),

dict(source_node='EMPIRICAL_FORMER_VH_VL_SURFACE_VHH', target_node='ANTIGEN_INTERFACE_ORGANIZATION', relationship='remains partially antigen-accessible in',
     evidence_level='LEVEL 1', supporting_stage=9, supporting_result='shieldable_engagement_fraction: 66.4% of VHH show nonzero 15-position engagement (mean 0.186)',
     adjusted_support='n/a (descriptive distributional fact)', robustness='confirmed via manual validation', collinearity_issue='none',
     causal_claim_allowed='observed structural/functional fact', figure_style='SOLID', final_include=True, notes='Module 6, supports Claim I -- not universally occluded'),

dict(source_node='SEQUENCE_ADAPTATION_PLUS_CDR3_ACCOMMODATION', target_node='VHH_vs_VH_GROUP_DIFFERENCES', relationship='co-occur with (group-level)',
     evidence_level='LEVEL 1', supporting_stage=11, supporting_result='7/7 primary Family-A group hypotheses significant after BH-FDR (q<1e-6 all)',
     adjusted_support='yes, formally tested', robustness='STABLE (PDB-cluster, bootstrap, permutation)', collinearity_issue='none (group indicator only)',
     causal_claim_allowed='"VHH and conventional VH/Fv formally differ across all seven predeclared structural/functional measures"',
     figure_style='SOLID', final_include=True, notes='Summarizes Module 7 group-comparison family'),

dict(source_node='STAGE6_SHIELDING', target_node='STAGE7_BSA_PACKING', relationship='HIGHLY COUPLED (not two mechanisms)',
     evidence_level='LEVEL 4', supporting_stage=11, supporting_result='r=0.969, VIF=15.5, condition number=26,287',
     adjusted_support='cannot be separated', robustness='confirmed severe collinearity', collinearity_issue='SEVERE',
     causal_claim_allowed='NONE -- must never be drawn as two separate arrows to the same downstream node',
     figure_style='DOUBLE_LINE_COUPLED', final_include=True, notes='Represented as ONE module (CDR3 geometric accommodation) containing both, per Section 8'),

dict(source_node='STAGE8_PPC', target_node='ANTIGEN_INTERFACE_ORGANIZATION', relationship='weak pairwise association (not FDR-significant)',
     evidence_level='SECONDARY_CONTEXTUAL', supporting_stage=11, supporting_result='SecondaryC_PPC: std-beta=0.087, q=0.10',
     adjusted_support='no (fails secondary FDR family)', robustness='NOT robust', collinearity_issue='none', causal_claim_allowed='none',
     figure_style='OMIT', final_include=False, notes='Excluded from final figure -- not FDR-significant'),

dict(source_node='STAGE8_PNC', target_node='ANTIGEN_INTERFACE_ORGANIZATION', relationship='weak pairwise association (not FDR-significant)',
     evidence_level='SECONDARY_CONTEXTUAL', supporting_stage=11, supporting_result='SecondaryC_PNC: std-beta=-0.080, q=0.10',
     adjusted_support='no', robustness='NOT robust', collinearity_issue='none', causal_claim_allowed='none', figure_style='OMIT', final_include=False, notes=''),

dict(source_node='STAGE8_APBS', target_node='ANTIGEN_INTERFACE_ORGANIZATION', relationship='weak pairwise association (not FDR-significant)',
     evidence_level='SECONDARY_CONTEXTUAL', supporting_stage=11, supporting_result='positive/negative APBS components q=0.27/0.80',
     adjusted_support='no', robustness='NOT robust', collinearity_issue='none', causal_claim_allowed='none', figure_style='OMIT', final_include=False, notes=''),

dict(source_node='GORDON_THETA_PHI', target_node='ANY', relationship='UNRESOLVED, not used',
     evidence_level='UNRESOLVED', supporting_stage=6, supporting_result='Text S2 unavailable, no corroborating code (Stage 6 audit)',
     adjusted_support='n/a', robustness='n/a', collinearity_issue='n/a', causal_claim_allowed='none', figure_style='OMIT', final_include=False, notes='Explicitly excluded from all Stage 7-12 work'),
]


def build():
    df = pd.DataFrame(EDGES)
    out = ROOT / "results/tables/stage12_mechanism_edges.csv"
    df.to_csv(out, index=False)
    print(f"Wrote {out}, {len(df)} edges ({df['final_include'].sum()} included, {(~df['final_include']).sum()} omitted)")
    return df


if __name__ == '__main__':
    build()
