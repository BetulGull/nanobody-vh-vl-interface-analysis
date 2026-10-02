"""
Final figures -- Section 6/27: predeclared CDR3-length x group comparative
extension. ONE model form, fit per outcome:

    Y = beta0 + beta1*CDR3_length + beta2*group_bin + beta3*(CDR3_length*group_bin)

group_bin = 1 for VHH, 0 for VH. beta3 (length_x_group) is the coefficient of
primary interest -- no causal interpretation, no stepwise selection.
slope_VH = beta1 (CDR3_length main effect, VH baseline)
slope_VHH = beta1 + beta3

Six predefined FDR families, corrected separately (never pooled):
  FAMILY1_regional_charge          -- sequence charge proxy at pH 7.4, by region
  FAMILY2_regional_hydropathy      -- mean Kyte-Doolittle, by region
  FAMILY3_surface_hydrophobic      -- H_native, PSH, largest hydrophobic patch
  FAMILY4_surface_electrostatics   -- PPC/PNC, APBS positive/negative mean
  FAMILY5_antigen_allocation       -- contact fractions / engagement fractions
  FAMILY6_interaction_chemistry    -- deferred until Stage-9 Arpeggio
                                       interaction-composition variables are
                                       aggregated for Figure 8 (not yet built)

Robust (HC3) SE throughout, matching stage11_stats.fit_ols_hc3 convention.
No new model fit beyond this one predeclared form; this is the sole allowed
new scientific extension per the task spec.
"""
import numpy as np
import pandas as pd
import statsmodels.api as sm
from pathlib import Path
import sys

ROOT = Path("/Users/betulgul/Desktop/staj")
sys.path.insert(0, str(ROOT / "src/qc"))
from stage11_stats import bh_fdr

FAMILIES = {
    'FAMILY1_regional_charge': [
        'FR1_charge', 'CDR1_charge', 'FR2_charge', 'CDR2_charge', 'FR3_charge',
        'CDR3_charge', 'FR4_charge', 'all_CDR_charge', 'full_domain_charge',
        'FULL_EMPIRICAL_INTERFACE_28_charge',
    ],
    'FAMILY2_regional_hydropathy': [
        'FR1_mean_hydropathy', 'CDR1_mean_hydropathy', 'FR2_mean_hydropathy',
        'CDR2_mean_hydropathy', 'FR3_mean_hydropathy', 'CDR3_mean_hydropathy',
        'FR4_mean_hydropathy', 'all_CDR_mean_hydropathy', 'full_domain_mean_hydropathy',
        'FULL_EMPIRICAL_INTERFACE_28_mean_hydropathy',
    ],
    'FAMILY3_surface_hydrophobic': [
        'H_native_full_interface_sasa_weighted_mean_KD', 'former_interface_core_PSH',
        'local_thr4.0_largest_patch_residue_count',
    ],
    'FAMILY4_surface_electrostatics': [
        'former_interface_core_PPC', 'former_interface_core_PNC',
        'interface_potential_positive_mean', 'interface_potential_negative_mean',
    ],
    'FAMILY5_antigen_allocation': [
        'framework_contact_fraction', 'cdr3_contact_fraction',
        'shieldable_engagement_fraction', 'full_interface_engagement_fraction',
    ],
}


def fit_interaction(d, outcome):
    sub = d[[outcome, 'CDR3_length', 'group_bin']].dropna().copy()
    sub['length_x_group'] = sub['CDR3_length'] * sub['group_bin']
    X = sm.add_constant(sub[['CDR3_length', 'group_bin', 'length_x_group']])
    y = sub[outcome]
    model = sm.OLS(y, X).fit(cov_type='HC3')
    n_vhh = int((sub['group_bin'] == 1).sum())
    n_vh = int((sub['group_bin'] == 0).sum())
    slope_vh = model.params['CDR3_length']
    slope_vhh = model.params['CDR3_length'] + model.params['length_x_group']
    beta_int = model.params['length_x_group']
    ci = model.conf_int().loc['length_x_group']
    p = model.pvalues['length_x_group']
    return {
        'outcome': outcome, 'N_VHH': n_vhh, 'N_VH': n_vh,
        'slope_VH': slope_vh, 'slope_VHH': slope_vhh,
        'beta_interaction': beta_int, 'CI_low': ci[0], 'CI_high': ci[1],
        'p_interaction': p,
    }


def main():
    d = pd.read_parquet(ROOT / "data/processed/final_figures_model_input.parquet").copy()
    d['group_bin'] = (d['group'] == 'VHH').astype(int)

    rows = []
    for family, outcomes in FAMILIES.items():
        fam_rows = [fit_interaction(d, o) for o in outcomes]
        pvals = [r['p_interaction'] for r in fam_rows]
        q = bh_fdr(pvals)
        for r, qi in zip(fam_rows, q):
            r2 = {'family': family}
            r2.update(r)
            r2['BH_FDR_q'] = qi
            rows.append(r2)

    out = pd.DataFrame(rows)[[
        'family', 'outcome', 'N_VHH', 'N_VH', 'slope_VH', 'slope_VHH',
        'beta_interaction', 'CI_low', 'CI_high', 'p_interaction', 'BH_FDR_q',
    ]]
    out.to_csv(ROOT / "results/tables/final_cdr3_length_group_interactions.csv", index=False)
    print(f"Wrote {len(out)} rows across {out['family'].nunique()} families")
    print(out.groupby('family').size())
    return out


if __name__ == '__main__':
    main()
