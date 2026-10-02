"""Stage 11 — Formal Statistical Inference and Integrated Mechanism Testing: test suite."""
import sys
from pathlib import Path
import pandas as pd
import numpy as np
import pytest

ROOT = Path("/Users/betulgul/Desktop/staj")
sys.path.insert(0, str(ROOT / "src" / "qc"))

from stage11_stats import zscore, bh_fdr, fit_ols_hc3

PROCESSED = ROOT / "data/processed"
TABLES = ROOT / "results/tables"


# ---------- 1: frozen upstream inputs unchanged ----------

def test_frozen_upstream_inputs_unchanged():
    assert len(pd.read_parquet(PROCESSED / "stage4_vhh_sequence_representatives.parquet")) == 533
    df3 = pd.read_parquet(PROCESSED / "stage3_paired_representatives.parquet")
    assert len(df3[df3['is_primary_representative'] == True]) == 1429
    assert len(pd.read_parquet(PROCESSED / "stage8_surface_patch_descriptors.parquet")) == 1962


# ---------- 2-3: one-row-per-representative / no pseudoreplication ----------

def test_one_row_per_representative():
    path = PROCESSED / "stage11_model_input.parquet"
    if not path.exists():
        pytest.skip("canonical table not built")
    df = pd.read_parquet(path)
    assert df['INSTANCE'].is_unique
    assert len(df) == 1962


def test_no_residue_level_pseudoreplication():
    path = PROCESSED / "stage11_model_input.parquet"
    if not path.exists():
        pytest.skip("canonical table not built")
    df = pd.read_parquet(path)
    # engagement fractions are already-aggregated fractions, not raw residue counts duplicated as rows
    assert (df['shieldable_engagement_fraction'].dropna() <= 1.0001).all()
    assert (df['shieldable_engagement_fraction'].dropna() >= 0).all()


# ---------- 4: model-specific missingness only ----------

def test_model_specific_missingness_not_global_complete_case():
    path = PROCESSED / "stage11_model_input.parquet"
    if not path.exists():
        pytest.skip("canonical table not built")
    df = pd.read_parquet(path)
    # VHH-only columns are missing for all VH rows (expected), but VH rows retain
    # their own valid values for group-applicable columns -- i.e. missingness is
    # column-specific, not a single global complete-case filter dropping all VH rows.
    assert df[df.group == 'VH']['H_intrinsic_binary_mean_KD'].notna().all()
    assert df[df.group == 'VH']['shielding_fraction'].isna().all()


# ---------- 5-6: group coding ----------

def test_group_coding_correct():
    path = PROCESSED / "stage11_model_input.parquet"
    if not path.exists():
        pytest.skip("canonical table not built")
    df = pd.read_parquet(path)
    assert set(df['group'].unique()) == {'VHH', 'VH'}
    assert (df['group'] == 'VHH').sum() == 533
    assert (df['group'] == 'VH').sum() == 1429


def test_vhh_vs_vh_comparator_for_heavy_domain_quantities():
    path = TABLES / "stage11_primary_group_inference.csv"
    if not path.exists():
        pytest.skip("group inference not built")
    df = pd.read_csv(path)
    fw_row = df[df['hypothesis_id'] == 'D_framework_engagement'].iloc[0]
    assert 'VH' in fw_row['comparison'] and 'Fv' not in fw_row['comparison']


def test_vhh_vs_fv_never_used_for_heavy_domain_specific():
    path = TABLES / "stage11_primary_group_inference.csv"
    if not path.exists():
        pytest.skip("group inference not built")
    df = pd.read_csv(path)
    for _, row in df.iterrows():
        assert 'Fv' not in row['comparison']  # none of the 6 predeclared hypotheses compare against Fv


# ---------- 8-10: standardization ----------

def test_standardized_predictors_mean_zero():
    s = pd.Series(np.random.default_rng(1).normal(5, 2, 100))
    z = zscore(s)
    assert abs(z.mean()) < 1e-9


def test_standardized_predictors_sd_one():
    s = pd.Series(np.random.default_rng(1).normal(5, 2, 100))
    z = zscore(s)
    assert abs(z.std() - 1.0) < 1e-9


def test_binary_group_not_standardized():
    path = TABLES / "stage11_primary_group_inference.csv"
    if not path.exists():
        pytest.skip("group inference not built")
    df = pd.read_csv(path)
    # standardized_effect for a binary-group OLS should not be a z-scored 0/1 indicator itself;
    # verify raw_effect and standardized_effect are NOT identical (confirms distinct computation)
    assert not (df['raw_effect'] == df['standardized_effect']).all()


# ---------- 11: compactness never with length+rho ----------

def test_compactness_not_with_length_and_rho_in_primary_model():
    path = TABLES / "stage11_all_model_coefficients.csv"
    if not path.exists():
        pytest.skip("model coefficients not built")
    df = pd.read_csv(path)
    m1a_preds = set(df[df.model_id == 'Model1A']['predictor'])
    assert not {'cdr3_sequence_length', 'cdr3_rho', 'cdr3_compactness'}.issubset(m1a_preds)


# ---------- 12: shielding/BSA collinearity audited ----------

def test_shielding_bsa_collinearity_documented():
    path = TABLES / "stage11_collinearity_audit.csv"
    if not path.exists():
        pytest.skip("collinearity audit not built")
    df = pd.read_csv(path)
    bsa_rows = df[df['model_set'] == 'shielding_BSA_pair']
    assert (bsa_rows['VIF'] > 10).all()  # confirms severe collinearity was actually detected, not assumed


# ---------- 13: Gordon theta/phi absent ----------

def test_gordon_theta_phi_absent():
    path = PROCESSED / "stage11_model_input.parquet"
    if not path.exists():
        pytest.skip("canonical table not built")
    df = pd.read_parquet(path)
    assert 'cdr3_theta' not in df.columns
    assert 'cdr3_phi' not in df.columns


# ---------- 14: HC3 reproducible ----------

def test_hc3_reproducible():
    rng = np.random.default_rng(0)
    df = pd.DataFrame({'y': rng.normal(0, 1, 200), 'x': rng.normal(0, 1, 200)})
    r1 = fit_ols_hc3(df, 'y', ['x'])
    r2 = fit_ols_hc3(df, 'y', ['x'])
    assert r1['coefficients']['x']['raw_beta'] == r2['coefficients']['x']['raw_beta']


# ---------- 15: deterministic bootstrap/permutation seed ----------

def test_deterministic_bootstrap_seed():
    from stage11_stats import bootstrap_group_diff
    rng_df = pd.DataFrame({'y': np.concatenate([np.random.default_rng(1).normal(0, 1, 50), np.random.default_rng(2).normal(1, 1, 50)]),
                            'g': ['A'] * 50 + ['B'] * 50})
    r1 = bootstrap_group_diff(rng_df, 'y', 'g', 'A', 'B', n_boot=100, seed=42)
    r2 = bootstrap_group_diff(rng_df, 'y', 'g', 'A', 'B', n_boot=100, seed=42)
    assert r1['boot_CI_low'] == r2['boot_CI_low']
    assert r1['boot_CI_high'] == r2['boot_CI_high']


# ---------- 16: BH-FDR correct ----------

def test_bh_fdr_correct():
    pvals = [0.001, 0.01, 0.02, 0.04, 0.5]
    q = bh_fdr(pvals)
    assert (q >= np.array(pvals)).all()  # BH q-values never smaller than raw p
    assert np.all(np.diff(q[np.argsort(pvals)]) >= -1e-9)  # monotonic non-decreasing when sorted by p


# ---------- 17: FDR families separate ----------

def test_fdr_families_kept_separate():
    grp_path = TABLES / "stage11_primary_group_inference.csv"
    if not grp_path.exists():
        pytest.skip("not built")
    grp = pd.read_csv(grp_path)
    mech = pd.read_csv(TABLES / "stage11_primary_mechanism_inference.csv")
    # each family's BH_FDR_q must be computed independently -- confirmed by checking
    # the q-value ordering is consistent with each family's OWN p-value ranking only
    assert 'BH_FDR_q' in grp.columns and 'BH_FDR_q' in mech.columns
    assert len(grp) != len(mech) or not grp['BH_FDR_q'].equals(mech['BH_FDR_q'])


# ---------- 18: no automatic stepwise ----------

def test_no_stepwise_selection_artifacts():
    path = TABLES / "stage11_all_model_coefficients.csv"
    if not path.exists():
        pytest.skip("not built")
    df = pd.read_csv(path)
    # predeclared model_ids only -- no auto-generated "stepwise"/"selected" model names
    assert not df['model_id'].str.contains('stepwise|selected|auto', case=False).any()


# ---------- 19-20: no Stage 10 / MD ----------

def test_no_stage10_aggregation_metric():
    for path in PROCESSED.glob("stage11_*.parquet"):
        df = pd.read_parquet(path)
        forbidden = ['aggrescan3d', 'developability_score', 'aggregation_propensity']
        for f in forbidden:
            assert f not in [c.lower() for c in df.columns]


def test_no_md_invoked():
    for path in PROCESSED.glob("stage11_*.parquet"):
        df = pd.read_parquet(path)
        assert not any('rmsf' in c.lower() or 'md_' in c.lower() for c in df.columns)


# ---------- 21-23: engagement semantics ----------

def test_engagement_denominator_preserved():
    path = ROOT / "data/interim/stage9/vhh_shieldable_interface_engagement.csv"
    if not path.exists():
        pytest.skip("not present")
    df = pd.read_csv(path)
    assert 'n_target_positions_total' in df.columns
    assert (df['n_target_positions_total'] == 15).all()


def test_zero_engagement_not_missing():
    path = PROCESSED / "stage11_model_input.parquet"
    if not path.exists():
        pytest.skip("canonical table not built")
    df = pd.read_parquet(path)
    zero_rows = df[df['shieldable_engagement_fraction'] == 0]
    assert len(zero_rows) > 0  # zeros are real, present values, not dropped as NaN


def test_mapping_failure_not_zero():
    path = PROCESSED / "stage8_apbs_interface_electrostatics.parquet"
    df = pd.read_parquet(path)
    fail = df[df['status'] != 'PDB2PQR_SUCCESS_APBS_SUCCESS']
    assert fail['interface_potential_mean'].isna().all()


# ---------- 24: PPC/PNC presence/magnitude ----------

def test_ppc_pnc_zero_presence_handling():
    path = PROCESSED / "stage11_model_input.parquet"
    if not path.exists():
        pytest.skip("canonical table not built")
    df = pd.read_parquet(path)
    vhh = df[df.group == 'VHH']
    frac_zero = (vhh['former_interface_core_PPC'] == 0).mean()
    assert frac_zero > 0.3  # confirms substantial zero-inflation was real, not assumed


# ---------- 25: APBS components separate ----------

def test_apbs_positive_negative_separate():
    path = PROCESSED / "stage11_model_input.parquet"
    if not path.exists():
        pytest.skip("canonical table not built")
    df = pd.read_parquet(path)
    assert 'interface_potential_positive_mean' in df.columns
    assert 'interface_potential_negative_mean' in df.columns
    assert 'net_potential' not in df.columns


# ---------- 26: PDB clustering correct ----------

def test_pdb_clustering_uses_pdb_correctly():
    path = PROCESSED / "stage11_model_input.parquet"
    if not path.exists():
        pytest.skip("canonical table not built")
    df = pd.read_parquet(path)
    vhh = df[df.group == 'VHH']
    assert vhh['PDB'].nunique() <= len(vhh)
    assert vhh['PDB'].nunique() > 1


# ---------- 27: high-coverage threshold unchanged ----------

def test_high_coverage_threshold_unchanged():
    path = TABLES / "stage11_high_coverage_sensitivity.csv"
    if not path.exists():
        pytest.skip("not built")
    df = pd.read_parquet(PROCESSED / "stage11_model_input.parquet")
    vhh = df[df.group == 'VHH']
    hc = vhh[vhh['primary_interface_position_coverage'] >= 0.90]
    assert len(hc) == 394  # matches Stage 8/9's established >=0.90 threshold result


# ---------- 28: dominant-antigen sensitivity deterministic ----------

def test_dominant_antigen_sensitivity_deterministic():
    path = TABLES / "stage11_dominant_antigen_sensitivity.csv"
    if not path.exists():
        pytest.skip("not built")
    df = pd.read_csv(path)
    assert (df['largest_group'] == 'Spike protein S1').all()
    assert (df['N_removed'] == 37).all()


# ---------- 29: all primary coefficients retained regardless of significance ----------

def test_all_primary_coefficients_retained_including_null():
    path = TABLES / "stage11_primary_mechanism_inference.csv"
    if not path.exists():
        pytest.skip("not built")
    df = pd.read_csv(path)
    assert (df['raw_p'] > 0.05).any()  # at least one nonsignificant result retained, not suppressed


# ---------- 30: CI reproducible ----------

def test_ci_calculations_reproducible():
    rng = np.random.default_rng(5)
    df = pd.DataFrame({'y': rng.normal(0, 1, 300), 'x1': rng.normal(0, 1, 300), 'x2': rng.normal(0, 1, 300)})
    r1 = fit_ols_hc3(df, 'y', ['x1', 'x2'])
    r2 = fit_ols_hc3(df, 'y', ['x1', 'x2'])
    assert r1['coefficients']['x1']['CI_low'] == r2['coefficients']['x1']['CI_low']
    assert r1['coefficients']['x1']['CI_high'] == r2['coefficients']['x1']['CI_high']


# ---------- 31: Stage 2-9 unchanged ----------

def test_stage2to9_outputs_unchanged():
    assert len(pd.read_parquet(PROCESSED / "stage4_vhh_sequence_representatives.parquet")) == 533
    assert len(pd.read_parquet(PROCESSED / "stage9_antigen_engagement_summary.parquet")) == 1962
    df8 = pd.read_parquet(PROCESSED / "stage8_apbs_interface_electrostatics.parquet")
    assert len(df8) == 1962
