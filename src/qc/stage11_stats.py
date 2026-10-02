"""Stage 11 — Core statistical helper functions (OLS/HC3, standardized coefficients,
bootstrap, permutation, PDB-cluster-robust SE). Predeclared model families only."""
import numpy as np
import pandas as pd
import statsmodels.api as sm

SEED = 42


def zscore(series):
    return (series - series.mean()) / series.std()


def fit_ols_hc3(data, outcome, predictors, group_col=None):
    """Fit OLS with HC3 robust SE. Returns dict with coefficients, CI, p, R2, N."""
    d = data[[outcome] + predictors + ([group_col] if group_col else [])].dropna()
    if group_col:
        # encode group as 0/1 if categorical passed as predictor name
        pass
    X = sm.add_constant(d[predictors])
    y = d[outcome]
    model = sm.OLS(y, X).fit(cov_type='HC3')
    result = {
        'N': len(d), 'R2': model.rsquared, 'adj_R2': model.rsquared_adj,
        'coefficients': {}
    }
    for p in predictors:
        result['coefficients'][p] = {
            'raw_beta': model.params[p], 'CI_low': model.conf_int().loc[p, 0],
            'CI_high': model.conf_int().loc[p, 1], 'p': model.pvalues[p],
        }
    result['model'] = model
    result['data'] = d
    return result


def fit_ols_hc3_standardized(data, outcome, predictors):
    """Same as fit_ols_hc3 but on z-scored outcome+predictors, for standardized betas."""
    d = data[[outcome] + predictors].dropna().copy()
    for c in [outcome] + predictors:
        d[c] = zscore(d[c])
    X = sm.add_constant(d[predictors])
    y = d[outcome]
    model = sm.OLS(y, X).fit(cov_type='HC3')
    return {p: model.params[p] for p in predictors}


def fit_ols_cluster(data, outcome, predictors, cluster_col):
    d = data[[outcome] + predictors + [cluster_col]].dropna()
    X = sm.add_constant(d[predictors])
    y = d[outcome]
    model = sm.OLS(y, X).fit(cov_type='cluster', cov_kwds={'groups': d[cluster_col]})
    return {p: {'raw_beta': model.params[p], 'CI_low': model.conf_int().loc[p, 0],
                'CI_high': model.conf_int().loc[p, 1], 'p': model.pvalues[p]} for p in predictors}


def bootstrap_group_diff(data, outcome, group_col, group_a, group_b, n_boot=2000, seed=SEED):
    rng = np.random.default_rng(seed)
    a = data[data[group_col] == group_a][outcome].dropna().values
    b = data[data[group_col] == group_b][outcome].dropna().values
    obs_diff = a.mean() - b.mean()
    diffs = np.empty(n_boot)
    for i in range(n_boot):
        ra = rng.choice(a, size=len(a), replace=True)
        rb = rng.choice(b, size=len(b), replace=True)
        diffs[i] = ra.mean() - rb.mean()
    ci_low, ci_high = np.percentile(diffs, [2.5, 97.5])
    return {'observed_diff': obs_diff, 'boot_CI_low': ci_low, 'boot_CI_high': ci_high, 'n_boot': n_boot}


def permutation_test_group_diff(data, outcome, group_col, group_a, group_b, n_perm=5000, seed=SEED):
    rng = np.random.default_rng(seed)
    sub = data[data[group_col].isin([group_a, group_b])][[outcome, group_col]].dropna()
    a = sub[sub[group_col] == group_a][outcome].values
    b = sub[sub[group_col] == group_b][outcome].values
    obs_diff = a.mean() - b.mean()
    combined = np.concatenate([a, b])
    n_a = len(a)
    perm_diffs = np.empty(n_perm)
    for i in range(n_perm):
        perm = rng.permutation(combined)
        perm_diffs[i] = perm[:n_a].mean() - perm[n_a:].mean()
    p_value = (np.sum(np.abs(perm_diffs) >= np.abs(obs_diff)) + 1) / (n_perm + 1)
    return {'observed_diff': obs_diff, 'perm_p': p_value, 'n_perm': n_perm}


def bh_fdr(pvalues):
    """Benjamini-Hochberg FDR correction. Returns array of q-values in original order."""
    pvals = np.asarray(pvalues, dtype=float)
    n = len(pvals)
    order = np.argsort(pvals)
    ranked = pvals[order]
    q = ranked * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    result = np.empty(n)
    result[order] = q
    return result


def compute_vif_cn(data, predictors):
    from statsmodels.stats.outliers_influence import variance_inflation_factor
    X = sm.add_constant(data[predictors].dropna())
    vifs = {c: variance_inflation_factor(X.values, i) for i, c in enumerate(X.columns) if c != 'const'}
    cn = np.linalg.cond(X.values)
    return vifs, cn
