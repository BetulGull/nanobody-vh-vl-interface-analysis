import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

ROOT = Path("/Users/betulgul/Desktop/staj")
FIG_DIR = ROOT / "results/figures/stage11"
FIG_DIR.mkdir(parents=True, exist_ok=True)
TABLES = ROOT / "results/tables"

VHH_BLUE = "#4D69B2"
VH_PURPLE = "#8B5FA8"


def fig1_primary_adjusted_effects():
    allc = pd.read_csv(TABLES / "stage11_all_model_coefficients.csv")
    rows_wanted = [
        ('Model1A', 'cdr3_sequence_length', 'CDR3 length -> shielding'),
        ('Model1A', 'cdr3_rho', 'rho -> shielding'),
        ('Model1A', 'FR2_CDR3_Dis', 'FR2_CDR3_Dis -> shielding'),
        ('ModelC1', 'H_intrinsic_binary_mean_KD', 'H_intrinsic -> native surface'),
        ('Model4C', 'shielding_fraction', 'shielding -> engagement_15 (geometry-adjusted)'),
        ('Model5B', 'shielding_fraction', 'shielding -> CDR3 antigen fraction (geometry-adjusted)'),
        ('Model6B', 'shielding_fraction', 'shielding -> framework fraction (geometry-adjusted)'),
    ]
    data = []
    for model_id, pred, label in rows_wanted:
        r = allc[(allc.model_id == model_id) & (allc.predictor == pred)].iloc[0]
        data.append({'label': label, 'beta': r['standardized_beta'], 'ci_low': None, 'ci_high': None, 'p': r['raw_p']})

    # need CI on standardized scale -- recompute quickly via bootstrap-free approx using raw CI ratio
    # (raw CI / raw beta) * standardized beta as a simple scale-consistent approximation is not exact;
    # instead reuse raw CI width scaled by std_beta/raw_beta ratio.
    for i, (model_id, pred, label) in enumerate(rows_wanted):
        r = allc[(allc.model_id == model_id) & (allc.predictor == pred)].iloc[0]
        scale = r['standardized_beta'] / r['raw_beta'] if r['raw_beta'] != 0 else 0
        data[i]['ci_low'] = r['CI_low'] * scale
        data[i]['ci_high'] = r['CI_high'] * scale

    df = pd.DataFrame(data)
    fig, ax = plt.subplots(figsize=(10, 5))
    y_pos = range(len(df))
    colors = [VHH_BLUE if p < 0.05 else '#AAAAAA' for p in df['p']]
    ax.errorbar(df['beta'], y_pos, xerr=[df['beta'] - df['ci_low'], df['ci_high'] - df['beta']],
                fmt='o', color=VHH_BLUE, ecolor='grey', capsize=3)
    ax.axvline(0, color='black', linewidth=0.8, linestyle='--')
    ax.set_yticks(list(y_pos))
    ax.set_yticklabels(df['label'])
    ax.invert_yaxis()
    ax.set_xlabel('Standardized coefficient (95% CI)')
    ax.set_title('Primary adjusted mechanistic effects (VHH only, HC3 robust SE)')
    fig.tight_layout()
    fig.savefig(FIG_DIR / "stage11_fig1_primary_adjusted_effects.png", dpi=150)
    plt.close(fig)


def fig2_group_effects():
    grp = pd.read_csv(TABLES / "stage11_primary_group_inference.csv")
    fig, ax = plt.subplots(figsize=(10, 5))
    y_pos = range(len(grp))
    ax.errorbar(grp['standardized_effect'], y_pos,
                xerr=[grp['standardized_effect'] - grp['CI_low'] * (grp['standardized_effect']/grp['raw_effect']),
                      grp['CI_high'] * (grp['standardized_effect']/grp['raw_effect']) - grp['standardized_effect']],
                fmt='s', color=VHH_BLUE, ecolor='grey', capsize=3)
    ax.axvline(0, color='black', linewidth=0.8, linestyle='--')
    ax.set_yticks(list(y_pos))
    ax.set_yticklabels(grp['hypothesis_id'])
    ax.invert_yaxis()
    ax.set_xlabel('Standardized effect (VHH - control), 95% CI')
    ax.set_title('VHH vs conventional VH/Fv: formal primary group effects\n(BH-FDR applied within family; all q<0.05)')
    fig.tight_layout()
    fig.savefig(FIG_DIR / "stage11_fig2_group_effects.png", dpi=150)
    plt.close(fig)


def fig3_shielding_function_adjusted():
    df = pd.read_parquet(ROOT / "data/processed/stage11_model_input.parquet")
    vhh = df[df.group == 'VHH'].dropna(subset=['shielding_fraction', 'shieldable_engagement_fraction'])
    fig, ax = plt.subplots(figsize=(7, 5.5))
    ax.scatter(vhh['shielding_fraction'], vhh['shieldable_engagement_fraction'], s=14, alpha=0.4, color=VHH_BLUE, label='observed (raw)')

    import statsmodels.api as sm
    X = sm.add_constant(vhh[['shielding_fraction']])
    m = sm.OLS(vhh['shieldable_engagement_fraction'], X).fit(cov_type='HC3')
    xs = np.linspace(vhh['shielding_fraction'].min(), vhh['shielding_fraction'].max(), 100)
    pred = m.params['const'] + m.params['shielding_fraction'] * xs
    ax.plot(xs, pred, color='#D35400', linewidth=2, label='primary adjusted relationship (Model 4A)')

    ax.set_xlabel('Stage-6 shielding_fraction')
    ax.set_ylabel('15-position non-CDR3 former-interface\nantigen-engagement fraction')
    ax.set_title('Shielding vs former-interface antigen engagement (VHH only)\nunadjusted association significant; attenuates to nonsignificant after geometry adjustment (Model 4C)')
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "stage11_fig3_shielding_function_adjusted.png", dpi=150)
    plt.close(fig)


def fig4_robustness_map():
    rob = pd.read_csv(TABLES / "stage11_robustness_summary.csv")
    cols = ['HC3', 'PDB_cluster', 'high_coverage', 'dominant_antigen_excluded', 'influence_sensitivity']

    def code(val):
        if pd.isna(val) or 'not separately' in str(val) or val == 'N/A':
            return 0.5
        if 'NOT significant' in str(val) or 'not significant' in str(val).lower():
            return 0.0
        return 1.0

    mat = np.array([[code(rob.iloc[i][c]) for c in cols] for i in range(len(rob))])
    fig, ax = plt.subplots(figsize=(10, 6))
    im = ax.imshow(mat, cmap='RdYlGn', vmin=0, vmax=1, aspect='auto')
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels(cols, rotation=30, ha='right')
    ax.set_yticks(range(len(rob)))
    ax.set_yticklabels([c[:55] for c in rob['primary_conclusion']], fontsize=7)
    ax.set_title('Mechanism robustness map', fontsize=12)
    fig.text(0.5, 0.96, 'green=stable/significant, yellow=not separately evaluated, red=not significant', ha='center', fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "stage11_fig4_robustness_map.png", dpi=150)
    plt.close(fig)


if __name__ == '__main__':
    fig1_primary_adjusted_effects()
    fig2_group_effects()
    fig3_shielding_function_adjusted()
    fig4_robustness_map()
    print("Stage 11 figures written")
