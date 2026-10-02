"""
Stage 7 — Required figures (visualization gate: all three are justified by
statistically strong, pre-computed exploratory correlations already on file
in stage7_exploratory_correlations.csv; no figure is generated speculatively).
"""
import sys
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path("/Users/betulgul/Desktop/staj")
FIG_DIR = ROOT / "results/figures/stage7"
FIG_DIR.mkdir(parents=True, exist_ok=True)
TABLES = ROOT / "results/tables"


def fig1_packing_vs_shielding():
    df = pd.read_csv(TABLES / "stage7_stage4to7_combined.csv")
    sub = df[['nonadj_BSA_pair_total', 'shielding_fraction']].dropna()
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.scatter(sub['nonadj_BSA_pair_total'], sub['shielding_fraction'], s=12, alpha=0.5, color='#2b6cb0')
    corr = sub['nonadj_BSA_pair_total'].corr(sub['shielding_fraction'])
    ax.set_xlabel('BSA_pair_total (pairwise mutual buried surface area, Å²)\nCDR3 vs NONADJACENT_PACKING_TARGET')
    ax.set_ylabel('Stage 6 shielding_fraction')
    ax.set_title(f'CDR3-target packing vs CDR3 shielding (n={len(sub)}, Pearson r={corr:.2f})')
    fig.tight_layout()
    fig.savefig(FIG_DIR / "stage7_fig1_packing_vs_shielding.png", dpi=150)
    plt.close(fig)


def fig2_position_pair_heatmap():
    pair_map = pd.read_csv(TABLES / "stage7_position_pair_interaction_map.csv")
    cdr3_positions = sorted(pair_map['cdr3_imgt_position'].unique())
    target_positions = sorted(pair_map['target_imgt_position'].unique())
    mat = np.zeros((len(cdr3_positions), len(target_positions)))
    ci = {p: i for i, p in enumerate(cdr3_positions)}
    ti = {p: i for i, p in enumerate(target_positions)}
    for _, row in pair_map.iterrows():
        mat[ci[row['cdr3_imgt_position']], ti[row['target_imgt_position']]] = row['interaction_frequency']

    fig, ax = plt.subplots(figsize=(7, 6))
    im = ax.imshow(mat, cmap='viridis', aspect='auto')
    ax.set_xticks(range(len(target_positions)))
    ax.set_xticklabels(target_positions, rotation=90)
    ax.set_yticks(range(len(cdr3_positions)))
    ax.set_yticklabels(cdr3_positions)
    ax.set_xlabel('Target IMGT position (NONADJACENT_PACKING_TARGET)')
    ax.set_ylabel('CDR3 IMGT position')
    ax.set_title('CDR3 x target position-pair interaction frequency\n(fraction of 533 VHH instances with >=1 specific noncovalent contact)')
    fig.colorbar(im, ax=ax, label='Interaction frequency')
    fig.tight_layout()
    fig.savefig(FIG_DIR / "stage7_fig2_position_pair_heatmap.png", dpi=150)
    plt.close(fig)


def fig3_composition_by_quartile():
    df = pd.read_csv(TABLES / "stage7_stage4to7_combined.csv")
    long_df = pd.read_csv(TABLES / "stage7_interaction_records_long.csv", low_memory=False)

    df = df.dropna(subset=['shielding_fraction']).copy()
    df['shielding_quartile'] = pd.qcut(df['shielding_fraction'], 4, labels=['Q1_low', 'Q2', 'Q3', 'Q4_high'])
    inst_to_q = dict(zip(df['INSTANCE'], df['shielding_quartile']))

    primary = long_df[long_df['include_primary_noncovalent']].copy()
    primary['shielding_quartile'] = primary['INSTANCE'].map(inst_to_q)
    primary = primary.dropna(subset=['shielding_quartile'])

    label_cols = ['hbond', 'weak_hbond', 'ionic', 'aromatic', 'hydrophobic', 'carbonyl', 'polar', 'weak_polar', 'vdw', 'vdw_clash']
    comp = {q: {lbl: 0 for lbl in label_cols} for q in df['shielding_quartile'].cat.categories}
    for _, row in primary.iterrows():
        labels = row['contact_labels'].split(';')
        for lbl in labels:
            if lbl in comp[row['shielding_quartile']]:
                comp[row['shielding_quartile']][lbl] += 1

    comp_df = pd.DataFrame(comp).T
    comp_frac = comp_df.div(comp_df.sum(axis=1), axis=0)

    fig, ax = plt.subplots(figsize=(8, 5))
    bottom = np.zeros(len(comp_frac))
    colors = plt.cm.tab10(np.linspace(0, 1, len(label_cols)))
    for lbl, color in zip(label_cols, colors):
        ax.bar(comp_frac.index, comp_frac[lbl], bottom=bottom, label=lbl, color=color)
        bottom += comp_frac[lbl].values
    ax.set_ylabel('Fraction of specific-interaction atom records')
    ax.set_xlabel('Shielding-fraction quartile (predeclared, descriptive grouping)')
    ax.set_title('Interaction-type composition by CDR3-shielding quartile')
    ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "stage7_fig3_composition_by_quartile.png", dpi=150)
    plt.close(fig)


if __name__ == '__main__':
    fig1_packing_vs_shielding()
    fig2_position_pair_heatmap()
    fig3_composition_by_quartile()
    print("Figures written to", FIG_DIR)
