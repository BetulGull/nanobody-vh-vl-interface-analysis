import sys
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path("/Users/betulgul/Desktop/staj")
FIG_DIR = ROOT / "results/figures/stage9"
FIG_DIR.mkdir(parents=True, exist_ok=True)
TABLES = ROOT / "results/tables"
PROCESSED = ROOT / "data/processed"

VHH_BLUE = "#4D69B2"
VH_PURPLE = "#8B5FA8"


def fig1_framework_contribution():
    df = pd.read_parquet(PROCESSED / "stage9_antigen_engagement_summary.parquet")
    vhh = df[df['group'] == 'VHH']['framework_contact_fraction'].dropna()
    vh = df[df['group'] == 'VH']['framework_contact_fraction'].dropna()
    fig, ax = plt.subplots(figsize=(6, 5))
    bp = ax.boxplot([vhh, vh], tick_labels=['VHH', 'conventional VH'], patch_artist=True, showfliers=False)
    for patch, color in zip(bp['boxes'], [VHH_BLUE, VH_PURPLE]):
        patch.set_facecolor(color)
    ax.set_ylabel('Framework antigen-contact fraction')
    ax.set_title(f'Framework contribution to antigen engagement\n(N={len(vhh)} VHH, N={len(vh)} conventional VH)')
    fig.tight_layout()
    fig.savefig(FIG_DIR / "stage9_fig1_framework_contribution.png", dpi=150)
    plt.close(fig)


def fig2_region_composition():
    df = pd.read_parquet(PROCESSED / "stage9_antigen_engagement_summary.parquet")
    vhh = df[df['group'] == 'VHH']
    vh = df[df['group'] == 'VH']

    vhh_comp = vhh[['CDR1_contact_residue_count', 'CDR2_contact_residue_count', 'CDR3_contact_residue_count', 'framework_contact_residue_count']].mean()
    vhh_comp.index = ['CDR1', 'CDR2', 'CDR3', 'Framework']
    vhh_frac = vhh_comp / vhh_comp.sum()

    vh_vh_comp = vh[['CDR1_contact_residue_count', 'CDR2_contact_residue_count', 'CDR3_contact_residue_count']].mean()
    vh_fw = vh['framework_contact_residue_count'].mean() - 0  # VH-only framework portion within total
    vh_vl_total = vh['VL_contact_count'].mean()
    # approximate VH framework vs VH CDR vs VL total split for a compact composition figure
    vh_labels = ['VH CDR1', 'VH CDR2', 'VH CDR3', 'VH FW', 'VL (total)']
    vh_vals = [vh_vh_comp['CDR1_contact_residue_count'], vh_vh_comp['CDR2_contact_residue_count'],
               vh_vh_comp['CDR3_contact_residue_count'], vh_fw, vh_vl_total]
    vh_frac2 = np.array(vh_vals) / sum(vh_vals)

    fig, axes = plt.subplots(1, 2, figsize=(11, 5))
    axes[0].pie(vhh_frac, labels=vhh_frac.index, autopct='%1.0f%%', colors=['#7FB3D5', '#5499C7', '#2874A6', VHH_BLUE])
    axes[0].set_title('VHH paratope composition')
    axes[1].pie(vh_frac2, labels=vh_labels, autopct='%1.0f%%', colors=['#D2B4DE', '#C39BD3', '#A569BD', VH_PURPLE, '#F5B7B1'])
    axes[1].set_title('Conventional Fv paratope composition')
    fig.suptitle('Region-level antigen-contact residue composition (mean fractions)')
    fig.tight_layout()
    fig.savefig(FIG_DIR / "stage9_fig2_region_composition.png", dpi=150)
    plt.close(fig)


def fig3_former_interface_position_map():
    df = pd.read_csv(PROCESSED / "stage9_former_interface_antigen_position_summary.csv")
    fig, ax1 = plt.subplots(figsize=(12, 5))
    x = range(len(df))
    ax1.bar([i - 0.2 for i in x], df['VHH_antigen_contact_frequency'], width=0.4, color=VHH_BLUE, label='VHH antigen-contact freq')
    ax1.bar([i + 0.2 for i in x], df['VH_antigen_contact_frequency'], width=0.4, color=VH_PURPLE, label='VH antigen-contact freq')
    ax1.set_xticks(list(x))
    ax1.set_xticklabels(df['imgt_position'].astype(str), rotation=90)
    ax1.set_ylabel('Antigen-contact frequency')
    ax1.set_xlabel('IMGT position (Stage-3 empirical former interface)')

    ax2 = ax1.twinx()
    ax2.plot(x, df['stage3_interface_frequency_5'], color='black', marker='o', markersize=3, linewidth=1, label='Stage-3 VL burial freq')
    ax2.set_ylabel('Stage-3 interface_frequency_5')

    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper left', fontsize=8)
    ax1.set_title('Former-interface position-level integration: VL burial vs antigen engagement (VHH vs VH)')
    fig.tight_layout()
    fig.savefig(FIG_DIR / "stage9_fig3_former_interface_position_map.png", dpi=150)
    plt.close(fig)


def fig4_shielding_vs_engagement():
    merged = pd.read_csv(ROOT / "data/interim/stage9/vhh_shielding_engagement_merged.csv")
    sub = merged[['shielding_fraction', 'shieldable_engagement_fraction']].dropna()
    corr = sub['shielding_fraction'].corr(sub['shieldable_engagement_fraction'])
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.scatter(sub['shielding_fraction'], sub['shieldable_engagement_fraction'], s=12, alpha=0.5, color=VHH_BLUE)
    ax.set_xlabel('Stage-6 shielding_fraction')
    ax.set_ylabel('15-position non-CDR3 former-interface\nantigen-engagement fraction')
    title = 'Shielding vs former-interface antigen engagement (VHH only)' if abs(corr) < 0.3 else \
        ('Inverse relationship: shielding vs antigen engagement (VHH only)' if corr < 0 else 'Positive relationship: shielding vs antigen engagement (VHH only)')
    ax.set_title(f'{title}\n(N={len(sub)}, Pearson r={corr:.2f})')
    fig.tight_layout()
    fig.savefig(FIG_DIR / "stage9_fig4_shielding_vs_engagement.png", dpi=150)
    plt.close(fig)


if __name__ == '__main__':
    fig1_framework_contribution()
    fig2_region_composition()
    fig3_former_interface_position_map()
    fig4_shielding_vs_engagement()
    print("Stage 9 figures written")
