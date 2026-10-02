import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
from pathlib import Path

ROOT = Path("/Users/betulgul/Desktop/staj")
FIG_DIR = ROOT / "results/figures/stage12"
FIG_DIR.mkdir(parents=True, exist_ok=True)
TABLES = ROOT / "results/tables"

VHH_BLUE = "#4D69B2"
VH_PURPLE = "#8B5FA8"
COUPLED_ORANGE = "#D35400"
NEUTRAL_GREY = "#888888"


def fig_final_mechanism():
    edges = pd.read_csv(TABLES / "stage12_mechanism_edges.csv")
    edges = edges[edges['final_include'] == True]

    fig, ax = plt.subplots(figsize=(11, 9))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 12)
    ax.axis('off')

    boxes = {
        'ABSENT_VL': (5, 11, 'Absent VL'),
        'FORMER_INTERFACE': (5, 9.7, 'Empirical former VH-VL surface\n(28 IMGT positions, Stage 3)'),
        'SEQ_ADAPT': (2.3, 8.2, 'Sequence-level\nchemical adaptation\n(Stage 4)'),
        'CDR3_MODULE': (7.7, 8.2, 'CDR3 geometric accommodation\n(shielding + packing/BSA)\n[COUPLED, Stages 6-7]'),
        'NATIVE_SURF': (2.3, 6.5, 'Native exposed\nsurface chemistry\n(Stage 5)'),
        'GEOM_AXIS': (7.7, 6.5, 'CDR3 geometry axis\n(length, rho, FR2_CDR3_Dis)'),
        'PATCH_PHENO': (5, 4.8, 'Surface patch /\nelectrostatic phenotype\n(Stage 8: PSH, patches, APBS)'),
        'ANTIGEN_ORG': (5, 3.0, 'Antigen-interface organization\n(Stage 9: framework/CDR3\nfraction, former-interface\nengagement)'),
        'GROUP_DIFF': (5, 1.3, 'Formal VHH vs conventional\ngroup differences\n(Stage 11, all 7 hypotheses q<1e-6)'),
    }

    for key, (x, y, label) in boxes.items():
        color = VHH_BLUE if key != 'CDR3_MODULE' else COUPLED_ORANGE
        box = FancyBboxPatch((x - 1.15, y - 0.55), 2.3, 1.1, boxstyle="round,pad=0.05",
                              edgecolor=color, facecolor='white', linewidth=1.8)
        ax.add_patch(box)
        ax.text(x, y, label, ha='center', va='center', fontsize=8)

    def arrow(a, b, style='solid', color='black', label=''):
        xa, ya, _ = boxes[a]
        xb, yb, _ = boxes[b]
        ls = {'solid': '-', 'dashed': '--', 'double': '-'}[style]
        lw = 2.5 if style == 'double' else 1.6
        arr = FancyArrowPatch((xa, ya - 0.55), (xb, yb + 0.55), arrowstyle='-|>', mutation_scale=15,
                               linewidth=lw, linestyle=ls, color=color, connectionstyle="arc3,rad=0.0")
        ax.add_patch(arr)
        if label:
            mx, my = (xa + xb) / 2, (ya + yb) / 2
            ax.text(mx, my, label, fontsize=6.5, ha='center', color=color, style='italic')

    arrow('ABSENT_VL', 'FORMER_INTERFACE', 'solid')
    arrow('FORMER_INTERFACE', 'SEQ_ADAPT', 'solid')
    arrow('FORMER_INTERFACE', 'CDR3_MODULE', 'solid')
    arrow('SEQ_ADAPT', 'NATIVE_SURF', 'solid', VHH_BLUE, 'H_intrinsic->surface\nR2=0.52, p=2e-128')
    arrow('CDR3_MODULE', 'GEOM_AXIS', 'double', COUPLED_ORANGE, 'shielding<->BSA\nr=0.97 (inseparable)')
    arrow('NATIVE_SURF', 'PATCH_PHENO', 'solid', VHH_BLUE)
    arrow('CDR3_MODULE', 'PATCH_PHENO', 'solid', COUPLED_ORANGE, 'shielding indep.\np=2e-17')
    arrow('GEOM_AXIS', 'ANTIGEN_ORG', 'solid', VHH_BLUE, 'CDR3 length indep.\np<1e-10')
    arrow('CDR3_MODULE', 'ANTIGEN_ORG', 'dashed', NEUTRAL_GREY, 'shielding: descriptive only\n(NOT indep. after adj., p=0.07-0.24)')
    arrow('PATCH_PHENO', 'ANTIGEN_ORG', 'solid', VHH_BLUE, 'PSH q=3.5e-10')
    arrow('ANTIGEN_ORG', 'GROUP_DIFF', 'solid')
    arrow('PATCH_PHENO', 'GROUP_DIFF', 'solid')

    legend_elems = [
        mpatches.Patch(facecolor='none', edgecolor='black', label='SOLID: robust adjusted / directly measured'),
        mpatches.Patch(facecolor='none', edgecolor=NEUTRAL_GREY, label='DASHED: descriptive only, not independent after adjustment'),
        mpatches.Patch(facecolor='none', edgecolor=COUPLED_ORANGE, label='DOUBLE LINE: highly coupled, not separable (VIF>10)'),
    ]
    ax.legend(handles=legend_elems, loc='lower left', fontsize=7, bbox_to_anchor=(-0.02, -0.02))
    ax.set_title('Integrated structural-to-functional mechanism (VHH viability without VL)\nEvery edge maps to stage12_mechanism_edges.csv', fontsize=11)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "stage12_fig_final_mechanism.png", dpi=150)
    plt.close(fig)


def fig_position_level_integration():
    df = pd.read_csv(ROOT / "data/processed/stage9_former_interface_antigen_position_summary.csv")
    df = df.sort_values('imgt_position')

    fig, ax1 = plt.subplots(figsize=(13, 5.5))
    x = range(len(df))
    ax1.bar([i - 0.2 for i in x], df['VHH_antigen_contact_frequency'], width=0.4, color=VHH_BLUE, label='VHH antigen-contact freq')
    ax1.bar([i + 0.2 for i in x], df['VH_antigen_contact_frequency'], width=0.4, color=VH_PURPLE, label='VH antigen-contact freq')

    tetrad_positions = {42, 49, 50, 52}
    xticklabels = []
    for _, row in df.iterrows():
        lbl = str(int(row['imgt_position']))
        if row['imgt_position'] in tetrad_positions:
            lbl = f"*{lbl}*"
        xticklabels.append(lbl)
    ax1.set_xticks(list(x))
    ax1.set_xticklabels(xticklabels, rotation=90, fontsize=8)
    ax1.set_ylabel('Antigen-contact frequency')
    ax1.set_xlabel('IMGT position (Stage-3 empirical former interface; *tetrad*)')

    ax2 = ax1.twinx()
    ax2.plot(x, df['stage3_interface_frequency_5'], color='black', marker='o', markersize=3, linewidth=1, label='Stage-3 VL burial freq')
    ax2.set_ylabel('Stage-3 interface_frequency_5 (VL burial)')

    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper left', fontsize=8)
    ax1.set_title('Final position-level integration: VL burial vs VHH/VH antigen engagement\n(28-position empirical former interface; tetrad positions marked)')
    fig.tight_layout()
    fig.savefig(FIG_DIR / "stage12_fig_position_level_integration.png", dpi=150)
    plt.close(fig)


def fig_final_effect_evidence():
    grp = pd.read_csv(TABLES / "stage11_primary_group_inference.csv")
    mech = pd.read_csv(TABLES / "stage11_primary_mechanism_inference.csv")

    rows = []
    for _, r in grp.iterrows():
        rows.append({'label': f"[GROUP] {r['hypothesis_id']}", 'beta': r['standardized_effect'], 'p': r['raw_p'],
                     'ci_low': r['CI_low'] * (r['standardized_effect'] / r['raw_effect']) if r['raw_effect'] != 0 else 0,
                     'ci_high': r['CI_high'] * (r['standardized_effect'] / r['raw_effect']) if r['raw_effect'] != 0 else 0})
    key_mech = ['ModelC2_shielding_fraction', 'ModelD2_shielding_fraction', 'Model4A_shielding_fraction',
                'Model4C_shielding_fraction', 'Model5A_shielding_fraction', 'Model5B_shielding_fraction',
                'Model6A_shielding_fraction', 'Model6B_shielding_fraction']
    for hid in key_mech:
        r = mech[mech['hypothesis_id'] == hid].iloc[0]
        scale = r['standardized_beta'] / r['raw_beta'] if r['raw_beta'] != 0 else 0
        rows.append({'label': f"[MECH] {hid}", 'beta': r['standardized_beta'], 'p': r['raw_p'],
                     'ci_low': r['CI_low'] * scale, 'ci_high': r['CI_high'] * scale})

    df = pd.DataFrame(rows)
    fig, ax = plt.subplots(figsize=(9, 8))
    y_pos = range(len(df))
    colors = [VHH_BLUE if p < 0.05 else '#CC3333' for p in df['p']]
    for i, (y, row) in enumerate(zip(y_pos, df.itertuples())):
        ax.errorbar(row.beta, y, xerr=[[row.beta - row.ci_low], [row.ci_high - row.beta]], fmt='o', color=colors[i], ecolor='grey', capsize=3)
    ax.axvline(0, color='black', linewidth=0.8, linestyle='--')
    ax.set_yticks(list(y_pos))
    ax.set_yticklabels(df['label'], fontsize=7)
    ax.invert_yaxis()
    ax.set_xlabel('Standardized effect (95% CI)')
    ax.set_title('Final adjusted-effect evidence figure\nBlue = significant (p<0.05); Red = not significant (attenuated/null results shown, not hidden)', fontsize=9)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "stage12_fig_final_effect_evidence.png", dpi=150)
    plt.close(fig)


def fig_vhh_vs_conventional_summary():
    grp = pd.read_csv(TABLES / "stage11_primary_group_inference.csv")
    fig, ax = plt.subplots(figsize=(8, 5.5))
    y_pos = range(len(grp))
    colors_bar = [VHH_BLUE if v > 0 else VH_PURPLE for v in grp['standardized_effect']]
    ax.barh(list(y_pos), grp['standardized_effect'], color=colors_bar)
    ax.axvline(0, color='black', linewidth=0.8)
    ax.set_yticks(list(y_pos))
    labels = ['H_intrinsic', 'Native surface', 'PSH', 'Largest hydrophobic patch',
              'Framework antigen fraction', '15-position engagement', 'CDR3 antigen fraction']
    ax.set_yticklabels(labels)
    ax.invert_yaxis()
    ax.set_xlabel('Standardized effect (VHH - conventional VH/Fv)')
    ax.set_title('Strongest validated VHH-vs-conventional group differences\n(all BH-FDR q<1e-6; blue=VHH higher, purple=VHH lower)')
    fig.tight_layout()
    fig.savefig(FIG_DIR / "stage12_fig_vhh_vs_conventional_summary.png", dpi=150)
    plt.close(fig)


if __name__ == '__main__':
    fig_final_mechanism()
    fig_position_level_integration()
    fig_final_effect_evidence()
    fig_vhh_vs_conventional_summary()
    print("Stage 12 figures written")
