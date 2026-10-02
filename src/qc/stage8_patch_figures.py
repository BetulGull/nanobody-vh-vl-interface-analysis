import sys
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path("/Users/betulgul/Desktop/staj")
FIG_DIR = ROOT / "results/figures/stage8"
FIG_DIR.mkdir(parents=True, exist_ok=True)
TABLES = ROOT / "results/tables"
PROCESSED = ROOT / "data/processed"

VHH_BLUE = "#4D69B2"
VH_PURPLE = "#8B5FA8"


def fig1_hydrophobic_patch_burden():
    df = pd.read_parquet(PROCESSED / "stage8_surface_patch_descriptors.parquet")
    df = df[df['status'] == 'OK']
    fig, axes = plt.subplots(1, 2, figsize=(10, 5))
    for ax, col, title in [(axes[0], 'former_interface_vicinity_PSH', 'Vicinity PSH'),
                             (axes[1], 'former_interface_core_PSH', 'Core PSH')]:
        vhh = df[df['group'] == 'VHH'][col].dropna()
        vh = df[df['group'] == 'VH'][col].dropna()
        bp = ax.boxplot([vhh, vh], tick_labels=['VHH', 'VH'], patch_artist=True, showfliers=False)
        for patch, color in zip(bp['boxes'], [VHH_BLUE, VH_PURPLE]):
            patch.set_facecolor(color)
        ax.set_title(title)
        ax.set_ylabel('PSH (adapted former-interface)')
    fig.suptitle('Former-interface hydrophobic patch burden: VHH vs VH (VL-removed geometric reference)')
    fig.tight_layout()
    fig.savefig(FIG_DIR / "stage8_fig1_hydrophobic_patch_burden.png", dpi=150)
    plt.close(fig)


def fig2_charge_patches():
    df = pd.read_parquet(PROCESSED / "stage8_surface_patch_descriptors.parquet")
    df = df[df['status'] == 'OK']
    fig, axes = plt.subplots(1, 2, figsize=(9, 5))
    for ax, col, title in [(axes[0], 'former_interface_core_PPC', 'PPC (positive patch)'),
                             (axes[1], 'former_interface_core_PNC', 'PNC (negative patch)')]:
        vhh = df[df['group'] == 'VHH'][col].dropna()
        vh = df[df['group'] == 'VH'][col].dropna()
        bp = ax.boxplot([vhh, vh], tick_labels=['VHH', 'VH'], patch_artist=True, showfliers=False)
        for patch, color in zip(bp['boxes'], [VHH_BLUE, VH_PURPLE]):
            patch.set_facecolor(color)
        ax.set_title(title)
        ax.set_ylabel('Patch score')
    fig.suptitle('Former-interface discrete charge patches: VHH vs VH (kept separate, never net)')
    fig.tight_layout()
    fig.savefig(FIG_DIR / "stage8_fig2_charge_patches.png", dpi=150)
    plt.close(fig)


def fig3_patch_topology():
    df = pd.read_parquet(PROCESSED / "stage8_patch_topology.parquet")
    df = df[df['status'] == 'OK']
    col = 'local_thr4.0_largest_patch_residue_count'
    fig, ax = plt.subplots(figsize=(6, 5))
    vhh = df[df['group'] == 'VHH'][col].dropna()
    vh = df[df['group'] == 'VH'][col].dropna()
    bp = ax.boxplot([vhh, vh], tick_labels=['VHH', 'VH'], patch_artist=True, showfliers=False)
    for patch, color in zip(bp['boxes'], [VHH_BLUE, VH_PURPLE]):
        patch.set_facecolor(color)
    ax.set_ylabel('Largest hydrophobic patch (residue count)')
    ax.set_title('PROJECT-SPECIFIC surface patch topology (local, 4.0A threshold)\nNot TAP/TNP PSH')
    fig.tight_layout()
    fig.savefig(FIG_DIR / "stage8_fig3_patch_topology.png", dpi=150)
    plt.close(fig)


if __name__ == '__main__':
    fig1_hydrophobic_patch_burden()
    fig2_charge_patches()
    fig3_patch_topology()
    print("Patch figures written")
