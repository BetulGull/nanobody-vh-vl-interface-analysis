"""
Stage 12 — Panel A (full-cohort PSH histogram, Raybould Fig 2A style) and
final multi-panel assembly (A-E).
"""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.gridspec import GridSpec
from pathlib import Path
from PIL import Image


def autocrop_white(png_path, pad=15):
    """Crop uniform white margins from a rendered PNG, keeping a small pad."""
    im = Image.open(png_path).convert("RGB")
    arr = np.array(im)
    non_white = np.any(arr < 250, axis=2)
    rows = np.any(non_white, axis=1)
    cols = np.any(non_white, axis=0)
    if not rows.any() or not cols.any():
        return
    rmin, rmax = np.where(rows)[0][[0, -1]]
    cmin, cmax = np.where(cols)[0][[0, -1]]
    rmin = max(0, rmin - pad); rmax = min(arr.shape[0], rmax + pad)
    cmin = max(0, cmin - pad); cmax = min(arr.shape[1], cmax + pad)
    im.crop((cmin, rmin, cmax, rmax)).save(png_path)

ROOT = Path("/Users/betulgul/Desktop/staj")
FIG_DIR = ROOT / "results/figures/stage12"

VHH_BLUE = "#4D69B2"
VH_PURPLE = "#8B5FA8"


def panel_A():
    df = pd.read_parquet(ROOT / "data/processed/stage8_surface_patch_descriptors.parquet")
    df = df[df['status'] == 'OK']
    vhh = df[df.group == 'VHH']['former_interface_vicinity_PSH']
    vh = df[df.group == 'VH']['former_interface_vicinity_PSH']
    n_vhh, n_vh = len(vhh), len(vh)

    bins = np.linspace(min(vhh.min(), vh.min()), max(vhh.max(), vh.max()), 31)

    fig, ax = plt.subplots(figsize=(5.2, 4.6))
    ax.hist(vhh, bins=bins, weights=np.ones(n_vhh) / n_vhh, color=VHH_BLUE, alpha=0.55, label='VHH', edgecolor='none')
    ax.hist(vh, bins=bins, weights=np.ones(n_vh) / n_vh, color=VH_PURPLE, alpha=0.55, label='conventional VH (VL-removed geometric reference)', edgecolor='none')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlabel('Former-interface vicinity PSH score\n(Kyte-Doolittle)')
    ax.set_ylabel('Proportion')
    ax.legend(frameon=False, fontsize=8, loc='upper right')
    ax.text(0.02, 0.98, 'A', transform=ax.transAxes, fontsize=16, fontweight='bold', va='top')
    fig.tight_layout()
    fig.savefig(FIG_DIR / "stage12_panel_A_psh_histogram.png", dpi=300)
    plt.close(fig)
    return n_vhh, n_vh


def colorbar_only():
    fig, ax = plt.subplots(figsize=(6.2, 1.5))
    fig.subplots_adjust(bottom=0.42, top=0.82, left=0.12, right=0.88)
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
        'raybould', [(0.02, 0.02, 0.55), (0.0, 0.55, 0.75), (0.0, 0.65, 0.55), (0.30, 0.85, 0.15)])
    norm = matplotlib.colors.Normalize(vmin=1.0, vmax=2.0)
    cb = matplotlib.colorbar.ColorbarBase(ax, cmap=cmap, norm=norm, orientation='horizontal',
                                            ticks=[1.0, 1.33, 1.66, 2.0])
    cb.ax.set_xticklabels(['1.00\nHydrophilic', '1.33', '1.66', '2.00\nHydrophobic'], fontsize=9)
    cb.set_label('Normalized Kyte-Doolittle hydropathy (frozen Stage-8 scale)', fontsize=8)
    fig.savefig(FIG_DIR / "stage12_panel_colorbar.png", dpi=300)
    plt.close(fig)


def assemble(n_vhh, n_vh):
    exemplars = pd.read_csv(ROOT / "results/structures/stage12_raybould_style_high_psh_exemplars.csv")
    vhh_ex = exemplars[exemplars.group == 'VHH'].iloc[0]
    vh_ex = exemplars[exemplars.group == 'VH'].iloc[0]

    for pname in ["stage12_panel_B_structural_reference.png", "stage12_panel_C_vhh_surface.png",
                  "stage12_panel_D_vh_surface.png", "stage12_panel_E_vhh_patch.png", "stage12_panel_E_vh_patch.png"]:
        autocrop_white(FIG_DIR / pname)

    fig = plt.figure(figsize=(14, 14.5))
    gs = GridSpec(4, 2, figure=fig, height_ratios=[1.15, 1.3, 1.3, 0.5], hspace=0.35, wspace=0.05)

    axA = fig.add_subplot(gs[0, 0])
    imgA = mpimg.imread(FIG_DIR / "stage12_panel_A_psh_histogram.png")
    axA.imshow(imgA)
    axA.axis('off')

    axB = fig.add_subplot(gs[0, 1])
    imgB = mpimg.imread(FIG_DIR / "stage12_panel_B_structural_reference.png")
    axB.imshow(imgB)
    axB.axis('off')
    axB.text(0.02, 0.98, 'B', transform=axB.transAxes, fontsize=16, fontweight='bold', va='top', color='black')
    axB.set_title('VH-VL native pairing\n(orange = frozen 28-position former interface;\npink = VL, translucent)', fontsize=9)

    axC = fig.add_subplot(gs[1, 0])
    imgC = mpimg.imread(FIG_DIR / "stage12_panel_C_vhh_surface.png")
    axC.imshow(imgC)
    axC.axis('off')
    axC.text(0.02, 0.98, 'C', transform=axC.transAxes, fontsize=16, fontweight='bold', va='top')
    axC.set_title(f"VHH exemplar (PSH 95th-percentile nearest match)\n{vhh_ex['PDB']} / {vhh_ex['INSTANCE']}, PSH={vhh_ex['PSH']:.1f}", fontsize=9)

    axD = fig.add_subplot(gs[1, 1])
    imgD = mpimg.imread(FIG_DIR / "stage12_panel_D_vh_surface.png")
    axD.imshow(imgD)
    axD.axis('off')
    axD.text(0.02, 0.98, 'D', transform=axD.transAxes, fontsize=16, fontweight='bold', va='top')
    axD.set_title(f"conventional VH exemplar (PSH 95th-percentile nearest match)\n{vh_ex['PDB']} / {vh_ex['INSTANCE']}, PSH={vh_ex['PSH']:.1f} (VL hidden for surface viz)", fontsize=9)

    axE1 = fig.add_subplot(gs[2, 0])
    imgE1 = mpimg.imread(FIG_DIR / "stage12_panel_E_vhh_patch.png")
    axE1.imshow(imgE1)
    axE1.axis('off')
    axE1.text(0.02, 0.98, 'E', transform=axE1.transAxes, fontsize=16, fontweight='bold', va='top')
    axE1.set_title('VHH: largest connected hydrophobic patch (black outline)', fontsize=9)

    axE2 = fig.add_subplot(gs[2, 1])
    imgE2 = mpimg.imread(FIG_DIR / "stage12_panel_E_vh_patch.png")
    axE2.imshow(imgE2)
    axE2.axis('off')
    axE2.set_title('VH: largest connected hydrophobic patch (black outline)', fontsize=9)

    axCB = fig.add_subplot(gs[3, :])
    imgCB = mpimg.imread(FIG_DIR / "stage12_panel_colorbar.png")
    axCB.imshow(imgCB)
    axCB.axis('off')

    fig.suptitle(f'Former-interface physicochemical surface phenotype: VHH (N={n_vhh}) vs conventional VH (N={n_vh})\n'
                 'Rendering logic adapted from Raybould et al. 2019 Fig. 2 (CDR-vicinity target) -> empirical former VH-VL interface target',
                 fontsize=11, y=0.995)

    fig.savefig(FIG_DIR / "stage12_fig_raybould_style_vhh_vs_vh_psh_surface.png", dpi=300, bbox_inches='tight')
    fig.savefig(FIG_DIR / "stage12_fig_raybould_style_vhh_vs_vh_psh_surface.pdf", bbox_inches='tight')
    try:
        fig.savefig(FIG_DIR / "stage12_fig_raybould_style_vhh_vs_vh_psh_surface.tiff", dpi=300, bbox_inches='tight')
    except Exception as e:
        print("TIFF export skipped:", e)
    plt.close(fig)


if __name__ == '__main__':
    n_vhh, n_vh = panel_A()
    colorbar_only()
    assemble(n_vhh, n_vh)
    print(f"Assembled final figure. VHH N={n_vhh}, VH N={n_vh}")
