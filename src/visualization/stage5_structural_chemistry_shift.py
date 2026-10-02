"""
Stage 5 — Figure 1: Structural Surface Chemistry (MAIN_CANDIDATE)

VHH_NATIVE vs VH_WITHOUT_VL_GEOMETRIC for interface_sasa_weighted_mean_KD
and interface_positiveKD_sasa_fraction. Violin + box design.
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
import pandas as pd

PROJECT_ROOT = Path(__file__).parent.parent.parent
OUT_DIR = PROJECT_ROOT / "results" / "figures" / "stage5"
OUT_DIR.mkdir(parents=True, exist_ok=True)

VHH_BLUE = "#4D69B2"
VH_PURPLE = "#8670B1"


def main():
    df = pd.read_parquet(PROJECT_ROOT / "data/processed/stage5_sequence_surface_chemistry.parquet")

    vhh = df[df['structural_state'] == 'VHH_NATIVE']
    vh = df[df['structural_state'] == 'VH_WITHOUT_VL_GEOMETRIC']

    fig, axes = plt.subplots(1, 2, figsize=(12, 6), facecolor='white')

    for ax, col, title in [
        (axes[0], 'interface_sasa_weighted_mean_KD', 'SASA-weighted mean hydropathy\n(empirical interface, evaluable residues)'),
        (axes[1], 'interface_positiveKD_sasa_fraction', 'Fraction of interface SASA\ncarried by positive-KD residues'),
    ]:
        vhh_vals = vhh[col].dropna()
        vh_vals = vh[col].dropna()

        parts = ax.violinplot([vhh_vals, vh_vals], positions=[0, 1], showmeans=False, showextrema=False)
        for pc, color in zip(parts['bodies'], [VHH_BLUE, VH_PURPLE]):
            pc.set_facecolor(color)
            pc.set_alpha(0.5)

        ax.boxplot([vhh_vals, vh_vals], positions=[0, 1], widths=0.15, patch_artist=True,
                   boxprops=dict(facecolor='white'), medianprops=dict(color='black', linewidth=2))

        ax.set_xticks([0, 1])
        ax.set_xticklabels([f'VHH_NATIVE\n(n={len(vhh_vals)})', f'VH_WITHOUT_VL\n_GEOMETRIC\n(n={len(vh_vals)})'], fontsize=9)
        ax.set_ylabel(col, fontsize=9)
        ax.set_title(title, fontsize=10)

    fig.suptitle('Stage 5 — Native 3D Solvent-Exposed Surface Chemistry at the Empirical Former VH-VL Interface\n'
                 '(no inferential p-values shown; descriptive comparison only)',
                 fontsize=11, fontweight='bold')
    plt.tight_layout()

    png_path = OUT_DIR / "stage5_structural_chemistry_shift.png"
    svg_path = OUT_DIR / "stage5_structural_chemistry_shift.svg"
    fig.savefig(png_path, dpi=300, bbox_inches='tight')
    fig.savefig(svg_path, bbox_inches='tight')
    plt.close(fig)

    print(f"Wrote {png_path}")
    print(f"Wrote {svg_path}")


if __name__ == "__main__":
    main()
