"""
Stage 5 — Figure 2: Intrinsic vs Structural Chemistry (MAIN_OR_SUPPLEMENTARY_MECHANISM)
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
import pandas as pd
import numpy as np

PROJECT_ROOT = Path(__file__).parent.parent.parent
OUT_DIR = PROJECT_ROOT / "results" / "figures" / "stage5"
OUT_DIR.mkdir(parents=True, exist_ok=True)

VHH_BLUE = "#4D69B2"
VH_PURPLE = "#8670B1"


def main():
    df = pd.read_parquet(PROJECT_ROOT / "data/processed/stage5_sequence_surface_chemistry.parquet")

    vhh = df[df['structural_state'] == 'VHH_NATIVE'].dropna(subset=['H_intrinsic_binary_mean_KD', 'interface_sasa_weighted_mean_KD'])
    vh = df[df['structural_state'] == 'VH_WITHOUT_VL_GEOMETRIC'].dropna(subset=['H_intrinsic_binary_mean_KD', 'interface_sasa_weighted_mean_KD'])

    fig, ax = plt.subplots(figsize=(8, 7), facecolor='white')

    ax.scatter(vh['H_intrinsic_binary_mean_KD'], vh['interface_sasa_weighted_mean_KD'],
               s=15, color=VH_PURPLE, alpha=0.4, label=f'Conventional VH (n={len(vh)})')
    ax.scatter(vhh['H_intrinsic_binary_mean_KD'], vhh['interface_sasa_weighted_mean_KD'],
               s=15, color=VHH_BLUE, alpha=0.4, label=f'VHH (n={len(vhh)})')

    for sub, color in [(vh, VH_PURPLE), (vhh, VHH_BLUE)]:
        z = np.polyfit(sub['H_intrinsic_binary_mean_KD'], sub['interface_sasa_weighted_mean_KD'], 1)
        xs = np.linspace(sub['H_intrinsic_binary_mean_KD'].min(), sub['H_intrinsic_binary_mean_KD'].max(), 50)
        ax.plot(xs, np.polyval(z, xs), color=color, linewidth=2)

    ax.set_xlabel('Stage 4: H_intrinsic_binary_mean_KD (sequence chemistry)', fontsize=10)
    ax.set_ylabel('Stage 5: interface_sasa_weighted_mean_KD (structural exposure)', fontsize=10)
    ax.set_title('Does Sequence Chemistry Translate into Structural Solvent-Exposed Chemistry?\n'
                 '(group trend lines shown; correlation, not causation)', fontsize=11, fontweight='bold')
    ax.legend(fontsize=9)

    plt.tight_layout()
    png_path = OUT_DIR / "stage5_intrinsic_vs_structural.png"
    svg_path = OUT_DIR / "stage5_intrinsic_vs_structural.svg"
    fig.savefig(png_path, dpi=300, bbox_inches='tight')
    fig.savefig(svg_path, bbox_inches='tight')
    plt.close(fig)

    print(f"Wrote {png_path}")
    print(f"Wrote {svg_path}")


if __name__ == "__main__":
    main()
