"""
Stage 4 — Optional Figure: H_intrinsic Distribution (VHH vs VH)

Generated because the distributions are cleanly separated with similar
spread, materially aiding interpretation of the population-level shift
(not just a few outlier sequences).
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
import pandas as pd

PROJECT_ROOT = Path(__file__).parent.parent.parent
OUT_DIR = PROJECT_ROOT / "results" / "figures" / "stage4"
OUT_DIR.mkdir(parents=True, exist_ok=True)

VHH_BLUE = "#4D69B2"
VH_PURPLE = "#8670B1"


def main():
    df = pd.read_parquet(PROJECT_ROOT / "data/processed/stage4_sequence_interface_chemistry.parquet")

    fig, axes = plt.subplots(1, 2, figsize=(13, 5), facecolor='white')

    for ax, col, title in [
        (axes[0], 'H_intrinsic_binary_mean_KD', 'H_intrinsic_binary_mean_KD\n(mean Kyte-Doolittle over primary binary interface)'),
        (axes[1], 'H_intrinsic_weighted_freq5', 'H_intrinsic_weighted_freq5\n(interface_frequency_5-weighted KD)'),
    ]:
        vhh_vals = df[df['group'] == 'VHH'][col].dropna()
        vh_vals = df[df['group'] == 'VH'][col].dropna()

        ax.hist(vh_vals, bins=30, alpha=0.6, color=VH_PURPLE, label=f'Conventional VH (n={len(vh_vals)})', density=True)
        ax.hist(vhh_vals, bins=30, alpha=0.6, color=VHH_BLUE, label=f'VHH (n={len(vhh_vals)})', density=True)
        ax.axvline(vhh_vals.mean(), color=VHH_BLUE, linestyle='--', linewidth=1.5)
        ax.axvline(vh_vals.mean(), color=VH_PURPLE, linestyle='--', linewidth=1.5)
        ax.set_xlabel(col, fontsize=9)
        ax.set_ylabel('Density', fontsize=10)
        ax.set_title(title, fontsize=10)
        ax.legend(fontsize=8)

    fig.suptitle('Stage 4 — Intrinsic Sequence Hydropathy at the Empirical Interface (VHH vs VH)',
                 fontsize=12, fontweight='bold')
    plt.tight_layout()

    png_path = OUT_DIR / "stage4_hintrinsic_distribution.png"
    svg_path = OUT_DIR / "stage4_hintrinsic_distribution.svg"
    fig.savefig(png_path, dpi=300, bbox_inches='tight')
    fig.savefig(svg_path, bbox_inches='tight')
    plt.close(fig)

    print(f"Wrote {png_path}")
    print(f"Wrote {svg_path}")


if __name__ == "__main__":
    main()
