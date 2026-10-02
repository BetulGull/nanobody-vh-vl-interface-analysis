"""Stage 6 — Figure 1: Direct CDR3 Shielding Distribution (MAIN_CANDIDATE)"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
import pandas as pd
import numpy as np

PROJECT_ROOT = Path(__file__).parent.parent.parent
OUT_DIR = PROJECT_ROOT / "results" / "figures" / "stage6"
OUT_DIR.mkdir(parents=True, exist_ok=True)

VHH_BLUE = "#4D69B2"


def main():
    df = pd.read_parquet(PROJECT_ROOT / "data/processed/stage6_vhh_cdr3_shielding.parquet")

    fig, axes = plt.subplots(1, 2, figsize=(12, 5.5), facecolor='white')

    for ax, col, title in [
        (axes[0], 'shielding_fraction', 'Shielding fraction\n(ΔSASA_shield / SASA_noCDR3)'),
        (axes[1], 'delta_sasa_shield', 'ΔSASA_shield (Å²)\n(absolute CDR3 occlusion)'),
    ]:
        vals = df[col].dropna()
        sorted_vals = np.sort(vals)
        ecdf = np.arange(1, len(sorted_vals) + 1) / len(sorted_vals)
        ax.plot(sorted_vals, ecdf, color=VHH_BLUE, linewidth=2)
        ax.fill_between(sorted_vals, 0, ecdf, alpha=0.15, color=VHH_BLUE)
        ax.set_xlabel(col, fontsize=10)
        ax.set_ylabel('Cumulative fraction of VHHs', fontsize=10)
        ax.set_title(f'{title}\n(N={len(vals)}, all values shown, none hidden)', fontsize=10)
        ax.axvline(vals.median(), color='black', linestyle='--', linewidth=1, label=f'median={vals.median():.2f}')
        ax.legend(fontsize=8)

    fig.suptitle('Stage 6 — Direct CDR3 Shielding of the Empirical Former VH-VL Surface (Shieldable Subset, N=533 VHH)',
                 fontsize=11, fontweight='bold')
    plt.tight_layout()

    png_path = OUT_DIR / "stage6_shielding_distribution.png"
    svg_path = OUT_DIR / "stage6_shielding_distribution.svg"
    fig.savefig(png_path, dpi=300, bbox_inches='tight')
    fig.savefig(svg_path, bbox_inches='tight')
    plt.close(fig)
    print(f"Wrote {png_path}")
    print(f"Wrote {svg_path}")


if __name__ == "__main__":
    main()
