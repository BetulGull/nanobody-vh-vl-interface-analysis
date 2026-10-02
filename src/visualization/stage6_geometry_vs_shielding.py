"""
Stage 6 — Figure 2: Geometry vs Shielding (MAIN_CANDIDATE_OR_MECHANISM)

Pre-specified biologically relevant descriptors: cdr3_sequence_length,
compactness, FR2_CDR3_Dis vs shielding_fraction. Not selected post-hoc for
largest correlation -- these are exactly the three geometry descriptors the
spec designates as biologically relevant.
"""

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
    df = pd.read_parquet(PROJECT_ROOT / "data/processed/stage6_cdr3_shielding_geometry_combined.parquet")

    fig, axes = plt.subplots(1, 3, figsize=(15, 5), facecolor='white')

    panels = [
        (axes[0], 'cdr3_sequence_length', 'CDR3 sequence length (residues)'),
        (axes[1], 'cdr3_compactness', 'Gordon compactness (length / ρ)'),
        (axes[2], 'FR2_CDR3_Dis', 'Zhang FR2_CDR3_Dis (Å)'),
    ]

    for ax, col, label in panels:
        sub = df.dropna(subset=[col, 'shielding_fraction'])
        ax.scatter(sub[col], sub['shielding_fraction'], s=12, color=VHH_BLUE, alpha=0.5)
        z = np.polyfit(sub[col], sub['shielding_fraction'], 1)
        xs = np.linspace(sub[col].min(), sub[col].max(), 50)
        ax.plot(xs, np.polyval(z, xs), color='black', linewidth=1.5)
        r = sub[col].corr(sub['shielding_fraction'])
        ax.set_xlabel(label, fontsize=9)
        ax.set_ylabel('shielding_fraction', fontsize=9)
        ax.set_title(f'r = {r:.2f} (N={len(sub)}, exploratory)', fontsize=9)

    fig.suptitle('Stage 6 — CDR3 Geometry vs Direct Shielding of the Former VH-VL Interface\n'
                 '(pre-specified descriptors; no p-value stars; formal inference deferred to Stage 11)',
                 fontsize=11, fontweight='bold')
    plt.tight_layout()

    png_path = OUT_DIR / "stage6_geometry_vs_shielding.png"
    svg_path = OUT_DIR / "stage6_geometry_vs_shielding.svg"
    fig.savefig(png_path, dpi=300, bbox_inches='tight')
    fig.savefig(svg_path, bbox_inches='tight')
    plt.close(fig)
    print(f"Wrote {png_path}")
    print(f"Wrote {svg_path}")


if __name__ == "__main__":
    main()
