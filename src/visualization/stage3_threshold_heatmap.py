"""
Stage 3 — Figure 2: Threshold Sensitivity Heatmap (SUPPLEMENTARY_SENSITIVITY)
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
import pandas as pd
import numpy as np

PROJECT_ROOT = Path(__file__).parent.parent.parent
OUT_DIR = PROJECT_ROOT / "results" / "figures" / "stage3"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def main():
    df = pd.read_csv(PROJECT_ROOT / "results/tables/stage3_interface_threshold_sensitivity.csv")

    delta_thresholds = sorted(df['delta_sasa_threshold'].unique())
    freq_thresholds = sorted(df['frequency_threshold'].unique())
    obs_thresholds = sorted(df['observation_fraction_threshold'].unique())

    fig, axes = plt.subplots(1, len(obs_thresholds), figsize=(15, 5), facecolor='white', sharey=True)

    for ax, obs_t in zip(axes, obs_thresholds):
        sub = df[df['observation_fraction_threshold'] == obs_t]
        pivot_n = sub.pivot(index='delta_sasa_threshold', columns='frequency_threshold', values='n_consensus_positions')
        pivot_n = pivot_n.reindex(index=delta_thresholds, columns=freq_thresholds)

        im = ax.imshow(pivot_n.values, cmap='viridis', aspect='auto')
        ax.set_xticks(range(len(freq_thresholds)))
        ax.set_xticklabels(freq_thresholds)
        ax.set_yticks(range(len(delta_thresholds)))
        ax.set_yticklabels(delta_thresholds)
        ax.set_xlabel('Frequency threshold')
        if ax is axes[0]:
            ax.set_ylabel('ΔSASA threshold (Å²)')
        ax.set_title(f'obs. fraction ≥ {obs_t}', fontsize=10)

        for i in range(len(delta_thresholds)):
            for j in range(len(freq_thresholds)):
                val = pivot_n.values[i, j]
                if not np.isnan(val):
                    ax.text(j, i, int(val), ha='center', va='center', color='white', fontsize=9)

        cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label('N consensus positions', fontsize=8)

    fig.suptitle('Stage 3 — Threshold Sensitivity: Consensus Interface Position Count', fontsize=13, fontweight='bold')
    plt.tight_layout()

    png_path = OUT_DIR / "stage3_threshold_sensitivity.png"
    svg_path = OUT_DIR / "stage3_threshold_sensitivity.svg"
    fig.savefig(png_path, dpi=300, bbox_inches='tight')
    fig.savefig(svg_path, bbox_inches='tight')
    plt.close(fig)

    print(f"Wrote {png_path}")
    print(f"Wrote {svg_path}")


if __name__ == "__main__":
    main()
