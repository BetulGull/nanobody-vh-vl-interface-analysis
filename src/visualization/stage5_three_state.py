"""
Stage 5 — Optional Figure: Three-State Structural Comparison

VH_WITH_VL_PAIRED / VH_WITHOUT_VL_GEOMETRIC / VHH_NATIVE, contextualizing
normal VL burial, the VL-removal geometric exposure, and native VHH exposure.
Generated because the ordering is non-obvious (VHH_NATIVE is even more
hydrophilic than VH_WITH_VL_PAIRED), materially aiding interpretation.
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
import pandas as pd

PROJECT_ROOT = Path(__file__).parent.parent.parent
OUT_DIR = PROJECT_ROOT / "results" / "figures" / "stage5"
OUT_DIR.mkdir(parents=True, exist_ok=True)

COLORS = {
    'VH_WITH_VL_PAIRED': "#EF7F7E",
    'VH_WITHOUT_VL_GEOMETRIC': "#8670B1",
    'VHH_NATIVE': "#4D69B2",
}
LABELS = {
    'VH_WITH_VL_PAIRED': 'VH + VL\n(paired, normal burial)',
    'VH_WITHOUT_VL_GEOMETRIC': 'VH, VL removed\n(geometric counterfactual)',
    'VHH_NATIVE': 'VHH\n(native, no VL ever)',
}


def main():
    df = pd.read_parquet(PROJECT_ROOT / "data/processed/stage5_sequence_surface_chemistry.parquet")

    states = ['VH_WITH_VL_PAIRED', 'VH_WITHOUT_VL_GEOMETRIC', 'VHH_NATIVE']
    fig, ax = plt.subplots(figsize=(8, 6), facecolor='white')

    data = [df[df['structural_state'] == s]['interface_sasa_weighted_mean_KD'].dropna() for s in states]
    positions = [0, 1, 2]

    parts = ax.violinplot(data, positions=positions, showextrema=False)
    for pc, state in zip(parts['bodies'], states):
        pc.set_facecolor(COLORS[state])
        pc.set_alpha(0.6)

    ax.boxplot(data, positions=positions, widths=0.15, patch_artist=True,
               boxprops=dict(facecolor='white'), medianprops=dict(color='black', linewidth=2))

    ax.set_xticks(positions)
    ax.set_xticklabels([f"{LABELS[s]}\n(n={len(d)})" for s, d in zip(states, data)], fontsize=9)
    ax.set_ylabel('interface_sasa_weighted_mean_KD', fontsize=10)
    ax.set_title('Three Structural States at the Empirical Former VH-VL Interface\n'
                 '(VH-without-VL is a same-coordinate geometric counterfactual, NOT an equilibrated isolated VH)',
                 fontsize=10, fontweight='bold')

    plt.tight_layout()
    png_path = OUT_DIR / "stage5_three_state_comparison.png"
    svg_path = OUT_DIR / "stage5_three_state_comparison.svg"
    fig.savefig(png_path, dpi=300, bbox_inches='tight')
    fig.savefig(svg_path, bbox_inches='tight')
    plt.close(fig)

    print(f"Wrote {png_path}")
    print(f"Wrote {svg_path}")


if __name__ == "__main__":
    main()
