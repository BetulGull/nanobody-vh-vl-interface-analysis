"""
Stage 4 — Figure 2: Positional Divergence (SUPPLEMENTARY_SEQUENCE_DIVERGENCE)

JSD across the 28 primary empirical interface positions, tetrad annotated,
framework vs CDR3 distinguished.
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
import pandas as pd

PROJECT_ROOT = Path(__file__).parent.parent.parent
OUT_DIR = PROJECT_ROOT / "results" / "figures" / "stage4"
OUT_DIR.mkdir(parents=True, exist_ok=True)

TETRAD_COLOR = "#00BCD4"
CDR_COLOR = "#EF7F7E"
FR_COLOR = "#8670B1"


def main():
    df = pd.read_parquet(PROJECT_ROOT / "data/processed/stage4_position_chemistry_comparison.parquet")
    interface = df[df['primary_binary_interface'] == True].copy()
    interface = interface.sort_values(['imgt_position', 'insertion_code']).reset_index(drop=True)
    interface['x_order'] = range(len(interface))
    interface['is_cdr'] = interface['region'].str.startswith('CDR')

    fig, ax = plt.subplots(figsize=(12, 5), facecolor='white')

    colors = [CDR_COLOR if c else FR_COLOR for c in interface['is_cdr']]
    ax.bar(interface['x_order'], interface['JSD'], color=colors, width=0.7)

    for _, row in interface[interface['vhh_tetrad_position']].iterrows():
        ax.plot(row['x_order'], row['JSD'] + 0.03, marker='v', color=TETRAD_COLOR, markersize=8)

    ax.set_xticks(interface['x_order'])
    ax.set_xticklabels(
        [f"{int(r['imgt_position'])}{r['insertion_code']}" for _, r in interface.iterrows()],
        rotation=90, fontsize=8,
    )
    ax.set_xlabel('Empirical VH-VL interface position (IMGT)', fontsize=10)
    ax.set_ylabel('Jensen-Shannon divergence\n(VHH vs VH amino-acid distribution)', fontsize=10)
    ax.set_title('Stage 4 — Positional Sequence Divergence at the Empirical Interface', fontsize=12, fontweight='bold')

    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D
    legend_elements = [
        Patch(facecolor=FR_COLOR, label='Framework position'),
        Patch(facecolor=CDR_COLOR, label='CDR position'),
        Line2D([0], [0], marker='v', color='w', markerfacecolor=TETRAD_COLOR, markersize=8, label='VHH tetrad position'),
    ]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=8)

    plt.tight_layout()
    png_path = OUT_DIR / "stage4_positional_divergence.png"
    svg_path = OUT_DIR / "stage4_positional_divergence.svg"
    fig.savefig(png_path, dpi=300, bbox_inches='tight')
    fig.savefig(svg_path, bbox_inches='tight')
    plt.close(fig)

    print(f"Wrote {png_path}")
    print(f"Wrote {svg_path}")


if __name__ == "__main__":
    main()
