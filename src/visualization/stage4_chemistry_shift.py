"""
Stage 4 — Figure 1: Positional Chemistry Shift (MAIN_CANDIDATE)

x-axis: ordered exact IMGT interface positions (28 primary empirical
interface positions only)
y-axis: delta_mean_KD = VHH mean KD - VH mean KD
Point size/opacity encodes Stage-3 interface_frequency_5 (subtly).
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
import pandas as pd
import numpy as np

PROJECT_ROOT = Path(__file__).parent.parent.parent
OUT_DIR = PROJECT_ROOT / "results" / "figures" / "stage4"
OUT_DIR.mkdir(parents=True, exist_ok=True)

VH_PURPLE = "#8670B1"
TETRAD_COLOR = "#00BCD4"
REGION_COLORS = {
    'FR1': '#E8E8E8', 'CDR1': '#FDE2C8', 'FR2': '#E8E8E8', 'CDR2': '#FDE2C8',
    'FR3': '#E8E8E8', 'CDR3': '#FDE2C8', 'FR4': '#E8E8E8',
}
IMGT_REGIONS = {
    'FR1': (1, 26), 'CDR1': (27, 38), 'FR2': (39, 55), 'CDR2': (56, 65),
    'FR3': (66, 104), 'CDR3': (105, 117), 'FR4': (118, 128),
}
TETRAD_POSITIONS = {42, 49, 50, 52}


def main():
    df = pd.read_parquet(PROJECT_ROOT / "data/processed/stage4_position_chemistry_comparison.parquet")
    interface = df[df['primary_binary_interface'] == True].copy()
    interface = interface.sort_values(['imgt_position', 'insertion_code']).reset_index(drop=True)
    interface['x_order'] = range(len(interface))

    fig, ax = plt.subplots(figsize=(14, 6.5), facecolor='white')

    for region, (start, end) in IMGT_REGIONS.items():
        mask = (interface['imgt_position'] >= start) & (interface['imgt_position'] <= end)
        if mask.any():
            xs = interface.loc[mask, 'x_order']
            ax.axvspan(xs.min() - 0.5, xs.max() + 0.5, color=REGION_COLORS[region], alpha=0.6, zorder=0)

    ax.axhline(0, color='black', linewidth=1, linestyle='-', zorder=1)

    sizes = 60 + 300 * interface['interface_frequency_5']
    colors = [TETRAD_COLOR if t else VH_PURPLE for t in interface['vhh_tetrad_position']]
    ax.scatter(interface['x_order'], interface['delta_mean_KD'], s=sizes, c=colors,
               alpha=0.75, edgecolors='black', linewidths=0.5, zorder=3)

    for _, row in interface[interface['vhh_tetrad_position']].iterrows():
        ax.annotate(f"{int(row['imgt_position'])}", (row['x_order'], row['delta_mean_KD']),
                    textcoords="offset points", xytext=(0, 10), ha='center', fontsize=8, fontweight='bold')

    ax.set_xticks(interface['x_order'])
    ax.set_xticklabels(
        [f"{int(r['imgt_position'])}{r['insertion_code']}" for _, r in interface.iterrows()],
        rotation=90, fontsize=8,
    )
    ax.set_xlabel('Empirical VH-VL interface position (IMGT)', fontsize=11)
    ax.set_ylabel(r'$\Delta$mean KD  =  VHH mean KD $-$ VH mean KD', fontsize=11)
    ax.set_title(
        'Stage 4 — Intrinsic Sequence Chemistry Shift at the Empirical VH-VL Interface\n'
        '(point size = Stage-3 interface_frequency_5; teal = VHH tetrad position, annotation only)',
        fontsize=12, fontweight='bold'
    )

    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor=VH_PURPLE, markersize=10, label='Non-tetrad interface position'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor=TETRAD_COLOR, markersize=10, label='VHH tetrad position (42/49/50/52)'),
    ]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=8, framealpha=0.9)

    for region, (start, end) in IMGT_REGIONS.items():
        mask = (interface['imgt_position'] >= start) & (interface['imgt_position'] <= end)
        if mask.any():
            xs = interface.loc[mask, 'x_order']
            ymax = ax.get_ylim()[1]
            ax.text((xs.min() + xs.max()) / 2, ymax * 0.95, region, ha='center', fontsize=8, style='italic')

    plt.tight_layout()
    png_path = OUT_DIR / "stage4_chemistry_shift.png"
    svg_path = OUT_DIR / "stage4_chemistry_shift.svg"
    fig.savefig(png_path, dpi=300, bbox_inches='tight')
    fig.savefig(svg_path, bbox_inches='tight')
    plt.close(fig)

    print(f"Wrote {png_path}")
    print(f"Wrote {svg_path}")


if __name__ == "__main__":
    main()
