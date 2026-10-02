"""
Stage 3 — Figure 1: Empirical VH-VL Interface Profile (MAIN_CANDIDATE)

x-axis: ordered VH IMGT positions
y-axis: interface_frequency_5 (PRIMARY nonredundant map)
Encodes: region shading, tetrad annotation, observation-fraction (opacity)
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path
import pandas as pd
import numpy as np

PROJECT_ROOT = Path(__file__).parent.parent.parent
OUT_DIR = PROJECT_ROOT / "results" / "figures" / "stage3"
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
PRIMARY_FREQ_CUTOFF = 0.25


def main():
    df = pd.read_parquet(PROJECT_ROOT / "data/processed/stage3_empirical_vh_vl_interface_primary.parquet")
    df = df.sort_values(['imgt_position', 'insertion_code']).reset_index(drop=True)
    df['x_order'] = range(len(df))
    n_structures_label = df['n_observed'].max()  # approx N for caption; exact denom varies per position

    fig, (ax_main, ax_obs) = plt.subplots(
        2, 1, figsize=(16, 7), facecolor='white',
        gridspec_kw={'height_ratios': [4, 1]}, sharex=True
    )

    # Region shading
    for region, (start, end) in IMGT_REGIONS.items():
        mask = (df['imgt_position'] >= start) & (df['imgt_position'] <= end)
        if mask.any():
            xs = df.loc[mask, 'x_order']
            ax_main.axvspan(xs.min() - 0.5, xs.max() + 0.5, color=REGION_COLORS[region], alpha=0.6, zorder=0)
            ax_obs.axvspan(xs.min() - 0.5, xs.max() + 0.5, color=REGION_COLORS[region], alpha=0.6, zorder=0)

    # Main bar plot: interface_frequency_5, opacity encodes observation_fraction
    for _, row in df.iterrows():
        alpha = max(0.15, min(1.0, row['observation_fraction']))
        ax_main.bar(row['x_order'], row['interface_frequency_5'], color=VH_PURPLE, alpha=alpha, width=0.8, zorder=2)

    ax_main.axhline(PRIMARY_FREQ_CUTOFF, color='black', linestyle='--', linewidth=1, label=f'Primary frequency cutoff = {PRIMARY_FREQ_CUTOFF}')

    # Tetrad annotation (markers above bars, not colored as predefined interface)
    tetrad_rows = df[df['vhh_tetrad_position']]
    for _, row in tetrad_rows.iterrows():
        ax_main.plot(row['x_order'], 1.05, marker='v', color=TETRAD_COLOR, markersize=8, zorder=3)
    if len(tetrad_rows) > 0:
        ax_main.plot([], [], marker='v', color=TETRAD_COLOR, linestyle='none', markersize=8, label='VHH tetrad position (42/49/50/52) — annotation only')

    ax_main.set_ylim(0, 1.15)
    ax_main.set_ylabel('Interface frequency\n(ΔSASA ≥ 5 Å²)', fontsize=10)
    ax_main.set_title(
        'Empirical VH-VL Interface Profile — Primary Nonredundant Map\n'
        f'({df["n_observed"].median():.0f} median structures per position observed)',
        fontsize=12, fontweight='bold'
    )
    ax_main.legend(loc='upper right', fontsize=8, framealpha=0.9)

    # Region labels on top
    for region, (start, end) in IMGT_REGIONS.items():
        mask = (df['imgt_position'] >= start) & (df['imgt_position'] <= end)
        if mask.any():
            xs = df.loc[mask, 'x_order']
            ax_main.text((xs.min() + xs.max()) / 2, 1.10, region, ha='center', fontsize=8, style='italic')

    # Observation-fraction secondary track
    ax_obs.bar(df['x_order'], df['observation_fraction'], color='gray', width=0.8, zorder=2)
    ax_obs.set_ylabel('Observation\nfraction', fontsize=9)
    ax_obs.set_ylim(0, 1.05)
    ax_obs.set_xlabel('VH IMGT position (ordered)', fontsize=10)

    # X tick labels: show every ~10th position
    tick_idx = list(range(0, len(df), max(1, len(df) // 20)))
    ax_obs.set_xticks(df['x_order'].iloc[tick_idx])
    ax_obs.set_xticklabels(
        [f"{r['imgt_position']}{r['insertion_code']}" for _, r in df.iloc[tick_idx].iterrows()],
        rotation=90, fontsize=7,
    )

    plt.tight_layout()
    png_path = OUT_DIR / "stage3_interface_profile.png"
    svg_path = OUT_DIR / "stage3_interface_profile.svg"
    fig.savefig(png_path, dpi=300, bbox_inches='tight')
    fig.savefig(svg_path, bbox_inches='tight')
    plt.close(fig)

    print(f"Wrote {png_path}")
    print(f"Wrote {svg_path}")


if __name__ == "__main__":
    main()
