"""Stage 6 — Figure 3: Position-Level Shielding (SUPPLEMENTARY_SPATIAL_MECHANISM)"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
import pandas as pd

PROJECT_ROOT = Path(__file__).parent.parent.parent
OUT_DIR = PROJECT_ROOT / "results" / "figures" / "stage6"
OUT_DIR.mkdir(parents=True, exist_ok=True)

VHH_BLUE = "#4D69B2"


def main():
    df = pd.read_parquet(PROJECT_ROOT / "data/processed/stage6_position_shielding_summary.parquet")
    df = df.sort_values(['imgt_position', 'insertion_code']).reset_index(drop=True)
    df['x_order'] = range(len(df))

    fig, ax1 = plt.subplots(figsize=(10, 5.5), facecolor='white')

    ax1.bar(df['x_order'], df['mean_delta_sasa_shield'], color=VHH_BLUE, width=0.6)
    ax1.set_ylabel('mean ΔSASA_shield (Å²)', fontsize=10, color=VHH_BLUE)
    ax1.set_xticks(df['x_order'])
    ax1.set_xticklabels([f"{int(r['imgt_position'])}{r['insertion_code']}" for _, r in df.iterrows()], rotation=0, fontsize=9)
    ax1.set_xlabel('Shieldable empirical-interface IMGT position', fontsize=10)

    ax2 = ax1.twinx()
    ax2.plot(df['x_order'], df['mean_interface_frequency_5'], color='gray', marker='o', markersize=4, linewidth=1, alpha=0.7)
    ax2.set_ylabel('Stage-3 interface_frequency_5', fontsize=9, color='gray')
    ax2.set_ylim(0, 1.05)

    ax1.set_title('Stage 6 — Which Former-Interface Positions Are Preferentially Occluded by CDR3?', fontsize=11, fontweight='bold')
    plt.tight_layout()

    png_path = OUT_DIR / "stage6_position_shielding.png"
    svg_path = OUT_DIR / "stage6_position_shielding.svg"
    fig.savefig(png_path, dpi=300, bbox_inches='tight')
    fig.savefig(svg_path, bbox_inches='tight')
    plt.close(fig)
    print(f"Wrote {png_path}")
    print(f"Wrote {svg_path}")


if __name__ == "__main__":
    main()
