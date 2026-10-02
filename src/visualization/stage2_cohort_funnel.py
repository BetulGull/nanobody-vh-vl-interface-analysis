"""
Stage-2 Cohort Funnel Figure

Justified per the visualization gate after the closure audit: structure-level
QC now materially reduces the dataset (VHH -16.5%, paired -6.9%).
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
import pandas as pd

PROJECT_ROOT = Path(__file__).parent.parent.parent
OUT_DIR = PROJECT_ROOT / "results" / "figures" / "stage2"
OUT_DIR.mkdir(parents=True, exist_ok=True)

VHH_BLUE = "#4D69B2"
CONV_PINK = "#EF7F7E"

funnel = pd.read_csv(PROJECT_ROOT / "results" / "tables" / "stage2_cohort_funnel.csv")

vhh_rows = funnel[funnel["stage"].str.contains("VHH")]
paired_rows = funnel[funnel["stage"].str.contains("paired")]

fig, axes = plt.subplots(1, 2, figsize=(12, 6), facecolor='white')

for ax, rows, color, title in [
    (axes[0], vhh_rows, VHH_BLUE, "VHH structural QC funnel"),
    (axes[1], paired_rows, CONV_PINK, "Paired VH-VL structural QC funnel"),
]:
    labels = [s.replace("Stage-1 ", "").replace("Stage-2 ", "").replace(" VHH", "").replace(" paired", "")
              for s in rows["stage"]]
    values = rows["n"].tolist()
    y_pos = range(len(labels))
    ax.barh(list(y_pos), values, color=color)
    ax.set_yticks(list(y_pos))
    ax.set_yticklabels(labels, fontsize=9)
    ax.invert_yaxis()
    for i, v in enumerate(values):
        ax.text(v + max(values) * 0.01, i, f"N={v}", va='center', fontsize=9)
    ax.set_xlabel("N structures")
    ax.set_title(title, fontsize=11, fontweight='bold')
    ax.set_xlim(0, max(values) * 1.15)

fig.suptitle("Stage 2 — Structural QC Cohort Funnel (post closure-audit fix)", fontsize=13, fontweight='bold')
plt.tight_layout()

png_path = OUT_DIR / "stage2_cohort_funnel.png"
svg_path = OUT_DIR / "stage2_cohort_funnel.svg"
fig.savefig(png_path, dpi=300, bbox_inches='tight')
fig.savefig(svg_path, bbox_inches='tight')
plt.close(fig)

print(f"Wrote {png_path}")
print(f"Wrote {svg_path}")
