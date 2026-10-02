"""
Stage 3 — Figure Manifest and README Generation
"""

from pathlib import Path
import pandas as pd

PROJECT_ROOT = Path(__file__).parent.parent.parent
FIG_DIR = PROJECT_ROOT / "results" / "figures" / "stage3"
FIG_DIR.mkdir(parents=True, exist_ok=True)


def main():
    primary_map = pd.read_parquet(PROJECT_ROOT / "data/processed/stage3_empirical_vh_vl_interface_primary.parquet")
    n_positions = len(primary_map)
    n_primary_interface = int(primary_map['primary_binary_interface'].sum())

    manifest_rows = [
        {
            'figure_id': 'FIG_S3_01', 'figure_role': 'MAIN_CANDIDATE',
            'filename': 'stage3_interface_profile.png',
            'title': 'Empirical VH-VL Interface Profile (Primary Nonredundant Map)',
            'scientific_question': 'Which VH IMGT positions are consistently buried by VL across independent paired structures?',
            'source_table': 'data/processed/stage3_empirical_vh_vl_interface_primary.parquet',
            'source_script': 'src/visualization/stage3_interface_profile.py',
            'variables_shown': 'interface_frequency_5 (y), imgt_position (x), region, observation_fraction, vhh_tetrad_position',
            'n_used': n_positions, 'n_excluded': 0,
            'filters': 'PRIMARY_NONREDUNDANT map only',
            'notes': f'{n_primary_interface} positions meet primary binary interface definition (freq>=0.25, obs_frac>=0.25, threshold=5A2)',
        },
        {
            'figure_id': 'FIG_S3_02', 'figure_role': 'SUPPLEMENTARY_SENSITIVITY',
            'filename': 'stage3_threshold_sensitivity.png',
            'title': 'Threshold Sensitivity: Consensus Interface Position Count',
            'scientific_question': 'How sensitive is the consensus interface definition to the choice of ΔSASA / frequency / observation thresholds?',
            'source_table': 'results/tables/stage3_interface_threshold_sensitivity.csv',
            'source_script': 'src/visualization/stage3_threshold_heatmap.py',
            'variables_shown': 'n_consensus_positions across delta_sasa_threshold x frequency_threshold, faceted by observation_fraction_threshold',
            'n_used': 27, 'n_excluded': 0,
            'filters': 'none (full threshold grid)',
            'notes': 'Primary definition (5A2, 0.25, 0.25) is one cell in this grid, not chosen after inspecting the grid',
        },
    ]

    manifest_df = pd.DataFrame(manifest_rows)
    manifest_out = FIG_DIR / "figure_manifest.csv"
    manifest_df.to_csv(manifest_out, index=False)
    print(f"Wrote {manifest_out}")

    readme = f"""# Stage 3 Figures

## FIG_S3_01 — Empirical VH-VL Interface Profile (MAIN_CANDIDATE)

**What it shows:** For every observed VH IMGT position (primary nonredundant
structural cohort), the fraction of structures in which that position had
ΔSASA ≥ 5 Å² due to VL presence (`interface_frequency_5`). Region shading
(FR/CDR), a dashed line at the primary consensus cutoff (0.25), and small
markers for the 4 VHH tetrad positions (42/49/50/52, annotation only — not
treated as predefined interface) are shown. A secondary track below shows
`observation_fraction` per position, since a position observed in very few
structures can otherwise look misleadingly "confident."

**What it does NOT show or prove:** binding energy, thermodynamic stability,
hydrophobicity, or anything about VHH. It is a purely geometric,
population-level empirical summary of {n_positions} VH IMGT positions across
the primary nonredundant paired-structure cohort. {n_primary_interface}
positions meet the primary binary consensus definition; this is a starting
definition, not a final biological claim — see the threshold-sensitivity
figure and tables for how it shifts under other reasonable thresholds.

## FIG_S3_02 — Threshold Sensitivity (SUPPLEMENTARY_SENSITIVITY)

**What it shows:** How the number of consensus interface positions changes
across the full grid of per-structure ΔSASA thresholds (1/5/10 Å²),
consensus frequency thresholds (0.10/0.25/0.50), and minimum
observation-fraction thresholds (0.10/0.25/0.50). The primary definition is
one specific cell in this grid, defined *before* looking at these results,
not chosen after inspecting them.

**What it does NOT show or prove:** which threshold is "correct" — no single
threshold combination is presented as uniquely biologically valid. This
figure exists purely to make the sensitivity of the definition transparent.
"""

    readme_path = FIG_DIR / "README.md"
    with open(readme_path, 'w') as f:
        f.write(readme)
    print(f"Wrote {readme_path}")


if __name__ == "__main__":
    main()
