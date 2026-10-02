"""Stage 6 — Figure Manifest and README"""

from pathlib import Path
import pandas as pd

PROJECT_ROOT = Path(__file__).parent.parent.parent
FIG_DIR = PROJECT_ROOT / "results" / "figures" / "stage6"
FIG_DIR.mkdir(parents=True, exist_ok=True)


def main():
    summaries = pd.read_parquet(PROJECT_ROOT / "data/processed/stage6_vhh_cdr3_shielding.parquet")
    n = len(summaries)

    rows = [
        {'figure_id': 'FIG_S6_01', 'figure_role': 'MAIN_CANDIDATE',
         'filename': 'stage6_shielding_distribution.png',
         'title': 'Direct CDR3 Shielding Distribution',
         'scientific_question': 'How much does CDR3 geometrically occlude the shieldable former VH-VL surface, across the full VHH cohort?',
         'source_table': 'data/processed/stage6_vhh_cdr3_shielding.parquet',
         'source_script': 'src/visualization/stage6_shielding_distribution.py',
         'variables_shown': 'shielding_fraction, delta_sasa_shield (ECDF)',
         'n_used': n, 'n_excluded': 0, 'filters': 'shieldable-interface target only (15 non-CDR3 positions)',
         'notes': 'Full distribution shown, no VHHs hidden. All 533 show positive shielding.'},
        {'figure_id': 'FIG_S6_02', 'figure_role': 'MAIN_CANDIDATE_OR_MECHANISM',
         'filename': 'stage6_geometry_vs_shielding.png',
         'title': 'CDR3 Geometry vs Direct Shielding',
         'scientific_question': 'Which CDR3 geometric descriptors (pre-specified: length, compactness, FR2_CDR3_Dis) associate with shielding magnitude?',
         'source_table': 'data/processed/stage6_cdr3_shielding_geometry_combined.parquet',
         'source_script': 'src/visualization/stage6_geometry_vs_shielding.py',
         'variables_shown': 'cdr3_sequence_length, cdr3_compactness, FR2_CDR3_Dis vs shielding_fraction',
         'n_used': n, 'n_excluded': 0, 'filters': 'none',
         'notes': 'Pre-specified descriptors, not selected for largest correlation. No p-value stars; formal inference deferred to Stage 11.'},
        {'figure_id': 'FIG_S6_03', 'figure_role': 'SUPPLEMENTARY_SPATIAL_MECHANISM',
         'filename': 'stage6_position_shielding.png',
         'title': 'Position-Level Shielding Map',
         'scientific_question': 'Which specific former-interface positions are most occluded by CDR3?',
         'source_table': 'data/processed/stage6_position_shielding_summary.parquet',
         'source_script': 'src/visualization/stage6_position_shielding.py',
         'variables_shown': 'mean_delta_sasa_shield per shieldable IMGT position, interface_frequency_5 (secondary track)',
         'n_used': 15, 'n_excluded': 13, 'filters': 'shieldable-interface positions only (13 CDR3 positions excluded by design)',
         'notes': 'Position 118 (FR4, immediately C-terminal to CDR3) and 42 (tetrad) show the largest mean shielding.'},
    ]

    df = pd.DataFrame(rows)
    df.to_csv(FIG_DIR / "figure_manifest.csv", index=False)
    print(f"Wrote {FIG_DIR / 'figure_manifest.csv'}")

    readme = """# Stage 6 Figures

## FIG_S6_01 — Direct CDR3 Shielding Distribution (MAIN_CANDIDATE)
Full ECDF of shielding_fraction and delta_sasa_shield across all 533 VHH
representatives, using the matched shieldable-interface target (28-position
empirical interface minus the 13 CDR3 positions). No VHHs are hidden,
including any with low shielding.

## FIG_S6_02 — CDR3 Geometry vs Shielding (MAIN_CANDIDATE_OR_MECHANISM)
Pre-specified geometric descriptors (CDR3 length, Gordon compactness, Zhang
FR2_CDR3_Dis) plotted against shielding_fraction. No p-value stars; these
are exploratory correlations, not confirmatory statistics (Stage 11).

## FIG_S6_03 — Position-Level Shielding (SUPPLEMENTARY_SPATIAL_MECHANISM)
Mean delta_sasa_shield at each of the 15 shieldable former-interface
positions, with Stage-3 interface_frequency_5 shown as a secondary
reference track.

**What none of these figures show or prove:** that CDR3 stabilizes VHH,
that shielding prevents aggregation, or that any specific atomistic
interaction holds CDR3 in place (Stage 7). These are purely geometric
solvent-accessible-area occlusion measurements on experimentally observed,
unrelaxed coordinates.
"""
    with open(FIG_DIR / "README.md", 'w') as f:
        f.write(readme)
    print(f"Wrote {FIG_DIR / 'README.md'}")


if __name__ == "__main__":
    main()
