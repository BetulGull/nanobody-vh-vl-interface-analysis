"""Stage 4 — Figure Manifest and README Generation"""

from pathlib import Path
import pandas as pd

PROJECT_ROOT = Path(__file__).parent.parent.parent
FIG_DIR = PROJECT_ROOT / "results" / "figures" / "stage4"
FIG_DIR.mkdir(parents=True, exist_ok=True)


def main():
    comparison = pd.read_parquet(PROJECT_ROOT / "data/processed/stage4_position_chemistry_comparison.parquet")
    n_interface = int(comparison['primary_binary_interface'].sum())

    manifest_rows = [
        {
            'figure_id': 'FIG_S4_01', 'figure_role': 'MAIN_CANDIDATE',
            'filename': 'stage4_chemistry_shift.png',
            'title': 'Intrinsic Sequence Chemistry Shift at the Empirical VH-VL Interface',
            'scientific_question': 'At positions the empirical Stage-3 map identifies as buried by VL, does VHH intrinsic hydropathy differ from conventional VH?',
            'source_table': 'data/processed/stage4_position_chemistry_comparison.parquet',
            'source_script': 'src/visualization/stage4_chemistry_shift.py',
            'variables_shown': 'delta_mean_KD (y), imgt_position (x), interface_frequency_5 (point size), vhh_tetrad_position (color), region (background shading)',
            'n_used': n_interface, 'n_excluded': 214 - n_interface,
            'filters': 'primary_binary_interface == True only (28 positions)',
            'notes': 'Sequence-level intrinsic hydropathy (Kyte-Doolittle), NOT solvent-exposed hydrophobicity. Tetrad annotated, not enforced.',
        },
        {
            'figure_id': 'FIG_S4_02', 'figure_role': 'SUPPLEMENTARY_SEQUENCE_DIVERGENCE',
            'filename': 'stage4_positional_divergence.png',
            'title': 'Positional Sequence Divergence (Jensen-Shannon) at the Empirical Interface',
            'scientific_question': 'Which empirical interface positions show the strongest VHH-vs-VH amino-acid distribution divergence?',
            'source_table': 'data/processed/stage4_position_chemistry_comparison.parquet',
            'source_script': 'src/visualization/stage4_divergence_supplementary.py',
            'variables_shown': 'JSD (y), imgt_position (x), region (CDR vs FR color), vhh_tetrad_position (marker)',
            'n_used': n_interface, 'n_excluded': 214 - n_interface,
            'filters': 'primary_binary_interface == True only',
            'notes': 'Descriptive ranking only; no significance testing performed (deferred to Stage 11).',
        },
        {
            'figure_id': 'FIG_S4_03', 'figure_role': 'OPTIONAL_DISTRIBUTION',
            'filename': 'stage4_hintrinsic_distribution.png',
            'title': 'H_intrinsic Distribution: VHH vs Conventional VH',
            'scientific_question': 'Is the VHH-vs-VH chemistry shift a population-wide effect or driven by a subset of sequences?',
            'source_table': 'data/processed/stage4_sequence_interface_chemistry.parquet',
            'source_script': 'src/visualization/stage4_hintrinsic_distribution.py',
            'variables_shown': 'H_intrinsic_binary_mean_KD, H_intrinsic_weighted_freq5 (histograms), group means (dashed lines)',
            'n_used': 533 + 1429, 'n_excluded': 0,
            'filters': 'primary nonredundant representatives only (533 VHH, 1429 VH)',
            'notes': 'Generated because distributions are cleanly separated with comparable spread, confirming a population-level shift rather than an outlier-driven one.',
        },
    ]

    manifest_df = pd.DataFrame(manifest_rows)
    manifest_out = FIG_DIR / "figure_manifest.csv"
    manifest_df.to_csv(manifest_out, index=False)
    print(f"Wrote {manifest_out}")

    readme = f"""# Stage 4 Figures

## FIG_S4_01 — Intrinsic Sequence Chemistry Shift (MAIN_CANDIDATE)

**What it shows:** For each of the {n_interface} Stage-3 primary empirical
VH-VL interface positions, the difference in mean Kyte-Doolittle intrinsic
hydropathy between VHH and conventional VH sequences at that exact IMGT
position. Point size encodes the Stage-3 `interface_frequency_5` (how
consistently that position is buried by VL across paired structures); the
four VHH tetrad positions (42/49/50/52) are annotated in teal — as
annotation only, not as a predefined or privileged subset.

**What it does NOT show or prove:** solvent-exposed hydrophobicity (this is
sequence-level Kyte-Doolittle chemistry only — no SASA calculated in this
stage), binding energy, aggregation propensity, or anything about VHH
stability/solubility directly. A negative Δmean KD means VHH is
intrinsically less hydrophobic (by sequence alone) than conventional VH at
that position — nothing more.

## FIG_S4_02 — Positional Divergence (SUPPLEMENTARY_SEQUENCE_DIVERGENCE)

**What it shows:** Jensen-Shannon divergence between VHH and VH amino-acid
distributions at each of the 28 primary interface positions, distinguishing
framework from CDR positions, with tetrad positions marked.

**What it does NOT show or prove:** statistical significance — no
multiple-testing-corrected inference is performed in Stage 4 (deferred to
Stage 11). This is a descriptive ranking only.

## FIG_S4_03 — H_intrinsic Distribution (OPTIONAL_DISTRIBUTION)

**What it shows:** Full distributions (not just means) of
`H_intrinsic_binary_mean_KD` and `H_intrinsic_weighted_freq5` for the
primary nonredundant VHH (n=533) and VH (n=1429) representative sets,
confirming the shift seen in FIG_S4_01 is a population-wide pattern.

**What it does NOT show or prove:** causality, mechanism, or that this
reflects actual solvent exposure in either group.
"""
    readme_path = FIG_DIR / "README.md"
    with open(readme_path, 'w') as f:
        f.write(readme)
    print(f"Wrote {readme_path}")


if __name__ == "__main__":
    main()
