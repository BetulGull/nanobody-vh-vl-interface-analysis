"""Stage 5 — Figure Manifest and README Generation"""

from pathlib import Path
import pandas as pd

PROJECT_ROOT = Path(__file__).parent.parent.parent
FIG_DIR = PROJECT_ROOT / "results" / "figures" / "stage5"
FIG_DIR.mkdir(parents=True, exist_ok=True)


def main():
    seq_chem = pd.read_parquet(PROJECT_ROOT / "data/processed/stage5_sequence_surface_chemistry.parquet")
    n_vhh = len(seq_chem[seq_chem['structural_state'] == 'VHH_NATIVE'])
    n_vh = len(seq_chem[seq_chem['structural_state'] == 'VH_WITHOUT_VL_GEOMETRIC'])

    manifest_rows = [
        {
            'figure_id': 'FIG_S5_01', 'figure_role': 'MAIN_CANDIDATE',
            'filename': 'stage5_structural_chemistry_shift.png',
            'title': 'Native 3D Solvent-Exposed Surface Chemistry: VHH vs VH-Without-VL',
            'scientific_question': 'Does the sequence-level VHH-vs-VH chemistry difference translate into actual solvent-exposed structural chemistry at the empirical former interface?',
            'source_table': 'data/processed/stage5_sequence_surface_chemistry.parquet',
            'source_script': 'src/visualization/stage5_structural_chemistry_shift.py',
            'variables_shown': 'interface_sasa_weighted_mean_KD, interface_positiveKD_sasa_fraction (violin+box, VHH_NATIVE vs VH_WITHOUT_VL_GEOMETRIC)',
            'n_used': n_vhh + n_vh, 'n_excluded': 0,
            'filters': 'primary nonredundant representatives only',
            'notes': 'No inferential p-values shown, per spec. VH_WITHOUT_VL_GEOMETRIC is a same-coordinate counterfactual, not an equilibrated isolated VH.',
        },
        {
            'figure_id': 'FIG_S5_02', 'figure_role': 'MAIN_OR_SUPPLEMENTARY_MECHANISM',
            'filename': 'stage5_intrinsic_vs_structural.png',
            'title': 'Sequence Chemistry vs Structural Solvent Exposure',
            'scientific_question': 'Does Stage-4 intrinsic sequence chemistry predict Stage-5 structural solvent-exposed chemistry?',
            'source_table': 'data/processed/stage5_sequence_surface_chemistry.parquet',
            'source_script': 'src/visualization/stage5_intrinsic_vs_structural.py',
            'variables_shown': 'H_intrinsic_binary_mean_KD (x) vs interface_sasa_weighted_mean_KD (y), per group trend lines',
            'n_used': n_vhh + n_vh, 'n_excluded': 0,
            'filters': 'primary nonredundant representatives only',
            'notes': 'Correlation only, no causal claim. Trend lines are simple linear fits for visual reference.',
        },
        {
            'figure_id': 'FIG_S5_03', 'figure_role': 'OPTIONAL_THREE_STATE',
            'filename': 'stage5_three_state_comparison.png',
            'title': 'Three Structural States at the Empirical Interface',
            'scientific_question': 'How does exposure change across VH+VL (normal burial), VH-without-VL (geometric counterfactual), and native VHH?',
            'source_table': 'data/processed/stage5_sequence_surface_chemistry.parquet',
            'source_script': 'src/visualization/stage5_three_state.py',
            'variables_shown': 'interface_sasa_weighted_mean_KD across all three structural_state values',
            'n_used': n_vhh + 2 * n_vh, 'n_excluded': 0,
            'filters': 'primary nonredundant representatives only',
            'notes': 'Generated because the ordering is non-obvious (VHH_NATIVE more hydrophilic than VH_WITH_VL_PAIRED), materially aiding interpretation. VH-without-VL explicitly labeled as not physically equilibrated.',
        },
    ]

    manifest_df = pd.DataFrame(manifest_rows)
    manifest_out = FIG_DIR / "figure_manifest.csv"
    manifest_df.to_csv(manifest_out, index=False)
    print(f"Wrote {manifest_out}")

    readme = """# Stage 5 Figures

## FIG_S5_01 — Native 3D Solvent-Exposed Surface Chemistry (MAIN_CANDIDATE)

**What it shows:** Distributions of SASA-weighted mean Kyte-Doolittle
hydropathy and positive-KD exposed-area fraction across the 28-position
empirical former VH-VL interface, comparing native VHH structures against
the same-coordinate VL-removed conventional VH geometric counterfactual.

**What it does NOT show or prove:** binding energy, thermodynamic stability,
solubility, or that CDR3 causes any observed difference (Stage 5 never
removes CDR3). No inferential p-values are shown -- this is a descriptive
comparison; formal significance testing is deferred to Stage 11.

## FIG_S5_02 — Intrinsic vs Structural Chemistry (MAIN_OR_SUPPLEMENTARY_MECHANISM)

**What it shows:** Whether Stage-4's sequence-level intrinsic hydropathy
correlates with Stage-5's actual solvent-exposed structural chemistry, for
VHH and VH separately.

**What it does NOT show or prove:** causality. A correlation confirms
sequence chemistry "shows up" structurally; it does not establish which
drives which, or rule out confounds (e.g. local packing, CDR3 geometry).

## FIG_S5_03 — Three-State Comparison (OPTIONAL_THREE_STATE)

**What it shows:** The full progression from normally-buried paired VH+VL,
through the VL-removal geometric counterfactual, to native intact VHH.

**What it does NOT show or prove:** that VH-without-VL is a physically
equilibrated, experimentally realizable state -- it is a same-coordinate
computational reference only.
"""
    readme_path = FIG_DIR / "README.md"
    with open(readme_path, 'w') as f:
        f.write(readme)
    print(f"Wrote {readme_path}")


if __name__ == "__main__":
    main()
