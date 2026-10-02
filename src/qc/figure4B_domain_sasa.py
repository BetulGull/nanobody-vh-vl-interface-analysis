"""
Figure 4B -- domain-level SASA and Zhang-style surface hydrophobicity.

Zhang et al. 2025: "Surface hydrophobicity was calculated as the sum of
the solvent-accessible surface areas (SASAs) of amino acids in the VHH...
each multiplied by its respective hydrophobicity index."

Reuses the frozen, already-validated Stage-5 per-residue FreeSASA table
(stage5_interface_surface_chemistry_residues.parquet) -- covers ALL 127
IMGT domain positions (not just the interface) for both VHH_NATIVE and
VH_WITHOUT_VL_GEOMETRIC structural states. No new FreeSASA run needed;
this is a pure aggregation of an already-frozen per-residue table.
"""
import pandas as pd
from pathlib import Path

ROOT = Path("/Users/betulgul/Desktop/staj")


def main():
    res = pd.read_parquet(ROOT / "data/processed/stage5_interface_surface_chemistry_residues.parquet")
    sub = res[res.structural_state.isin(['VHH_NATIVE', 'VH_WITHOUT_VL_GEOMETRIC']) & (res.coordinate_observed == True)].copy()

    agg = sub.groupby(['INSTANCE', 'group', 'structural_state']).apply(
        lambda g: pd.Series({
            'domain_SASA': g['sasa'].sum(),
            'domain_surface_hydrophobicity_ZhangSASAxKD': (g['sasa'] * g['kyte_doolittle']).sum(),
            'n_residues_observed': len(g),
        }), include_groups=False
    ).reset_index()

    out_path = ROOT / "data/processed/figure4B_domain_sasa_hydrophobicity.parquet"
    agg.to_parquet(out_path, index=False)
    agg.to_csv(str(out_path).replace('.parquet', '.csv'), index=False)
    print(f"Wrote {len(agg)} rows")
    print(agg.groupby('group')[['domain_SASA', 'domain_surface_hydrophobicity_ZhangSASAxKD']].describe().T)
    return agg


if __name__ == '__main__':
    main()
