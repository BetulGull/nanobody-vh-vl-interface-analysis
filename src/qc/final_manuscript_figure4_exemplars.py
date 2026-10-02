"""
Final manuscript figures -- Figure 4 exemplar selection (corrected).

The FR2_CDR3_Dis multimodality validation gate
(results/final_manuscript_figures/audits/figure4_fr2_cdr3_distance_distribution_validation.md)
found the KDE bimodal appearance is NOT bandwidth-robust and the formal
Hartigan dip test is non-significant (p=0.081). Per the gate's fallback
rule, low/high exemplars use deterministic Q25/Q75 quantiles of
FR2_CDR3_Dis, not KDE mode locations. Labels are "LOW/HIGH FR2-CDR3
DISTANCE", never bent/extended.
"""
import pandas as pd
from pathlib import Path

ROOT = Path("/Users/betulgul/Desktop/staj")


def main():
    d = pd.read_parquet(ROOT / "data/processed/final_figures_model_input.parquet")
    d['resolution'] = pd.to_numeric(d['resolution'], errors='coerce')
    vhh = d[d.group == 'VHH'].copy()
    elig = vhh[(vhh['primary_interface_position_coverage'] == 1.0) & (vhh['resolution'].notna())]

    q25 = elig['FR2_CDR3_Dis'].quantile(0.25)
    q75 = elig['FR2_CDR3_Dis'].quantile(0.75)

    rows = []
    for label, panel, target in [('low', 'Figure4-D', q25), ('high', 'Figure4-E', q75)]:
        e = elig.copy()
        e['dist'] = (e['FR2_CDR3_Dis'] - target).abs()
        e = e.sort_values(['dist', 'resolution', 'SABDAB_ID'])
        sel = e.iloc[0]
        rows.append({
            'panel': panel, 'label': f'{label.upper()} FR2-CDR3 DISTANCE', 'group': 'VHH',
            'INSTANCE': sel['INSTANCE'], 'SABDAB_ID': sel.get('SABDAB_ID', None),
            'selection_metric': 'FR2_CDR3_Dis', 'target_quantile': f'Q{25 if label=="low" else 75}',
            'target_value': round(target, 4), 'actual_value': round(sel['FR2_CDR3_Dis'], 4),
            'coverage': sel['primary_interface_position_coverage'], 'resolution': sel['resolution'],
            'tie_break': 'lowest resolution then SABDAB_ID alphabetical',
            'selection_reason': (f'nearest-to-{"Q25" if label=="low" else "Q75"} FR2_CDR3_Dis among '
                                  f'coverage==1.0 eligible VHH representatives (N={len(elig)} eligible); '
                                  f'quantile-based per figure4 multimodality gate, KDE bimodality not robust'),
        })

    out = pd.DataFrame(rows)
    out_path = ROOT / "results/final_manuscript_figures/tables/figure4_fr2_cdr3_exemplars_quantile.csv"
    out.to_csv(out_path, index=False)
    print(out.to_string())
    return out


if __name__ == '__main__':
    main()
