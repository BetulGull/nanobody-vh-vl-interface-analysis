"""
Final figures -- Section 20/Step 6: freeze ALL deterministic 3D exemplar
selections needed across Figures 1, 3, 4, 6, 7 BEFORE any rendering occurs.

Every rule below reads only frozen upstream tables (final_figures_model_input,
stage9 antigen-contact long table) -- no new structural computation, no
visual inspection. Ties are broken deterministically (lowest resolution,
then SABDAB_ID alphabetical) and recorded.

Writes: results/structures/final/final_figure_3d_exemplar_selections.csv
Columns: figure, panel, group, PDB, INSTANCE, sequence_group_id,
selection_metric, target_quantile_or_value, actual_value, coverage,
resolution, tie_break, selection_reason
"""
import pandas as pd
import numpy as np
from pathlib import Path
from scipy.stats import gaussian_kde
from scipy.signal import argrelextrema

ROOT = Path("/Users/betulgul/Desktop/staj")


def eligible(d, group, coverage_min=1.0):
    g = d[d.group == group].copy()
    g['resolution'] = pd.to_numeric(g['resolution'], errors='coerce')
    return g[(g['primary_interface_position_coverage'] >= coverage_min) & (g['resolution'].notna())]


def nearest(elig, metric_col, target, extra_tiebreak_cols=('resolution', 'SABDAB_ID')):
    e = elig.copy()
    e['dist'] = (e[metric_col] - target).abs()
    e = e.sort_values(['dist', *extra_tiebreak_cols])
    return e.iloc[0], e['dist'].iloc[0]


def row_from(sel, dist, figure, panel, group, metric, target, tie_break, reason):
    return {
        'figure': figure, 'panel': panel, 'group': group,
        'PDB': sel['PDB'] if 'PDB' in sel else sel['SABDAB_ID'].split('_')[0],
        'INSTANCE': sel['INSTANCE'], 'sequence_group_id': sel.get('sequence_group_id', np.nan),
        'selection_metric': metric, 'target_quantile_or_value': target,
        'actual_value': sel[metric] if metric in sel else np.nan,
        'coverage': sel['primary_interface_position_coverage'], 'resolution': sel['resolution'],
        'tie_break': tie_break, 'selection_reason': reason,
    }


def main():
    d = pd.read_parquet(ROOT / "data/processed/final_figures_model_input.parquet")
    reg = pd.read_parquet(ROOT / "data/processed/final_regional_sequence_physicochemistry.parquet")[['INSTANCE', 'sequence_group_id']]
    d = d.merge(reg, on='INSTANCE', how='left')
    d['PDB'] = d['INSTANCE'].str.replace('pdb_', '', regex=False).str.split('-').str[0]

    rows = []

    # --- Figure 1 Panel A/B: native Fv exemplar (nearest-to-median resolution, coverage==1.0, VH group) ---
    elig = eligible(d, 'VH')
    med_res = elig['resolution'].median()
    sel, dist = nearest(elig, 'resolution', med_res)
    rows.append(row_from(sel, dist, 'Figure1', 'A_B', 'VH', 'resolution', med_res,
                          'lowest resolution then SABDAB_ID alphabetical',
                          'nearest-to-median resolution among coverage==1.0 eligible VH representatives (typical native Fv)'))

    # --- Figure 4 Panels D/E: low/high FR2_CDR3_Dis VHH exemplars (objectively-detected KDE modes) ---
    x = d[d.group == 'VHH']['FR2_CDR3_Dis'].dropna().values
    kde = gaussian_kde(x)
    grid = np.linspace(x.min(), x.max(), 500)
    dens = kde(grid)
    maxima_idx = argrelextrema(dens, np.greater)[0]
    minima_idx = argrelextrema(dens, np.less)[0]
    peaks = sorted(grid[maxima_idx])
    trough = grid[minima_idx][0] if len(minima_idx) else float(np.median(x))
    low_target, high_target = peaks[0], peaks[-1]
    elig_vhh = eligible(d, 'VHH')
    sel_low, dist_low = nearest(elig_vhh, 'FR2_CDR3_Dis', low_target)
    sel_high, dist_high = nearest(elig_vhh, 'FR2_CDR3_Dis', high_target)
    reason_common = (f'KDE-detected 2 local maxima at {low_target:.2f}/{high_target:.2f} Angstrom, '
                      f'1 local minimum (trough) at {trough:.2f} Angstrom, on the VHH FR2_CDR3_Dis '
                      f'distribution -- objectively detected in this cohort, not assumed from Zhang 2025; '
                      f'labeled Low/High distance, NOT bent/extended, absent further formal validation')
    rows.append(row_from(sel_low, dist_low, 'Figure4', 'D', 'VHH', 'FR2_CDR3_Dis', low_target,
                          'lowest resolution then SABDAB_ID alphabetical', reason_common + ' (low-distance mode)'))
    rows.append(row_from(sel_high, dist_high, 'Figure4', 'E', 'VHH', 'FR2_CDR3_Dis', high_target,
                          'lowest resolution then SABDAB_ID alphabetical', reason_common + ' (high-distance mode)'))

    # --- Figure 6 bottom row: 2x2 CDR3-length Q25/Q75 exemplars, VHH and VH separately ---
    for grp in ['VHH', 'VH']:
        elig_g = eligible(d, grp)
        q25 = elig_g['CDR3_length'].quantile(0.25)
        q75 = elig_g['CDR3_length'].quantile(0.75)
        sel_s, dist_s = nearest(elig_g, 'CDR3_length', q25)
        sel_l, dist_l = nearest(elig_g, 'CDR3_length', q75)
        rows.append(row_from(sel_s, dist_s, 'Figure6', f'{grp}_short', grp, 'CDR3_length', q25,
                              'lowest resolution then SABDAB_ID alphabetical',
                              f'nearest-to-Q25 CDR3_length among coverage==1.0 eligible {grp} representatives'))
        rows.append(row_from(sel_l, dist_l, 'Figure6', f'{grp}_long', grp, 'CDR3_length', q75,
                              'lowest resolution then SABDAB_ID alphabetical',
                              f'nearest-to-Q75 CDR3_length among coverage==1.0 eligible {grp} representatives'))

    # --- Figure 7 G/H: antigen-complex exemplars, matched-antigen validation ---
    contacts = pd.read_parquet(ROOT / "data/processed/stage9_antibody_antigen_contacts_long.parquet")
    contacts = contacts.merge(d[['INSTANCE', 'antigen_name']], on='INSTANCE', how='left')
    candidates = sorted(set(d[d.group == 'VHH']['antigen_name'].dropna()) & set(d[d.group == 'VH']['antigen_name'].dropna()))
    best = None
    for ag in candidates:
        sub = contacts[(contacts.antigen_name == ag) & (contacts.contact_at_4A == True)]
        vhh_res = set(sub[sub.group == 'VHH']['antigen_residue_id'])
        vh_res = set(sub[sub.group == 'VH']['antigen_residue_id'])
        n_vhh_inst = sub[sub.group == 'VHH']['INSTANCE'].nunique()
        n_vh_inst = sub[sub.group == 'VH']['INSTANCE'].nunique()
        if n_vhh_inst < 1 or n_vh_inst < 1 or len(vhh_res) == 0 or len(vh_res) == 0:
            continue
        inter = vhh_res & vh_res
        union = vhh_res | vh_res
        jac = len(inter) / len(union) if union else 0.0
        if best is None or jac > best[1]:
            best = (ag, jac, n_vhh_inst, n_vh_inst, len(inter), len(union))

    matched_validated = best is not None and best[1] >= 0.30 and best[4] >= 10
    if matched_validated:
        ag, jac, n_vhh_inst, n_vh_inst, n_shared, n_union = best
        cand_vhh = d[(d.group == 'VHH') & (d.antigen_name == ag)]
        cand_vhh = eligible(cand_vhh.assign(primary_interface_position_coverage=cand_vhh['primary_interface_position_coverage']), 'VHH', coverage_min=0.0)
        cand_vhh['resolution'] = pd.to_numeric(cand_vhh['resolution'], errors='coerce')
        sel_vhh = cand_vhh.sort_values(['primary_interface_position_coverage', 'resolution'], ascending=[False, True]).iloc[0]
        cand_vh = d[(d.group == 'VH') & (d.antigen_name == ag)]
        cand_vh['resolution'] = pd.to_numeric(cand_vh['resolution'], errors='coerce')
        sel_vh = cand_vh.sort_values(['primary_interface_position_coverage', 'resolution'], ascending=[False, True]).iloc[0]
        reason = (f'VALIDATED matched-antigen pair: antigen_name="{ag}", Jaccard antigen-contact-residue '
                  f'overlap at 4A = {jac:.3f} ({n_shared}/{n_union} shared residues), VHH_complexes={n_vhh_inst}, '
                  f'VH_complexes={n_vh_inst}. Selected the highest-coverage/lowest-resolution structure per group. '
                  f'NOTE: this is a matched-antigen, epitope-OVERLAPPING pair (not confirmed identical epitope like '
                  f'Gordon Fig 14 6WAQ/6YLA) -- labeled precisely as such, never claimed as fully matched.')
        rows.append(row_from(sel_vhh, 0.0, 'Figure7', 'G', 'VHH', 'coverage_then_resolution', np.nan, 'coverage desc, resolution asc', reason))
        rows.append(row_from(sel_vh, 0.0, 'Figure7', 'H', 'VH', 'coverage_then_resolution', np.nan, 'coverage desc, resolution asc', reason))
    else:
        elig_vhh_ag = d[(d.group == 'VHH') & d.antigen_name.notna()].copy()
        elig_vhh_ag['resolution'] = pd.to_numeric(elig_vhh_ag['resolution'], errors='coerce')
        elig_vhh_ag = elig_vhh_ag[elig_vhh_ag.resolution.notna()]
        med = elig_vhh_ag['resolution'].median()
        sel_vhh, dist_vhh = nearest(elig_vhh_ag, 'resolution', med)
        elig_vh_ag = d[(d.group == 'VH') & d.antigen_name.notna()].copy()
        elig_vh_ag['resolution'] = pd.to_numeric(elig_vh_ag['resolution'], errors='coerce')
        elig_vh_ag = elig_vh_ag[elig_vh_ag.resolution.notna()]
        med2 = elig_vh_ag['resolution'].median()
        sel_vh, dist_vh = nearest(elig_vh_ag, 'resolution', med2)
        reason = ('No validated matched-antigen pair found (best candidate Jaccard < 0.30 or < 10 shared '
                   'contact residues) -- using INDEPENDENT exemplars (different antigens), explicitly labeled '
                   'as independent, NOT a Gordon-Fig-14-style matched-epitope comparison.')
        rows.append(row_from(sel_vhh, dist_vhh, 'Figure7', 'G', 'VHH', 'resolution', med, 'lowest resolution then SABDAB_ID alphabetical', reason))
        rows.append(row_from(sel_vh, dist_vh, 'Figure7', 'H', 'VH', 'resolution', med2, 'lowest resolution then SABDAB_ID alphabetical', reason))

    out = pd.DataFrame(rows)
    out_path = ROOT / "results/structures/final/final_figure_3d_exemplar_selections.csv"
    out.to_csv(out_path, index=False)
    print(f"Wrote {len(out)} exemplar selections to {out_path}")
    print(out[['figure', 'panel', 'group', 'INSTANCE', 'selection_metric', 'actual_value', 'coverage', 'resolution']].to_string())
    print(f"\nMatched-antigen validation: {'VALIDATED' if matched_validated else 'NOT validated'} -- {best}")
    return out


if __name__ == '__main__':
    main()
