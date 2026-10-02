"""
Stage 5 — Interactive 3D Native Surface Chemistry Viewers

Deterministic representative selection: closest to cohort medians of
interface_sasa_weighted_mean_KD and interface_total_sasa, high coordinate
coverage, resolution as secondary tie-break. Residue-level exposed
chemistry, NOT a hydrophobic patch map.
"""

import sys
from pathlib import Path
import pandas as pd
import numpy as np
import py3Dmol
import logging

PROJECT_ROOT = Path(__file__).parent.parent.parent
STRUCTURES_DIR = PROJECT_ROOT / "results" / "structures" / "interactive"
STRUCTURES_DIR.mkdir(parents=True, exist_ok=True)

logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger(__name__)

VHH_BLUE = "#4D69B2"
VH_PURPLE = "#8670B1"


def select_representative(seq_chem: pd.DataFrame, structural_state: str, meta: pd.DataFrame) -> pd.Series:
    sub = seq_chem[seq_chem['structural_state'] == structural_state].copy()
    sub = sub[sub['interface_coordinate_coverage_given_present'] >= 0.90]

    med_kd = sub['interface_sasa_weighted_mean_KD'].median()
    med_sasa = sub['interface_total_sasa'].median()

    sub['dist'] = np.sqrt(
        ((sub['interface_sasa_weighted_mean_KD'] - med_kd) / sub['interface_sasa_weighted_mean_KD'].std()) ** 2 +
        ((sub['interface_total_sasa'] - med_sasa) / sub['interface_total_sasa'].std()) ** 2
    )

    sub = sub.merge(meta[['INSTANCE', 'resolution']], on='INSTANCE', how='left')
    sub['res_num'] = pd.to_numeric(sub['resolution'], errors='coerce')
    sub = sub.sort_values(['dist', 'res_num', 'INSTANCE'])
    return sub.iloc[0]


def render_vhh(row: pd.Series, vhh_meta: pd.DataFrame, combined: pd.DataFrame) -> Path:
    meta_row = vhh_meta[vhh_meta['INSTANCE'] == row['INSTANCE']].iloc[0]
    path = PROJECT_ROOT / meta_row['extracted_structure_path']
    with open(path) as f:
        cif_data = f.read()
    hchain = meta_row['Hchain']

    pos_data = combined[(combined['INSTANCE'] == row['INSTANCE']) & (combined['structural_state'] == 'VHH_NATIVE') &
                          combined['sasa'].notna()]

    view = py3Dmol.view(width=900, height=700)
    view.addModel(cif_data, 'cif')
    view.setStyle({}, {'cartoon': {'color': 'lightgray'}})
    view.setStyle({'chain': hchain}, {'cartoon': {'color': VHH_BLUE}})

    max_sasa_kd = (pos_data['sasa'] * pos_data['kyte_doolittle'].fillna(0)).abs().max() or 1.0
    for _, prow in pos_data.iterrows():
        if pd.isna(prow['kyte_doolittle']):
            continue
        contribution = prow['sasa'] * prow['kyte_doolittle']
        intensity = min(1.0, abs(contribution) / max_sasa_kd)
        if contribution > 0:
            color = f"#{255:02x}{int(255 - 180*intensity):02x}{0:02x}"  # yellow->orange (hydrophobic)
        else:
            color = f"#{int(255-180*intensity):02x}{int(255-100*intensity):02x}{255:02x}"  # blue-ish (hydrophilic)
        resi = f"{int(prow['imgt_position'])}{prow['insertion_code']}" if prow['insertion_code'] else int(prow['imgt_position'])
        try:
            view.addStyle({'chain': hchain, 'resi': str(resi)}, {'cartoon': {'color': color}})
        except Exception:
            pass

    view.zoomTo()
    html_body = view._make_html()

    full_html = f"""<!DOCTYPE html>
<html>
<head><title>Stage 5 VHH Native Surface — {meta_row['PDB']}</title></head>
<body>
<h2>Stage 5 — VHH Native Structural Surface Chemistry</h2>
<p><b>PDB:</b> {meta_row['PDB']} | <b>INSTANCE:</b> {meta_row['INSTANCE']}</p>
<p>Base VHH variable domain (blue). Empirical former-interface residues colored by
residue-level SASA-weighted hydropathy contribution (SASA_i &times; KD_i):
orange/yellow = hydrophobic exposed area, blue = hydrophilic exposed area.
This is residue-level exposed chemistry, NOT a hydrophobic patch/cluster map.</p>
<p>Selected deterministically as the structure closest to the cohort medians of
interface_sasa_weighted_mean_KD ({row['interface_sasa_weighted_mean_KD']:.2f}) and
interface_total_sasa ({row['interface_total_sasa']:.1f} &Aring;&sup2;), among
&ge;90% coordinate-coverage structures, resolution as tie-break.</p>
{html_body}
</body>
</html>"""

    out_path = STRUCTURES_DIR / f"stage5_vhh_native_surface_{meta_row['PDB']}_{meta_row['INSTANCE']}.html"
    with open(out_path, 'w') as f:
        f.write(full_html)
    logger.info(f"Wrote {out_path}")
    return out_path, meta_row


def render_vh(row: pd.Series, vh_meta: pd.DataFrame, combined: pd.DataFrame) -> Path:
    meta_row = vh_meta[vh_meta['INSTANCE'] == row['INSTANCE']].iloc[0]
    path = PROJECT_ROOT / meta_row['extracted_structure_path']
    with open(path) as f:
        cif_data = f.read()
    hchain = meta_row['Hchain']
    lchain = meta_row['Lchain']

    pos_data = combined[(combined['INSTANCE'] == row['INSTANCE']) & (combined['structural_state'] == 'VH_WITHOUT_VL_GEOMETRIC') &
                          combined['sasa'].notna()]

    # Render only the VH chain (VL physically present in file but hidden/not styled to represent the geometric-removal reference)
    view = py3Dmol.view(width=900, height=700)
    view.addModel(cif_data, 'cif')
    view.setStyle({'chain': hchain}, {'cartoon': {'color': 'lightgray'}})
    view.setStyle({'chain': lchain}, {'cartoon': {'hidden': True}})

    max_sasa_kd = (pos_data['sasa'] * pos_data['kyte_doolittle'].fillna(0)).abs().max() or 1.0
    for _, prow in pos_data.iterrows():
        if pd.isna(prow['kyte_doolittle']):
            continue
        contribution = prow['sasa'] * prow['kyte_doolittle']
        intensity = min(1.0, abs(contribution) / max_sasa_kd)
        if contribution > 0:
            color = f"#{255:02x}{int(255 - 180*intensity):02x}{0:02x}"
        else:
            color = f"#{int(255-180*intensity):02x}{int(255-100*intensity):02x}{255:02x}"
        resi = f"{int(prow['imgt_position'])}{prow['insertion_code']}" if prow['insertion_code'] else int(prow['imgt_position'])
        try:
            view.addStyle({'chain': hchain, 'resi': str(resi)}, {'cartoon': {'color': color}})
        except Exception:
            pass

    view.zoomTo({'chain': hchain})
    html_body = view._make_html()

    full_html = f"""<!DOCTYPE html>
<html>
<head><title>Stage 5 VH-Without-VL Surface — {meta_row['PDB']}</title></head>
<body>
<h2>Stage 5 — Conventional VH Geometric-Counterfactual Surface Chemistry</h2>
<p><b>PDB:</b> {meta_row['PDB']} | <b>INSTANCE:</b> {meta_row['INSTANCE']}</p>
<p><b>IMPORTANT:</b> VL is hidden/removed ONLY for this geometric SASA reference
calculation (VH_WITHOUT_VL_GEOMETRIC). This is a same-coordinate computational
counterfactual, NOT an experimentally isolated or equilibrated VH structure.
The physical VL chain remains in the deposited coordinates but is not displayed.</p>
<p>VH variable domain (gray). Empirical former-interface residues colored by
residue-level SASA-weighted hydropathy contribution in this VL-removed geometric
state: orange/yellow = hydrophobic exposed area, blue = hydrophilic exposed area.</p>
<p>Selected deterministically as the structure closest to the cohort medians of
interface_sasa_weighted_mean_KD ({row['interface_sasa_weighted_mean_KD']:.2f}) and
interface_total_sasa ({row['interface_total_sasa']:.1f} &Aring;&sup2;).</p>
{html_body}
</body>
</html>"""

    out_path = STRUCTURES_DIR / f"stage5_vh_without_vl_surface_{meta_row['PDB']}_{meta_row['INSTANCE']}.html"
    with open(out_path, 'w') as f:
        f.write(full_html)
    logger.info(f"Wrote {out_path}")
    return out_path, meta_row


def main():
    seq_chem = pd.read_parquet(PROJECT_ROOT / "data/processed/stage5_sequence_surface_chemistry.parquet")
    combined = pd.read_parquet(PROJECT_ROOT / "data/processed/stage5_interface_surface_chemistry_residues.parquet")
    vhh_meta = pd.read_parquet(PROJECT_ROOT / "data/processed/stage4_vhh_sequence_representatives.parquet")
    vh_meta = pd.read_parquet(PROJECT_ROOT / "data/processed/stage3_paired_representatives.parquet")

    vhh_rep = select_representative(seq_chem, 'VHH_NATIVE', vhh_meta)
    vh_rep = select_representative(seq_chem, 'VH_WITHOUT_VL_GEOMETRIC', vh_meta)

    logger.info(f"VHH representative: {vhh_rep['INSTANCE']}")
    logger.info(f"VH representative: {vh_rep['INSTANCE']}")

    vhh_path, vhh_meta_row = render_vhh(vhh_rep, vhh_meta, combined)
    vh_path, vh_meta_row = render_vh(vh_rep, vh_meta, combined)

    selection_df = pd.DataFrame([
        {'category': 'VHH_NATIVE', 'INSTANCE': vhh_rep['INSTANCE'], 'PDB': vhh_meta_row['PDB'],
         'interface_sasa_weighted_mean_KD': vhh_rep['interface_sasa_weighted_mean_KD'],
         'interface_total_sasa': vhh_rep['interface_total_sasa'],
         'selection_rule': 'closest to cohort medians (KD, total SASA), >=90% coord coverage, resolution tie-break'},
        {'category': 'VH_WITHOUT_VL_GEOMETRIC', 'INSTANCE': vh_rep['INSTANCE'], 'PDB': vh_meta_row['PDB'],
         'interface_sasa_weighted_mean_KD': vh_rep['interface_sasa_weighted_mean_KD'],
         'interface_total_sasa': vh_rep['interface_total_sasa'],
         'selection_rule': 'closest to cohort medians (KD, total SASA), >=90% coord coverage, resolution tie-break'},
    ])
    selection_out = PROJECT_ROOT / "results" / "structures" / "stage5_representative_selection.csv"
    selection_df.to_csv(selection_out, index=False)
    logger.info(f"Wrote {selection_out}")

    viewer_index_path = PROJECT_ROOT / "results" / "structures" / "viewer_index.csv"
    new_entries = pd.DataFrame([
        {'instance': vhh_rep['INSTANCE'], 'pdb': vhh_meta_row['PDB'], 'hchain': vhh_meta_row['Hchain'],
         'lchain': 'NA', 'antigen_chains': '', 'interactive_html_path': str(vhh_path.relative_to(PROJECT_ROOT)),
         'static_figure_path': ''},
        {'instance': vh_rep['INSTANCE'], 'pdb': vh_meta_row['PDB'], 'hchain': vh_meta_row['Hchain'],
         'lchain': vh_meta_row['Lchain'], 'antigen_chains': '', 'interactive_html_path': str(vh_path.relative_to(PROJECT_ROOT)),
         'static_figure_path': ''},
    ])
    if viewer_index_path.exists():
        existing = pd.read_csv(viewer_index_path)
        combined_idx = pd.concat([existing, new_entries], ignore_index=True)
    else:
        combined_idx = new_entries
    combined_idx.to_csv(viewer_index_path, index=False)
    logger.info(f"Updated {viewer_index_path}")

    return vhh_path, vh_path


if __name__ == "__main__":
    main()
