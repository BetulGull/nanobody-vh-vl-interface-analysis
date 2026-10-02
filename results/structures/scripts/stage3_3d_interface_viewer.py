"""
Stage 3 — Interactive 3D Empirical VH-VL Interface Viewer

Deterministic representative selection: closest to cohort medians of
total_vh_delta_sasa_due_to_vl and n_interface_residues_5 among PRIMARY
representative structures, resolution as secondary tie-break.
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

VH_PURPLE = "#8670B1"
VL_PINK = "#EF7F7E"


def select_representative_structure(summary_df: pd.DataFrame) -> pd.Series:
    primary = summary_df[summary_df['is_primary_representative'] == True].copy()
    primary = primary[primary['sasa_qc_pass'] == True]

    med_total = primary['total_vh_delta_sasa_due_to_vl'].median()
    med_n5 = primary['n_interface_residues_5'].median()

    primary['dist'] = np.sqrt(
        (primary['total_vh_delta_sasa_due_to_vl'] - med_total) ** 2 +
        (primary['n_interface_residues_5'] - med_n5) ** 2
    )
    primary = primary.sort_values(['dist', 'resolution', 'INSTANCE'])
    return primary.iloc[0]


def render(row: pd.Series, position_map: pd.DataFrame, out_name: str):
    master = pd.read_parquet(PROJECT_ROOT / "data/processed/stage2_structure_qc_master.parquet")
    meta = master[master['INSTANCE'] == row['INSTANCE']].iloc[0]
    extracted_path = PROJECT_ROOT / meta['extracted_structure_path']

    with open(extracted_path) as f:
        cif_data = f.read()

    hchain = meta['Hchain']
    lchain = meta['Lchain']

    freq_by_pos = dict(zip(
        zip(position_map['imgt_position'], position_map['insertion_code'].fillna('')),
        position_map['interface_frequency_5']
    ))

    view = py3Dmol.view(width=950, height=750)
    view.addModel(cif_data, 'cif')
    view.setStyle({}, {'cartoon': {'color': 'lightgray'}})
    view.setStyle({'chain': lchain}, {'cartoon': {'color': VL_PINK}})
    view.setStyle({'chain': hchain}, {'cartoon': {'color': VH_PURPLE}})

    # Color VH residues by population interface_frequency_5 (yellow->orange)
    for (imgt_pos, ins), freq in freq_by_pos.items():
        if freq is None or freq <= 0:
            continue
        resi = f"{imgt_pos}{ins}" if ins else int(imgt_pos)
        intensity = min(1.0, freq)
        r = 255
        g = int(255 - 140 * intensity)
        b = 0
        color = f"#{r:02x}{g:02x}{b:02x}"
        try:
            view.addStyle({'chain': hchain, 'resi': str(resi)}, {'cartoon': {'color': color}})
        except Exception:
            pass

    view.zoomTo()
    html_body = view._make_html()

    full_html = f"""<!DOCTYPE html>
<html>
<head><title>Stage 3 Empirical VH-VL Interface — {meta['PDB']}</title></head>
<body>
<h2>Stage 3 — Empirical VH-VL Interface (Representative Structure)</h2>
<p><b>PDB:</b> {meta['PDB']} | <b>INSTANCE:</b> {meta['INSTANCE']} | <b>SABDAB_ID:</b> {meta.get('SABDAB_ID','N/A')}</p>
<p><b>VH</b> (Hchain {hchain}) = Gordon purple. <b>VL</b> (Lchain {lchain}) = Gordon pink.
VH residues additionally colored yellow&rarr;orange by <b>POPULATION interface_frequency_5</b>
(fraction of structures across the primary nonredundant cohort in which that exact VH IMGT
position had &Delta;SASA &ge; 5 &Aring;&sup2; due to VL).</p>
<p><b>IMPORTANT:</b> the coloring shown is a <u>population-level empirical interface frequency</u>
projected onto this one representative VH, NOT a measurement unique to this specific structure.
This structure's own values: total ΔSASA due to VL = {row['total_vh_delta_sasa_due_to_vl']:.1f} &Aring;&sup2;,
n interface residues (&ge;5 &Aring;&sup2;) = {int(row['n_interface_residues_5'])}. This structure was
selected deterministically as the one closest to the cohort medians of these two structure-specific
metrics — it did not generate the population frequency shown in the coloring.</p>
{html_body}
</body>
</html>"""

    out_path = STRUCTURES_DIR / out_name
    with open(out_path, 'w') as f:
        f.write(full_html)
    logger.info(f"Wrote {out_path}")
    return out_path, meta


def main():
    summary_df = pd.read_parquet(PROJECT_ROOT / "data/processed/stage3_structure_delta_sasa_summary.parquet")
    position_map = pd.read_parquet(PROJECT_ROOT / "data/processed/stage3_empirical_vh_vl_interface_primary.parquet")

    rep_row = select_representative_structure(summary_df)
    logger.info(f"Selected representative: {rep_row['INSTANCE']} (PDB {rep_row['PDB']})")
    logger.info(f"  total_vh_delta_sasa_due_to_vl = {rep_row['total_vh_delta_sasa_due_to_vl']:.2f}")
    logger.info(f"  n_interface_residues_5 = {rep_row['n_interface_residues_5']}")

    out_name = f"stage3_empirical_vh_vl_interface_{rep_row['PDB']}_{rep_row['INSTANCE']}.html"
    html_path, meta = render(rep_row, position_map, out_name)

    # Update viewer_index.csv
    viewer_index_path = PROJECT_ROOT / "results" / "structures" / "viewer_index.csv"
    new_entry = pd.DataFrame([{
        'instance': rep_row['INSTANCE'], 'pdb': rep_row['PDB'],
        'hchain': meta['Hchain'], 'lchain': meta['Lchain'],
        'antigen_chains': '', 'interactive_html_path': str(html_path.relative_to(PROJECT_ROOT)),
        'static_figure_path': '',
    }])
    if viewer_index_path.exists():
        existing = pd.read_csv(viewer_index_path)
        combined = pd.concat([existing, new_entry], ignore_index=True)
    else:
        combined = new_entry
    combined.to_csv(viewer_index_path, index=False)
    logger.info(f"Updated {viewer_index_path}")

    return html_path, rep_row


if __name__ == "__main__":
    main()
