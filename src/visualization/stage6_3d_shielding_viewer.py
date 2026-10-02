"""
Stage 6 — Interactive 3D Native-vs-No-CDR3 Viewer

Deterministic selection: closest to cohort medians of shielding_fraction and
delta_sasa_shield, high coordinate coverage, resolution as tie-break.
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
CDR3_TEAL = "#00BCD4"
INTERFACE_ORANGE = "#FFB300"
CDR3_START, CDR3_END = 105, 117


def select_representative(summaries: pd.DataFrame, vhh_meta: pd.DataFrame) -> pd.Series:
    sub = summaries[summaries['shieldable_interface_coordinate_coverage'] >= 0.90].copy()
    med_frac = sub['shielding_fraction'].median()
    med_delta = sub['delta_sasa_shield'].median()

    sub['dist'] = np.sqrt(
        ((sub['shielding_fraction'] - med_frac) / sub['shielding_fraction'].std()) ** 2 +
        ((sub['delta_sasa_shield'] - med_delta) / sub['delta_sasa_shield'].std()) ** 2
    )
    sub = sub.merge(vhh_meta[['INSTANCE', 'resolution']], on='INSTANCE', how='left')
    sub['res_num'] = pd.to_numeric(sub['resolution'], errors='coerce')
    sub = sub.sort_values(['dist', 'res_num', 'INSTANCE'])
    return sub.iloc[0]


def main():
    summaries = pd.read_parquet(PROJECT_ROOT / "data/processed/stage6_vhh_cdr3_shielding.parquet")
    residues = pd.read_parquet(PROJECT_ROOT / "data/processed/stage6_cdr3_shielding_residues.parquet")
    vhh_meta = pd.read_parquet(PROJECT_ROOT / "data/processed/stage4_vhh_sequence_representatives.parquet")

    rep = select_representative(summaries, vhh_meta)
    meta_row = vhh_meta[vhh_meta['INSTANCE'] == rep['INSTANCE']].iloc[0]
    hchain = meta_row['Hchain']
    path = PROJECT_ROOT / meta_row['extracted_structure_path']

    with open(path) as f:
        cif_data = f.read()

    pos_data = residues[(residues['INSTANCE'] == rep['INSTANCE']) & residues['delta_sasa_shield'].notna()]
    max_delta = pos_data['delta_sasa_shield'].abs().max() or 1.0

    # Native view: full VHH with CDR3 highlighted and shieldable interface colored by delta_sasa_shield
    view_native = py3Dmol.view(width=900, height=700)
    view_native.addModel(cif_data, 'cif')
    view_native.setStyle({'chain': hchain}, {'cartoon': {'color': VHH_BLUE}})
    view_native.setStyle({'chain': hchain, 'resi': [str(i) for i in range(CDR3_START, CDR3_END + 1)]}, {'cartoon': {'color': CDR3_TEAL}})
    for _, prow in pos_data.iterrows():
        intensity = min(1.0, abs(prow['delta_sasa_shield']) / max_delta)
        r, g, b = 255, int(255 - 140 * intensity), 0
        color = f"#{r:02x}{g:02x}{b:02x}"
        resi = f"{int(prow['imgt_position'])}{prow['insertion_code']}" if prow['insertion_code'] else int(prow['imgt_position'])
        try:
            view_native.addStyle({'chain': hchain, 'resi': str(resi)}, {'cartoon': {'color': color}})
        except Exception:
            pass
    view_native.zoomTo({'chain': hchain})
    html_native = view_native._make_html()

    # No-CDR3 view: same structure, CDR3 hidden
    view_nocdr3 = py3Dmol.view(width=900, height=700)
    view_nocdr3.addModel(cif_data, 'cif')
    view_nocdr3.setStyle({'chain': hchain}, {'cartoon': {'color': VHH_BLUE}})
    view_nocdr3.setStyle({'chain': hchain, 'resi': [str(i) for i in range(CDR3_START, CDR3_END + 1)]}, {'cartoon': {'hidden': True}})
    for _, prow in pos_data.iterrows():
        intensity = min(1.0, abs(prow['delta_sasa_shield']) / max_delta)
        r, g, b = 255, int(255 - 140 * intensity), 0
        color = f"#{r:02x}{g:02x}{b:02x}"
        resi = f"{int(prow['imgt_position'])}{prow['insertion_code']}" if prow['insertion_code'] else int(prow['imgt_position'])
        try:
            view_nocdr3.addStyle({'chain': hchain, 'resi': str(resi)}, {'cartoon': {'color': color}})
        except Exception:
            pass
    view_nocdr3.zoomTo({'chain': hchain})
    html_nocdr3 = view_nocdr3._make_html()

    full_html = f"""<!DOCTYPE html>
<html>
<head><title>Stage 6 CDR3 Shielding — {meta_row['PDB']}</title>
<style>
.tab {{ overflow: hidden; border-bottom: 1px solid #ccc; }}
.tab button {{ background: #eee; border: none; padding: 10px 20px; cursor: pointer; font-size: 14px; }}
.tab button.active {{ background: {VHH_BLUE}; color: white; }}
.tabcontent {{ display: none; }}
.tabcontent.active {{ display: block; }}
</style>
</head>
<body>
<h2>Stage 6 — Native VHH vs Same-Coordinate No-CDR3 SASA Counterfactual</h2>
<p><b>PDB:</b> {meta_row['PDB']} | <b>INSTANCE:</b> {rep['INSTANCE']}</p>
<p><b>IMPORTANT:</b> The "CDR3 hidden" view is a <u>geometric SASA counterfactual</u> —
the exact same experimentally observed coordinates with CDR3 atoms omitted from the
solvent-accessibility calculation only. It is NOT a physical CDR3-deletion mutant,
not relaxed, not repacked, not molecular-dynamics-simulated.</p>
<p>VHH base = blue. CDR3 (IMGT 105-117) = teal (native view only). Shieldable
former-interface residues (the 15 non-CDR3 positions of the 28-position empirical
interface) colored orange by this structure's own delta_sasa_shield magnitude.</p>
<p>This structure's values: delta_sasa_shield = {rep['delta_sasa_shield']:.1f} &Aring;&sup2;,
shielding_fraction = {rep['shielding_fraction']:.3f} — selected deterministically as
closest to the cohort medians (not cherry-picked), among &ge;90% coordinate-coverage
structures, resolution as tie-break.</p>
<div class="tab">
  <button class="tablinks active" onclick="showTab('native')">A. Native VHH (CDR3 visible)</button>
  <button class="tablinks" onclick="showTab('nocdr3')">B. No-CDR3 SASA counterfactual (CDR3 hidden)</button>
</div>
<div id="native" class="tabcontent active">{html_native}</div>
<div id="nocdr3" class="tabcontent">{html_nocdr3}</div>
<script>
function showTab(id) {{
  document.querySelectorAll('.tabcontent').forEach(el => el.classList.remove('active'));
  document.querySelectorAll('.tablinks').forEach(el => el.classList.remove('active'));
  document.getElementById(id).classList.add('active');
  event.currentTarget.classList.add('active');
}}
</script>
</body>
</html>"""

    out_path = STRUCTURES_DIR / f"stage6_cdr3_shielding_{meta_row['PDB']}_{rep['INSTANCE']}.html"
    with open(out_path, 'w') as f:
        f.write(full_html)
    logger.info(f"Wrote {out_path}")

    # Selection audit
    selection_df = pd.DataFrame([{
        'INSTANCE': rep['INSTANCE'], 'PDB': meta_row['PDB'],
        'delta_sasa_shield': rep['delta_sasa_shield'], 'shielding_fraction': rep['shielding_fraction'],
        'selection_rule': 'closest to cohort medians (shielding_fraction, delta_sasa_shield), >=90% coord coverage, resolution tie-break',
    }])
    selection_out = PROJECT_ROOT / "results" / "structures" / "stage6_representative_selection.csv"
    selection_df.to_csv(selection_out, index=False)
    logger.info(f"Wrote {selection_out}")

    viewer_index_path = PROJECT_ROOT / "results" / "structures" / "viewer_index.csv"
    new_entry = pd.DataFrame([{
        'instance': rep['INSTANCE'], 'pdb': meta_row['PDB'], 'hchain': hchain, 'lchain': 'NA',
        'antigen_chains': '', 'interactive_html_path': str(out_path.relative_to(PROJECT_ROOT)), 'static_figure_path': '',
    }])
    if viewer_index_path.exists():
        existing = pd.read_csv(viewer_index_path)
        combined = pd.concat([existing, new_entry], ignore_index=True)
    else:
        combined = new_entry
    combined.to_csv(viewer_index_path, index=False)
    logger.info(f"Updated {viewer_index_path}")

    return out_path, rep


if __name__ == "__main__":
    main()
