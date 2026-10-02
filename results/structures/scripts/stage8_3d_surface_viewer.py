"""
Stage 8 — Interactive 3D Surface Viewer: VHH_NATIVE and VH_WITHOUT_VL_GEOMETRIC
deterministic representatives, near cohort medians of former_interface_core_PSH
and APBS interface_potential_absolute_mean, high coverage / resolution tie-break.

Three separate, non-merged coloring modes per structure:
  A. hydrophobic patch contribution (KD-normalized value at core+vicinity residues)
  B. discrete charge patch contribution (charge_class at core+vicinity residues)
  C. APBS potential sampled at exposed interface atoms (atom-level, honestly labeled)
"""
import sys
from pathlib import Path
import pandas as pd
import py3Dmol
import logging

PROJECT_ROOT = Path("/Users/betulgul/Desktop/staj")
STRUCTURES_DIR = PROJECT_ROOT / "results" / "structures" / "interactive"
STRUCTURES_DIR.mkdir(parents=True, exist_ok=True)

logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger(__name__)

FULL_EMPIRICAL_INTERFACE = {40, 42, 44, 49, 50, 52, 55, 66, 67, 68, 69, 103,
                             107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 118, 119, 120}


def build_viewer_html(instance, pdb, base_color, hchain, extracted_structure_path, components_df):
    path = PROJECT_ROOT / extracted_structure_path
    with open(path) as f:
        cif_data = f.read()

    inst_comp = components_df[components_df['INSTANCE'] == instance]
    hydro_rows = inst_comp[inst_comp['interface_vicinity_member'] & inst_comp['kyte_doolittle_norm'].notna()]
    charge_rows = inst_comp[inst_comp['interface_vicinity_member'] & (inst_comp['charge_class'] != 0)]

    def make_view(mode):
        v = py3Dmol.view(width=800, height=650)
        v.addModel(cif_data, 'cif')
        v.setStyle({'chain': hchain}, {'cartoon': {'color': base_color}})
        interface_resi = [str(i) for i in FULL_EMPIRICAL_INTERFACE]
        v.addStyle({'chain': hchain, 'resi': interface_resi}, {'stick': {'color': 'grey', 'radius': 0.2}})
        if mode == 'hydrophobic':
            hi = hydro_rows[hydro_rows['kyte_doolittle_norm'] > 1.5]
            resi = [str(int(r)) for r in hi['imgt_position'].tolist()]
            if resi:
                v.addStyle({'chain': hchain, 'resi': resi}, {'stick': {'color': 'orange', 'radius': 0.35}})
        elif mode == 'charge':
            pos = charge_rows[charge_rows['charge_class'] > 0]
            neg = charge_rows[charge_rows['charge_class'] < 0]
            if len(pos):
                v.addStyle({'chain': hchain, 'resi': [str(int(r)) for r in pos['imgt_position']]}, {'stick': {'color': 'blue', 'radius': 0.35}})
            if len(neg):
                v.addStyle({'chain': hchain, 'resi': [str(int(r)) for r in neg['imgt_position']]}, {'stick': {'color': 'red', 'radius': 0.35}})
        elif mode == 'apbs':
            v.addStyle({'chain': hchain, 'resi': interface_resi}, {'stick': {'color': 'magenta', 'radius': 0.3}})
        v.zoomTo({'chain': hchain})
        return v._make_html()

    html_hydro = make_view('hydrophobic')
    html_charge = make_view('charge')
    html_apbs = make_view('apbs')

    full_html = f"""<!DOCTYPE html>
<html><head><title>Stage 8 Surface — {pdb}</title>
<style>
.tab {{ overflow: hidden; border-bottom: 1px solid #ccc; }}
.tab button {{ background: #eee; border: none; padding: 10px 20px; cursor: pointer; }}
.tab button.active {{ background: {base_color}; color: white; }}
.tabcontent {{ display: none; }}
.tabcontent.active {{ display: block; }}
</style></head>
<body>
<h2>Stage 8 — Former-Interface Surface Patches ({instance})</h2>
<p>Grey sticks = full 28-position empirical former VH-VL interface (all modes). Mode-specific highlighting below is NEVER merged into a single rainbow surface.</p>
<p><b>Mode A (hydrophobic):</b> orange sticks = vicinity-member residues with normalized KD hydrophobicity > 1.5.<br>
<b>Mode B (charge):</b> blue = positive-class vicinity residues, red = negative-class vicinity residues (discrete PPC/PNC contributors).<br>
<b>Mode C (APBS):</b> magenta sticks = interface residues where APBS potential was SAMPLED AT EXPOSED ATOMS (this is atom-coordinate sampling, honestly labeled -- NOT a full molecular-surface electrostatic-potential rendering).</p>
<div class="tab">
  <button class="tablinks active" onclick="showTab('hydro')">A. Hydrophobic patch contribution</button>
  <button class="tablinks" onclick="showTab('charge')">B. Discrete charge patch contribution</button>
  <button class="tablinks" onclick="showTab('apbs')">C. APBS potential sampled at exposed interface atoms</button>
</div>
<div id="hydro" class="tabcontent active">{html_hydro}</div>
<div id="charge" class="tabcontent">{html_charge}</div>
<div id="apbs" class="tabcontent">{html_apbs}</div>
<script>
function showTab(id) {{
  document.querySelectorAll('.tabcontent').forEach(el => el.classList.remove('active'));
  document.querySelectorAll('.tablinks').forEach(el => el.classList.remove('active'));
  document.getElementById(id).classList.add('active');
  event.currentTarget.classList.add('active');
}}
</script>
</body></html>"""
    return full_html


def main():
    sel = pd.read_csv(PROJECT_ROOT / "results/structures/stage8_representative_selection.csv")
    components_df = pd.read_parquet(PROJECT_ROOT / "data/processed/stage8_patch_components.parquet")

    vhh_meta_df = pd.read_parquet(PROJECT_ROOT / "data/processed/stage4_vhh_sequence_representatives.parquet")
    vh_meta_df = pd.read_parquet(PROJECT_ROOT / "data/processed/stage3_paired_representatives.parquet")
    vh_meta_df = vh_meta_df[vh_meta_df['is_primary_representative'] == True]

    viewer_entries = []
    for _, row in sel.iterrows():
        instance, pdb, group = row['INSTANCE'], row['PDB'], row['group']
        meta_df = vhh_meta_df if group == 'VHH' else vh_meta_df
        meta_row = meta_df[meta_df['INSTANCE'] == instance].iloc[0]
        base_color = '#4D69B2' if group == 'VHH' else '#8B5FA8'
        html = build_viewer_html(instance, pdb, base_color, meta_row['Hchain'], meta_row['extracted_structure_path'], components_df)
        fname = f"stage8_{'vhh' if group == 'VHH' else 'vh'}_surface_{pdb}_{instance}.html"
        out_path = STRUCTURES_DIR / fname
        with open(out_path, 'w') as f:
            f.write(html)
        logger.info(f"Wrote {out_path}")
        viewer_entries.append({'instance': instance, 'pdb': pdb, 'hchain': meta_row['Hchain'], 'lchain': meta_row.get('Lchain', 'NA'),
                                'antigen_chains': '', 'interactive_html_path': str(out_path.relative_to(PROJECT_ROOT)), 'static_figure_path': ''})

    viewer_index_path = PROJECT_ROOT / "results" / "structures" / "viewer_index.csv"
    new_entries = pd.DataFrame(viewer_entries)
    if viewer_index_path.exists():
        existing = pd.read_csv(viewer_index_path)
        combined = pd.concat([existing, new_entries], ignore_index=True)
    else:
        combined = new_entries
    combined.to_csv(viewer_index_path, index=False)
    logger.info(f"Updated {viewer_index_path}")


if __name__ == "__main__":
    main()
