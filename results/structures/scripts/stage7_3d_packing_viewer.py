"""
Stage 7 — Interactive 3D CDR3-Target Packing Viewer

Deterministic selection: pdb_00004grw-E, closest simultaneously to cohort
medians of shielding_fraction and BSA_pair_total (see
stage7_representative_selection.csv). Distinguishes native structure / CDR3 /
NONADJACENT_PACKING_TARGET / positions actually involved in a specific
noncovalent interaction from mere spatial proximity.
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

VHH_GREY = "#B0B0B0"
CDR3_TEAL = "#00BCD4"
TARGET_ORANGE = "#FFB300"
SPECIFIC_RED = "#E53935"
CDR3_START, CDR3_END = 105, 117
NONADJACENT_PACKING_TARGET = {40, 42, 44, 49, 50, 52, 55, 66, 67, 68, 69, 103, 119, 120}


def main():
    rep_sel = pd.read_csv(PROJECT_ROOT / "results/structures/stage7_representative_selection.csv")
    rep = rep_sel.iloc[0]
    instance = rep['INSTANCE']

    vhh_meta = pd.read_parquet(PROJECT_ROOT / "data/processed/stage4_vhh_sequence_representatives.parquet")
    meta_row = vhh_meta[vhh_meta['INSTANCE'] == instance].iloc[0]
    hchain = meta_row['Hchain']
    path = PROJECT_ROOT / meta_row['extracted_structure_path']
    with open(path) as f:
        cif_data = f.read()

    long_df = pd.read_csv(PROJECT_ROOT / "results/tables/stage7_interaction_records_long.csv", low_memory=False)
    inst_records = long_df[long_df['INSTANCE'] == instance]
    specific = inst_records[inst_records['include_primary_noncovalent']]
    specific_target_positions = sorted(set(specific['target_imgt_position']))
    specific_cdr3_positions = sorted(set(specific['cdr3_imgt_position']))

    view = py3Dmol.view(width=900, height=700)
    view.addModel(cif_data, 'cif')
    view.setStyle({'chain': hchain}, {'cartoon': {'color': VHH_GREY}})
    view.setStyle({'chain': hchain, 'resi': [str(i) for i in range(CDR3_START, CDR3_END + 1)]}, {'cartoon': {'color': CDR3_TEAL}})
    view.addStyle({'chain': hchain, 'resi': [str(i) for i in NONADJACENT_PACKING_TARGET]}, {'stick': {'color': TARGET_ORANGE, 'radius': 0.25}})
    if specific_target_positions:
        view.addStyle({'chain': hchain, 'resi': [str(i) for i in specific_target_positions]}, {'stick': {'color': SPECIFIC_RED, 'radius': 0.35}})
    if specific_cdr3_positions:
        view.addStyle({'chain': hchain, 'resi': [str(i) for i in specific_cdr3_positions]}, {'stick': {'color': SPECIFIC_RED, 'radius': 0.35}})
    view.zoomTo({'chain': hchain})
    html_view = view._make_html()

    bsa_val = rep['nonadj_BSA_pair_total']
    shield_val = rep['shielding_fraction']

    full_html = f"""<!DOCTYPE html>
<html>
<head><title>Stage 7 CDR3-Target Packing — {meta_row['PDB']}</title></head>
<body>
<h2>Stage 7 — CDR3-Target Intradomain Packing (Arpeggio-typed)</h2>
<p><b>PDB:</b> {meta_row['PDB']} | <b>INSTANCE:</b> {instance}</p>
<p>Grey = VHH backbone. Teal = CDR3 (IMGT 105-117). Orange sticks = full
NONADJACENT_PACKING_TARGET set (14 former-interface positions, IMGT 118
excluded as directly covalently bonded to CDR3 -- see
stage7_covalent_adjacent_policy.md) regardless of whether a specific
interaction was detected there in THIS structure. Red sticks = residues
(CDR3 and/or target) that participate in at least one Arpeggio-typed specific
noncovalent interaction (hbond/ionic/aromatic/hydrophobic/polar/vdw/etc,
excluding 'proximal'-only and covalent records) in this specific structure --
this is the distinction between mere spatial proximity (orange alone) and an
actual typed noncovalent contact (red).</p>
<p>This structure's values: BSA_pair_total = {bsa_val:.1f} &Aring;&sup2;,
shielding_fraction = {shield_val:.3f} — selected deterministically as
simultaneously closest to both cohort medians (not cherry-picked).</p>
<p>n specific noncovalent CDR3-target atom records in this structure: {len(specific)}</p>
{html_view}
</body>
</html>"""

    out_path = STRUCTURES_DIR / f"stage7_cdr3_packing_{meta_row['PDB']}_{instance}.html"
    with open(out_path, 'w') as f:
        f.write(full_html)
    logger.info(f"Wrote {out_path}")

    viewer_index_path = PROJECT_ROOT / "results" / "structures" / "viewer_index.csv"
    new_entry = pd.DataFrame([{
        'instance': instance, 'pdb': meta_row['PDB'], 'hchain': hchain, 'lchain': 'NA',
        'antigen_chains': '', 'interactive_html_path': str(out_path.relative_to(PROJECT_ROOT)), 'static_figure_path': '',
    }])
    if viewer_index_path.exists():
        existing = pd.read_csv(viewer_index_path)
        combined = pd.concat([existing, new_entry], ignore_index=True)
    else:
        combined = new_entry
    combined.to_csv(viewer_index_path, index=False)
    logger.info(f"Updated {viewer_index_path}")

    return out_path


if __name__ == "__main__":
    main()
