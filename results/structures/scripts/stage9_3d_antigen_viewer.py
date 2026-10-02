"""
Stage 9 — Interactive 3D Antigen-Interface Viewer.
Deterministic representatives near joint cohort medians of
framework_contact_fraction and shieldable_engagement_fraction.
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
CDR3_START, CDR3_END = 105, 117


def build_vhh_viewer(instance, pdb, hchain, extracted_structure_path, long_df):
    path = PROJECT_ROOT / extracted_structure_path
    with open(path) as f:
        cif_data = f.read()

    inst_contacts = long_df[(long_df['INSTANCE'] == instance) & long_df['contact_at_4p5A']]
    fw_contact_positions = sorted(set(inst_contacts[(~inst_contacts['antibody_imgt_position'].between(CDR3_START, CDR3_END)) &
                                                       (inst_contacts['antibody_region'] != 'UNKNOWN') &
                                                       (~inst_contacts['antibody_region'].isin(['CDR1', 'CDR2']))]['antibody_imgt_position']))

    v = py3Dmol.view(width=900, height=700)
    v.addModel(cif_data, 'cif')
    v.setStyle({'chain': hchain}, {'cartoon': {'color': '#4D69B2'}})
    v.setStyle({'chain': hchain, 'resi': [str(i) for i in range(CDR3_START, CDR3_END + 1)]}, {'cartoon': {'color': 'cyan'}})
    v.addStyle({'chain': hchain, 'resi': [str(i) for i in FULL_EMPIRICAL_INTERFACE]}, {'stick': {'color': 'orange', 'radius': 0.25}})
    if fw_contact_positions:
        v.addStyle({'chain': hchain, 'resi': [str(i) for i in fw_contact_positions]}, {'stick': {'color': 'red', 'radius': 0.4}})
    antigen_chains = sorted(inst_contacts['antigen_chain'].unique())
    for ac in antigen_chains:
        v.setStyle({'chain': ac}, {'cartoon': {'color': '#66BB6A', 'opacity': 0.7}})
    v.zoomTo({'chain': hchain})
    html = v._make_html()

    full_html = f"""<!DOCTYPE html>
<html><head><title>Stage 9 VHH-Antigen — {pdb}</title></head><body>
<h2>Stage 9 — VHH-Antigen Interface ({instance})</h2>
<p>Blue = VHH framework/CDR1/CDR2 backbone. Cyan = CDR3. Orange sticks = full 28-position empirical former VH-VL interface (regardless of antigen contact). Red sticks = framework residues actually contacting antigen (&lt;=4.5&Aring;) in this structure. Green (restrained) = native antigen chain(s).</p>
<p>Native antigen-bound structure (antigen NOT computationally removed, unlike Stages 5/8).</p>
{html}
</body></html>"""
    out_path = STRUCTURES_DIR / f"stage9_vhh_antigen_{pdb}_{instance}.html"
    with open(out_path, 'w') as f:
        f.write(full_html)
    logger.info(f"Wrote {out_path}")
    return out_path


def build_fv_viewer(instance, pdb, hchain, lchain, extracted_structure_path, long_df):
    path = PROJECT_ROOT / extracted_structure_path
    with open(path) as f:
        cif_data = f.read()

    inst_contacts = long_df[(long_df['INSTANCE'] == instance) & long_df['contact_at_4p5A']]
    vh_fw_positions = sorted(set(inst_contacts[(inst_contacts['antibody_domain'] == 'VH') &
                                                  (~inst_contacts['antibody_region'].isin(['CDR1', 'CDR2', 'CDR3']))]['antibody_imgt_position']))

    v = py3Dmol.view(width=900, height=700)
    v.addModel(cif_data, 'cif')
    v.setStyle({'chain': hchain}, {'cartoon': {'color': '#8B5FA8'}})
    v.setStyle({'chain': lchain}, {'cartoon': {'color': '#F48FB1'}})
    if vh_fw_positions:
        v.addStyle({'chain': hchain, 'resi': [str(i) for i in vh_fw_positions]}, {'stick': {'color': 'red', 'radius': 0.4}})
    antigen_chains = sorted(inst_contacts['antigen_chain'].unique())
    for ac in antigen_chains:
        v.setStyle({'chain': ac}, {'cartoon': {'color': '#66BB6A', 'opacity': 0.7}})
    v.zoomTo({'chain': hchain})
    html = v._make_html()

    full_html = f"""<!DOCTYPE html>
<html><head><title>Stage 9 Fv-Antigen — {pdb}</title></head><body>
<h2>Stage 9 — Conventional Fv-Antigen Interface ({instance})</h2>
<p>Purple = VH. Pink = VL. Red sticks = VH framework residues contacting antigen (&lt;=4.5&Aring;). Green (restrained) = native antigen chain(s).</p>
{html}
</body></html>"""
    out_path = STRUCTURES_DIR / f"stage9_fv_antigen_{pdb}_{instance}.html"
    with open(out_path, 'w') as f:
        f.write(full_html)
    logger.info(f"Wrote {out_path}")
    return out_path


def main():
    sel = pd.read_csv(PROJECT_ROOT / "results/structures/stage9_representative_selection.csv")
    long_df = pd.read_parquet(PROJECT_ROOT / "data/processed/stage9_antibody_antigen_contacts_long.parquet")

    vhh_meta = pd.read_parquet(PROJECT_ROOT / "data/processed/stage4_vhh_sequence_representatives.parquet")
    vh_meta = pd.read_parquet(PROJECT_ROOT / "data/processed/stage3_paired_representatives.parquet")
    vh_meta = vh_meta[vh_meta['is_primary_representative'] == True]

    vhh_row = sel[sel['group'] == 'VHH'].iloc[0]
    vhh_meta_row = vhh_meta[vhh_meta['INSTANCE'] == vhh_row['INSTANCE']].iloc[0]
    out1 = build_vhh_viewer(vhh_row['INSTANCE'], vhh_row['PDB'], vhh_meta_row['Hchain'], vhh_meta_row['extracted_structure_path'], long_df)

    vh_row = sel[sel['group'] == 'VH'].iloc[0]
    vh_meta_row = vh_meta[vh_meta['INSTANCE'] == vh_row['INSTANCE']].iloc[0]
    out2 = build_fv_viewer(vh_row['INSTANCE'], vh_row['PDB'], vh_meta_row['Hchain'], vh_meta_row['Lchain'], vh_meta_row['extracted_structure_path'], long_df)

    entries = []
    for out_path, row, meta_row in [(out1, vhh_row, vhh_meta_row), (out2, vh_row, vh_meta_row)]:
        entries.append({'instance': row['INSTANCE'], 'pdb': row['PDB'], 'hchain': meta_row['Hchain'],
                         'lchain': meta_row.get('Lchain', 'NA'), 'antigen_chains': '',
                         'interactive_html_path': str(out_path.relative_to(PROJECT_ROOT)), 'static_figure_path': ''})
    new_entries = pd.DataFrame(entries)
    viewer_index_path = PROJECT_ROOT / "results" / "structures" / "viewer_index.csv"
    if viewer_index_path.exists():
        existing = pd.read_csv(viewer_index_path)
        combined = pd.concat([existing, new_entries], ignore_index=True)
    else:
        combined = new_entries
    combined.to_csv(viewer_index_path, index=False)
    logger.info(f"Updated {viewer_index_path}")


if __name__ == "__main__":
    main()
