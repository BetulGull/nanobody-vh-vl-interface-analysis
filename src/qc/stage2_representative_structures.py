"""
Stage-2 Representative Structure Selection and 3D Visualization

Deterministic selection (lexicographically smallest INSTANCE among fully
Gordon-passing structures in each final cohort — no aesthetic cherry-picking)
of one VHH-antigen and one paired VH-VL-antigen complex, rendered as
standalone interactive HTML via py3Dmol.
"""

import sys
from pathlib import Path
import pandas as pd
import py3Dmol
import logging

sys.path.insert(0, str(Path(__file__).parent.parent))

PROJECT_ROOT = Path(__file__).parent.parent.parent
STRUCTURES_DIR = PROJECT_ROOT / "results" / "structures" / "interactive"
STRUCTURES_DIR.mkdir(parents=True, exist_ok=True)

logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger(__name__)


def select_representative(df: pd.DataFrame) -> pd.Series:
    """Deterministic selection: lexicographically smallest INSTANCE. No aesthetic choice."""
    sorted_df = df.sort_values("INSTANCE").reset_index(drop=True)
    return sorted_df.iloc[0]


def render_structure(row: pd.Series, out_name: str, title: str):
    """Render a structure with H(/L)chain and antigen chains highlighted."""
    extracted_path = PROJECT_ROOT / row['extracted_structure_path']
    with open(extracted_path) as f:
        cif_data = f.read()

    hchain = row['Hchain']
    lchain = row['Lchain']
    lchain = None if (pd.isna(lchain) or str(lchain).strip() == 'NA') else lchain
    antigen_chains = str(row.get('antigen_chains', '')).split(';')
    antigen_chains = [c for c in antigen_chains if c]

    view = py3Dmol.view(width=900, height=700)
    view.addModel(cif_data, 'cif')
    view.setStyle({}, {'cartoon': {'color': 'lightgray'}})

    view.setStyle({'chain': hchain}, {'cartoon': {'color': '#4D69B2'}})  # VHH/VH blue
    if lchain:
        view.setStyle({'chain': lchain}, {'cartoon': {'color': '#8670B1'}})  # VL purple
    for ac in antigen_chains:
        view.setStyle({'chain': ac}, {'cartoon': {'color': '#6BAA63'}})  # antigen green

    view.zoomTo()

    html = view._make_html()
    full_html = f"""<!DOCTYPE html>
<html>
<head><title>{title}</title></head>
<body>
<h2>{title}</h2>
<p>PDB: {row['PDB']} | INSTANCE: {row['INSTANCE']} | SABDAB_ID: {row.get('SABDAB_ID', 'N/A')}</p>
<p>Hchain: {hchain} (blue) | Lchain: {lchain or 'N/A (VHH)'} (purple) | Antigen: {','.join(antigen_chains)} (green)</p>
{html}
</body>
</html>"""

    out_path = STRUCTURES_DIR / out_name
    with open(out_path, 'w') as f:
        f.write(full_html)
    logger.info(f"Wrote {out_path}")
    return out_path


def main():
    master = pd.read_parquet(PROJECT_ROOT / "data" / "processed" / "stage2_structure_qc_master.parquet")

    vhh_strict_gordon = master[
        (master["primary_vhh_candidate_strict"] == True) &
        (master["qc_stage2_gordon_complex_pass"] == True)
    ]
    vhh_rep = select_representative(vhh_strict_gordon)

    paired_gordon = master[
        (master["cohort"] == "PAIRED") &
        (master["qc_stage2_gordon_complex_pass"] == True)
    ]
    paired_rep = select_representative(paired_gordon)

    vhh_path = render_structure(
        vhh_rep, f"vhh_example_{vhh_rep['PDB']}_{vhh_rep['INSTANCE']}.html",
        "Representative Strict-Primary VHH-Antigen Complex (Stage-2 QC-Passing)"
    )
    paired_path = render_structure(
        paired_rep, f"paired_example_{paired_rep['PDB']}_{paired_rep['INSTANCE']}.html",
        "Representative Paired VH-VL-Antigen Complex (Stage-2 QC-Passing)"
    )

    selection_df = pd.DataFrame([
        {
            'category': 'VHH_representative',
            'instance': vhh_rep['INSTANCE'],
            'pdb': vhh_rep['PDB'],
            'selection_rule': 'lexicographically smallest INSTANCE among primary_vhh_candidate_strict + qc_stage2_gordon_complex_pass',
            'html_path': str(vhh_path.relative_to(PROJECT_ROOT)),
        },
        {
            'category': 'PAIRED_representative',
            'instance': paired_rep['INSTANCE'],
            'pdb': paired_rep['PDB'],
            'selection_rule': 'lexicographically smallest INSTANCE among PAIRED cohort + qc_stage2_gordon_complex_pass',
            'html_path': str(paired_path.relative_to(PROJECT_ROOT)),
        },
    ])
    selection_out = PROJECT_ROOT / "results" / "structures" / "representative_selection.csv"
    selection_df.to_csv(selection_out, index=False)
    logger.info(f"Wrote {selection_out}")

    # viewer_index.csv
    viewer_index = pd.DataFrame([
        {
            'instance': vhh_rep['INSTANCE'], 'pdb': vhh_rep['PDB'],
            'hchain': vhh_rep['Hchain'], 'lchain': 'NA',
            'antigen_chains': vhh_rep.get('antigen_chains', ''),
            'interactive_html_path': str(vhh_path.relative_to(PROJECT_ROOT)),
            'static_figure_path': '',
        },
        {
            'instance': paired_rep['INSTANCE'], 'pdb': paired_rep['PDB'],
            'hchain': paired_rep['Hchain'], 'lchain': paired_rep['Lchain'],
            'antigen_chains': paired_rep.get('antigen_chains', ''),
            'interactive_html_path': str(paired_path.relative_to(PROJECT_ROOT)),
            'static_figure_path': '',
        },
    ])
    viewer_index_out = PROJECT_ROOT / "results" / "structures" / "viewer_index.csv"
    viewer_index.to_csv(viewer_index_out, index=False)
    logger.info(f"Wrote {viewer_index_out}")

    return vhh_path, paired_path


if __name__ == "__main__":
    main()
