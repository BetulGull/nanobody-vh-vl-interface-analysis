"""
Final manuscript figures -- PyMOL structural panel renders for Figures 1, 4, 6, 7.
Run with: micromamba run -n vhh-interface-viz pymol -cq scripts/final_manuscript_render_structures.py

Colors follow the project-wide convention: VHH=blue, conventional VH=purple,
VL=pink, antigen=restrained green, CDR3=teal, former-interface(28-pos)=orange,
shieldable-15=yellow, tetrad=magenta spheres, everything else=light grey.
White background, orthographic projection, ray-traced, consistent per-panel
camera captured once and reused within each comparison pair (Section 6.C).
"""
from pymol import cmd
from pathlib import Path

ROOT = Path("/Users/betulgul/Desktop/staj")
STRUCT_DIR = ROOT / "data/interim/final_manuscript_structures"
OUT_DIR = ROOT / "results/final_manuscript_figures/main"
OUT_DIR.mkdir(parents=True, exist_ok=True)

FULL_28 = [40, 42, 44, 49, 50, 52, 55, 66, 67, 68, 69, 103, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 118, 119, 120]
CDR3_RANGE = list(range(105, 118))
SHIELDABLE_15 = [p for p in FULL_28 if p not in CDR3_RANGE]
TETRAD = [42, 49, 50, 52]

COL = {
    'vhh_blue': [0.20, 0.40, 0.80], 'vh_purple': [0.50, 0.30, 0.70], 'vl_pink': [0.95, 0.60, 0.75],
    'antigen_green': [0.55, 0.75, 0.55], 'cdr3_teal': [0.0, 0.60, 0.60],
    'interface_orange': [1.0, 0.60, 0.0], 'shieldable_yellow': [1.0, 0.85, 0.2],
    'tetrad_magenta': [0.85, 0.1, 0.65], 'grey': [0.82, 0.82, 0.82],
}
for name, rgb in COL.items():
    cmd.set_color(name, rgb)


def base_settings():
    cmd.bg_color('white')
    cmd.set('ray_opaque_background', 1)
    cmd.set('orthoscopic', 1)
    cmd.set('ray_shadows', 0)
    cmd.set('ambient', 0.4)
    cmd.set('specular', 0.15)
    cmd.set('cartoon_fancy_helices', 1)
    cmd.set('antialias', 2)


def ray_png(path, width=1400, height=1400, dpi=600):
    cmd.ray(width, height)
    cmd.png(str(path), dpi=dpi)


# ============ FIGURE 1: native Fv anchor ============
def figure1():
    cmd.delete('all')
    base_settings()
    cmd.load(str(STRUCT_DIR / 'figure1_native_fv.pdb'), 'fv')
    cmd.hide('everything')
    cmd.show('cartoon', 'fv')
    cmd.color('vh_purple', 'chain H')
    cmd.color('vl_pink', 'chain L')
    cmd.set('cartoon_transparency', 0.0)
    cmd.orient('fv')
    cmd.zoom('fv', buffer=3)
    cmd.turn('y', 20)
    view = cmd.get_view()

    # Panel A: CDRs annotated on VH+VL
    cmd.color('cdr3_teal', f'chain H and resi {"+".join(map(str, CDR3_RANGE))}')
    cmd.set_view(view)
    ray_png(OUT_DIR / 'Figure1_panelA_native_fv.png')

    # Panel B: 28-position empirical interface highlighted on VH, VL semi-transparent
    cmd.color('vh_purple', 'chain H')
    cmd.color('interface_orange', f'chain H and resi {"+".join(map(str, FULL_28))}')
    cmd.color('shieldable_yellow', f'chain H and resi {"+".join(map(str, SHIELDABLE_15))}')
    cmd.show('spheres', f'chain H and resi {"+".join(map(str, TETRAD))} and name CA')
    cmd.color('tetrad_magenta', f'chain H and resi {"+".join(map(str, TETRAD))} and name CA')
    cmd.set('sphere_scale', 0.6)
    cmd.set('cartoon_transparency', 0.65, 'chain L')
    cmd.set_view(view)
    ray_png(OUT_DIR / 'Figure1_panelB_former_interface.png')
    print('Figure 1 done')


# ============ FIGURE 4: FR2-CDR3 distance low/high exemplars ============
def figure4():
    cmd.delete('all')
    base_settings()
    cmd.load(str(STRUCT_DIR / 'fig4_low.pdb'), 'low')
    cmd.load(str(STRUCT_DIR / 'fig4_high_aligned.pdb'), 'high')
    for obj in ['low', 'high']:
        cmd.hide('everything', obj)
        cmd.show('cartoon', obj)
        cmd.color('vhh_blue', obj)
        cmd.color('cdr3_teal', f'{obj} and resi {"+".join(map(str, CDR3_RANGE))}')
        cmd.show('spheres', f'{obj} and resi 42 and name CA')
        cmd.color('tetrad_magenta', f'{obj} and resi 42 and name CA')
        cmd.set('sphere_scale', 0.7, obj)

    cmd.orient('low')
    cmd.zoom('low', buffer=3)
    view = cmd.get_view()

    cmd.disable('high')
    cmd.enable('low')
    cmd.set_view(view)
    ray_png(OUT_DIR / 'Figure4_panelD_low_distance.png')

    cmd.disable('low')
    cmd.enable('high')
    cmd.set_view(view)
    ray_png(OUT_DIR / 'Figure4_panelE_high_distance.png')
    print('Figure 4 done')


# ============ FIGURE 6: CDR3-length exemplars (VHH short vs long) ============
def figure6():
    cmd.delete('all')
    base_settings()
    cmd.load(str(STRUCT_DIR / 'fig6_vhh_short.pdb'), 'short')
    cmd.load(str(STRUCT_DIR / 'fig6_vhh_long_aligned.pdb'), 'long')
    for obj in ['short', 'long']:
        cmd.hide('everything', obj)
        cmd.show('cartoon', obj)
        cmd.color('vhh_blue', obj)
        cmd.color('cdr3_teal', f'{obj} and resi {"+".join(map(str, CDR3_RANGE))}')
        cmd.color('interface_orange', f'{obj} and resi {"+".join(map(str, FULL_28))} and not resi {"+".join(map(str, CDR3_RANGE))}')

    cmd.orient('short')
    cmd.zoom('short', buffer=3)
    view = cmd.get_view()

    cmd.disable('long')
    cmd.enable('short')
    cmd.set_view(view)
    ray_png(OUT_DIR / 'Figure6_panel_vhh_short_cdr3.png')

    cmd.disable('short')
    cmd.enable('long')
    cmd.set_view(view)
    ray_png(OUT_DIR / 'Figure6_panel_vhh_long_cdr3.png')
    print('Figure 6 done')


# ============ FIGURE 7: antigen-engagement complexes, matched epitope ============
def figure7():
    cmd.delete('all')
    base_settings()
    cmd.load(str(STRUCT_DIR / 'fig7_vhh_complex.pdb'), 'vhh_cplx')
    cmd.load(str(STRUCT_DIR / 'fig7_vh_complex_aligned.pdb'), 'vh_cplx')
    matched_epitope_resi = [41, 42, 43, 44, 45, 47, 66, 67, 69, 71, 72]  # 8k4q-frame numbering, validated matched set
    epitope_sel = '+'.join(str(r) for r in matched_epitope_resi)

    for obj, ab_chains in [('vhh_cplx', ['A']), ('vh_cplx', ['H', 'L'])]:
        cmd.hide('everything', obj)
        cmd.show('cartoon', f'{obj} and chain {"+".join(ab_chains)}')
        cmd.show('surface', f'{obj} and chain G')
        cmd.color('antigen_green', f'{obj} and chain G')
        cmd.color('interface_orange', f'{obj} and chain G and resi {epitope_sel}')
        if obj == 'vhh_cplx':
            cmd.color('vhh_blue', f'{obj} and chain A')
        else:
            cmd.color('vh_purple', f'{obj} and chain H')
            cmd.color('vl_pink', f'{obj} and chain L')
        cmd.set('transparency', 0.1, f'{obj} and chain G')

    cmd.orient('vhh_cplx')
    cmd.zoom('vhh_cplx', buffer=5)
    view = cmd.get_view()

    cmd.disable('vh_cplx')
    cmd.enable('vhh_cplx')
    cmd.set_view(view)
    ray_png(OUT_DIR / 'Figure7_panelG_vhh_complex.png')

    cmd.disable('vhh_cplx')
    cmd.enable('vh_cplx')
    cmd.set_view(view)
    ray_png(OUT_DIR / 'Figure7_panelH_vh_complex.png')
    print('Figure 7 done')


figure1()
figure4()
figure6()
figure7()
print('ALL STRUCTURAL PANELS RENDERED')
