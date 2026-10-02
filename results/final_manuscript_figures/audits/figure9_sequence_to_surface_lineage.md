# Figure 9 — Sequence-to-Surface Lineage Audit

Every quantity below traces to an existing frozen file, verified by direct load in this session (not memory).

| metric | source file | exact definition | cohort | position subset | frozen/new | direct residue attribution? | limitations |
|---|---|---|---|---|---|---|---|
| A. Empirical former-interface | `data/processed/stage3_empirical_vh_vl_interface_primary.csv` | `primary_binary_interface==True` | both | 28 entries (25 unique integers + 111A/112A/112B) | frozen | n/a | see Figure-1 25-vs-28 audit; all 28 used here, none dropped |
| B. VHH AA frequency | `data/processed/stage4_position_aa_frequencies.csv` (already aggregated in `Figure2_panelA_VHH_sequence_logo.csv`) | count/N_observed per (position, AA), X dropped+renormalized | VHH N=533 | 28 | frozen (reused from Figure 2) | yes | none |
| C. VH AA frequency | same file, `Figure2_panelB_VH_sequence_logo.csv` | same | VH N=1429 | 28 | frozen (reused) | yes | none |
| D. JSD by position | `data/processed/stage4_position_chemistry_comparison.csv` (reused as `Figure2_panelD_JSD.csv`) | Jensen-Shannon divergence between complete VHH/VH AA distributions | both | 28 | frozen (reused) | position-level only, unsigned | none |
| E. Mean hydropathy by position | `stage4_position_chemistry_comparison.csv`: `mean_KD_VHH`, `mean_KD_VH`, `delta_mean_KD` | mean Kyte-Doolittle over observed residues at that position | both | 28 | frozen | yes | none |
| F. Sequence charge by position | **NEW**, derived directly from `stage4_interface_projection_long_all_instances.parquet` (same frozen per-residue table used for the Figure-1/2 charge work), mean of {D:-1,E:-1,R:+1,K:+1,H:+0.1} per position per group | both | 28 | new (allowed: "mean sequence-charge difference from existing residue assignments") | yes | none |
| G. Native exposure (SASA) | `data/processed/stage5_interface_surface_chemistry_residues.parquet`, mean `sasa` per position, `structural_state` VHH_NATIVE / VH_WITHOUT_VL_GEOMETRIC | both | 28 (127 domain positions available, restricted to 28) | frozen | yes | VH uses VL-omitted-from-calculation convention, not isolation |
| H. CDR3 shielding by position | `data/processed/stage6_position_shielding_summary.csv`: `mean_delta_sasa_shield` | VHH only (shielding is a VHH-specific geometric masking concept — CDR3 does not exist in the same sense to "shield" in a paired VH-VL structure where VL already covers the surface) | VHH only | 15 (SHIELDABLE_NON_CDR3_INTERFACE subset; 13 of the 28 CDR3-range positions have no shielding value by definition) | frozen | yes | VH has no equivalent value; NA for VH and for the 13 CDR3-range positions |
| I. PSH | `data/processed/stage8_surface_patch_descriptors.parquet`: `former_interface_core_PSH` | frozen Stage-8 patch score | both | aggregate (not position-resolved) | frozen | no (whole-structure scalar) | reused unmodified from Figure 5/earlier Figure 9 |
| J. PPC | same file: `former_interface_core_PPC` | frozen Stage-8 patch score, zero-inflated | both | aggregate | frozen | no | reused unmodified |
| K. PNC | same file: `former_interface_core_PNC` | frozen Stage-8 patch score, zero-inflated | both | aggregate | frozen | no | reused unmodified |
| L. PPC patch residue membership | **NOT AVAILABLE FROM VALIDATED OUTPUTS** | PPC is a Raybould-style pairwise-sum score (`stage8_patch_scoring.py::compute_patch_scores`, sum over residue PAIRS of \|Q1\|\|Q2\|/d²) — there is no connected-component "largest positive patch" concept analogous to PSH's `stage8_patch_topology.py` (which is explicitly a hydrophobic-only module, `HYDROPHOBIC_AA` residue set, no charge equivalent exists in the frozen codebase) | n/a | n/a | n/a | n/a | confirmed by reading `stage8_patch_topology.py` in full: no PPC/PNC connectivity function exists anywhere in the project |
| M. PNC patch residue membership | **NOT AVAILABLE FROM VALIDATED OUTPUTS** | same reasoning as L | n/a | n/a | n/a | n/a | n/a |
| N. APBS positive/negative components | `final_figures_model_input.parquet`: `interface_potential_positive_mean`, `interface_potential_negative_mean` | frozen Stage-8 APBS, per-structure means, kept separate (never netted) | VHH N=521, VH N=1386 (APBS-valid subset) | aggregate | frozen | no | reused unmodified |
| O. Local APBS residue-level values | `data/interim/stage8/apbs_cohort/{instance}/{instance}.pqr.dx` (raw APBS potential grid, frozen, per-structure) exists and is reused for 3D surface rendering only (grid-to-surface coloring, not a residue-additive decomposition) | both (where .dx exists for the chosen exemplar) | continuum grid, not per-residue | frozen | **no** — a 3D surface color is a spatial visualization of the existing grid, not a residue-attributed value; never described as a per-residue APBS contribution | grid-based rendering only |
| P. Former-interface antigen contact frequency | `data/processed/stage9_framework_antigen_contact_frequency.csv`: `VHH_antigen_contact_frequency`, `VH_antigen_contact_frequency` | frozen Stage-9 per-position contact frequency across the full cohort | both | subset of 28 (this table is framework-scoped; CDR3-range positions are absent, matching the Former-Interface-Atlas finding from the earlier session) | frozen | yes | NA for CDR3-range positions (107-116); used only as a small secondary annotation track, per instruction |
| Existing 3D structures/scripts | `data/interim/final_manuscript_structures/`, `data/interim/figure5_median_structures/`, `scripts/figure5_render_median_panels.py`, `results/final_manuscript_figures/main/Figure5_panel{F,G,H}*.png` | Figure-5 median-PSH exemplars (VHH `pdb_00006jb2-A`, VH `pdb_00004k2u-I-M`) reused directly per explicit instruction not to select new exemplars unnecessarily | both | n/a | frozen (reused) | n/a | reused for Layer C (3D position mapping) and the PSH module |

## 25-vs-28 handling

All 28 frozen interface entries are used throughout (confirmed via `figure1_empirical_interface_position_audit.csv`, all `available_in_stage4_sequence_table==True`) — no silent reduction to 25. Positions without a value for a given metric (e.g., CDR3-range positions for CDR3 shielding, which is defined only over the 15-position non-CDR3 subset) are shown as NA in the position-level bridge table, never zero-filled.

## Summary of what is genuinely new in this revision

1. Per-position mean sequence charge (VHH, VH, delta) — a direct mean of existing per-residue charge assignments, explicitly an allowed derived summary.
2. Per-position native exposure aggregation from the frozen Stage-5 residue table (a mean of existing SASA values, not a new SASA calculation).
3. PPC/PNC prevalence (proportion with patch>0) and conditional-magnitude (among patch>0) statistics — standard statistical summaries of already-existing frozen scalar values, not new biological metrics.
4. Deterministic exemplar selection for PPC-positive, PNC-positive, and APBS multivariate-median representative structures (selection only, no new calculation of the underlying values).
5. 3D surface rendering of per-residue `charge_class` (an existing Stage-8 per-residue classification) and of the frozen APBS potential grid — visualization only.

Everything else in this figure reuses frozen values without modification.
