# Figure 10 — Feature-by-Position Audit

## Selected union positions (14, IMGT order)

40, 42, 49, 50, 52, 55 (FR2) | 66, 69 (FR3) | 107, 112, 112B, 115, 116 (CDR3) | 118 (FR4)

Each position is the union of the top-5 (predeclared, non-manual `TOPN=5` ranking rule, identical to the rule already used in Figure 9) positions for at least one of the six included features. No position was added or removed by hand.

## Metric semantics (verified from source tables, not memory)

| Feature | Semantics | Sign convention (verified) | Source file | Source column |
|---|---|---|---|---|
| Sequence divergence | Unsigned | JSD ≥ 0 always | `stage4_position_chemistry_comparison.csv` | `JSD` |
| Hydropathy shift | Signed, VHH−VH | Verified `delta_mean_KD == mean_KD_VHH - mean_KD_VH` by direct recomputation (`np.allclose`, exact match) | `stage4_position_chemistry_comparison.csv` | `delta_mean_KD` |
| Charge shift | Signed, VHH−VH | Verified `delta_charge == VHH_mean_sequence_charge - VH_mean_sequence_charge` by direct recomputation | `Figure9_position_level_bridge.csv` (itself derived from `stage4_interface_projection_long_all_instances.parquet`) | `delta_charge` |
| Native exposure Δ | Signed, VHH−VH | Computed as `VHH_native_exposure - VH_native_exposure`, both already-existing per-group means from Stage 5 | `Figure9_position_level_bridge.csv` | `VHH_native_exposure`, `VH_native_exposure` |
| CDR3 shielding | Unsigned, **VHH-only** | No conventional-VH equivalent exists by definition (shielding is a VHH-specific geometric masking concept) | `stage6_position_shielding_summary.csv` | `mean_delta_sasa_shield` |
| Antigen engagement | Signed, VHH−VH | **Existing frozen column**, not derived in this session — verified `difference_VHHminusVH == VHH_antigen_contact_frequency - VH_antigen_contact_frequency` by direct recomputation | `stage9_framework_antigen_contact_frequency.csv` | `difference_VHHminusVH` |

## Provenance of each row's underlying values

- Sequence divergence, hydropathy shift: existing frozen Stage-4 columns, reused verbatim (same as Figure 1/2/9).
- Charge shift: reused directly from `Figure9_position_level_bridge.csv::delta_charge`, itself an already-audited transparent derivation (mean of an existing per-residue charge-score lookup, VHH minus VH) — not recomputed in this session.
- Native exposure Δ: a simple, transparent difference of two already-existing per-group means (same derivation class explicitly allowed in the Figure 9 no-invention protocol).
- CDR3 shielding: existing frozen Stage-6 column, reused verbatim, shown only for VHH and visually distinguished (separate sequential violet color family, explicit "(VHH only)" row label) so it can never be mistaken for a signed VHH-vs-VH comparison.
- Antigen engagement: existing frozen Stage-9 column (`difference_VHHminusVH`) — this is not something newly derived in this session; the frozen table already contains the signed difference.

## No unsupported residue-level metrics added

Per the explicit constraint, no residue-level PPC, PNC, PSH, or APBS attribution was added anywhere in Figure 10 — these quantities do not exist at position-level resolution in any validated output (see the Figure 9 audit, which already established this), and Figure 10 does not attempt to construct them.

## NA handling

CDR3 shielding and antigen engagement are undefined at CDR3-range positions (107, 112, 112B, 115, 116) by construction (shielding cannot apply to CDR3 itself; the antigen-contact-frequency table is framework-scoped). These cells are rendered as explicit gray hatched "N/A" cells, never zero-filled or silently omitted from the matrix.

## False-adjacency prevention

Positions are laid out with proportional integer spacing within each region block and an explicit visual gap (0.9 cell-widths) between region blocks, with a labeled gray region band above the matrix (FR2 / FR3 / CDR3 / FR4). Position 52 (FR2) and 66 (FR3) are visually separated by a clear gap and a region-band boundary rather than appearing adjacent, and likewise 69 (FR3) → 107 (CDR3) → 118 (FR4).

## Tetrad annotation

IMGT 42/49/50/52 are outlined with a thin magenta rectangle in every row (secondary annotation only — the cell fill color always encodes that row's own feature value, never overridden by the tetrad marker) and bolded/gold on the x-axis position labels.

## Verdict

All 12 mandatory visual-QA checks pass (see Figure10.png). No scientific values were computed beyond simple, already-precedented transparent differences of existing per-group values; all raw per-group inputs and the one existing signed-difference column are frozen, unmodified outputs.
