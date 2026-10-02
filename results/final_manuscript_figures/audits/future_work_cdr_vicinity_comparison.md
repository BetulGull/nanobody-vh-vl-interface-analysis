# Future Work (recorded only — not performed, not part of Figure 5)

## Question

Is VHH physicochemical remodeling restricted primarily to the surface exposed by VL loss (the former VH-VL interface, Figure 5's subject), or does the entire antigen-binding surface also exhibit a distinct physicochemical organization relative to conventional VH-VL?

## Proposed comparison

VHH CDR vicinity (Raybould-style: every surface-exposed CDR + anchor residue, plus other exposed residues within 4 A) vs the COMPLETE conventional VH-VL CDR vicinity (both chains' CDRs, since the conventional antigen-binding surface is built from 6 CDR loops, not 3).

## Three possible outcomes and their interpretation

**Scenario 1** — former-interface PSH differs (as shown in Figure 5), but CDR-vicinity PSH is similar between groups → the physicochemical remodeling is preferentially localized to the VL-loss surface, not a general shift over the whole antigen-binding surface.

**Scenario 2** — both former-interface and CDR-vicinity PSH are reduced in VHH → VHH exhibits broader surface physicochemical remodeling extending beyond the former VH-VL interface.

**Scenario 3** — former-interface PSH is reduced in VHH but CDR-vicinity PSH is equal or higher → VHH may spatially redistribute hydrophobic character away from the VL-loss surface while retaining or concentrating it in the antigen-binding region.

## Required methodological care if this is performed later

1. Use the exact Raybould-style CDR-vicinity definition (already source-verified in this project, see `results/audits/stage12_raybould_fig2_visualization_source_audit.md`), not a substitute.
2. Report the exact raw Raybould metric for reproducibility.
3. VHH has 3 CDR loops, conventional VH-VL has 6 — region size/exposed-residue-count/exposed-surface-area differ structurally between groups, so appropriate sensitivity analyses (e.g. per-residue or per-area normalization) must accompany the raw comparison, not replace it.
4. If bound (antigen-complex) structures are used, the antigen should be removed geometrically while retaining native antibody coordinates, to isolate the antibody's own intrinsic exposed CDR-vicinity surface.

## Status

Not calculated. Not added as a Figure 5 panel or any other figure. Recorded here only, per explicit instruction, for a future task.
