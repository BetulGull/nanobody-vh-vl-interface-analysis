# Figure 7 — Visual Engagement Audit

## 1-2. 28-position former-interface and contact definition (confirmed identical to Panel A)

`full_interface_engagement_fraction` (Panel A's frozen metric, `data/interim/stage9/{vhh,vh}_full_interface_engagement.csv`, denominator = 25 unique-integer `FULL_EMPIRICAL_INTERFACE` positions) was cross-validated against `stage9_antibody_antigen_contacts_long.parquet` restricted to `antibody_domain in {VHH, VH}` (VL excluded) at `contact_at_4p5A` — this threshold and domain filter exactly reproduces both frozen structure-specific values: `pdb_00008k4q-B` = 10/25 = 0.40, `pdb_00006wgl-A-B` = 7/25 = 0.28. Confirms the Panel B/C residue coloring uses the identical numerator/denominator logic as Panel A, not a re-derived approximation.

## 3-7. Mutually exclusive classification (verified programmatically)

For each structure: `F` = the 25-position `FULL_EMPIRICAL_INTERFACE` set; `A` = all antibody (VHH/VH domain only) positions with `contact_at_4p5A==True` (any IMGT position, not restricted to F). Rendered:
- `F ∩ A` → red (`former_interface_engaged`)
- `F \ A` → yellow (`former_interface_nonengaged`)
- `A \ F` → orange (`other_antigen_contact`)

Verified in `Figure7_structure_B/C_residue_classes.csv`: every position appears in exactly one `render_class`, `len(F∩A)` matches the frozen numerator exactly (10 for B, 7 for C) — no former-interface residue is colored orange, no non-former-interface residue is colored red.

## 9-10. Panel B/C engagement fractions match cohort source data

Panel B: 10/25 = 0.40 (matches `full_interface_engagement_fraction` for `pdb_00008k4q-B` in `final_figures_model_input.parquet` exactly). Panel C: 7/25 = 0.28 (matches `pdb_00006wgl-A-B` exactly).

## 11-14. Colors (verified against `Figure7_color_mapping.csv`)

VHH box fill `#7B9BE0` (blue family, same hue as VHH structural blue `#3060CC`); conventional-VH box fill `#B48FBF` (lavender, blended from VH purple `#7D4FA8` + VL pink `#F2A6C4`); antigen `#C8C8C8` (light gray, changed from the prior green); boxplot median/whiskers/outline `#2B2B2B` (dark charcoal, changed from the prior orange — orange is now reserved exclusively for `other_antigen_contact`).

## 15-16. Consistent semantics; no VL misclassification

Both panels use the identical color-to-class mapping. VL (`pdb_00006wgl` chain L) is rendered only as pink structural context (cartoon, translucent) — it is never included in the `F` (former-interface) or classification logic, which is built entirely from `antibody_domain=='VH'` (heavy-chain) rows.

## 17. Matched-epitope labeling

Rendered (dark gray, `#707070`) using the already-validated residue set from `Figure7_matched_pair_summary.csv` / the prior matched-epitope validation gate (`VALID_MATCHED_EPITOPE`, Jaccard=0.92, structural superposition RMSD=2.34 Å corroborated) — not a new or invented mapping. The antigen-based structural alignment used to render Panels B/C in a comparable frame reuses the exact same Kabsch transform (RMSD=2.3426 Å, N=176 common antigen Cα, re-verified identical to the frozen validation record) rather than re-deriving a new superposition.

## 18. No underlying frozen result changed

`full_interface_engagement_fraction`, the matched-epitope validation verdict, and the antigen sequence-identity/Jaccard values are all read-only in this revision — confirmed unchanged (`Figure7_cohort_engagement.csv` values identical to `final_figures_model_input.parquet`; `Figure7_matched_pair_summary.csv` identical to the frozen `figure7_matched_epitope_validation.csv`).

## Verdict

All 18 required checks pass.
