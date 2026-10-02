# Stage 12 — Raybould et al. 2019 Figure 2 Visualization Source Audit

**Source:** Raybould MIJ, Marks C, Krawczyk K, Taddese B, Nowak J, Lewis AP, Bujotzek A, Shi J, Deane CM. "Five computational developability guidelines for therapeutic antibody profiling." *PNAS* 116(10):4025-4030 (2019). Full text retrieved directly from PubMed Central: PMC6410772 (https://pmc.ncbi.nlm.nih.gov/articles/PMC6410772/). PNAS's own site (pnas.org) returned HTTP 403 to automated fetch; PMC was used as the authoritative open-access mirror of the identical published text — not reproduced from memory.

## A. What Figure 2A plots (verbatim caption)

> "(A) CDR vicinity PSH scores across the 137 CST (blue) and human VdH Ig-seq (red) models. The CSTs are underrepresented at higher PSH values."

A histogram/distribution comparison of the **CDR-vicinity PSH score** (a single scalar per antibody model) across two cohorts (clinical-stage therapeutics, "CST," vs a human repertoire reference, "Ig-seq"), using distinct colors per cohort.

## B. What Figure 2B renders (verbatim caption)

> "(B) Galiximab (Kyte and Doolittle CDR vicinity PSH score of 167.89) has a large surface-exposed patch of hydrophobicity in its CDRH3 loop. Heavy- and light-chain surfaces outside the CDR vicinity are colored in white."

A single real antibody's molecular surface, physicochemically colored **only within the CDR vicinity**; everything outside that region is explicitly white.

## C. What region is colored

Only the **CDR vicinity** — defined in Methods (below) — is colored by hydrophobicity. Not the whole Fv.

## D. What remains white

Explicitly stated: **"Heavy- and light-chain surfaces outside the CDR vicinity are colored in white."** This is a direct quote, not an inference.

## E. Exposure definition (verbatim Methods)

> "Residues defined as 'surface-exposed' have ≥7.5% relative exposure across side-chain atoms, compared with the open-chain form alanine-R-alanine, as calculated with the Shrake and Rupley algorithm."

**Confirms Stage 8's already-adopted threshold exactly** (7.5% side-chain relative ASA) — independently verified here from the primary source text itself, not merely inferred from the TNP secondary-evidence chain used in Stage 8's original source gate.

## F. The hydrophobicity scale

**Kyte & Doolittle** — confirmed by the Figure 2B caption itself ("Kyte and Doolittle CDR vicinity PSH score") and by the PSH equation's use of "H(R,S)" for "scheme S."

## G. The numerical hydrophobicity normalization

> "H(R,S) is the normalized hydrophobicity score (**between 1 and 2**) for residue R in scheme S."

Confirms the exact **[1, 2]** normalized range already used in Stage 8's implementation (`KD_NORM`, min-max rescaled). No different normalization is used anywhere in the source.

## H. How PSH differs from simple residue coloring

PSH is **not** a per-residue quantity. The full equation (verbatim):

> "∑<sub>R1R2</sub> H(R1,S)·H(R2,S) / r12²" where "R1 and R2 are two surface-exposed residues with a closest heavy-atom distance, r12, <7.5 Å."

PSH is a **pairwise, spatially-weighted sum** over all exposed-residue pairs within 7.5 Å of each other — it aggregates spatial clustering of hydrophobicity, not merely how hydrophobic any single residue is. The molecular-surface coloring in Figure 2B therefore shows the **per-residue hydrophobicity component** that feeds into PSH, not the PSH score itself rendered onto individual residues. This distinction is preserved explicitly in this project's own figure (Section 15 of the task spec) — color = residue hydropathy; spatial aggregation of colored residues is what produces a high PSH score, and the two must never be conflated in the caption.

## I. Why galiximab was used in Figure 2B

Stated directly: to illustrate a **large, real, surface-exposed hydrophobic patch in its CDRH3 loop** (PSH=167.89, presumably a high value within their cohort, though the exact percentile is not stated in the retrieved text) — i.e., Raybould selected an antibody that visibly demonstrates what a high-PSH surface looks like, analogous in purpose (illustrating the phenotype) to this project's own predeclared 95th-percentile exemplar-selection rule (Section 6 of the task spec) — but the source text does not state that Raybould used a 95th-percentile selection rule specifically; this project's adaptation of *how* to deterministically pick an illustrative high-PSH structure is our own addition, not reproduced from the paper (correctly flagged as such per Section 6's explicit instruction).

## J. Quantitative vs illustrative elements

- **Quantitative (Panel A):** the full-cohort PSH histogram — a genuine distributional/statistical comparison across many antibodies.
- **Illustrative (Panel B):** a single real molecular structure, chosen to *demonstrate* the concept of a hydrophobic patch, not to represent a statistical comparison by itself. The paper does not claim Panel B alone proves anything about CST-vs-Ig-seq differences — that claim rests on Panel A.

## Additional details retrieved (not explicitly requested by the numbered items above, but relevant)

- **Colour scale/legend for the molecular surface:** **NOT STATED IN TEXT.** No RGB values, no colorbar anchors, no explicit "1.00/1.33/1.66/2.00" gradient description was found in the retrieved Methods or caption text. This project's four-anchor color scale (deep blue→bright green, per the task spec's Section 13) is therefore explicitly a **graphical reconstruction**, not a reproduction of stated RGB values — documented honestly as such in this project's figure caption, per the task instruction.
- **Rendering software:** **NOT STATED IN TEXT** (no PyMOL/Chimera/ChimeraX mention found). This project uses PyMOL (available locally, `vhh-interface-viz` environment, confirmed functional via headless smoke test) as a reasonable, commonly-used, publication-capable choice — not because the source paper specifies it.
- **"Anchor residue":** **NOT explicitly defined** as a standalone term in the retrieved text; used only in the CDR-vicinity definition ("every surface-exposed IMGT-defined CDR and anchor residue"), implying framework residues immediately flanking each CDR loop that are also treated as part of the "core" target set (consistent with this project's own Stage-8 `TNP_STANDARD_CORE` definition, which used the TNP secondary-evidence chain rather than this primary text).

## Summary: what is being adapted vs reproduced exactly

| Element | Raybould 2019 | This project's adaptation |
|---|---|---|
| Target physicochemical region | CDR vicinity | Frozen Stage-3 28-position FULL_EMPIRICAL_INTERFACE (never substituted) |
| Exposure threshold | ≥7.5% relative side-chain ASA (Shrake-Rupley) | Identical (Stage 8, independently confirmed here) |
| Hydrophobicity scale | Kyte-Doolittle, normalized [1,2] | Identical (Stage 8 `KD_NORM`) |
| PSH pairwise cutoff | <7.5 Å | Identical (Stage 8) |
| Vicinity radius | 4 Å | Identical (Stage 8, used for the *vicinity* target, distinct from the *core* 28-position set) |
| What's colored | CDR vicinity only; rest white | Former-interface only; rest white (direct structural analogue) |
| Exemplar selection | Not stated as a formal rule (galiximab shown as illustration) | Predeclared 95th-percentile-of-PSH rule (this project's own addition, explicitly not attributed to Raybould) |
| Color scale RGB values | Not stated in text | Graphical reconstruction, documented as such |
| Rendering software | Not stated in text | PyMOL (this project's choice) |

**This project never claims Raybould studied VHH, used a 95th-percentile selection rule, or specified exact RGB values — all such claims are explicitly attributed to this project's own adaptation.**
