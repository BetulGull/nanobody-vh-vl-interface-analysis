# Gordon et al. 2023 — CDR-H3 Geometry Method Definition

**Source paper:** Gordon GL, Capel HL, Guluglu B, Richardson E, Stafford RL, Deane CM (2023). "A comparison of the binding sites of antibodies and single-domain antibodies." *Front Immunol* 14:1231623. `references/Gordon_2023.pdf` (confirmed present).

**Supplement status:** `references/Gordon_2023_supplement.pdf` — **NOT PRESENT** in this repository despite being an expected location. The paper's own Methods section (2.9, p.4) explicitly states the exact geometric construction is given in "Text S2" (a supplementary document). This document was not supplied.

**Author code status:** `external/gordon_2023_code/` — **NOT PRESENT**. However, `external/TNP/` (the "Therapeutic Nanobody Profiler," Gordon GL, Deane CM et al., same OPIG lab, a related but textually distinct later work) contains `lib/python/theraprofnano/CDR_Profiler/CDR3_Conf_Assigner.py`, which implements a CDR3 "reach"/compactness calculation. This is **not confirmed to be the official Gordon 2023 code release** (no explicit citation link found), but is strong corroborating evidence from the same authors implementing the same qualitative concept — used here as priority-2 evidence (supplementary-equivalent), never as a silent substitute for the actual Text S2.

---

## Priority-of-authority note (per instruction)

1. Author source code (if clearly corresponding to the paper) — **partially available** (TNP code, not confirmed as the exact Gordon 2023 release, covers ρ only).
2. Supplementary methods — **not available**.
3. Main paper methods — **available**, qualitative description only, defers exact formula to Text S2.

---

## 1. ρ (rho) — RESOLVED, implementable

**Main paper text (p.4, Section 2.9):**
> "Using the spherical coordinates method, ρ describes the reach of the CDR-H3 loop away from the rest of the VH domain. A CDR-H3 loop in an extended conformation will have a high ρ value whereas a loop of identical length that is folded against the VH domain will have a lower ρ value."

**Corroborating author code** (`CDR3_Conf_Assigner.py`, function `get_H3_anchor_line` + `main_compactness`):
```python
# Anchor set for CDR-H3 (IMGT numbering), per code:
anchors = {3: [102, 103, 118, 119]}          # 2 residues flanking each side
loop_boundaries = {3: list(range(105, 117+1))}  # CDR3 = IMGT 105-117

# c = centroid of ALL backbone atoms (N, CA, C) of the 4 anchor residues
c = mean(anchor_backbone_coords)   # anchor_backbone_coords includes N,CA,C for each of 102,103,118,119

# h = centroid of ALL backbone atoms (N, CA, C) of every CDR3 residue (105-117)
h = mean(cdr3_backbone_coords)

rho = ||h - c||   # Euclidean distance between the two centroids
```

**Discrepancy documented, not silently resolved:** the main paper's QC/exclusion criterion (Section 2.9) removes structures with missing backbone atoms in CDR-H loops "or anchors (**three** residues on either side of each loop)" — a 3-residue-each-side window used for *data curation*. The actual geometry-construction anchor set in the corroborating code uses only **2** residues each side (IMGT 102,103 and 118,119). These are not necessarily the same set for the same purpose: the paper's "3 each side" appears to be the QC/anchor-completeness filter, while the code's 4-point anchor set (2+2) is what is actually used to build the reference centroid for ρ. **This implementation uses the code's exact 2+2 anchor set for the geometric formula** (priority: code > paper text, per instruction), and separately applies a 3-residue-each-side backbone-completeness QC check (matching the paper's stated exclusion criterion) before trusting a ρ value.

**Implementation used in this project:**
- Anchor centroid `c` = mean of backbone atoms (N, CA, C) over IMGT 102, 103, 118, 119.
- CDR3 centroid `h` = mean of backbone atoms (N, CA, C) over all real, coordinate-observed CDR3 residues (IMGT 105–117, including insertion-coded positions).
- `rho = ||h - c||` (Å).
- QC gate: structure excluded from ρ if any of the 3-residue-each-side anchor window (IMGT 102–104, 118–120) or any CDR3 residue lacks complete backbone.

## 2. Compactness — RESOLVED, implementable

**Main paper text (p.4):**
> "Lastly, we divide the loop length by ρ to give a measure of compaction. A loop with low compactness uses its entire length to reach away from the VH domain, whereas high compactness corresponds to a loop that is packed against the VH."

**Corroborating code comment** (`main_compactness` docstring): *"return loop reach, used in calculation in main script to get compactness (cdr3 length / rho)"*.

**Formula:** `compactness = cdr3_length / rho`

**"Loop length" interpretation:** the paper does not specify sequence-residue-count vs. summed inter-residue chain length in Å for "loop length" at this specific point in the text; the code comment "cdr3 length" alongside the sequence-residue-count context of the paper's Figure 1 (CDR-H3 loop length measured in residues throughout the whole paper) supports **residue count** as the intended "length" unit here. Used: `cdr3_sequence_length` (real-residue count, from `_entity_poly_seq`, independent of resolution — see Stage 6 CDR3-length section) as the numerator.

## 3. θ (theta) — **UNRESOLVED, NOT IMPLEMENTED**

**Main paper text (p.4):**
> "θ gives a measure of the elevation of the loop. A low value corresponds to a CDR-H3 that extends directly up and away from the rest of the VH domain, whereas a high value indicates that the loop is 'folding' down. In the case of Ab structures, a high θ value corresponds to a loop that is packed into the groove created by the VH-VL interface."

This is a **qualitative** description only. The exact coordinate-frame construction, the reference axis/plane defining "elevation," and the precise trigonometric formula are stated to be in Text S2, which is not available. No corroborating code for θ was found anywhere in `external/TNP`. **This project does NOT invent a coordinate frame or formula for θ.** `cdr3_theta` is reported as `NaN` / `geometry_missing_reason = "GORDON_THETA_UNDEFINED_NO_SUPPLEMENT"` for every structure.

## 4. φ (phi) — **UNRESOLVED, NOT IMPLEMENTED**

**Main paper text (p.4):**
> "φ gives an indication of whether the CDR-H3 loop is horizontally oriented towards the rest of the VH domain or away from it. In the case of Ab structures, a high φ value indicates packing against the VL domain."

Same situation as θ: qualitative only, exact reference-vector/plane construction deferred to Text S2, no corroborating code found. **Not implemented.** `cdr3_phi` reported as `NaN` / `geometry_missing_reason = "GORDON_PHI_UNDEFINED_NO_SUPPLEMENT"`.

## 5. R3 coordinate system

Referenced in the main text ("examining their centers of geometry in reference to an R3 coordinate system (see Text S2 for method...)") as the overall framework within which ρ, θ, φ are all defined together. Since θ and φ construction is unavailable, the full R3 frame is not reconstructed; only the isolated centroid-distance operation needed for ρ (which does not require a full oriented coordinate frame — it is a simple Euclidean distance between two centroids) is implemented.

---

## Summary of what this project computes under the "Gordon" name

| Descriptor | Status | Formula used |
|---|---|---|
| `cdr3_rho` | **Implemented** | `‖centroid(CDR3 backbone) − centroid(anchor backbone: IMGT 102,103,118,119)‖` |
| `cdr3_compactness` | **Implemented** | `cdr3_sequence_length / cdr3_rho` |
| `cdr3_theta` | **Not implemented** | N/A — Text S2 unavailable, no corroborating code |
| `cdr3_phi` | **Not implemented** | N/A — Text S2 unavailable, no corroborating code |

No plausible substitute was fabricated for θ or φ under Gordon's name, per instruction.

---

## θ/φ Structured Dependency Record (updated after literature dependency re-check)

Re-checked per `results/audits/stage6_literature_dependency_audit.md` — no new local source material found. Structured per the required template; every field that cannot be sourced is marked UNRESOLVED rather than filled with an assumption.

**SOURCE:** Gordon et al. 2023, main paper only (`references/Gordon_2023.pdf`). Text S2 (the paper's own stated source of the exact method) is hosted externally at Frontiers' website and is **not present in this repository**.

**EXACT LOCATION:** Main paper, Section 2.9 "Orientation of CDR-H3 loops," p.4, paragraph beginning "We analyzed the general orientation..." The paper explicitly cites "(see Text S2 for method and Supplementary Figure S4)" for the actual construction — i.e., the main paper itself declares the location of the authoritative definition to be a document this project does not have.

**ORIGIN:** UNRESOLVED. The paper states θ/φ/ρ are computed "in reference to an R3 coordinate system" and "by examining their centers of geometry," but does not state, in the main text, what point is the coordinate-system origin (e.g., an anchor centroid, a fixed VH-domain reference point, or something else).

**REFERENCE FRAME:** UNRESOLVED. No axis construction, no stated reference vector/plane, no handedness given in the main text.

**ATOM / RESIDUE SELECTION:** Partially inferable by analogy to ρ (anchor residues IMGT 102,103,118,119 and CDR3 105–117, per the corroborating TNP code), but this is an *inference from ρ's construction*, not a confirmed statement that θ/φ use the identical atom selection — the paper describes all three as computed "using the spherical coordinates method" together, which is suggestive but not an explicit confirmation for θ/φ specifically. Not treated as sufficient to implement.

**THETA DEFINITION:** UNRESOLVED. Main text (qualitative only): "θ gives a measure of the elevation of the loop. A low value corresponds to a CDR-H3 that extends directly up and away from the rest of the VH domain, whereas a high value indicates that the loop is 'folding' down." No formula, no reference plane for "elevation," no equation.

**PHI DEFINITION:** UNRESOLVED. Main text (qualitative only): "φ gives an indication of whether the CDR-H3 loop is horizontally oriented towards the rest of the VH domain or away from it." No formula, no reference vector for "horizontal," no equation.

**ANGLE UNITS:** UNRESOLVED (not stated in main text; would be in Text S2).

**ANGLE RANGE:** UNRESOLVED (not stated; a "spherical coordinates" framing suggests conventional physics/math conventions such as θ∈[0,π] and φ∈[0,2π) or φ∈(−π,π], but this is a plausible guess from the general term "spherical coordinates," not a confirmed convention — explicitly not used as the basis for implementation).

**HANDEDNESS / SIGN CONVENTION:** UNRESOLVED.

**WRAP-AROUND CONVENTION:** UNRESOLVED (relevant only if φ is confirmed periodic, which itself is not confirmed here).

**IMPLEMENTATION INTERPRETATION:** None attempted. Every field above marked UNRESOLVED is a required component per the task's own checklist (§5: "Do NOT begin implementation if any required component remains unresolved"), and the majority of the required components are unresolved. Implementation is correctly withheld.

**RESIDUAL AMBIGUITIES:** Coordinate origin, reference frame/axes, exact θ formula, exact φ formula, angle units, angle range, handedness, and wrap-around convention are all unresolved. Only the atom/residue selection is even partially inferable, and only by analogy, not confirmation.

**GORDON_THETA_PHI_STATUS = UNRESOLVED_METHOD_DEPENDENCY**
