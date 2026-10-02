# Zhang et al. 2025 — FR2_CDR3_Dis Method Definition

**Source paper:** Zhang L, Yang T, Zhang Y, Yan J, Li J, Xu W, Zhu W, Wang X (2025). "Comparative analysis of CDR3 length-dependent patterns in VHHs." *Front Immunol* 16:1647230. `references/Zhang_2025.pdf` (confirmed present).

**Supplement status:** `references/Zhang_2025_supplement.pdf` — not present, **not needed**: the exact computational definition is given in full in the main paper's Materials and Methods section (Section 2, p.2), quoted verbatim below.

**Author code status:** `external/zhang_2025_code/` — not present, not needed for the same reason.

**Priority-of-authority resolution:** Main paper Methods text is sufficient and unambiguous; no discrepancy to resolve.

---

## Exact quoted definition (Materials and Methods, p.2)

> "CDRs were defined according to the IMGT numbering scheme... To calculate the minimum distance between position 42 and CDR3 residues (denoted as FR2_CDR3_Dis) for CDR3 conformation analysis, **the central coordinates of the residue at position 42 and each residue in CDR3 (excluding the first and last two residues) were obtained. The distance between two central coordinates was calculated as the distance between two residues.**"

**Corroborating description (Results, p.6, Figure 4 legend):**
> "FR2_CDR3_Dis: distance between CDR3 and the FR2 region, calculated as the minimum distance between the residue at IMGT 42 and any residue within CDR3 after excluding first and last two residues."

Both passages agree exactly: FR2_CDR3_Dis is the **minimum**, over all "internal" CDR3 residues, of the distance between the **centroid** ("central coordinate") of IMGT position 42 and the **centroid** of that internal CDR3 residue.

---

## Extracted definition components

| Component | Definition used |
|---|---|
| **Position 42 representation** | IMGT position 42 (FR2), single residue, exact IMGT identity (project's IMGT auth_seq_id/insertion-code source, per Stage 2-5 validated convention). |
| **CDR3 residues included** | IMGT 105–117 (project's standard CDR3 definition), **excluding the first two and last two residues in true sequence order** (i.e. real-residue order via `label_seq_id`, not merely the base IMGT integers 105/106/116/117 — see note below). |
| **"Central coordinate" of a residue** | Literal reading: centroid (arithmetic mean) of **all non-hydrogen atom coordinates** of that residue. The paper does not qualify this as CA-only, backbone-only, or side-chain-only; "central coordinate" is the natural term for a whole-residue geometric center, and this is the interpretation used. Documented explicitly as an interpretation, not an invented method, since the paper gives no further qualification. |
| **Distance operation** | Euclidean distance between the two centroids (position-42 centroid, internal-CDR3-residue centroid). |
| **Aggregation rule** | **Minimum** over all internal CDR3 residues (i.e., over the full internal set, take the smallest centroid-to-centroid distance). |
| **First/last two exclusion** | Confirmed exactly as specified in the original Stage-6 prompt's prior-project-review summary — this is the ONE piece of prior summary that is directly and exactly corroborated by the primary source text, so it is retained with confidence. |

**Note on "first and last two residues":** IMGT numbering places any CDR3-length insertions in the *middle* of the loop (around IMGT 111/112), never at the true sequence ends. Therefore the true first-two and last-two residues in real chain order are reliably IMGT 105, 106 (first two) and 116, 117 (last two) regardless of loop length/insertion depth — confirmed by inspecting the project's own residue-scheme ordering (`label_seq_id` order) rather than assumed. The "internal" set used for FR2_CDR3_Dis is therefore: every real, coordinate-observed CDR3 residue with true sequence order strictly between the 106-residue and the 116-residue (this includes all insertion-coded middle residues, e.g. 107–115 plus any 111X/112X insertions).

**Edge case, documented not silently handled:** a CDR3 with fewer than 5 real residues has an empty "internal" set after excluding 2+2. `FR2_CDR3_Dis` is reported as `NaN` with `geometry_missing_reason = "CDR3_TOO_SHORT_FOR_INTERNAL_SET"` for such structures (none silently dropped from the output table).

---

## Summary of what this project computes under the "Zhang" name

| Descriptor | Status | Formula used |
|---|---|---|
| `FR2_CDR3_Dis` | **Fully implemented, unambiguous** | `min over internal CDR3 residues r of ‖centroid(residue 42) − centroid(r)‖`, where "internal" excludes the true first two and last two real CDR3 residues in sequence order. |

No ambiguity remains for this metric; the full computational protocol was recoverable directly from the main paper's Methods section without requiring the supplement or author code.
