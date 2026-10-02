# Stage 6 — Gordon θ/φ Literature Dependency Audit

**Generated:** 2026-09-02
**Purpose:** Determine whether any genuinely new source material in this repository resolves the Gordon et al. 2023 θ/φ geometric definitions, which the prior Stage-6 audit found unresolvable from the previously available sources.

---

## Files/locations searched

| Location | Method | Result |
|---|---|---|
| `references/` | full directory listing | Only `Gordon_2023.pdf`, `Zhang_2025.pdf` present. No `Gordon_2023_supplement.pdf`, no `Zhang_2025_supplement.pdf`. Unchanged since prior audit. |
| `docs/` | full directory listing | `STAGE2_PROGRESS.md`, `visualization_standard.md`, `STAGE2_IMPLEMENTATION_PLAN.md` — none relevant to Gordon geometry. |
| Project root | filename search (`*gordon*`, `*supplement*`) | No matches. |
| `external/` | full subdirectory listing | Only `external/TNP` exists (already known, already inspected). No new external code directories. |
| `Gordon_2023.pdf` itself | full-text extraction (`pdftotext -layout`), all 18 pages | Confirmed 18 pages, entirely main-text. No embedded supplement pages. Explicit text (p.15, "Supplementary material" section): *"The Supplementary Material for this article can be found online at: https://www.frontiersin.org/articles/10.3389/fimmu.2023.1231623/full#supplementary-material"* — i.e., the supplement is hosted externally at the publisher's website, not bundled in this PDF, and not present anywhere in this repository. |
| `Gordon_2023.pdf` — code pointer | full-text search | *"The code generated and datasets analyzed for this study can be found at github.com/oxpig."* — a pointer to the OPIG GitHub **organization** generally, not a specific repository or file. No such repository is cloned or present anywhere in this project (only `external/TNP`, a different, later, related-but-distinct OPIG publication, is present). |
| `external/TNP` | full-text grep for `elevation`, `azimuth`, `orientation.*angle`, `spherical`, `R3` (case-insensitive, all `.py` files) | Only match: `CDR_Profiler/CDR3_Conf_Assigner.py` — the same file already fully inspected in the prior audit. Confirmed (re-read in full) to implement only `rho` (via `get_spherical_coordinates`, which literally returns only a scalar norm) and nothing resembling θ or φ. No new functions, no coordinate-frame construction, no angular computation anywhere in this file or elsewhere in TNP. |

## New relevant files discovered

**None.**

## Does Text S2 exist locally?

**No.** Confirmed absent from `references/`, `external/`, `docs/`, and the project root. The paper's own text confirms it is hosted only at the publisher's external website, which is outside this project's local source materials.

## Are exact θ/φ definitions now recoverable?

**No.** No new source material was found. The main paper (already inspected in the prior audit) gives only the qualitative description on p.4 (Section 2.9) and explicitly defers the exact coordinate-frame construction and formulas to Text S2, which remains unavailable. No corroborating author code for θ/φ exists anywhere in this repository.

---

## Conclusion

This search is definitive given the constraint to local repository material: no new evidence changes the prior audit's finding. **Gordon θ and φ remain an unresolved method dependency.** Per instruction, no substitute, approximation, or homemade coordinate-frame convention is implemented under the Gordon name.
