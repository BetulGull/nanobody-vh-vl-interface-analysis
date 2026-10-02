# Stage 8 — TAP Source Gate (TAP vs TNP, explicit)

**Search performed:** `external/`, `references/`, `docs/`, project root — for Raybould et al. 2019 (Therapeutic Antibody Profiler, TAP), PSH/PPC/PNC/SFvCSP, and any TAP-specific supplementary material or author code.

## 1. Is the TAP paper locally available?

**NO.** No Raybould 2019 PDF, no TAP-specific supplement, anywhere in the repository. `references/` contains only `Gordon_2023.pdf` and `Zhang_2025.pdf`.

## 2. Is TAP supplementary material locally available?

**NO.**

## 3. Is original TAP author code locally available?

**NO** direct Raybould/TAP repository. However:

## 4. What does `external/TNP` actually contain?

`external/TNP` is the Therapeutic **Nanobody** Profiler (Gordon, Deane et al., OPIG lab — same group as Stage 6's Gordon 2023). Direct inspection of its manuscript (`external/TNP/paper/manuscript/TNP_MAIN.pdf`, extracted via `pdftotext`) contains this **explicit, unambiguous statement**:

> "Patches of surface hydrophobicity and charge were calculated following methods from TAP [2]."

and in the Acknowledgements:

> "We extend thanks to Matthew Raybould and the co-authors of the original Therapeutic Antibody Profiler [2]."

Reference [2] is confirmed in the bibliography extraction to be exactly: *Raybould, M. I. J. et al. Five computational developability guidelines for therapeutic antibody profiling. PNAS 116, 4025-4030 (2019).*

This is **not** merely a name coincidence (the risk this gate was designed to catch) — it is the TNP authors' own explicit statement that their PSH/PPC/PNC implementation (`Hydrophobicity_and_Charge_Assigner.py`, output field names `Patch_Hydrophob_CDR`, `Patch_Pos_Charge_CDR`, `Patch_Neg_Charge_CDR`; plotting templates literally named `PSH_D3_Template.py`, `PPC_D3_Template.py`, `PNC_D3_Template.py`) is a direct reproduction of the TAP methodology, restricted to the H-chain/CDR-vicinity context appropriate for nanobodies (TAP's own primary anatomical context for these three metrics is already the CDR-vicinity paratope region, per the TNP paper's own framing: "patches of surface hydrophobicity, positive charge, and negative charge across the CDR vicinity" as one of "the TAP metrics that can be directly applied to nanobodies").

**Authority classification: this is priority-2 evidence (supplementary-methods-equivalent) — an explicit, citation-confirmed reproduction of the TAP mathematical definition, not an independent reinterpretation.** Used accordingly (never treated as literal priority-1 Raybould author code, but not dismissed as "merely similarly named" either).

## 5. Per-metric evidence and status

### PSH / PPC / PNC — **EXACTLY_RESOLVED** (formula); **PARTIALLY_RESOLVED** (one denominator sub-component, see below)

Full formula recovered directly from `Hydrophobicity_and_Charge_Assigner.py` (quoted/described exactly, not from memory):

1. **Hydrophobicity scale**: Kyte & Doolittle (1982) raw values, confirmed byte-for-byte against `external/TNP/lib/python/theraprofnano/Hydrophobicity_and_Charge_Profiler/dat/Hydrophobics.txt` column 0 (Ile 4.5, Val 4.2, ..., Arg -4.5 — exact KD table).
2. **Normalization**: each residue's KD value is min-max rescaled across the 20 canonical amino acids to the range [1, 2]: `annotation[aa] = (KD[aa] - min(KD)) / (max(KD) - min(KD)) + 1.0`. This is a **fixed, precomputed per-residue-type table**, independent of any specific structure.
3. **Exposure classification**: a residue is exposed if its **relative side-chain accessible surface area (as computed by the `psa` program's "Total Side" `Per.` column) is >= 7.5%**; otherwise it is excluded entirely from patch scoring (`if res_asa < 7.5: continue`). Confirmed exact threshold value from source, not from memory.
4. **Charge classes**: ASP/GLU = -1 (negative); ARG/LYS = +1 (positive); HIS = +0.1 (positive, "appropriate at pH 7.4 according to Henderson-Hasselbalch," explicit in source comment); all other residues = 0 (no charge class). This exact pH-7.4 His partial-charge convention is used as-is.
5. **Vicinity/neighborhood construction**: starting from a defined "core" residue set (in TNP: CDR ± 2 per IMGT/Chothia definition), any other **exposed** residue is added to the vicinity set if the **minimum pairwise distance between ANY atoms (all heavy atoms, not side-chain-only) of the two residues is < 4.0 Å**. This is the exact "4 Å exposed neighborhood" rule anticipated in the original Stage-8 prompt (Section 5) — confirmed from source, not assumed.
6. **Patch score** (for each of the three classes: hydrophobic, positive, negative): for every ordered pair of DISTINCT residues (r1, r2) **both belonging to the vicinity set**, both carrying a defined value for that class, and separated by minimum all-heavy-atom distance `d <= 7.5 Å`: accumulate `score[r1] * score[r2] / d^2` (a Coulombic-like 1/d² spatial decay weighting) into that class's total. `PSH` = hydrophobic-class total; `PPC` = positive-class total; `PNC` = negative-class total. Confirmed exact formula and both cutoffs (4.0 Å for vicinity membership, 7.5 Å for pairwise patch-scoring inclusion) from source.
7. **Salt-bridge neutralization** (a real, source-confirmed detail not anticipated in the original prompt): before patch scoring, Lys/Arg-Asp/Glu pairs with a side-chain donor-acceptor atom distance < 3.2 Å have BOTH residues' charge and hydrophobicity contributions zeroed out (assumed to be in a permanent, physiologically-neutral salt bridge). This project reproduces this rule exactly since it is part of the verified source formula, not an invented addition.

**Sub-component requiring a documented, honest substitution — the relative-ASA denominator:**

The `psa` program itself (author binary, `external/TNP/bin/psa_mac`, Mach-O x86_64, version 2.0 Feb 2000) is **present and directly executable** on this machine (confirmed: ran successfully via Rosetta on a real project structure, `pdb_00001jtp`). Its output format was inspected directly (not assumed) and matches the parser exactly: fixed-width columns, with byte offsets `[55:61]` = absolute side-chain Å² ("Total Side Sum") and `[61:67]` = relative side-chain percentage ("Total Side Per.") — confirmed to be exactly the two values `Hydrophobicity_and_Charge_Assigner.py` reads (`asa_area` and `asa` respectively).

`psa`'s internal per-residue-type maximum reference ASA table (the denominator that turns absolute Å² into the "Per." percentage) is **compiled into the binary** — no external, human-readable radius/reference-library file is bundled (`-R file` accepts an optional radius library, but none is shipped, so the built-in default is used; confirmed via `ls external/TNP/bin/` — only the four binaries/scripts, no `.lib`/`.dat`/`.rad` reference file). This one specific sub-component is therefore **PARTIALLY_RESOLVED**: the threshold value (7.5%) and the underlying quantity (side-chain relative ASA) are exactly source-confirmed, but the exact numeric reference-max-ASA table used internally by the compiled binary is not independently inspectable as text.

**Resolution adopted (source-faithful, not an approximation):** rather than attempting to reverse-engineer or substitute a different published max-ASA table (which would introduce an undocumented numerical discrepancy against the verified formula), **this project runs the actual `psa_mac` binary directly** on the same minimal single-chain PDB input format the original TNP code itself constructs (see `stage8_patch_scoring.py`), for both the standard CDR-vicinity descriptors AND the project's adapted former-interface descriptors. This uses the author-validated tool as the literal source of the relative-ASA percentage, sidestepping the need to reconstruct its internal table at all. FreeSASA (Stages 3-6's validated engine) is retained in parallel for the `total_sasa`/`sidechain_sasa` audit columns in the patch-component table (Section 14), so both are visible and never silently conflated — but the specific TAP/TNP-formula relative-ASA percentage and exposure classification used for PSH/PPC/PNC comes directly from `psa`, not from a FreeSASA-based approximation.

### SFvCSP — **UNRESOLVED_METHOD_DEPENDENCY**

No occurrence of "SFvCSP," "charge symmetry," or "Fv charge symmetry parameter" found anywhere in `external/TNP` (code or paper text) or elsewhere in the repository. TAP's SFvCSP is inherently a paired-Fv (VH+VL) descriptor; even if it were resolved, it would still be `NA` for VHH per the standing rule (Section 8 of the original Stage-8 spec). For conventional paired Fv, SFvCSP is **not implemented** in this project — genuinely unavailable, not fabricated. Reported as `NA` / `UNRESOLVED_METHOD_DEPENDENCY` for all structures, VHH and conventional VH alike.

## 6. CDR-vicinity rule (standard TNP context)

Confirmed exact: IMGT CDR definitions used by TNP are H1=27-38, H2=54-67 (**note: not the project's own Stage-6-frozen CDR2 56-65** — TNP's own internal CDR2 boundary differs slightly, used only for the *standard* TNP-context descriptor, never substituted into this project's frozen empirical interface or CDR3 definitions), H3=103-119 (also differs slightly from this project's frozen CDR3 105-117 used everywhere else — again, used only inside the standard-TNP-context calculation, kept strictly separate). These TNP-internal region boundaries are used **only** when computing the *standard* (non-adapted) TAP-style descriptor for direct comparability with TNP's own published numbers; every *adapted* former-interface descriptor in this project uses the project's own frozen Stage-3/4/6 definitions (28-position empirical interface, IMGT CDR3 105-117) exclusively.

## 7. Summary table

| Metric | Status | Source |
|---|---|---|
| PSH (formula) | EXACTLY_RESOLVED | `Hydrophobicity_and_Charge_Assigner.py`, citation-confirmed TAP reproduction |
| PPC (formula) | EXACTLY_RESOLVED | same |
| PNC (formula) | EXACTLY_RESOLVED | same |
| Relative side-chain ASA denominator | PARTIALLY_RESOLVED → resolved operationally by running the real `psa` binary directly | `external/TNP/bin/psa_mac`, empirically verified output format |
| Exposure threshold (7.5%) | EXACTLY_RESOLVED | source code, confirmed exact |
| CDR-vicinity distance rule (4.0 Å) | EXACTLY_RESOLVED | source code, confirmed exact |
| Patch pairwise cutoff (7.5 Å) | EXACTLY_RESOLVED | source code, confirmed exact |
| Hydrophobicity normalization | EXACTLY_RESOLVED | source code + data file, confirmed exact |
| Charge classes / pH-7.4 His convention | EXACTLY_RESOLVED | source code, confirmed exact |
| SFvCSP | UNRESOLVED_METHOD_DEPENDENCY | not present anywhere locally; not invented |

**Per Section D/E of the amendment: PSH, PPC, and PNC are all sufficiently resolved (formula + operational resolution of the one denominator ambiguity via the real author binary) to implement under their published names, both as standard TNP-context descriptors and as adapted former-interface descriptors (same mathematical rule, project's own anatomical target). SFvCSP is not implemented for any structure and is not fabricated for VHH.**
