# Stage 8 — APBS / PDB2PQR Environment and Parameters

## Versions (verified, not trusted from Stage-0 records)

```
micromamba run -n vhh-elec apbs --version     -> APBS 3.4.1
micromamba run -n vhh-elec pdb2pqr30 --version -> pdb2pqr 3.6.1
```
Matches Stage-0 records exactly.

## PDB2PQR command

```
pdb2pqr30 --ff PARSE --with-ph 7.4 --titration-state-method propka --keep-chain \
  --apbs-input <instance>.in <instance>.pdb <instance>.pqr
```

- **Force field:** PARSE (the standard/default APBS-compatible force field; explicitly chosen rather than left as an undocumented default — PDB2PQR requires one of `--ff`/`--userff`/`--clean` to be specified, so no run proceeds on an implicit choice).
- **Protonation/titration method:** `propka` (pKa-based titration), at **pH 7.4** (physiological, chosen and documented BEFORE inspecting any group results, per the amendment's protonation policy). A pH 7.0 sensitivity subset is planned (Section 23) alongside the primary 7.4 run, not a scan for the largest group difference.
- **Histidine handling:** propka assigns each histidine's protonation state individually based on its computed pKa at the target pH (not a fixed pKa=6.0 default) — this is visible in the PDB2PQR log for every structure (`INFO:Applying pKa values at a pH of 7.40`) and is retained per-structure as an audit trail via the log file.
- **`--keep-chain`**: required — without it PDB2PQR silently renames/collapses chain IDs; with it, the surrogate single-character chain ID assigned by this project's own PDB writer (reusing the Stage-7 `write_pdb_with_surrogate_mapping` pattern) is preserved unchanged into the PQR file, confirmed by direct inspection (see `stage8_pdb2pqr_mapping_validation.csv`).
- **`--apbs-input`**: auto-generates a molecule-centered APBS input file from the PQR, used as-is (not hand-tuned) for reproducible, algorithmic (not ad hoc) grid construction — see below.

## APBS input parameters (from the PDB2PQR-generated template, inspected directly, not assumed)

```
mg-auto
dime 129 129 129
cglen/fglen: molecule-extent-derived (per-structure, see below)
cgcent/fgcent: mol 1 (grid centered on the molecule)
lpbe (linearized Poisson-Boltzmann equation)
bcfl sdh (single Debye-Huckel boundary condition)
pdie 2.0 (solute/protein dielectric)
sdie 78.54 (solvent dielectric, water at ~298K)
srfm smol (smoothed molecular surface)
chgm spl2 (cubic B-spline charge discretization)
sdens 10.00
srad 1.40 (solvent probe radius, Å -- matches the FreeSASA probe radius used throughout Stages 3-7)
swin 0.30
temp 298.15 K
calcenergy total
```

No ionic-strength (`ion` keyword) was added by PDB2PQR's default template — confirmed by direct inspection of the generated `.in` file (no `ion charge ... conc ... radius ...` line present), meaning the smoke-test runs use **zero explicit ionic strength** (pure Poisson-Boltzmann with only `sdie`/`pdie`, no Debye screening from added salt). This is documented explicitly here rather than silently assumed; a project decision on whether to add explicit physiological ionic strength (e.g., 150 mM) for the production run will be made and documented in `stage8_apbs_runtime_benchmark.csv`'s companion decision (Section K), not silently left as PDB2PQR's un-stated default.

## Grid construction policy

`dime` (129³) is fixed by PDB2PQR's `psize.py` grid-sizing heuristic for `mg-auto` (multigrid-compatible grid dimension). `cglen`/`fglen`/`cgcent`/`fgcent` are computed **per-structure** from the actual molecular extent with an automatic solvent margin — this is a reproducible, algorithmic policy (same rule applied identically to every structure), not a hand-tuned one-off. Verified directly (§ smoke test below): the resulting grid **fully encloses** the molecule with substantial margin in all three dimensions for both a VHH and a VH test structure. Per-structure grid dimensions, spacing, and physical box size are recorded for every production structure (Section 27 requirement) rather than assumed comparable.

## Output format

OpenDX (`.dx`), potential in units of **kT/e** (dimensionless electrostatic potential, thermal-energy-per-unit-charge) — confirmed directly from the `.dx` file header (`# POTENTIAL (kT/e)`), not assumed from memory.
