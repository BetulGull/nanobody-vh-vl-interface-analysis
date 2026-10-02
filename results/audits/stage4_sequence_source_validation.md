
AUTHORITATIVE SEQUENCE SOURCE DOCUMENTATION (Stage 4, spec section 5)

Source field/table used:
  _entity_poly_seq (per-chain, label_seq_id-ordered, contiguous real-residue
  list) is the authoritative source for "does a real residue exist at all."
  _atom_site (cross-referenced via auth_seq_id/insertion code/model) supplies
  resolution status (observed vs unresolved) and coordinates when present.
  This is the same mechanism already validated in Stage 2/3
  (MMCIFStructure.build_residue_scheme), reused here without modification.

How IMGT positions are recovered:
  auth_seq_id + pdbx_PDB_ins_code in these SAbDab2-derived structures IS the
  IMGT position (validated extensively in Stage 2). For residues genuinely
  unresolved in density, exact IMGT identity is recovered by unambiguous
  anchor-bracket interpolation (see stage4_sequence_source.py) wherever the
  bracket is not itself ambiguous; genuinely ambiguous cases are flagged
  AMBIGUOUS_UNRESOLVED rather than guessed (see
  stage4_mapping_bias_audit.py / stage4_interface_mapping_status_by_position.csv
  for prevalence).

Agreement with observed structural residues:
  100% by construction and verified above: every residue reported as
  PRESENT_OBSERVED by SequenceSource.project() is required to match the
  amino acid identity found by direct atom_site lookup at that exact
  auth_seq_id/insertion-code -- these are not independent sources that could
  disagree, they are the same underlying data accessed through one audited
  code path.

How unresolved residues are represented:
  PRESENT_UNRESOLVED: real residue (from _entity_poly_seq), unambiguous
  interpolated IMGT position, amino acid identity known, coordinate_observed
  = False. Never conflated with ABSENT_CONFIRMED (a real residue does not
  exist at all at that slot -- e.g. an unused IMGT numbering gap for a
  shorter loop).

Positions that cannot be mapped confidently:
  AMBIGUOUS_UNRESOLVED -- see stage4_interface_mapping_status_by_position.csv
  for exact prevalence within the 28 primary interface positions (concentrated
  at the 3 insertion-coded ones: 111A, 112A, 112B).
