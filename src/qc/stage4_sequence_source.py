"""
Stage 4 — Authoritative Sequence Source and Position Projection

CRITICAL DESIGN NOTE (documented per spec §5 and the user's safeguard #1):

A residue can be:
  (a) present in the true polymer sequence AND resolved in coordinates
  (b) present in the true polymer sequence but NOT resolved in coordinates
  (c) genuinely absent from the polymer (e.g. an unused IMGT slot for a
      shorter loop, or truly beyond the end of the construct)

These states are NEVER collapsed. The authoritative source for "does a real
residue exist at all" is `_entity_poly_seq` (via
MMCIFStructure.build_residue_scheme, already validated in Stage 2/3);
`_atom_site` only tells us which of those real residues were resolved.

GENUINE AMBIGUITY (documented, not silently resolved): assigning an exact
IMGT auth_seq_id/insertion-code identity to an UNRESOLVED real residue
requires knowing IMGT's canonical gap-placement rule for that loop length,
which is NOT re-derived here (that would be equivalent to re-running ANARCI,
explicitly out of scope). This is only a problem for "ambiguous" base IMGT
positions — those that are observed carrying a non-empty insertion code
ANYWHERE in the Stage-3 cohort (empirically: {32, 33, 60, 61, 67, 68, 79,
82, 85, 111, 112, 122} — insertion codes are NOT limited to CDR3 in this
dataset, checked directly, not assumed). For all OTHER (non-ambiguous) base
positions, label_seq_id is contiguous and no insertion sibling ever
competes for that slot, so unresolved real residues between two observed
anchors can be assigned unambiguously by simple count-matched interpolation.

Four resulting states per (structure, target IMGT position):
  PRESENT_OBSERVED             - direct coordinate match
  PRESENT_UNRESOLVED           - real residue, unambiguously interpolated
  ABSENT_CONFIRMED             - provably no real residue occupies this slot
  AMBIGUOUS_UNRESOLVED         - real residue may exist nearby but exact
                                  slot cannot be determined without ANARCI;
                                  excluded from chemistry denominators,
                                  reported separately for audit (never
                                  silently merged into ABSENT or PRESENT).
"""

import sys
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from mmcif_qc import MMCIFStructure

VARIABLE_DOMAIN_START = 1
VARIABLE_DOMAIN_END = 128


def compute_ambiguous_base_positions(residue_df: pd.DataFrame) -> Set[int]:
    """
    Base IMGT auth positions that are observed carrying a non-empty
    insertion code ANYWHERE in the Stage-3 residue table. Empirically
    determined, never assumed to be limited to CDR3.
    """
    ins_rows = residue_df[residue_df['insertion_code'] != '']
    return set(int(p) for p in ins_rows['imgt_position'].unique())


class SequenceSource:
    """Authoritative per-structure real-residue scheme, IMGT-position projectable."""

    def __init__(self, structure: MMCIFStructure, chain_id: str, model_id: str, ambiguous_base_positions: Set[int]):
        self.structure = structure
        self.chain_id = chain_id
        self.model_id = model_id
        self.ambiguous_base_positions = ambiguous_base_positions

        full_scheme = structure.build_residue_scheme(chain_id, model_id)
        # Restrict to variable-domain-plausible residues: everything up to
        # the first residue whose OBSERVED auth_seq_id exceeds 128 by a
        # clear margin is excluded downstream by scoping target positions to
        # 1-128 anyway; keep the full scheme here since unresolved residues
        # near the boundary need context from both sides.
        self.scheme = full_scheme

        # Index of observed residues in label_seq_id order (already the scheme order)
        self.observed_indexed = [
            (i, r['auth_seq_id'], r['ins_code']) for i, r in enumerate(self.scheme)
            if r['observed'] and r['auth_seq_id'] is not None
        ]
        # Direct lookup: (auth_seq_id, ins_code) -> scheme index, for observed residues only
        self.observed_lookup = {
            (r['auth_seq_id'], r['ins_code']): i
            for i, r in enumerate(self.scheme) if r['observed'] and r['auth_seq_id'] is not None
        }
    
    def project(self, target_pos: int, target_ins: str) -> Dict:
        """
        Determine the state of one target (imgt_position, insertion_code)
        for this sequence.

        Returns dict: {state, amino_acid, coordinate_observed}
        """
        key = (target_pos, target_ins)

        # --- Case: direct observation ---
        if key in self.observed_lookup:
            idx = self.observed_lookup[key]
            return {
                'state': 'PRESENT_OBSERVED',
                'amino_acid': self.scheme[idx]['mon_id'],
                'coordinate_observed': True,
            }

        # --- Case: ambiguous base position, not directly observed ---
        if target_pos in self.ambiguous_base_positions:
            # Can we prove ABSENCE? If the bracket of real residues between the
            # nearest observed anchors strictly below/above this base position
            # contains FEWER real residues than would be needed to reach an
            # insertion this deep, it's confirmably absent. Otherwise ambiguous.
            lo_idx, hi_idx = self._find_anchors_by_base(target_pos)
            if hi_idx is None:
                return {'state': 'AMBIGUOUS_UNRESOLVED', 'amino_acid': None, 'coordinate_observed': False}

            bracket = self.scheme[lo_idx + 1:hi_idx]
            if len(bracket) == 0:
                return {'state': 'ABSENT_CONFIRMED', 'amino_acid': None, 'coordinate_observed': False}
            # Real residues exist in the bracket but we cannot prove which
            # exact insertion slot(s) they occupy without IMGT's placement
            # rule -> ambiguous, not silently assumed present or absent.
            return {'state': 'AMBIGUOUS_UNRESOLVED', 'amino_acid': None, 'coordinate_observed': False}

        # --- Case: non-ambiguous base position, not directly observed ---
        # No insertion sibling ever competes for this slot in this cohort, so
        # label_seq_id contiguity gives an unambiguous count-matched mapping.
        lo_idx, hi_idx = self._find_anchors_by_base(target_pos)
        if hi_idx is None:
            return {'state': 'AMBIGUOUS_UNRESOLVED', 'amino_acid': None, 'coordinate_observed': False}

        # lo_idx == -1 is the virtual "start of chain" boundary (see
        # _find_anchors_by_base docstring): treat it as auth position 0 so
        # range(lo_auth+1, hi_auth) correctly starts at IMGT position 1.
        lo_auth = self.scheme[lo_idx]['auth_seq_id'] if lo_idx >= 0 else 0
        hi_auth = self.scheme[hi_idx]['auth_seq_id']
        bracket = self.scheme[lo_idx + 1:hi_idx]

        # Expected number of non-ambiguous integer slots strictly between the
        # two anchors' base auth numbers (excluding any ambiguous base
        # numbers, which are handled by the branch above and must not be
        # double-counted here).
        expected_slots = [
            p for p in range(lo_auth + 1, hi_auth) if p not in self.ambiguous_base_positions
        ]

        if target_pos not in expected_slots:
            # target isn't even a candidate slot in this bracket at all
            # (e.g. it's an ambiguous base number handled elsewhere, or
            # outside the bracket range) -- should not normally reach here,
            # but treat conservatively as ambiguous rather than asserting
            # absence for a slot we didn't actually reason about.
            return {'state': 'AMBIGUOUS_UNRESOLVED', 'amino_acid': None, 'coordinate_observed': False}

        if len(bracket) == 0:
            # No real residues at all occupy this bracket -> every candidate
            # slot, including target, is unambiguously absent (a genuinely
            # shorter loop, IMGT's normal gapping behavior -- NOT missing
            # density, since there is nothing here to be missing).
            return {'state': 'ABSENT_CONFIRMED', 'amino_acid': None, 'coordinate_observed': False}

        if len(bracket) == len(expected_slots):
            # Exactly enough real residues to fill every candidate slot ->
            # unambiguous direct 1:1 mapping in label_seq_id order.
            slot_index = expected_slots.index(target_pos)
            residue = bracket[slot_index]
            return {
                'state': 'PRESENT_UNRESOLVED',
                'amino_acid': residue['mon_id'],
                'coordinate_observed': False,
            }

        # Partial fill: some but not all candidate slots are occupied by real
        # residues, and WHICH ones cannot be determined without IMGT's
        # canonical gap-placement rule (not re-derived here). Genuinely
        # ambiguous -- reported as such, never guessed.
        return {'state': 'AMBIGUOUS_UNRESOLVED', 'amino_acid': None, 'coordinate_observed': False}

    def _find_anchors_by_base(self, target_pos: int) -> Tuple[int, Optional[int]]:
        """
        Nearest observed scheme indices with auth_seq_id < target_pos (lo)
        and > target_pos (hi).

        lo defaults to -1 (a virtual "start of chain" boundary) when no
        observed anchor exists below the target -- this is always safe
        because label_seq_id=1 is, by definition, the first real residue of
        the entity's true polymer sequence (_entity_poly_seq); nothing can
        precede it, so an unresolved N-terminal stretch can be bracketed
        against the literal start of the scheme without any ambiguity.

        hi has NO such safe default: an unresolved run near the C-terminus
        of a chain that continues into a constant domain (e.g. FAB Hchain)
        cannot be safely bounded by "end of scheme" (that would risk
        reaching into CH1 territory), so hi_idx stays None (unresolvable
        without an observed anchor) when none is found -- correctly kept as
        AMBIGUOUS/unresolvable by the caller in that case.
        """
        lo_idx = -1
        for i, auth_pos, ins in self.observed_indexed:
            if auth_pos < target_pos:
                lo_idx = i
            elif auth_pos == target_pos:
                continue  # shouldn't reach here (direct-match branch handles ins=='' case already)
            else:
                break

        hi_idx = None
        for i, auth_pos, ins in reversed(self.observed_indexed):
            if auth_pos > target_pos:
                hi_idx = i
            elif auth_pos == target_pos:
                continue
            else:
                break

        return lo_idx, hi_idx
