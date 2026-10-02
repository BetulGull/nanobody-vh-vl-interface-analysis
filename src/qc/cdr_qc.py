"""
CDR QC — CDR membership and backbone completeness checking.

Uses IMGT numbering definitions (via mmcif_structure.build_full_imgt_mapping_with_insertions)
to identify and validate CDR residues, correctly preserving insertion codes.
"""

from typing import Dict, List, Tuple, Optional
import logging
import re

logger = logging.getLogger(__name__)

# IMGT CDR definitions (from stage2.yaml)
IMGT_CDR_RANGES = {
    'H1': (27, 38),
    'H2': (56, 65),
    'H3': (105, 117),
    'L1': (27, 38),
    'L2': (56, 65),
    'L3': (105, 117),
}

BACKBONE_ATOMS = {'N', 'CA', 'C', 'O'}


_RESIDUE_STR_RE = re.compile(r'^(-?\d+)([A-Za-z]*)$')


def _parse_residue_str(residue_str: str) -> Tuple[int, str]:
    """Split '112A' -> (112, 'A'); '112' -> (112, '')."""
    m = _RESIDUE_STR_RE.match(residue_str)
    if not m:
        raise ValueError(f"Cannot parse residue string: {residue_str!r}")
    return int(m.group(1)), m.group(2)


class CDRQuality:
    """
    Validate CDR completeness for a given chain, using actual mmCIF residue
    observations (auth_seq_id + insertion code = IMGT position, per SAbDab2
    convention; no ANARCI re-run).
    """

    def __init__(self):
        pass

    def is_cdr_position(self, chain_type: str, imgt_pos: int) -> Tuple[bool, Optional[str]]:
        """
        Check if an IMGT position falls within a CDR for a chain type.

        Args:
            chain_type: 'H' or 'L'
            imgt_pos: IMGT position (1-128)

        Returns:
            (is_cdr, cdr_name)
        """
        for cdr_name, (start, end) in IMGT_CDR_RANGES.items():
            if cdr_name[0] != chain_type:
                continue
            if start <= imgt_pos <= end:
                return True, cdr_name
        return False, None

    def get_cdr_name(self, chain_type: str, imgt_pos: int) -> Optional[str]:
        is_cdr, cdr_name = self.is_cdr_position(chain_type, imgt_pos)
        return cdr_name

    def get_cdr_range(self, cdr_name: str) -> Tuple[int, int]:
        return IMGT_CDR_RANGES.get(cdr_name, (None, None))

    def _find_cdr_real_residues(
        self,
        scheme: List[Dict],
        cdr_start: int,
        cdr_end: int,
    ) -> Optional[List[Dict]]:
        """
        Given a full residue scheme (ordered by label_seq_id — the true,
        gap-free polymer order; see MMCIFStructure.build_residue_scheme),
        find the REAL residues that structurally belong to a CDR range by
        anchoring on the nearest OBSERVED framework residues immediately
        flanking the range on each side, then taking every scheme entry
        strictly between those two anchors (by label_seq_id order).

        This is what distinguishes "no real residue occupies this IMGT slot"
        (fine — the anchors are adjacent to it in the true sequence, nothing
        is missing) from "a real residue occupies this position/slot but was
        not resolved in the density" (the residue appears as an entry between
        the anchors with observed=False).

        Returns None if no anchor could be found on either side (should not
        happen for real antibody structures with resolved framework — flagged
        as an edge case for audit rather than silently guessed at).
        """
        observed_indexed = [
            (i, r['auth_seq_id']) for i, r in enumerate(scheme)
            if r['observed'] and r['auth_seq_id'] is not None
        ]

        lo_idx = None
        for i, auth_pos in observed_indexed:
            if auth_pos <= cdr_start - 1:
                lo_idx = i  # keep advancing to the LATEST such anchor (closest before CDR)
            else:
                break

        hi_idx = None
        for i, auth_pos in reversed(observed_indexed):
            if auth_pos >= cdr_end + 1:
                hi_idx = i  # keep advancing to the EARLIEST such anchor (closest after CDR)
            else:
                break

        if lo_idx is None or hi_idx is None or hi_idx <= lo_idx:
            return None

        return scheme[lo_idx + 1:hi_idx]

    def check_cdr_completeness(
        self,
        mmcif_structure,
        chain_id: str,
        chain_type: str,
        model_id: str,
    ) -> Dict[str, bool]:
        """
        Check which CDRs have every REAL residue backbone-complete.

        Correctness rule (distinguishes two different kinds of "absence"):
          - An IMGT numbering slot with NO real residue (loop shorter than
            the scheme's maximum length, e.g. CDR-H1 27-38 gapped in the
            middle for an 8-residue loop) does NOT count against
            completeness — there was never a residue there to observe.
          - A REAL residue (present in the polymer sequence per
            _entity_poly_seq) that has zero atoms in the density, or an
            incomplete N/CA/C/O backbone, DOES count against completeness.

        This is determined via MMCIFStructure.build_residue_scheme, which
        uses label_seq_id (contiguous 1..N over the true polymer, independent
        of IMGT auth_seq_id) cross-referenced against _entity_poly_seq (the
        full real-residue list) to tell real-but-unobserved residues apart
        from IMGT placeholder gaps. The CDR's real-residue span is located by
        anchoring on the nearest observed framework residues immediately
        before/after the CDR range (see _find_cdr_real_residues).

        Returns:
            Dict mapping CDR name (e.g., 'H1') to bool. A CDR whose anchors
            could not be located (should not happen for well-resolved
            antibody framework) is reported as False (unknown treated
            conservatively as incomplete) — such cases are surfaced via
            get_missing_cdr_positions for audit.
        """
        scheme = mmcif_structure.build_residue_scheme(chain_id, model_id)

        result = {}
        for cdr_num in ['1', '2', '3']:
            cdr_name = f'{chain_type}{cdr_num}'
            cdr_start, cdr_end = self.get_cdr_range(cdr_name)
            if cdr_start is None:
                continue

            if not scheme:
                result[cdr_name] = False
                continue

            real_residues = self._find_cdr_real_residues(scheme, cdr_start, cdr_end)
            if real_residues is None:
                result[cdr_name] = False  # anchors not found — conservative
                continue

            if len(real_residues) == 0:
                # Genuinely zero real residues occupy this CDR span (extremely
                # short loop) — nothing to be incomplete.
                result[cdr_name] = True
                continue

            result[cdr_name] = all(
                r['observed'] and r['backbone_complete'] for r in real_residues
            )

        return result

    def check_h3_geometry_anchor(
        self,
        mmcif_structure,
        chain_id: str,
        model_id: str,
    ) -> bool:
        """
        Check if CDR-H3 anchor residues (3 before and after CDR-H3, IMGT
        102-104 and 118-120) are present with complete backbone.

        Stage-2 only; no geometry calculated — flags readiness for Stage-6.
        """
        full_mapping = mmcif_structure.build_full_imgt_mapping_with_insertions(chain_id, model_id)
        flanking_positions = [102, 103, 104, 118, 119, 120]

        for pos in flanking_positions:
            residue_variants = full_mapping.get(pos, [])
            if not residue_variants:
                return False

            for residue_str in residue_variants:
                base_pos, ins_code = _parse_residue_str(residue_str)
                atoms = mmcif_structure.get_atoms_for_chain_residue(
                    chain_id, base_pos, model_id, ins_code
                )
                if atoms:
                    atoms = mmcif_structure.select_altloc(atoms)
                if not atoms or not mmcif_structure.check_backbone_complete(atoms):
                    return False

        return True

    def get_missing_cdr_positions(
        self,
        mmcif_structure,
        chain_id: str,
        chain_type: str,
        model_id: str,
    ) -> Dict[str, List[Dict]]:
        """
        Identify REAL CDR residues that are missing or backbone-incomplete —
        i.e. entries in the true polymer sequence (_entity_poly_seq) within
        the CDR's real-residue span that are unobserved or backbone-
        incomplete. IMGT numbering placeholder gaps are never reported here,
        since they are not real residues (see check_cdr_completeness).

        Returns:
            Dict mapping CDR name to a list of dicts describing each missing
            real residue: {label_seq_id, mon_id, auth_seq_id, ins_code,
            reason} where reason is 'UNOBSERVED' or 'BACKBONE_INCOMPLETE'.
        """
        scheme = mmcif_structure.build_residue_scheme(chain_id, model_id)
        result = {}

        for cdr_num in ['1', '2', '3']:
            cdr_name = f'{chain_type}{cdr_num}'
            cdr_start, cdr_end = self.get_cdr_range(cdr_name)
            if cdr_start is None:
                continue

            if not scheme:
                result[cdr_name] = []
                continue

            real_residues = self._find_cdr_real_residues(scheme, cdr_start, cdr_end)
            if real_residues is None:
                result[cdr_name] = []
                continue

            missing = []
            for r in real_residues:
                if not r['observed']:
                    missing.append({
                        'label_seq_id': r['label_seq_id'], 'mon_id': r['mon_id'],
                        'auth_seq_id': None, 'ins_code': '', 'reason': 'UNOBSERVED',
                    })
                elif not r['backbone_complete']:
                    missing.append({
                        'label_seq_id': r['label_seq_id'], 'mon_id': r['mon_id'],
                        'auth_seq_id': r['auth_seq_id'], 'ins_code': r['ins_code'],
                        'reason': 'BACKBONE_INCOMPLETE',
                    })

            result[cdr_name] = missing

        return result

    def get_cdr_atom_coords(
        self,
        mmcif_structure,
        chain_id: str,
        chain_type: str,
        cdr_name: str,
        model_id: str,
    ) -> List[Tuple[float, float, float]]:
        """
        Get all atom coordinates (all atoms, not just backbone) for a given CDR.

        Used for Gordon contact-search input.
        """
        full_mapping = mmcif_structure.build_full_imgt_mapping_with_insertions(chain_id, model_id)
        cdr_start, cdr_end = self.get_cdr_range(cdr_name)
        if cdr_start is None:
            return []

        coords = []
        for pos in range(cdr_start, cdr_end + 1):
            for residue_str in full_mapping.get(pos, []):
                base_pos, ins_code = _parse_residue_str(residue_str)
                atoms = mmcif_structure.get_atoms_for_chain_residue(
                    chain_id, base_pos, model_id, ins_code
                )
                if atoms:
                    atoms = mmcif_structure.select_altloc(atoms)
                for atom in atoms:
                    coords.append((atom['x'], atom['y'], atom['z']))

        return coords
