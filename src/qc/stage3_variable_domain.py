"""
Stage 3 — Variable-Domain Atom Extraction

Extracts VH/VL variable-domain-only atom subsets for ΔSASA calculation.

CRITICAL: this is residue-IDENTITY based (uses the validated Stage-2 IMGT
auth_seq_id/insertion-code mapping — see src/qc/mmcif_qc.py
build_residue_scheme, itself cross-referenced against _entity_poly_seq),
NOT a chain-position or "first N residues" cutoff. A FAB Hchain that
physically contains both VH and CH1 in one chain is handled correctly
because CH1 residues carry IMGT/auth positions > 128 (validated below on a
real FAB structure, not assumed) — the boundary is drawn on IMGT identity,
not on counting residues from the chain's start.
"""

from typing import Dict, List, Any
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from mmcif_qc import MMCIFStructure

VARIABLE_DOMAIN_IMGT_START = 1
VARIABLE_DOMAIN_IMGT_END = 128

# Non-hydrogen element filter (hydrogens are almost never present in X-ray
# structures here, but exclude explicitly if ever deposited)
HYDROGEN_ELEMENTS = {'H', 'D'}


def extract_variable_domain_atoms(
    structure: MMCIFStructure,
    chain_id: str,
    model_id: str,
) -> List[Dict[str, Any]]:
    """
    Extract non-hydrogen atoms belonging ONLY to the IMGT variable domain
    (1-128) of a chain, using real-residue identity from build_residue_scheme.

    Returns:
        List of atom dicts, each carrying full residue-identity provenance:
        auth_asym_id, label_seq_id, auth_seq_id, ins_code, residue_name,
        atom_name, element, x, y, z, occupancy, altloc.
    """
    scheme = structure.build_residue_scheme(chain_id, model_id)

    atoms_out = []
    for r in scheme:
        if not r['observed'] or r['auth_seq_id'] is None:
            continue
        if not (VARIABLE_DOMAIN_IMGT_START <= r['auth_seq_id'] <= VARIABLE_DOMAIN_IMGT_END):
            continue

        residue_atoms = structure.get_atoms_for_chain_residue(
            chain_id, r['auth_seq_id'], model_id, r['ins_code']
        )
        if not residue_atoms:
            continue
        residue_atoms = structure.select_altloc(residue_atoms)

        for atom in residue_atoms:
            if atom['element'].strip().upper() in HYDROGEN_ELEMENTS:
                continue
            atoms_out.append({
                'auth_asym_id': chain_id,
                'label_seq_id': r['label_seq_id'],
                'auth_seq_id': r['auth_seq_id'],
                'ins_code': r['ins_code'],
                'residue_name': r['mon_id'],
                'atom_name': atom['name'],
                'element': atom['element'],
                'x': atom['x'],
                'y': atom['y'],
                'z': atom['z'],
                'occupancy': atom['occupancy'],
                'altloc': atom['altloc'],
            })

    return atoms_out


def get_variable_domain_imgt_span(atoms: List[Dict[str, Any]]) -> tuple:
    """Return (first_imgt_position_str, last_imgt_position_str) e.g. ('1', '128')."""
    if not atoms:
        return (None, None)
    positions = sorted(set((a['auth_seq_id'], a['ins_code']) for a in atoms), key=lambda p: (p[0], p[1]))
    first = f"{positions[0][0]}{positions[0][1]}"
    last = f"{positions[-1][0]}{positions[-1][1]}"
    return (first, last)
