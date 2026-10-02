"""
Antigen QC — Protein-antigen entity identification from mmCIF.

Identifies candidate protein-antigen chains based on mmCIF entity/polymer information.
Never infers from compound text; uses actual mmCIF mappings only.
"""

from typing import Dict, List, Set, Tuple, Optional
import logging

logger = logging.getLogger(__name__)


class AntigenQuality:
    """
    Identify protein-antigen entities in a structure.
    """

    # Valid protein polymer types (from stage2.yaml)
    PROTEIN_POLYMER_TYPES = {
        'polypeptide(L)',
        'polypeptide(D)',
    }

    # Exclude these from antigen consideration
    EXCLUDE_TYPES = {
        'water',
        'ion',
        'non-polymer',
        'nucleotide',
        'nucleic-acid',
        'branched',
    }

    def __init__(self):
        pass

    def _normalize_cif_null(self, value: str) -> str:
        """CIF uses '?' (unknown) and '.' (inapplicable) as null sentinels."""
        if value in ('?', '.', ''):
            return ''
        return value

    def identify_antigen_chains(
        self,
        mmcif_structure,
        model_id: str,
        antibody_chains: Set[str],  # Hchain + Lchain (if present)
    ) -> Tuple[List[str], Dict[str, str]]:
        """
        Identify candidate protein-antigen chains.

        Args:
            mmcif_structure: MMCIFStructure instance
            model_id: Model number
            antibody_chains: Set of chain IDs belonging to the antibody (to exclude)

        Returns:
            (antigen_chain_ids, chain_to_type_map)
            antigen_chain_ids: List of chains that are protein antigens
            chain_to_type_map: Dict mapping each antigen chain to its polymer type
        """
        antigen_chains = []
        chain_to_type = {}

        # Get all chains in this model
        all_chains = mmcif_structure.get_chains_for_model(model_id)

        for chain_id in all_chains:
            # Skip antibody chains
            if chain_id in antibody_chains:
                continue

            # Get entity ID for this chain
            entity_id = mmcif_structure.get_entity_id_for_chain(chain_id)
            if not entity_id:
                logger.debug(f"Chain {chain_id}: no entity ID mapping")
                continue

            # Get entity info
            if entity_id not in mmcif_structure.entity_info:
                logger.debug(f"Chain {chain_id} (entity {entity_id}): no entity info")
                continue

            entity_type = mmcif_structure.entity_info[entity_id].get('type', 'unknown')

            # Check if it's a protein polymer (types are case-sensitive: 'polypeptide(L)')
            if entity_type in self.PROTEIN_POLYMER_TYPES:
                antigen_chains.append(chain_id)
                chain_to_type[chain_id] = entity_type
                logger.debug(f"Chain {chain_id}: identified as {entity_type}")

            elif any(exclude in entity_type.lower() for exclude in self.EXCLUDE_TYPES):
                logger.debug(f"Chain {chain_id}: excluded type {entity_type}")

            else:
                logger.debug(f"Chain {chain_id}: unknown type {entity_type}")

        return antigen_chains, chain_to_type

    def get_antigen_atoms(
        self,
        mmcif_structure,
        antigen_chains: List[str],
        model_id: str,
    ) -> Dict[str, list]:
        """
        Extract atom coordinates for antigen chains.

        Args:
            mmcif_structure: MMCIFStructure instance
            antigen_chains: List of antigen chain IDs
            model_id: Model number

        Returns:
            Dict mapping chain_id to list of atom dicts
        """
        result = {}

        if '_atom_site' not in mmcif_structure.data:
            return result

        for chain_id in antigen_chains:
            chain_atoms = []
            for row in mmcif_structure.data['_atom_site']:
                if (row.get('_atom_site.auth_asym_id') == chain_id and
                    row.get('_atom_site.pdbx_PDB_model_num') == model_id):

                    atom_dict = {
                        'name': row.get('_atom_site.label_atom_id', '?').strip(),
                        'element': row.get('_atom_site.type_symbol', '?').strip(),
                        'x': float(row.get('_atom_site.Cartn_x', '0')),
                        'y': float(row.get('_atom_site.Cartn_y', '0')),
                        'z': float(row.get('_atom_site.Cartn_z', '0')),
                        'occupancy': float(row.get('_atom_site.occupancy', '1.0')),
                        'bfactor': float(row.get('_atom_site.B_iso_or_equiv', '0.0')),
                        'altloc': self._normalize_cif_null(row.get('_atom_site.label_alt_id', '').strip()) or 'A',
                        'residue_number': int(row.get('_atom_site.auth_seq_id', '0')),
                        'insertion_code': self._normalize_cif_null(row.get('_atom_site.pdbx_PDB_ins_code', '').strip()),
                    }
                    chain_atoms.append(atom_dict)

            result[chain_id] = chain_atoms

        return result

    def has_protein_antigen(
        self,
        mmcif_structure,
        model_id: str,
        antibody_chains: Set[str],
    ) -> bool:
        """
        Quick check: does this structure have at least one protein-antigen candidate?
        """
        antigen_chains, _ = self.identify_antigen_chains(mmcif_structure, model_id, antibody_chains)
        return len(antigen_chains) > 0
