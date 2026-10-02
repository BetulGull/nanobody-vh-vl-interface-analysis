"""
mmCIF QC — mmCIF file parsing and structural mapping.

Extracts chain/entity/polymer information, maps variable domains via IMGT numbering,
preserves insertion codes, and validates model selection.
"""

from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any
import logging
import gemmi

logger = logging.getLogger(__name__)


class mmCIFParseError(Exception):
    """Raised when mmCIF parsing fails."""
    pass


class MMCIFStructure:
    """
    Represents parsed mmCIF structure with chain/entity/polymer mappings.
    """

    def __init__(self, mmcif_path: str):
        """
        Parse mmCIF file and extract structural information.

        Args:
            mmcif_path: Path to .cif file
        """
        self.mmcif_path = Path(mmcif_path)
        if not self.mmcif_path.exists():
            raise FileNotFoundError(f"mmCIF file not found: {mmcif_path}")

        self.data = {}  # raw parsed mmCIF blocks
        self.models = []  # list of model numbers present
        self.chains_per_model = {}  # model_id -> [chain_ids]
        self.entity_info = {}  # entity_id -> {type, polymers_yes_no, polymer_type, ...}
        self.asym_to_entity = {}  # chain_id -> entity_id
        self.sequences = {}  # chain_id -> sequence string (1-letter code)

        self._parse_mmcif()

    def _parse_mmcif(self):
        """
        Parse mmCIF file using gemmi's CIF reader.

        A hand-rolled line-based tokenizer was tried initially and found to
        silently misalign columns on files containing multi-line semicolon
        ( ;...; ) text fields — a valid and common CIF construct. gemmi
        implements the full CIF grammar and is used here instead; only the
        categories this project needs (_atom_site, _entity, _entity_poly)
        are materialized into self.data, in the same
        {tag_name: value_string} row-dict shape the rest of this module
        (and cdr_qc.py / antigen_qc.py / pairing_qc.py) already expects.
        """
        try:
            doc = gemmi.cif.read(str(self.mmcif_path))
        except Exception as e:
            raise mmCIFParseError(f"gemmi failed to parse {self.mmcif_path}: {e}")

        if len(doc) == 0:
            raise mmCIFParseError(f"No data blocks found in {self.mmcif_path}")

        block = doc.sole_block() if len(doc) == 1 else doc[0]

        for category in ('_atom_site', '_entity', '_entity_poly', '_struct_asym', '_entity_poly_seq'):
            self.data[category] = self._read_category_as_dicts(block, category)

        if not self.data.get('_atom_site'):
            raise mmCIFParseError(f"No _atom_site records found in {self.mmcif_path}")

        self._extract_models_and_chains()
        self._extract_entity_info()
        self._extract_asym_mapping()

    def _read_category_as_dicts(self, block: "gemmi.cif.Block", category: str) -> List[Dict[str, str]]:
        """
        Read all columns of a mmCIF category (loop_ or single-row) into a
        list of {full_tag_name: value} dicts, e.g. {'_atom_site.auth_asym_id': 'A', ...}.

        Uses gemmi's find_mmcif_category, which correctly handles quoting,
        multi-line values, and column order as actually present in the file.
        """
        table = block.find_mmcif_category(category + '.')
        if not table or table.width() == 0:
            return []

        tags = list(table.tags)  # gemmi returns fully-qualified tags, e.g. '_atom_site.group_PDB'
        rows = []
        for row in table:
            d = {}
            for i, tag in enumerate(tags):
                d[tag] = row[i]
            rows.append(d)
        return rows

    def _extract_models_and_chains(self):
        """Extract model numbers and chain IDs."""
        if '_atom_site' not in self.data:
            logger.warning("No _atom_site block found")
            return

        atom_site = self.data['_atom_site']
        models_set = set()
        chains_set = set()

        for row in atom_site:
            model_id = row.get('_atom_site.pdbx_PDB_model_num', '1')
            chain_id = row.get('_atom_site.auth_asym_id', '?')
            models_set.add(model_id)
            chains_set.add(chain_id)

        self.models = sorted(list(models_set))
        if not self.models:
            self.models = ['1']

        # Build chains per model
        for model_id in self.models:
            model_chains = set()
            for row in atom_site:
                m = row.get('_atom_site.pdbx_PDB_model_num', '1')
                if m == model_id:
                    chain_id = row.get('_atom_site.auth_asym_id', '?')
                    if chain_id and chain_id != '?':
                        model_chains.add(chain_id)
            self.chains_per_model[model_id] = sorted(list(model_chains))

    def _extract_entity_info(self):
        """
        Extract entity type and polymer information.

        _entity.type only distinguishes 'polymer'/'non-polymer'/'water'.
        The actual polymer chemistry (polypeptide(L) vs polyribonucleotide etc.)
        is in _entity_poly.type — required to distinguish protein antigens
        from nucleic-acid chains.
        """
        if '_entity' not in self.data:
            return

        entity_data = self.data['_entity']
        for row in entity_data:
            entity_id = row.get('_entity.id', '?')
            entity_type = row.get('_entity.type', 'unknown')
            self.entity_info[entity_id] = {
                'type': entity_type,
            }

        # Overlay with more specific polymer type from _entity_poly, when present
        if '_entity_poly' in self.data:
            for row in self.data['_entity_poly']:
                entity_id = row.get('_entity_poly.entity_id', '?')
                poly_type = row.get('_entity_poly.type', None)
                if entity_id in self.entity_info and poly_type:
                    self.entity_info[entity_id]['type'] = poly_type

    def _extract_asym_mapping(self):
        """
        Map author chain IDs (auth_asym_id) to entity IDs.

        IMPORTANT: label_asym_id (used by _struct_asym) frequently differs from
        auth_asym_id (used everywhere else in this project, incl. SAbDab2
        Hchain/Lchain fields and antigen_chain). Mapping via _struct_asym alone
        silently mismatches. Instead, derive auth_asym_id -> entity_id directly
        from _atom_site, which carries label_entity_id per atom alongside
        auth_asym_id.
        """
        if '_atom_site' not in self.data:
            return

        for row in self.data['_atom_site']:
            auth_chain = row.get('_atom_site.auth_asym_id', '?')
            entity_id = row.get('_atom_site.label_entity_id', '?')
            if auth_chain and auth_chain != '?' and entity_id and entity_id != '?':
                if auth_chain not in self.asym_to_entity:
                    self.asym_to_entity[auth_chain] = entity_id

    def _normalize_cif_null(self, value: str) -> str:
        """CIF uses '?' (unknown) and '.' (inapplicable) as null sentinels."""
        if value in ('?', '.', ''):
            return ''
        return value

    def get_model(self, model_id: Optional[str] = None) -> str:
        """Get the model ID to use (default: first model)."""
        if model_id is None:
            return self.models[0] if self.models else '1'
        return model_id

    def get_chains_for_model(self, model_id: str) -> List[str]:
        """Get chain IDs present in a specific model."""
        return self.chains_per_model.get(model_id, [])

    def get_entity_id_for_chain(self, chain_id: str) -> Optional[str]:
        """Get entity ID for a chain."""
        return self.asym_to_entity.get(chain_id)

    def get_atoms_for_chain_residue(
        self,
        chain_id: str,
        imgt_position: int,
        model_id: str,
        insertion_code: str = '',
    ) -> List[Dict[str, Any]]:
        """
        Get all atoms for a specific residue in a chain.

        Args:
            chain_id: Chain identifier (e.g., 'A')
            imgt_position: IMGT position (e.g., 27)
            model_id: Model number
            insertion_code: Insertion code if present (e.g., 'A', '', 'B')

        Returns:
            List of atom dicts with fields: name, element, coords, occupancy, bfactor, altloc
        """
        if '_atom_site' not in self.data:
            return []

        atoms = []
        for row in self.data['_atom_site']:
            if (row.get('_atom_site.auth_asym_id') == chain_id and
                row.get('_atom_site.pdbx_PDB_model_num') == model_id):

                # Check residue number and insertion code
                res_num_str = row.get('_atom_site.auth_seq_id', '?')
                ins_code = self._normalize_cif_null(row.get('_atom_site.pdbx_PDB_ins_code', '').strip())

                try:
                    res_num = int(res_num_str)
                except (ValueError, TypeError):
                    continue

                if res_num == imgt_position and ins_code == insertion_code:
                    # NOTE: these mmCIF files carry label_atom_id / label_alt_id,
                    # not auth_atom_id / auth_alt_id (those fields do not exist here).
                    atom_dict = {
                        'name': row.get('_atom_site.label_atom_id', '?').strip(),
                        'element': row.get('_atom_site.type_symbol', '?').strip(),
                        'x': float(row.get('_atom_site.Cartn_x', '0')),
                        'y': float(row.get('_atom_site.Cartn_y', '0')),
                        'z': float(row.get('_atom_site.Cartn_z', '0')),
                        'occupancy': float(row.get('_atom_site.occupancy', '1.0')),
                        'bfactor': float(row.get('_atom_site.B_iso_or_equiv', '0.0')),
                        'altloc': self._normalize_cif_null(row.get('_atom_site.label_alt_id', '').strip()) or 'A',
                    }
                    atoms.append(atom_dict)

        return atoms

    def select_altloc(self, atoms: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """
        Select one conformer per atom using highest-occupancy policy.

        Args:
            atoms: List of atom dicts from a single residue

        Returns:
            List with one conformer per atom name (highest occupancy, alphabetical tie-break)
        """
        if not atoms:
            return []

        # Group by atom name
        by_name = {}
        for atom in atoms:
            name = atom['name']
            if name not in by_name:
                by_name[name] = []
            by_name[name].append(atom)

        # Select best conformer per atom
        result = []
        for name, conformers in by_name.items():
            # Sort by occupancy (descending), then altloc (ascending for alphabetical)
            conformers.sort(key=lambda x: (-x['occupancy'], x['altloc']))
            result.append(conformers[0])

        return result

    def get_backbone_atoms(self, atoms: List[Dict[str, Any]]) -> Dict[str, Optional[Tuple[float, float, float]]]:
        """
        Extract backbone atoms (N, CA, C, O) from atom list.

        Returns:
            Dict mapping atom name ('N', 'CA', 'C', 'O') to (x, y, z) or None if missing
        """
        backbone = {'N': None, 'CA': None, 'C': None, 'O': None}
        for atom in atoms:
            name = atom['name']
            if name in backbone:
                backbone[name] = (atom['x'], atom['y'], atom['z'])
        return backbone

    def check_backbone_complete(self, atoms: List[Dict[str, Any]]) -> bool:
        """Check if all 4 backbone atoms (N, CA, C, O) are present."""
        backbone = self.get_backbone_atoms(atoms)
        return all(v is not None for v in backbone.values())

    def build_residue_scheme(self, chain_id: str, model_id: str) -> List[Dict[str, Any]]:
        """
        Build the full, real polymer-residue scheme for a chain, distinguishing
        residues that truly exist in the construct sequence from IMGT numbering
        placeholders (which are not real residues at all).

        Source of truth for "which residues really exist": _entity_poly_seq,
        keyed by label_seq_id, which is a CONTIGUOUS 1..N index over the actual
        polymer chain (no gaps — unlike auth_seq_id/IMGT numbering, which
        deliberately skips numbers for shorter loops). Cross-referencing
        _atom_site's label_seq_id (present on every atom row, observed or not
        — well, only observed atoms appear in _atom_site at all) against this
        full list tells us, for every REAL residue, whether it was resolved in
        the density and, if so, its auth_seq_id/insertion code (its IMGT
        position) and backbone completeness.

        Returns:
            Ordered list (by label_seq_id) of dicts:
                {label_seq_id, mon_id, observed (bool),
                 auth_seq_id (int or None), ins_code (str),
                 backbone_complete (bool or None; None if unobserved)}
        """
        entity_id = self.get_entity_id_for_chain(chain_id)
        if entity_id is None:
            return []

        poly_seq_rows = [
            r for r in self.data.get('_entity_poly_seq', [])
            if r.get('_entity_poly_seq.entity_id') == entity_id
        ]
        if not poly_seq_rows:
            return []

        # label_seq_id -> observed atoms (any atom rows) for this chain/model
        observed_by_label: Dict[int, List[Dict[str, Any]]] = {}
        for row in self.data.get('_atom_site', []):
            if (row.get('_atom_site.auth_asym_id') != chain_id or
                row.get('_atom_site.pdbx_PDB_model_num') != model_id):
                continue
            label_seq_str = row.get('_atom_site.label_seq_id', '?')
            if label_seq_str in ('?', '.', ''):
                continue
            try:
                label_seq = int(label_seq_str)
            except (ValueError, TypeError):
                continue

            atom_dict = {
                'name': row.get('_atom_site.label_atom_id', '?').strip(),
                'x': float(row.get('_atom_site.Cartn_x', '0')),
                'y': float(row.get('_atom_site.Cartn_y', '0')),
                'z': float(row.get('_atom_site.Cartn_z', '0')),
                'occupancy': float(row.get('_atom_site.occupancy', '1.0')),
                'altloc': self._normalize_cif_null(row.get('_atom_site.label_alt_id', '').strip()) or 'A',
                'auth_seq_id': row.get('_atom_site.auth_seq_id', '?'),
                'ins_code': self._normalize_cif_null(row.get('_atom_site.pdbx_PDB_ins_code', '').strip()),
            }
            observed_by_label.setdefault(label_seq, []).append(atom_dict)

        scheme = []
        for row in sorted(poly_seq_rows, key=lambda r: int(r['_entity_poly_seq.num'])):
            label_seq = int(row['_entity_poly_seq.num'])
            mon_id = row.get('_entity_poly_seq.mon_id', '?')

            atoms_here = observed_by_label.get(label_seq, [])
            if atoms_here:
                selected = self.select_altloc(atoms_here)
                auth_seq_id = int(atoms_here[0]['auth_seq_id'])
                ins_code = atoms_here[0]['ins_code']
                backbone_complete = self.check_backbone_complete(selected)
                observed = True
            else:
                auth_seq_id = None
                ins_code = ''
                backbone_complete = None
                observed = False

            scheme.append({
                'label_seq_id': label_seq,
                'mon_id': mon_id,
                'observed': observed,
                'auth_seq_id': auth_seq_id,
                'ins_code': ins_code,
                'backbone_complete': backbone_complete,
            })

        return scheme

    def build_imgt_mapping(self, chain_id: str, model_id: str) -> Dict[int, str]:
        """
        Build IMGT-position -> observed-residue-string mapping for a chain.

        SAbDab2-derived structures carry IMGT numbering directly in auth_seq_id
        (+ pdbx_PDB_ins_code for insertion-coded positions, e.g. 111 vs 111A vs
        111B). This project does NOT rerun ANARCI; auth_seq_id/ins_code is the
        source of truth for IMGT position.

        Args:
            chain_id: auth_asym_id of the chain
            model_id: model number to restrict to

        Returns:
            Dict mapping integer IMGT position -> residue string, e.g.
            {105: '105', 112: '112', 113: '112A'} preserving insertion-coded
            residues distinctly (111 != 111A != 111B) via the string suffix,
            while the dict key remains the base integer position for CDR
            range lookups.
        """
        if '_atom_site' not in self.data:
            return {}

        seen = {}
        for row in self.data['_atom_site']:
            if (row.get('_atom_site.auth_asym_id') != chain_id or
                row.get('_atom_site.pdbx_PDB_model_num') != model_id):
                continue

            res_num_str = row.get('_atom_site.auth_seq_id', '?')
            ins_code = self._normalize_cif_null(row.get('_atom_site.pdbx_PDB_ins_code', '').strip())

            try:
                res_num = int(res_num_str)
            except (ValueError, TypeError):
                continue

            residue_str = f"{res_num}{ins_code}"

            # An IMGT position with insertion codes has multiple physical residues
            # (e.g., 111, 111A, 111B). Preserve all distinctly; keep base position
            # mapped to the FIRST (lowest insertion-code) residue found, and track
            # the full set for CDR-completeness checks that need every inserted
            # residue, not just the base one.
            if res_num not in seen:
                seen[res_num] = residue_str
            else:
                # Keep alphabetically-first insertion code as canonical if duplicate
                if residue_str < seen[res_num]:
                    seen[res_num] = residue_str

        return seen

    def build_full_imgt_mapping_with_insertions(self, chain_id: str, model_id: str) -> Dict[int, List[str]]:
        """
        Like build_imgt_mapping but returns ALL residue strings per base IMGT
        position, including every insertion-coded variant (111, 111A, 111B, ...).

        Needed for CDR-completeness checks so an insertion-coded loop residue is
        not silently dropped.
        """
        if '_atom_site' not in self.data:
            return {}

        result: Dict[int, List[str]] = {}
        for row in self.data['_atom_site']:
            if (row.get('_atom_site.auth_asym_id') != chain_id or
                row.get('_atom_site.pdbx_PDB_model_num') != model_id):
                continue

            res_num_str = row.get('_atom_site.auth_seq_id', '?')
            ins_code = self._normalize_cif_null(row.get('_atom_site.pdbx_PDB_ins_code', '').strip())

            try:
                res_num = int(res_num_str)
            except (ValueError, TypeError):
                continue

            residue_str = f"{res_num}{ins_code}"
            result.setdefault(res_num, [])
            if residue_str not in result[res_num]:
                result[res_num].append(residue_str)

        for res_num in result:
            result[res_num].sort()

        return result
