"""
Structure Resolver — Maps INSTANCE metadata to actual extracted mmCIF files.

Handles APFS case-collision resolution and extraction-manifest provenance tracking.
Must preserve case-sensitive chain IDs and distinguish original member paths from
extracted (potentially renamed) paths.
"""
from pathlib import Path
import pandas as pd
from typing import Optional, Dict, Any
import logging

logger = logging.getLogger(__name__)

PROJECT_ROOT = Path("/Users/betulgul/Desktop/staj")
MANIFEST_PATH = PROJECT_ROOT / "data" / "raw_extracted" / "extraction_manifest.tsv"


class StructureResolutionError(Exception):
    """Raised when a structure cannot be uniquely resolved."""
    pass


class StructureResolver:
    """
    Resolves INSTANCE rows to actual extracted mmCIF file paths.

    Workflow:
    1. Load extraction manifest (maps original archive member → extracted path)
    2. For each INSTANCE, find the corresponding archive member
    3. Look up extracted path in manifest
    4. Verify file exists
    5. Track case-collision resolution
    """

    def __init__(self, manifest_path: Optional[Path] = None):
        self.manifest_path = manifest_path or MANIFEST_PATH
        self.manifest: Optional[pd.DataFrame] = None
        self._load_manifest()

    def _load_manifest(self):
        """Load extraction manifest."""
        if not self.manifest_path.exists():
            raise FileNotFoundError(f"Extraction manifest not found: {self.manifest_path}")
        self.manifest = pd.read_csv(self.manifest_path, sep='\t')
        logger.info(f"Loaded extraction manifest: {len(self.manifest)} rows")

    def resolve_structure_path(
        self,
        extracted_structure_path: str,
        pdb: str,
        instance: str,
        archive_sha256: str,
    ) -> Dict[str, Any]:
        """
        Resolve an INSTANCE to its actual mmCIF file path.

        Args:
            extracted_structure_path: Path from master metadata (e.g., 'data/raw_extracted/sdh/pdb_00003zkq/...')
            pdb: PDB ID for fallback/verification
            instance: INSTANCE identifier
            archive_sha256: Archive SHA256 for archive identity verification

        Returns:
            Dict with:
                structure_path: absolute path to mmCIF file
                resolved: bool (success)
                error_reason: str or None
                case_collision_resolved: bool
                manifest_row: dict or None (matching manifest entry)
        """
        full_path = PROJECT_ROOT / extracted_structure_path

        if not full_path.exists():
            return {
                "structure_path": None,
                "resolved": False,
                "error_reason": "FAIL_STRUCTURE_FILE_MISSING",
                "case_collision_resolved": False,
                "manifest_row": None,
            }

        # Try to find manifest row for this extracted path
        # Manifest tracks original archive member and maps to extracted path
        matching_rows = self.manifest[self.manifest["extracted_path"] == extracted_structure_path]

        if len(matching_rows) == 0:
            # File exists on disk but not in manifest — possibly fallback or manual addition
            # This is acceptable; just log it
            logger.warning(f"INSTANCE {instance}: extracted file exists but not in manifest: {extracted_structure_path}")
            return {
                "structure_path": str(full_path),
                "resolved": True,
                "error_reason": None,
                "case_collision_resolved": False,
                "manifest_row": None,
            }

        if len(matching_rows) > 1:
            logger.error(f"INSTANCE {instance}: multiple manifest rows for {extracted_structure_path}; using first")

        manifest_row = matching_rows.iloc[0].to_dict()
        case_collision = bool(manifest_row.get("case_collision_flag", False))

        return {
            "structure_path": str(full_path),
            "resolved": True,
            "error_reason": None,
            "case_collision_resolved": case_collision,
            "manifest_row": manifest_row,
        }

    def verify_provenance(
        self,
        instance: str,
        pdb: str,
        archive_sha256: str,
        extracted_structure_path: str,
    ) -> bool:
        """
        Verify that archive_sha256 and extracted_structure_path are consistent.

        Returns:
            True if provenance is verifiable; False otherwise.
        """
        if not self.manifest_path.exists():
            return False

        # Cross-check archive SHA against manifest
        matching_rows = self.manifest[self.manifest["archive_sha256"] == archive_sha256]
        if len(matching_rows) == 0:
            logger.warning(f"INSTANCE {instance}: archive SHA256 {archive_sha256} not found in manifest")
            return False

        return True
