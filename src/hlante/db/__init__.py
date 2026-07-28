"""
hlante.db
========

Database interfaces used by HLAnte.

Modules
-------
gwas
    GWAS Catalog access (bulk TSV download with local indexed lookup).
pharmgkb
    PharmGKB bulk-dump parsing and querying.
curated
    Built-in curated HLA–disease association table.
imgt
    Local IPD-IMGT/HLA allele/G-group/P-group database management.
afnd
    Allele Frequency Net Database (AFND) lookups for input-quality
    scoring.

Each submodule exposes a client class that follows a common
protocol: update (``update``), query (``query_allele``), and report
local status (``version``).
"""

from __future__ import annotations

import hashlib
import logging
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Protocol

logger: logging.Logger = logging.getLogger(__name__)


class DatabaseError(Exception):
    """
    Base class for database-layer errors.
    """

    pass


class DatabaseIntegrityError(DatabaseError):
    """
    Raised when an installed database file fails checksum verification.
    """

    pass


def sha256_file(path: Path) -> str:
    """
    Return the SHA-256 hex digest of a file, streamed in chunks.
    """
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 16), b""):
            digest.update(chunk)
    return digest.hexdigest()


def atomic_install(tmp: Path, dest: Path, *, keep_backup: bool = True) -> str:
    """
    Atomically install a staged download and verify its integrity.

    Steps: (1) checksum the staged ``tmp`` file; (2) back up any existing
    ``dest`` to ``dest.bak``; (3) atomically replace ``dest`` via
    :meth:`Path.replace` (POSIX-atomic); (4) re-read ``dest`` and confirm its
    SHA-256 matches the staged checksum — detecting any corruption during the
    swap. On mismatch the previous file is restored from the ``.bak`` and a
    :class:`DatabaseIntegrityError` is raised.

    Returns
    -------
    str
        The verified SHA-256 hex digest of the installed file.
    """
    tmp_path, dest_path = Path(tmp), Path(dest)
    expected = sha256_file(tmp_path)
    backup = dest_path.with_suffix(dest_path.suffix + ".bak")
    had_existing = dest_path.exists()
    if had_existing and keep_backup:
        shutil.copy2(dest_path, backup)
    tmp_path.replace(dest_path)
    actual = sha256_file(dest_path)
    if actual != expected:
        if had_existing and keep_backup and backup.exists():
            backup.replace(dest_path)  # Roll back
        raise DatabaseIntegrityError(
            f"Checksum mismatch installing {dest_path}: "
            f"staged {expected[:12]}… != installed {actual[:12]}… "
            "(previous copy restored from .bak)"
        )
    return expected


class DatabaseDownloadError(DatabaseError):
    """
    Raised when a remote resource cannot be downloaded or refreshed.
    """

    pass


class DatabaseQueryError(DatabaseError):
    """
    Raised for errors encountered while executing a query.
    """

    pass


@dataclass(frozen=True)
class DBRecord:
    """
    Common raw record returned by every database module.

    Attributes
    ----------
    allele : str
        HLA allele in IPD-IMGT/HLA canonical form.
    phenotype : str
        Phenotype / disease / drug reaction text.
    source_db : str
        Originating database name.
    source_id : str
        Identifier within the source database.
    pmid : list of str
        Related PubMed IDs.
    attributes : dict
        Source-specific raw fields.
    """

    allele: str
    phenotype: str
    source_db: str
    source_id: str
    pmid: List[str]
    attributes: Dict[str, Any]


class BaseDB(Protocol):
    """
    Protocol implemented by every database interface.
    """

    name: str

    def update(self, target_dir: Path) -> Path:
        """
        Download or refresh the local copy.

        Parameters
        ----------
        target_dir : Path
            Directory in which to store the database files.

        Returns
        -------
        Path
            Path to the primary artefact written / refreshed.
        """
        ...

    def version(self) -> Optional[str]:
        """
        Return the version of the local copy.

        Returns
        -------
        str or None
            ``None`` when unknown.
        """
        ...

    def query_allele(self, allele: str) -> List[DBRecord]:
        """
        Return all records associated with an allele.

        Parameters
        ----------
        allele : str
            Allele in IPD-IMGT/HLA canonical form.

        Returns
        -------
        list of DBRecord
            Associated records.
        """
        ...


__all__ = [
    "BaseDB",
    "DBRecord",
    "DatabaseError",
    "DatabaseDownloadError",
    "DatabaseQueryError",
]
