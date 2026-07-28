"""
hlante.db.imgt
=============

IPD-IMGT/HLA local database management.

Downloads the relevant files from the IPD-IMGT/HLA GitHub mirror
(``ANHIG/IMGTHLA``) and serves them locally. This is the authoritative
source for allele-name validation, G-group / P-group mappings, and
resolution conversions within HLAnte.

Default download directory: ``~/.hlante/imgt_hla/``.
"""

from __future__ import annotations

import json
import logging
import re
import urllib.error
import urllib.request
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Set

from hlante.db import (
    DBRecord,
    DatabaseDownloadError,
    DatabaseQueryError,
    atomic_install,
    sha256_file,
)

logger: logging.Logger = logging.getLogger(__name__)


#: Default Git ref for the ANHIG/IMGTHLA mirror. ``"Latest"`` tracks the
#: moving release branch; pass an explicit release tag to
#: :func:`download_imgt_db` (``ref=...``) / ``db-update --ref`` to pin a
#: reproducible snapshot. The chosen ref is recorded in ``version.json``
#: alongside per-file SHA-256 checksums.
IMGT_DEFAULT_REF: str = "Latest"
IMGT_GITHUB_RAW_TEMPLATE: str = "https://raw.githubusercontent.com/ANHIG/IMGTHLA/{ref}"


IMGT_GITHUB_RAW_BASE: str = IMGT_GITHUB_RAW_TEMPLATE.format(ref=IMGT_DEFAULT_REF)
IMGT_ALLELE_LIST_URL: str = f"{IMGT_GITHUB_RAW_BASE}/Allelelist.txt"
IMGT_G_GROUP_URL: str = f"{IMGT_GITHUB_RAW_BASE}/wmda/hla_nom_g.txt"
IMGT_P_GROUP_URL: str = f"{IMGT_GITHUB_RAW_BASE}/wmda/hla_nom_p.txt"

DEFAULT_LOCAL_DIR: Path = Path.home() / ".hlante" / "imgt_hla"

ALLELE_LIST_FILENAME: str = "Allelelist.txt"
G_GROUP_FILENAME: str = "hla_nom_g.txt"
P_GROUP_FILENAME: str = "hla_nom_p.txt"
VERSION_FILENAME: str = "version.json"


def normalize_imgt_ref(ref: str) -> str:
    """
    Map a release number onto the branch name used by the ANHIG mirror.

    The mirror names its release branches without separators (``3640`` for
    IPD-IMGT/HLA 3.64.0), so the human-readable form that appears in the
    literature and in this project's documentation — ``3.64.0`` — has to be
    translated before it can be fetched. Anything that is not a dotted release
    number (``Latest``, a commit SHA, an explicit branch or tag) is passed
    through untouched.
    """
    candidate = ref.strip()
    if re.fullmatch(r"\d+\.\d+\.\d+", candidate):
        return candidate.replace(".", "")
    return candidate


class IMGTDatabaseError(DatabaseQueryError):
    """
    Raised for errors against the IPD-IMGT/HLA local database.
    """

    pass


@dataclass
class IMGTAllele:
    """
    Domain model for an IPD-IMGT/HLA allele.

    Attributes
    ----------
    accession : str
        IMGT accession (e.g., ``"HLA00001"``).
    name : str
        Canonical allele name (e.g., ``"A*01:01:01:01"``).
    g_group : str, optional
        G-group the allele belongs to.
    p_group : str, optional
        P-group the allele belongs to.
    """

    accession: str
    name: str
    g_group: Optional[str] = None
    p_group: Optional[str] = None


# ---------------------------------------------------------------------------
# File download
# ---------------------------------------------------------------------------


def _http_fetch(url: str, dest: Path, timeout: float = 60.0) -> Path:
    """
    Download a URL to a destination file.

    Parameters
    ----------
    url : str
        Source URL.
    dest : Path
        Target file path.
    timeout : float, optional
        Download timeout in seconds.

    Returns
    -------
    Path
        The written file path.

    Raises
    ------
    DatabaseDownloadError
        On HTTP or network error.
    """
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    try:
        logger.info("Downloading IMGT file: %s → %s", url, dest)
        with urllib.request.urlopen(url, timeout=timeout) as response:  # noqa: S310
            with tmp.open("wb") as out:
                while True:
                    chunk = response.read(1 << 15)
                    if not chunk:
                        break
                    out.write(chunk)
    except (urllib.error.URLError, TimeoutError, OSError) as exc:
        if tmp.exists():
            tmp.unlink(missing_ok=True)
        raise DatabaseDownloadError(f"IMGT file download failed ({url}): {exc}") from exc
    # Atomic swap with a .bak rollback and post-replace checksum verification
    # rather than a bare rename.
    atomic_install(tmp, dest)
    return dest


def _parse_allelelist_version(path: Path) -> Optional[str]:
    """
    Read the version string from the ``# version:`` header line of
    ``Allelelist.txt``.

    Parameters
    ----------
    path : Path
        Path to Allelelist.txt.

    Returns
    -------
    str or None
        Version string when present.
    """
    version_re = re.compile(r"#\s*version:\s*(.+)", re.IGNORECASE)
    try:
        with path.open("r", encoding="utf-8") as handle:
            for _ in range(20):
                line = handle.readline()
                if not line:
                    break
                match = version_re.search(line)
                if match:
                    return match.group(1).strip()
    except OSError:
        return None
    return None


def download_imgt_db(
    target_dir: Optional[Path] = None,
    *,
    include_groups: bool = True,
    force: bool = False,
    ref: str = IMGT_DEFAULT_REF,
) -> Path:
    """
    Download the IPD-IMGT/HLA source files from GitHub.

    Files downloaded:

    - ``Allelelist.txt`` — official allele list (name → accession).
    - ``wmda/hla_nom_g.txt`` — G-group definitions (optional).
    - ``wmda/hla_nom_p.txt`` — P-group definitions (optional).
    - ``version.json`` — version and download-timestamp metadata.

    Parameters
    ----------
    target_dir : Path, optional
        Download directory. Defaults to :data:`DEFAULT_LOCAL_DIR`.
    include_groups : bool, optional
        Download the G/P-group files too (default ``True``).
    force : bool, optional
        Overwrite existing files if ``True``.

    Returns
    -------
    Path
        The local database root directory.

    Raises
    ------
    DatabaseDownloadError
        When any download step fails.
    """
    root = Path(target_dir) if target_dir is not None else DEFAULT_LOCAL_DIR
    root.mkdir(parents=True, exist_ok=True)

    ref = normalize_imgt_ref(ref)
    base = IMGT_GITHUB_RAW_TEMPLATE.format(ref=ref)
    allele_url = f"{base}/Allelelist.txt"
    g_url = f"{base}/wmda/hla_nom_g.txt"
    p_url = f"{base}/wmda/hla_nom_p.txt"

    allele_dest = root / ALLELE_LIST_FILENAME
    if force or not allele_dest.exists():
        _http_fetch(allele_url, allele_dest)
    else:
        logger.info("Allelelist.txt is already present, skipping: %s", allele_dest)

    if include_groups:
        g_dest = root / G_GROUP_FILENAME
        p_dest = root / P_GROUP_FILENAME
        for url, dest in ((g_url, g_dest), (p_url, p_dest)):
            if force or not dest.exists():
                _http_fetch(url, dest)
            else:
                logger.info("%s is already present, skipping: %s", dest.name, dest)

    # Record per-file SHA-256 so the exact snapshot can be verified /
    # Reproduced, and the moving "Latest" ref cannot silently change the data
    # Underneath a cached install without the checksum revealing it.
    checksums: Dict[str, str] = {
        name: sha256_file(root / name)
        for name in (ALLELE_LIST_FILENAME, G_GROUP_FILENAME, P_GROUP_FILENAME)
        if (root / name).is_file()
    }
    version = _parse_allelelist_version(allele_dest)
    meta = {
        "version": version,
        "downloaded_at": datetime.now(timezone.utc).isoformat(),
        "source_base": base,
        "ref": ref,
        "sha256": checksums,
    }
    (root / VERSION_FILENAME).write_text(
        json.dumps(meta, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    logger.info("IPD-IMGT/HLA download complete (version=%s) → %s", version, root)
    return root


# ---------------------------------------------------------------------------
# Parse helpers
# ---------------------------------------------------------------------------


def parse_allelelist(path: Path) -> Dict[str, str]:
    """
    Parse ``Allelelist.txt`` into an ``{allele_name: accession}`` dict.

    Parameters
    ----------
    path : Path
        Path to Allelelist.txt.

    Returns
    -------
    dict of str to str
        Mapping from allele name to IMGT accession.

    Raises
    ------
    IMGTDatabaseError
        When the file cannot be read or its format is unrecognised.
    """
    if not path.is_file():
        raise IMGTDatabaseError(f"Allelelist.txt not found: {path}")

    alleles: Dict[str, str] = {}
    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            # Skip a header row like "AlleleID,Allele"
            if line.lower().startswith("alleleid"):
                continue
            parts = [p.strip() for p in line.split(",")]
            if len(parts) < 2:
                parts = [p.strip() for p in line.split()]
            if len(parts) < 2:
                continue
            accession, allele_name = parts[0], parts[1]
            if not accession or not allele_name:
                continue
            alleles[allele_name] = accession
    if not alleles:
        raise IMGTDatabaseError(f"No allele rows found in Allelelist.txt: {path}")
    return alleles


def parse_group_file(path: Path) -> Dict[str, List[str]]:
    """
    Parse a WMDA group file (``hla_nom_g.txt`` / ``hla_nom_p.txt``).

    Actual row format: ``<gene>*;<member1>/<member2>/...;<group_name_suffix>``
    Example: ``A*;01:01:01:01/01:01:01:02N/...;01:01:01G``

    Parameters
    ----------
    path : Path
        Path to the group file.

    Returns
    -------
    dict of str to list of str
        Group name (e.g., ``"A*01:01:01G"``) → list of member alleles.
    """
    groups: Dict[str, List[str]] = {}
    if not path.is_file():
        return groups

    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split(";")
            if len(parts) < 3:
                continue
            gene_prefix = parts[0].rstrip("*")
            members_field = parts[1].strip()  # Slash-separated allele members
            group_fields = parts[2].strip()  # G/P-group name suffix (e.g. "01:01:01G")
            if not group_fields:
                continue
            members = [f"{gene_prefix}*{m.strip()}" for m in members_field.split("/") if m.strip()]
            group_name = f"{gene_prefix}*{group_fields}"
            groups[group_name] = members
    return groups


# ---------------------------------------------------------------------------
# IMGTDB — high-level class interface (BaseDB-compatible)
# ---------------------------------------------------------------------------


class IMGTDB:
    """
    IPD-IMGT/HLA local-database client (BaseDB-compatible).

    Parameters
    ----------
    local_dir : Path, optional
        Directory containing Allelelist + group files.
    """

    name: str = "imgt"

    def __init__(self, local_dir: Optional[Path] = None) -> None:
        self.local_dir: Path = Path(local_dir) if local_dir else DEFAULT_LOCAL_DIR
        self._alleles: Dict[str, str] = {}
        self._g_groups: Dict[str, List[str]] = {}
        self._p_groups: Dict[str, List[str]] = {}
        self._loaded: bool = False

    def update(self, target_dir: Optional[Path] = None) -> Path:
        """
        Download IPD-IMGT/HLA source files from GitHub.

        Parameters
        ----------
        target_dir : Path, optional
            Write directory; defaults to ``self.local_dir``.

        Returns
        -------
        Path
            The local database root.

        Raises
        ------
        DatabaseDownloadError
            If the download fails.
        """
        dest = Path(target_dir) if target_dir else self.local_dir
        root = download_imgt_db(dest)
        self.local_dir = root
        self._loaded = False
        return root

    def version(self) -> Optional[str]:
        """
        Return the local IPD-IMGT/HLA version.

        Returns
        -------
        str or None
            Version string when known.
        """
        meta_path = self.local_dir / VERSION_FILENAME
        if meta_path.is_file():
            try:
                meta = json.loads(meta_path.read_text(encoding="utf-8"))
            except json.JSONDecodeError:
                meta = {}
            if meta.get("version"):
                return str(meta["version"])
        allele_path = self.local_dir / ALLELE_LIST_FILENAME
        return _parse_allelelist_version(allele_path) if allele_path.is_file() else None

    def _ensure_loaded(self) -> None:
        if self._loaded:
            return
        self._alleles = parse_allelelist(self.local_dir / ALLELE_LIST_FILENAME)
        self._g_groups = parse_group_file(self.local_dir / G_GROUP_FILENAME)
        self._p_groups = parse_group_file(self.local_dir / P_GROUP_FILENAME)
        self._loaded = True

    def query_allele(self, allele: str) -> List[DBRecord]:
        """
        Check whether an allele is known and return a structured record.

        Parameters
        ----------
        allele : str
            Allele name.

        Returns
        -------
        list of DBRecord
            One-element list when known, empty otherwise.
        """
        self._ensure_loaded()
        accession = self._alleles.get(allele)
        if accession is None:
            return []
        return [
            DBRecord(
                allele=allele,
                phenotype="",
                source_db=self.name,
                source_id=accession,
                pmid=[],
                attributes={"accession": accession},
            )
        ]

    def is_known(self, allele: str) -> bool:
        """
        Return whether the allele is present in IPD-IMGT/HLA.
        """
        self._ensure_loaded()
        return allele in self._alleles

    def g_group_of(self, allele: str) -> Optional[str]:
        """
        Return the G-group of a given allele, if any.
        """
        self._ensure_loaded()
        for group_name, members in self._g_groups.items():
            if allele in members or allele == group_name:
                return group_name
        return None

    def p_group_of(self, allele: str) -> Optional[str]:
        """
        Return the P-group of a given allele, if any.
        """
        self._ensure_loaded()
        for group_name, members in self._p_groups.items():
            if allele in members or allele == group_name:
                return group_name
        return None

    def members_of_g_group(self, g_group: str) -> Set[str]:
        """
        Return all alleles that belong to a G-group.
        """
        self._ensure_loaded()
        return set(self._g_groups.get(g_group, []))

    def members_of_p_group(self, p_group: str) -> Set[str]:
        """
        Return all alleles that belong to a P-group.
        """
        self._ensure_loaded()
        return set(self._p_groups.get(p_group, []))


__all__ = [
    "normalize_imgt_ref",
    "IMGTDB",
    "IMGTAllele",
    "IMGTDatabaseError",
    "IMGT_ALLELE_LIST_URL",
    "IMGT_G_GROUP_URL",
    "IMGT_P_GROUP_URL",
    "IMGT_GITHUB_RAW_BASE",
    "DEFAULT_LOCAL_DIR",
    "ALLELE_LIST_FILENAME",
    "G_GROUP_FILENAME",
    "P_GROUP_FILENAME",
    "VERSION_FILENAME",
    "download_imgt_db",
    "parse_allelelist",
    "parse_group_file",
]
