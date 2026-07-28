"""
hlante.normalizer
================

Allele normalization against IPD-IMGT/HLA.

This module maps allele calls in :class:`hlante.parser.HLAGenotype` to
their IPD-IMGT/HLA records and produces a :class:`NormalizedAllele`
for each call. Missing resolution (ambiguity), G-group / P-group
notation, and novel alleles are flagged as metadata on the output.

Public API
----------
- :func:`load_imgt_db` — load a local IPD-IMGT/HLA copy.
- :func:`normalize_allele` — normalize a single allele.
- :func:`resolve_ambiguity` — expand a truncated or group allele into
  the list of possible full alleles.
- :func:`batch_normalize` — parallel normalization of an iterable of
  :class:`HLAGenotype` records.
"""

from __future__ import annotations

import concurrent.futures
import json
import logging
import re
from dataclasses import dataclass
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

from hlante.db.imgt import (
    ALLELE_LIST_FILENAME,
    DEFAULT_LOCAL_DIR,
    G_GROUP_FILENAME,
    P_GROUP_FILENAME,
    VERSION_FILENAME,
    IMGTDatabaseError,
    parse_allelelist,
    parse_group_file,
)
from hlante.parser import HLAGenotype

logger: logging.Logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

#: Age threshold (days) after which a warning is emitted about a stale
#: local copy (~6 months).
STALE_THRESHOLD_DAYS: int = 180

#: Tokens representing absent or undetermined allele calls.
NULL_ALLELE_TOKENS: frozenset[str] = frozenset(
    {"", "*", "-", ".", "N/A", "NA", "NONE", "NOT TYPED", "NOTTYPED"}
)

#: HLA Class II gene prefixes.
CLASS_II_GENE_PREFIXES: Tuple[str, ...] = ("DR", "DQ", "DP", "DM", "DO", "DN")

_HLA_PREFIX: str = "HLA-"

_ALLELE_SUFFIX_RE: re.Pattern[str] = re.compile(r"([GP])$")
_NOMENCLATURE_SUFFIX_RE: re.Pattern[str] = re.compile(r"[NLSQCA]$")
_ALLELE_BASIC_RE: re.Pattern[str] = re.compile(r"^[A-Z]+\d*\*\d{2,3}(:\d{2,3})*[A-Z]?$")


# ---------------------------------------------------------------------------
# Exceptions
# ---------------------------------------------------------------------------


class HLANormalizerError(Exception):
    """
    Base class for errors raised during normalization.
    """

    pass


class IMGTDatabaseMissingError(HLANormalizerError):
    """
    Raised when the local IPD-IMGT/HLA copy cannot be found.
    """

    pass


class InvalidAlleleError(HLANormalizerError):
    """
    Raised for allele expressions that do not conform to IMGT nomenclature.
    """

    pass


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------


@dataclass
class NormalizedAllele:
    """
    Normalized allele record.

    Attributes
    ----------
    allele_name : str
        Input allele (normalized, without ``HLA-`` prefix).
    imgt_accession : str or None
        IPD-IMGT/HLA accession (e.g., ``"HLA00001"``). ``None`` if no
        exact match.
    protein_group : str or None
        G-group corresponding to the allele (e.g., ``"A*01:01:01G"``).
        ``None`` when not available.
    hla_class : str
        ``"I"`` or ``"II"``.
    gene : str
        Gene name with the ``HLA-`` prefix (e.g., ``"HLA-A"``,
        ``"HLA-DRB1"``).
    resolution_level : int
        Number of colon-separated fields (1-4).
    is_ambiguous : bool
        Ambiguous call (low resolution, G/P group, or multiple IMGT
        matches).
    is_novel : bool
        No exact or prefix match in IPD-IMGT/HLA.
    sample_id : str or None
        Source sample ID (populated by :func:`batch_normalize`).
    source_tool : str or None
        Tool that produced the call.
    source_locus : str or None
        Source locus name (e.g., ``"HLA-A"``).
    source_resolution : str or None
        Resolution label as reported by the parser.
    allele_index : int or None
        Index within the genotype — ``0`` for the first allele, ``1``
        for the second.
    expression_suffix : str or None
        IPD-IMGT/HLA expression suffix when present: ``N`` (null —
        antigen *not* expressed on the cell surface), ``L`` (low
        surface expression), ``S`` (soluble / secreted only), ``C``
        (cytoplasm only), ``A`` (aberrant expression), ``Q``
        (questionable expression). ``None`` for a normally expressed
        allele. Surface-expression-dependent HLA disease/drug
        associations do not apply to ``N`` alleles (see
        :func:`hlante.annotator.annotate_genotype`).
    """

    allele_name: str
    imgt_accession: Optional[str]
    protein_group: Optional[str]
    hla_class: str
    gene: str
    resolution_level: int
    is_ambiguous: bool
    is_novel: bool = False
    sample_id: Optional[str] = None
    source_tool: Optional[str] = None
    source_locus: Optional[str] = None
    source_resolution: Optional[str] = None
    allele_index: Optional[int] = None
    expression_suffix: Optional[str] = None
    imgt_match_category: str = ""
    imgt_match_candidates: int = 0
    #: Per-allele quality reported by the typing tool, when it provides one.
    #: Only T1K's native layout carries this; see :class:`hlante.parser.HLAGenotype`.
    caller_quality: Optional[float] = None

    @property
    def is_null(self) -> bool:
        """
        Whether this is a *null* allele (IPD-IMGT/HLA ``N`` suffix).

        A null allele carries a nucleotide change that abolishes
        surface expression of the antigen; HLA peptide-presentation —
        and therefore the disease- and drug-hypersensitivity
        associations that depend on it — does not apply.
        """
        return self.expression_suffix == "N"

    @property
    def is_low_or_aberrant_expression(self) -> bool:
        """
        Whether this allele has reduced or aberrant cell-surface
        expression (IPD-IMGT/HLA ``L`` / ``S`` / ``C`` / ``A`` / ``Q``
        suffix). Such alleles are annotated but flagged, because the
        clinical relevance of their (reduced) expression is uncertain.
        """
        return self.expression_suffix in {"L", "S", "C", "A", "Q"}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _strip_hla_prefix(allele: str) -> str:
    """
    Strip the leading ``HLA-`` prefix from an allele name.
    """
    if allele.upper().startswith(_HLA_PREFIX):
        return allele[len(_HLA_PREFIX) :]
    return allele


def _is_null_token(value: Optional[str]) -> bool:
    """
    Return whether the value is a null / untyped allele token.
    """
    if value is None:
        return True
    return value.strip().upper() in NULL_ALLELE_TOKENS


def _gene_of(allele_clean: str) -> str:
    """
    Return the ``HLA-<GENE>`` form of the allele's gene.
    """
    if "*" in allele_clean:
        gene = allele_clean.split("*", 1)[0]
    else:
        gene = allele_clean
    return f"{_HLA_PREFIX}{gene.upper()}"


def _hla_class_of(gene: str) -> str:
    """
    Return the HLA class (``"I"`` / ``"II"``) for a gene name.
    """
    bare = _strip_hla_prefix(gene).upper()
    if bare.startswith(CLASS_II_GENE_PREFIXES):
        return "II"
    return "I"


def _resolution_of(allele_clean: str) -> int:
    """
    Return the field-level resolution of an HLA allele.

    Resolution is the number of colon-separated fields after the
    asterisk, per IPD-IMGT/HLA nomenclature (Marsh et al. 2010).
    Current nomenclature admits at most four fields, so the returned
    integer is in 1..4. Earlier releases reported this on a digit
    scale (2 / 4 / 6 / 8), which conflated fields with digits.

    Examples
    --------
    - ``A*02``              → 1  (one field)
    - ``A*02:01``           → 2  (two fields)
    - ``A*02:01:01``        → 3  (three fields)
    - ``A*02:01:01:01``     → 4  (four fields)
    - ``DPB1*104:01:01``    → 3  (three fields, seven digits)
    - ``B*57:01G`` / ``B*57:01P`` → 2  (G/P suffix stripped)

    The G/P-group suffix (``G``/``P``) and nomenclature suffix
    letters (``N``/``L``/``S``/``Q``/``C``/``A``) are stripped before
    counting.
    """
    if "*" in allele_clean:
        _, rest = allele_clean.split("*", 1)
    else:
        rest = allele_clean
    # Strip G/P group suffix first, then any nomenclature suffix.
    rest = _ALLELE_SUFFIX_RE.sub("", rest)
    rest = _NOMENCLATURE_SUFFIX_RE.sub("", rest)
    if not rest:
        return 1
    return min(len(rest.split(":")), 4)


def _is_group_allele(allele_clean: str) -> Optional[str]:
    """
    Return ``"G"`` / ``"P"`` when the allele ends with that suffix; else
    ``None``.
    """
    match = _ALLELE_SUFFIX_RE.search(allele_clean)
    if match:
        return match.group(1)
    return None


def _expression_suffix_of(allele_clean: str) -> Optional[str]:
    """
    Return the IPD-IMGT/HLA expression suffix (``N``/``L``/``S``/``C``/
    ``A``/``Q``) of an allele, or ``None``.

    The expression suffix is the trailing letter of the allele name and
    encodes the surface-expression status of the protein (Marsh et al.
    2010): ``N`` null (not expressed), ``L`` low, ``S`` secreted/soluble,
    ``C`` cytoplasm-only, ``A`` aberrant, ``Q`` questionable. The
    G/P-group designators (``G``/``P``) are *not* expression suffixes and
    are excluded.
    """
    if _is_group_allele(allele_clean) is not None:
        return None
    match = _NOMENCLATURE_SUFFIX_RE.search(allele_clean)
    return match.group() if match else None


def _is_not_present_placeholder(allele_clean: str) -> bool:
    """
    Detect the all-zero "gene not present / not typed" placeholder some
    typing tools (notably HLA-HD) emit for the present-or-absent loci
    ``HLA-DRB3/4/5`` — e.g. ``DRB3*00:00`` or ``DRB3*00:00N``.

    Allele field ``00`` does not denote a real IPD-IMGT/HLA allele
    (numbering begins at ``01``); an all-``00`` field set therefore means
    "no allele at this locus" and must be treated as a null/absent call
    rather than annotated as a novel allele.
    """
    if "*" not in allele_clean:
        return False
    _, rest = allele_clean.split("*", 1)
    rest = _ALLELE_SUFFIX_RE.sub("", rest)  # Strip any G/P group suffix
    rest = _NOMENCLATURE_SUFFIX_RE.sub("", rest)  # Strip any expression suffix
    fields = [f for f in rest.split(":") if f]
    return bool(fields) and all(set(f) == {"0"} for f in fields)


# ---------------------------------------------------------------------------
# IPD-IMGT/HLA loader
# ---------------------------------------------------------------------------


def load_imgt_db(db_path: Optional[Path] = None) -> Dict[str, Any]:
    """
    Load the local IPD-IMGT/HLA copy.

    When ``db_path`` is a directory, :data:`Allelelist.txt` and
    (optionally) :data:`hla_nom_g.txt` / :data:`hla_nom_p.txt` are
    expected inside it. When it is a file, the file is treated as
    ``Allelelist.txt``.

    Returned dictionary
    -------------------
    ``{"alleles", "g_groups", "p_groups", "allele_to_g", "allele_to_p",
    "version", "downloaded_at", "path", "is_stale"}``

    Parameters
    ----------
    db_path : Path, optional
        Local IPD-IMGT/HLA directory. Falls back to
        :data:`hlante.db.imgt.DEFAULT_LOCAL_DIR`.

    Returns
    -------
    dict
        Parsed database content.

    Raises
    ------
    IMGTDatabaseMissingError
        When required files are missing.
    """
    root = Path(db_path) if db_path is not None else DEFAULT_LOCAL_DIR

    if root.is_file():
        allele_path = root
        root_dir = root.parent
    else:
        root_dir = root
        allele_path = root_dir / ALLELE_LIST_FILENAME

    if not allele_path.is_file():
        raise IMGTDatabaseMissingError(
            f"IPD-IMGT/HLA Allelelist not found: {allele_path}.\n"
            f"The IPD-IMGT/HLA release is a required, one-time download "
            f"(~10 MB). Install it with:\n"
            f"    hlante db-update --db imgt\n"
            f"To reproduce a specific release, pin it with "
            f"`--imgt-ref` (for example: "
            f"`hlante db-update --db imgt --imgt-ref 3.64.0`)."
        )

    try:
        alleles = parse_allelelist(allele_path)
    except IMGTDatabaseError as exc:  # pragma: no cover - defensive
        raise IMGTDatabaseMissingError(str(exc)) from exc

    g_groups = parse_group_file(root_dir / G_GROUP_FILENAME)
    p_groups = parse_group_file(root_dir / P_GROUP_FILENAME)

    allele_to_g: Dict[str, str] = {}
    for group, members in g_groups.items():
        for member in members:
            allele_to_g[member] = group

    allele_to_p: Dict[str, str] = {}
    for group, members in p_groups.items():
        for member in members:
            allele_to_p[member] = group

    version, downloaded_at = _read_version_meta(root_dir, allele_path)
    is_stale = _is_stale(downloaded_at)
    if is_stale and downloaded_at is not None:
        age_days = (datetime.now(timezone.utc) - downloaded_at).days
        logger.warning(
            "Local IPD-IMGT/HLA copy is %d day(s) old (> %d-day threshold). "
            "Run `download_imgt_db(force=True)` to refresh.",
            age_days,
            STALE_THRESHOLD_DAYS,
        )

    logger.info(
        "IPD-IMGT/HLA loaded: %d allele(s), %d G-group(s), %d P-group(s) (version=%s)",
        len(alleles),
        len(g_groups),
        len(p_groups),
        version,
    )

    return {
        "alleles": alleles,
        "g_groups": g_groups,
        "p_groups": p_groups,
        "allele_to_g": allele_to_g,
        "allele_to_p": allele_to_p,
        "version": version,
        "downloaded_at": downloaded_at,
        "path": root_dir,
        "is_stale": is_stale,
    }


def _read_version_meta(root: Path, allele_path: Path) -> Tuple[Optional[str], Optional[datetime]]:
    """
    Read the version and download timestamp from ``version.json`` (or
    fall back to the Allelelist header + mtime).

    Returns
    -------
    tuple of (str or None, datetime or None)
    """
    meta_path = root / VERSION_FILENAME
    version: Optional[str] = None
    downloaded_at: Optional[datetime] = None

    if meta_path.is_file():
        try:
            meta = json.loads(meta_path.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            meta = {}
        version = meta.get("version")
        dl = meta.get("downloaded_at")
        if isinstance(dl, str):
            try:
                downloaded_at = datetime.fromisoformat(dl)
                if downloaded_at.tzinfo is None:
                    downloaded_at = downloaded_at.replace(tzinfo=timezone.utc)
            except ValueError:
                downloaded_at = None

    if version is None:
        from hlante.db.imgt import _parse_allelelist_version  # Local import

        version = _parse_allelelist_version(allele_path)

    if downloaded_at is None:
        try:
            downloaded_at = datetime.fromtimestamp(allele_path.stat().st_mtime, tz=timezone.utc)
        except OSError:
            downloaded_at = None

    return version, downloaded_at


def _is_stale(downloaded_at: Optional[datetime]) -> bool:
    """
    Return whether ``downloaded_at`` is older than
    :data:`STALE_THRESHOLD_DAYS`.
    """
    if downloaded_at is None:
        return False
    threshold = timedelta(days=STALE_THRESHOLD_DAYS)
    now = datetime.now(timezone.utc)
    if downloaded_at.tzinfo is None:
        downloaded_at = downloaded_at.replace(tzinfo=timezone.utc)
    return now - downloaded_at > threshold


# ---------------------------------------------------------------------------
# Single-allele normalization
# ---------------------------------------------------------------------------


#: How a submitted allele name resolves against the loaded IPD-IMGT/HLA
#: release. The categories are mutually exclusive and exhaustive:
#:
#: - ``exact`` — the name is listed verbatim in ``Allelelist.txt``.
#: - ``prefix_unique`` — not listed, but exactly one listed allele extends it.
#:   This is a property of the current release, not evidence that the call
#:   identifies a single allele in general.
#: - ``prefix_multiple`` — several listed alleles extend the name; the call
#:   denotes a set, and ``imgt_match_candidates`` gives its size.
#: - ``g_group`` / ``p_group`` — a G- or P-suffixed name found in the
#:   corresponding WHO nomenclature file; the count is the group's membership.
#: - ``unmatched`` — nothing in the release matches, including a G/P-suffixed
#:   name whose group is absent. Reported as putatively novel.
MATCH_CATEGORIES: Tuple[str, ...] = (
    "exact",
    "prefix_unique",
    "prefix_multiple",
    "g_group",
    "p_group",
    "unmatched",
)


def classify_imgt_match(
    clean: str,
    alleles_map: Dict[str, str],
    g_groups: Dict[str, List[str]],
    p_groups: Dict[str, List[str]],
) -> Tuple[str, int, List[str]]:
    """
    Classify how *clean* resolves against the loaded release.

    Returns the category, the number of IPD-IMGT/HLA alleles the name
    denotes, and the matching allele names for a prefix match (empty for
    the other categories, which do not need the list).

    Reporting these separately answers a question the aggregate
    "normalisation rate" cannot: a rate of 100% is compatible with almost
    every call denoting a large set of alleles rather than one.
    """
    suffix = _is_group_allele(clean)
    if suffix == "G":
        members = g_groups.get(clean) or []
        return ("g_group", len(members), []) if members else ("unmatched", 0, [])
    if suffix == "P":
        members = p_groups.get(clean) or []
        return ("p_group", len(members), []) if members else ("unmatched", 0, [])

    if clean in alleles_map:
        return ("exact", 1, [])

    prefix = clean + ":"
    matches = [a for a in alleles_map if a.startswith(prefix)]
    if not matches:
        return ("unmatched", 0, [])
    return ("prefix_unique" if len(matches) == 1 else "prefix_multiple", len(matches), matches)


def normalize_allele(
    allele: Optional[str],
    imgt_db: Dict[str, Any],
) -> Optional[NormalizedAllele]:
    """
    Turn a single allele expression into a :class:`NormalizedAllele`.

    Behaviour
    ---------
    - Null tokens (``"*"``, ``"-"``, ``"Not typed"``, etc.) return
      ``None`` and emit a debug log entry.
    - G/P-group inputs (e.g., ``A*01:01:01G``) populate
      ``protein_group`` and are marked as ambiguous.
    - Exact IPD-IMGT/HLA match populates ``imgt_accession``.
    - Prefix matches (no exact match) yield ``is_ambiguous=True``;
      otherwise the allele is flagged ``is_novel=True``.
    - Resolution below 4 digits (``A*02``) is always ambiguous.

    Parameters
    ----------
    allele : str or None
        Allele to normalize. ``HLA-`` prefix is stripped automatically.
    imgt_db : dict
        Dictionary produced by :func:`load_imgt_db`.

    Returns
    -------
    NormalizedAllele or None
        ``None`` for null tokens; a populated record otherwise.

    Raises
    ------
    InvalidAlleleError
        When the allele does not conform to IMGT nomenclature.
    """
    if _is_null_token(allele):
        logger.debug("Null allele token skipped: %r", allele)
        return None

    assert allele is not None  # For mypy/ruff
    clean = _strip_hla_prefix(allele.strip())
    if not _ALLELE_BASIC_RE.match(clean):
        raise InvalidAlleleError(f"Does not conform to IMGT nomenclature: {allele!r}")

    # HLA-HD-style "gene not present" placeholder (all-zero fields, e.g.
    # DRB3*00:00 / DRB3*00:00N) is an absent call, not a novel allele.
    if _is_not_present_placeholder(clean):
        logger.debug("Gene-not-present placeholder skipped: %r", allele)
        return None

    expression_suffix = _expression_suffix_of(clean)

    alleles_map: Dict[str, str] = imgt_db["alleles"]
    g_groups: Dict[str, List[str]] = imgt_db["g_groups"]
    p_groups: Dict[str, List[str]] = imgt_db["p_groups"]
    allele_to_g: Dict[str, str] = imgt_db["allele_to_g"]

    gene = _gene_of(clean)
    hla_class = _hla_class_of(gene)
    resolution_level = _resolution_of(clean)
    group_suffix = _is_group_allele(clean)

    imgt_accession: Optional[str] = None
    protein_group: Optional[str] = None
    is_ambiguous: bool = False
    is_novel: bool = False

    category, candidates, prefix_matches = classify_imgt_match(
        clean, alleles_map, g_groups, p_groups
    )

    if group_suffix == "G":
        members = g_groups.get(clean)
        if members:
            imgt_accession = alleles_map.get(members[0])
            # Only assert the group when the release actually lists it; a
            # G-suffixed name absent from hla_nom_g.txt must not be echoed
            # Back as though it were a recognised group.
            protein_group = clean
        else:
            is_novel = True
        is_ambiguous = True
    elif group_suffix == "P":
        members = p_groups.get(clean)
        if members:
            imgt_accession = alleles_map.get(members[0])
            # Map to the corresponding G-group if available.
            protein_group = allele_to_g.get(members[0])
        else:
            is_novel = True
        is_ambiguous = True
    elif category == "exact":
        imgt_accession = alleles_map.get(clean)
        protein_group = allele_to_g.get(clean)
        if protein_group is None and resolution_level < 4:
            # Try to infer a G-group from longer equivalents.
            for longer, group in allele_to_g.items():
                if longer.startswith(clean + ":"):
                    protein_group = group
                    break
        is_ambiguous = resolution_level < 2
    elif prefix_matches:
        is_ambiguous = True
        g_of_first = allele_to_g.get(prefix_matches[0])
        if g_of_first and all(allele_to_g.get(a) == g_of_first for a in prefix_matches):
            protein_group = g_of_first
    else:
        is_novel = True
        is_ambiguous = True

    return NormalizedAllele(
        allele_name=clean,
        imgt_accession=imgt_accession,
        protein_group=protein_group,
        hla_class=hla_class,
        gene=gene,
        resolution_level=resolution_level,
        is_ambiguous=is_ambiguous,
        is_novel=is_novel,
        expression_suffix=expression_suffix,
        imgt_match_category=category,
        imgt_match_candidates=candidates,
    )


# ---------------------------------------------------------------------------
# Ambiguity resolution
# ---------------------------------------------------------------------------


def resolve_ambiguity(
    allele: str,
    imgt_db: Dict[str, Any],
) -> List[str]:
    """
    Expand a truncated or group allele into the list of possible full
    allele names, sorted alphabetically.

    - ``"A*02:01:01G"`` → members of the G-group.
    - ``"A*02:01P"`` → members of the P-group.
    - ``"A*02"`` → every full allele whose name begins with ``A*02:``.
    - A fully specified (8-digit) allele with no group suffix is
      returned as a single-item list.

    Parameters
    ----------
    allele : str
        Allele expression.
    imgt_db : dict
        Dictionary produced by :func:`load_imgt_db`.

    Returns
    -------
    list of str
        Sorted list of candidate full allele names.
    """
    if _is_null_token(allele):
        return []

    clean = _strip_hla_prefix(allele.strip())
    if not _ALLELE_BASIC_RE.match(clean):
        raise InvalidAlleleError(f"Does not conform to IMGT nomenclature: {allele!r}")

    g_groups: Dict[str, List[str]] = imgt_db["g_groups"]
    p_groups: Dict[str, List[str]] = imgt_db["p_groups"]
    alleles_map: Dict[str, str] = imgt_db["alleles"]

    group_suffix = _is_group_allele(clean)
    if group_suffix == "G":
        return sorted(g_groups.get(clean, []))
    if group_suffix == "P":
        return sorted(p_groups.get(clean, []))

    if clean in alleles_map:
        return [clean]

    prefix = clean + ":"
    matches = [a for a in alleles_map if a.startswith(prefix)]
    matches.sort()
    return matches


# ---------------------------------------------------------------------------
# Batch normalization
# ---------------------------------------------------------------------------


def batch_normalize(
    genotypes: Sequence[HLAGenotype],
    db_path: Optional[Path] = None,
    *,
    max_workers: Optional[int] = None,
    imgt_db: Optional[Dict[str, Any]] = None,
) -> List[NormalizedAllele]:
    """
    Parallel-normalize the allele1/allele2 values in a sequence of
    :class:`HLAGenotype` records.

    Output order is deterministic: for each input genotype, ``allele1``
    comes first, then ``allele2`` when present. Null / invalid inputs
    that :func:`normalize_allele` returns ``None`` for are skipped.

    Parameters
    ----------
    genotypes : sequence of HLAGenotype
        Inputs to normalize.
    db_path : Path, optional
        Local IPD-IMGT/HLA directory. Used when ``imgt_db`` is not
        provided.
    max_workers : int, optional
        Thread-pool size. Defaults to ``min(32, len * 2)``.
    imgt_db : dict, optional
        Pre-loaded IPD-IMGT/HLA dictionary (for performance).

    Returns
    -------
    list of NormalizedAllele
        Ordered normalized alleles.
    """
    if imgt_db is None:
        imgt_db = load_imgt_db(db_path)

    tasks: List[Tuple[Any, ...]] = []
    for genotype in genotypes:
        tasks.append((genotype, 0, genotype.allele1))
        if genotype.allele2 is not None:
            tasks.append((genotype, 1, genotype.allele2))

    if not tasks:
        return []

    if max_workers is None:
        max_workers = max(1, min(32, len(tasks) * 2))

    results: List[Optional[NormalizedAllele]] = [None] * len(tasks)
    errors: List[BaseException] = []

    def _worker(idx: int, genotype: HLAGenotype, allele_idx: int, allele: str) -> None:
        try:
            norm = normalize_allele(allele, imgt_db)
            if norm is not None:
                norm.sample_id = genotype.sample_id
                norm.source_tool = genotype.tool
                norm.source_locus = genotype.locus
                norm.source_resolution = genotype.resolution
                norm.allele_index = allele_idx
                norm.caller_quality = (
                    genotype.caller_quality1
                    if allele_idx == 0
                    else genotype.caller_quality2
                )
            results[idx] = norm
        except BaseException as exc:  # noqa: BLE001
            errors.append(exc)

    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as pool:
        futures = [
            pool.submit(_worker, idx, genotype, allele_idx, allele)
            for idx, (genotype, allele_idx, allele) in enumerate(tasks)
        ]
        for fut in concurrent.futures.as_completed(futures):
            # Propagate worker errors.
            fut.result()

    if errors:
        # Deterministic: raise the first error seen.
        raise errors[0]

    return [norm for norm in results if norm is not None]


__all__ = [
    "NormalizedAllele",
    "MATCH_CATEGORIES",
    "classify_imgt_match",
    "HLANormalizerError",
    "IMGTDatabaseMissingError",
    "InvalidAlleleError",
    "STALE_THRESHOLD_DAYS",
    "NULL_ALLELE_TOKENS",
    "load_imgt_db",
    "normalize_allele",
    "resolve_ambiguity",
    "batch_normalize",
]
