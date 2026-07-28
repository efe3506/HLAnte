"""
hlante.db.gwas
=============

GWAS Catalog interface — **bulk TSV dump** approach.

Background
----------
The first iteration of this module attempted to filter associations
through the REST API using query parameters such as ``?hla=`` or
``?strongestAllele=``. In the GWAS Catalog's Spring Data REST setup
these parameters are **silently ignored** (the server returns the
unfiltered top-N), and no official endpoint supports allele-level
filtering (only ``findByRsId``, ``findByEfoTrait``, ``findByPubmedId``,
and ``findByStudyAccessionId`` are available). This module therefore
downloads the bulk ``gwas-catalog-associations_ontology-annotated-
full.zip`` dump and indexes the HLA alleles in the
``STRONGEST SNP-RISK ALLELE`` column for local lookup.

Operating modes
---------------
- Normal: :meth:`GWASClient.update` downloads the zip and extracts it
  under ``local_dir``; :meth:`query_allele` performs in-memory hash
  lookups on the index.
- Offline: requires a local copy; returns an empty list when missing.
- Test: point ``local_dir`` at a fixture directory or inject a
  ``fetcher`` that returns a synthetic zip payload.
"""

from __future__ import annotations

import csv
import dataclasses
import io
import logging
import re
import time
import urllib.error
import urllib.request
import zipfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Tuple, cast

logger: logging.Logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

GWAS_CATALOG_BASE: str = "https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest"
GWAS_CATALOG_FULL_ZIP_URL: str = (
    f"{GWAS_CATALOG_BASE}/gwas-catalog-associations_ontology-annotated-full.zip"
)

# Backward compatibility: the REST URL constant is kept so that imports
# From older code keep working, but it is no longer used.
GWAS_CATALOG_API_BASE: str = "https://www.ebi.ac.uk/gwas/rest/api"

DEFAULT_LOCAL_DIR: Path = Path.home() / ".hlante" / "gwas"
DEFAULT_CACHE_DIR: Path = DEFAULT_LOCAL_DIR  # Backward-compat alias

DEFAULT_P_VALUE_THRESHOLD: float = 5e-8
DEFAULT_MAX_RETRIES: int = 3
DEFAULT_TIMEOUT: float = 300.0

# Backward compatibility shims for older signatures.
DEFAULT_CACHE_TTL = None
DEFAULT_RATE_PER_SEC: float = 0.0

GWAS_TSV_FILENAME_DEFAULT: str = "gwas-catalog-download-associations-alt-full.tsv"
GWAS_HLA_SUBSET_FILENAME: str = "gwas-hla-subset.tsv"

# Required TSV columns
_COL_STRONGEST: str = "STRONGEST SNP-RISK ALLELE"
_COL_DISEASE: str = "DISEASE/TRAIT"
_COL_MAPPED_TRAIT: str = "MAPPED_TRAIT"
_COL_MAPPED_GENE: str = "MAPPED_GENE"
_COL_PVALUE: str = "P-VALUE"
_COL_OR_BETA: str = "OR or BETA"
_COL_PUBMED: str = "PUBMEDID"
_COL_STUDY_ACC: str = "STUDY ACCESSION"

# Captures any of "HLA-B*57:01", "HLA-B*57:01-?", "B*57:01",
# "B*57:01-?"; the trailing single-letter risk-allele suffix
# (``-?``, ``-A``, ``-G`` etc.) is ignored.
_ALLELE_RE: re.Pattern[str] = re.compile(
    r"^(?:HLA-)?([A-Z]+\d*\*\d{2,3}(?::\d{2,3})*)(?:-[A-Z?])?$"
)


# ---------------------------------------------------------------------------
# Exceptions
# ---------------------------------------------------------------------------


class GWASDatabaseError(Exception):
    """
    Raised for GWAS Catalog download/parse errors.
    """

    pass


class GWASDownloadError(GWASDatabaseError):
    """
    Raised when the bulk dump cannot be downloaded.
    """

    pass


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------


@dataclass
class GWASHit:
    """
    A single GWAS Catalog hit.

    Attributes
    ----------
    trait : str
        Disease / phenotype (``DISEASE/TRAIT`` or ``MAPPED_TRAIT``).
    p_value : float or None
        Statistical significance.
    odds_ratio : float or None
        OR / beta value (``OR or BETA`` column).
    pmid : str or None
        PubMed identifier.
    study_accession : str or None
        GWAS Catalog study accession (``GCSTxxxxxxx``).
    allele : str
        Allele used during the query (gene-prefixed, HLA-less form).
    url : str, optional
        Source URL (``https://www.ebi.ac.uk/gwas/studies/<acc>``).
    """

    trait: str
    p_value: Optional[float]
    odds_ratio: Optional[float]
    pmid: Optional[str]
    study_accession: Optional[str]
    allele: str
    url: Optional[str] = None
    trait_was_deprecated: bool = False
    effect_size_warning: str = ""
    #: Resolution of the catalogue record that matched, as a statement about
    #: the record itself: ``allele`` (three or four fields), ``subtype`` (two
    #: fields) or ``allele_group`` (one field). It does NOT say whether the
    #: query was broadened — see :attr:`match_was_broadened` for that.
    annotation_scope: str = "allele"
    #: Number of keys in the GWAS index that sit under the matched prefix.
    index_siblings: int = 1
    #: The catalogue key that actually matched, after any truncation.
    matched_allele: str = ""
    #: True when the query allele had to be truncated before a record was
    #: found, i.e. the association is reported for a less specific name than
    #: the one submitted.
    match_was_broadened: bool = False


# ---------------------------------------------------------------------------
# Obsolete EFO remapping
# ---------------------------------------------------------------------------


#: Curated static mapping of GWAS Catalog "obsolete_*" traits to
#: their current EFO replacement terms (or a best-effort label).
#: Used by :func:`_remap_trait`. Planned: replace with dynamic
#: ``efo.obo`` ingestion in v0.2.0 via :func:`_load_efo_obsolete_map`.
OBSOLETE_EFO_MAP: Dict[str, str] = {
    "obsolete_myositis": "inflammatory myopathy (EFO:0000574)",
    "obsolete_uveitis": "uveitis (EFO:0004284)",
    "obsolete_juvenile idiopathic arthritis": "juvenile idiopathic arthritis (EFO:0000685)",
    "obsolete_sclerosing cholangitis": "primary sclerosing cholangitis (EFO:0004268)",
    "obsolete_neuromyelitis optica": "neuromyelitis optica spectrum disorder (EFO:0009800)",
    "obsolete_Autoimmune Hepatitis": "autoimmune hepatitis (EFO:0000537)",
    "obsolete_late-onset myasthenia gravis": "myasthenia gravis (EFO:0000712)",
    "obsolete_juvenile dermatomyositis": "juvenile dermatomyositis (EFO:0009795)",
    "obsolete_drug-induced liver injury": "drug-induced liver injury (EFO:0004228)",
    "obsolete_cervical carcinoma": "cervical carcinoma (EFO:0001061)",
    "obsolete_toxic epidermal necrolysis": "toxic epidermal necrolysis (EFO:0004197)",
    "obsolete_Stevens-Johnson syndrome": "Stevens-Johnson syndrome (EFO:0004190)",
    "obsolete_Behcet's syndrome": "Behcet syndrome (EFO:0000246)",
}


def _remap_trait(trait: str) -> Tuple[str, bool]:
    """
    Remap a GWAS Catalog trait string.

    Parameters
    ----------
    trait : str
        Raw trait label as provided by the GWAS Catalog.

    Returns
    -------
    (remapped_trait, was_deprecated) : tuple
        When the raw label carries at least one ``obsolete_…`` token:

        * If the full label is listed in :data:`OBSOLETE_EFO_MAP`, the
          clean current EFO term is returned (e.g. ``"inflammatory
          myopathy (EFO:0000574)"``).
        * The GWAS Catalog sometimes packs several traits into one
          cell separated by ``", "`` (e.g.
          ``"obsolete_uveitis, obsolete_juvenile idiopathic
          arthritis"``). In that case the label is split, each part
          is remapped individually, and the parts are rejoined. Parts
          that start with ``obsolete_`` but are not in
          :data:`OBSOLETE_EFO_MAP` pass through unchanged — the
          ``obsolete_`` prefix already signals the status.

        In every branch where at least one ``obsolete_`` token is
        present, ``was_deprecated`` is ``True`` so downstream
        consumers (TSV ``trait_was_deprecated`` metadata, Markdown
        provenance note, JSON field) can emit their own provenance
        message without the trait cell carrying a trailing ``[...]``
        suffix.

        Non-obsolete traits are returned unchanged with ``False``.
    """
    if not trait:
        return trait, False
    # Fast path: exact full-string match against the curated map.
    if trait in OBSOLETE_EFO_MAP:
        return OBSOLETE_EFO_MAP[trait], True
    if "obsolete_" not in trait:
        return trait, False
    # Compound path: split the Catalog cell on ", " and remap each
    # Part independently. A single obsolete token anywhere in the
    # Cell is enough to set ``was_deprecated``.
    parts = trait.split(", ")
    if len(parts) == 1:
        return trait, True
    remapped: List[str] = []
    was_deprecated = False
    for part in parts:
        if part.startswith("obsolete_"):
            was_deprecated = True
            remapped.append(OBSOLETE_EFO_MAP.get(part, part))
        else:
            remapped.append(part)
    return ", ".join(remapped), was_deprecated


def _load_efo_obsolete_map(efo_obo_path: Path) -> Dict[str, str]:
    """
    Parse EBI EFO ``.obo`` to build obsolete → replacement map.

    Intended as a future hook for ``db-update --db gwas --efo-obo
    path/to/efo.obo``. Stub implementation returns an empty dict.

    Notes
    -----
    Full implementation planned for v0.2.0. Download
    ``efo.obo`` from https://www.ebi.ac.uk/efo/efo.obo and parse
    ``[Term]`` blocks with ``is_obsolete: true`` plus ``replaced_by``
    / ``consider`` tags.
    """
    return {}


# ---------------------------------------------------------------------------
# Extreme effect-size classification
# ---------------------------------------------------------------------------


#: Heuristic keywords that flag a GWAS trait as a continuous /
#: quantitative measure. Effect sizes reported for these traits are
#: often regression betas or per-SD effects that are NOT comparable
#: to binary-disease odds ratios.
_QUANTITATIVE_TRAIT_KEYWORDS: Tuple[str, ...] = (
    "amount",
    "level",
    "measurement",
    "count",
    "concentration",
    "seropositivity",
    "response",
    "ratio",
)


def _classify_effect_size(
    or_value: Optional[float],
    trait: str,
) -> str:
    """
    Classify a GWAS effect size for user warnings.

    Returns a pipe-joined warning code:

    - ``""`` — normal binary-disease effect size.
    - ``"extreme_value"`` — ``or_value > 10``; unusual for a
      binary-disease odds ratio.
    - ``"extreme_value|quantitative_trait_effect"`` — ``or_value >
      10`` *and* the trait text matches quantitative keywords; the
      reported value may be a β coefficient, not an OR.
    """
    if or_value is None or or_value <= 10:
        return ""
    trait_lower = (trait or "").lower()
    is_quant = any(kw in trait_lower for kw in _QUANTITATIVE_TRAIT_KEYWORDS)
    if is_quant:
        return "extreme_value|quantitative_trait_effect"
    return "extreme_value"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _strip_hla(allele: str) -> str:
    if allele.upper().startswith("HLA-"):
        return allele[4:]
    return allele


def _coerce_float(value: Any) -> Optional[float]:
    if value is None:
        return None
    stripped = str(value).strip()
    if not stripped or stripped.lower() in {"na", "nan", "none", "-"}:
        return None
    try:
        return float(stripped)
    except ValueError:
        return None


def _extract_allele_name(raw: str) -> Optional[str]:
    """
    Extract the canonical allele name from a ``STRONGEST SNP-RISK
    ALLELE`` cell.

    Examples:
    - ``"HLA-B*57:01"`` → ``"B*57:01"``
    - ``"B*57:01-?"`` → ``"B*57:01"``
    - ``"rs12345-T"`` → ``None`` (rs-based, not HLA)
    """
    if not raw:
        return None
    cell = raw.strip()
    match = _ALLELE_RE.match(cell)
    if match:
        return match.group(1)
    return None


_NOMENCLATURE_SUFFIX_RE: re.Pattern[str] = re.compile(r"[A-Z]$")

#: Colon-group count → user-facing TSV label (digit-based).
_FIELD_LABELS: Dict[int, str] = {
    1: "one-field",
    2: "two-field",
    3: "three-field",
    4: "four-field",
}

#: Placeholder returned when no match was found at any fallback level.
RESOLUTION_LABEL_NONE: str = "none"


def _colon_group_count(allele: str) -> int:
    """
    Count colon-separated groups in an allele.

    Examples: ``"A*02:01:01:01"`` → 4, ``"A*02"`` → 1, ``"A"`` → 0.
    """
    core = _strip_hla(allele)
    if "*" not in core:
        return 0
    _, rest = core.split("*", 1)
    rest = _NOMENCLATURE_SUFFIX_RE.sub("", rest)
    if not rest:
        return 0
    return len(rest.split(":"))


def _truncate_to_fields(allele: str, n_fields: int) -> str:
    """
    Truncate an allele to at most ``n_fields`` colon-groups.

    ``n_fields`` here denotes the *colon-group count*, which is distinct
    from the field-count label used in the user-facing report.
    Examples:

    - ``_truncate_to_fields("A*02:01:01:01", 2)`` → ``"A*02:01"``
    - ``_truncate_to_fields("DRB1*03:01:01", 4)`` → ``"DRB1*03:01:01"``
      (already at 3 groups, no truncation needed)
    - ``_truncate_to_fields("B*57:01G", 2)`` → ``"B*57:01G"`` (G suffix kept)
    - ``_truncate_to_fields("HLA-A*02:01", 1)`` → ``"A*02"`` (HLA- prefix
      stripped)

    Parameters
    ----------
    allele : str
        Source allele.
    n_fields : int
        Colon-group count to keep (≥1). Values ≤ 0 yield no truncation.

    Returns
    -------
    str
        Truncated allele.
    """
    core = _strip_hla(allele)
    if "*" not in core or n_fields < 1:
        return core
    gene, rest = core.split("*", 1)
    suffix_match = _NOMENCLATURE_SUFFIX_RE.search(rest)
    suffix = suffix_match.group() if suffix_match else ""
    body = rest[: -len(suffix)] if suffix else rest
    parts = body.split(":") if body else []
    kept = parts[:n_fields]
    return f"{gene}*{':'.join(kept)}{suffix}"


def _get_resolution_levels(allele: str) -> List[int]:
    """
    Return colon-group counts from the current resolution down to ``1``.

    Examples: ``"A*02:01:01:01"`` → ``[4, 3, 2, 1]``,
    ``"A*02:01"`` → ``[2, 1]``, ``"A*02"`` → ``[1]``.
    """
    current = _colon_group_count(allele)
    if current <= 0:
        return []
    return list(range(current, 0, -1))


def _field_label(colon_groups: int) -> str:
    """
    Map a colon-group count to the digit-based label used by the report.
    """
    return _FIELD_LABELS.get(colon_groups, f"{colon_groups}-field")


# ---------------------------------------------------------------------------
# GWAS client
# ---------------------------------------------------------------------------


class GWASClient:
    """
    GWAS Catalog bulk-TSV client.

    Parameters
    ----------
    local_dir : Path, optional
        Directory of the downloaded TSV (or the HLA subset).
        Defaults to :data:`DEFAULT_LOCAL_DIR`.
    p_value_threshold : float, optional
        Upper bound on p-values for returned hits.
    offline : bool, optional
        If ``True`` the client never downloads; a local copy is required.
    fetcher : callable, optional
        URL → bytes callable (returning zip content). Overridable for
        tests.
    max_retries : int, optional
        Download retry limit.
    timeout : float, optional
        HTTP timeout (seconds).
    sleep : callable, optional
        Back-off sleep function (``lambda _: None`` in tests).
    """

    def __init__(
        self,
        local_dir: Optional[Path] = None,
        *,
        p_value_threshold: float = DEFAULT_P_VALUE_THRESHOLD,
        offline: bool = False,
        fetcher: Optional[Callable[[str], bytes]] = None,
        max_retries: int = DEFAULT_MAX_RETRIES,
        timeout: float = DEFAULT_TIMEOUT,
        sleep: Callable[[float], None] = time.sleep,
        # Backward-compat parameters (now inert)
        cache_dir: Optional[Path] = None,
        rate_per_sec: float = DEFAULT_RATE_PER_SEC,  # noqa: ARG002
        cache_ttl: Any = None,  # noqa: ARG002
    ) -> None:
        if local_dir is None and cache_dir is not None:
            local_dir = cache_dir  # Support the legacy keyword name
        self.local_dir: Path = Path(local_dir) if local_dir else DEFAULT_LOCAL_DIR
        self.p_value_threshold: float = p_value_threshold
        self.offline: bool = offline
        self.max_retries: int = max_retries
        self.timeout: float = timeout
        self._fetcher: Callable[[str], bytes] = fetcher or self._default_fetcher
        self._sleep: Callable[[float], None] = sleep
        self._by_allele: Dict[str, List[GWASHit]] = {}
        self._loaded: bool = False

    # ---- Public API ----

    def update(self) -> Path:
        """
        Download the bulk GWAS zip, extract it, and emit an HLA subset.

        Returns
        -------
        Path
            Root download directory.

        Raises
        ------
        GWASDownloadError
            If downloading or extraction fails.
        """
        if self.offline:
            raise GWASDownloadError("Cannot download GWAS data in offline mode.")
        self.local_dir.mkdir(parents=True, exist_ok=True)
        raw = self._download_with_retry(GWAS_CATALOG_FULL_ZIP_URL)
        try:
            with zipfile.ZipFile(io.BytesIO(raw)) as zf:
                zf.extractall(self.local_dir)
        except (zipfile.BadZipFile, OSError) as exc:
            raise GWASDownloadError(f"GWAS zip extraction failed: {exc}") from exc
        # Emit an HLA subset so subsequent loads avoid parsing the full TSV.
        self._write_hla_subset()
        self._loaded = False
        logger.info("GWAS Catalog updated: %s", self.local_dir)
        return self.local_dir

    def load(self) -> None:
        """
        Load the local TSV into memory (idempotent).

        Raises
        ------
        GWASDatabaseError
            If the expected TSV cannot be located.
        """
        if self._loaded:
            return
        tsv_path = self._locate_tsv()
        if tsv_path is None:
            raise GWASDatabaseError(
                f"GWAS TSV not found in {self.local_dir}. "
                "Call update() first (or `hlante db-update --db gwas`)."
            )
        self._by_allele = self._parse_tsv(tsv_path)
        self._loaded = True
        logger.info(
            "GWAS Catalog loaded: %d unique HLA allele(s), %d hit(s)",
            len(self._by_allele),
            sum(len(v) for v in self._by_allele.values()),
        )

    def query_allele(self, allele: str) -> List[GWASHit]:
        """
        Look up hits for an HLA allele in the in-memory index.

        Parameters
        ----------
        allele : str
            Allele in IPD-IMGT/HLA format (``HLA-`` prefix optional).

        Returns
        -------
        list of GWASHit
            Hits whose p-value passes :attr:`p_value_threshold`.
        """
        try:
            self.load()
        except GWASDatabaseError as exc:
            if self.offline:
                logger.warning(
                    "GWAS: offline mode and local TSV missing (%s); returning empty result.",
                    exc,
                )
                return []
            raise

        bare = _strip_hla(allele)
        hits = self._by_allele.get(bare, [])
        # A genome-wide-significance filter must require an actual p-value at or
        # Below threshold. A hit whose p-value cell did not parse (None) has not
        # Demonstrably met the bar and is excluded, rather than passed through
        # as if significant.
        return [h for h in hits if h.p_value is not None and h.p_value <= self.p_value_threshold]

    def query_allele_with_fallback(
        self,
        allele: str,
        min_resolution: int = 2,
    ) -> Tuple[List["GWASHit"], str]:
        """
        Query with stepwise truncation when a full-resolution hit is
        unavailable.

        Example: if ``A*03:02:01:01`` has no record in GWAS, the client
        tries ``A*03:02:01`` → ``A*03:02`` → ``A*03`` in order.

        Parameters
        ----------
        allele : str
            Allele to query (HLA- prefix optional).
        min_resolution : int, optional
            Lowest *digit-based* level to descend to (2, 4, 6, or 8).
            Default ``2`` — allows descent down to ``A*02``.

        Returns
        -------
        (hits, resolution_label) : tuple
            - ``hits``: hits passing the p-value threshold (possibly empty).
            - ``resolution_label``: the resolution tier at which a match
              was found (``"one-field"`` … ``"four-field"``) or ``"none"``
              when no tier matched.
        """
        # Translate the digit-count limit to a colon-group floor.
        min_colon_groups = max(1, min_resolution // 2)
        for n_colon in _get_resolution_levels(allele):
            if n_colon < min_colon_groups:
                break
            truncated = _truncate_to_fields(allele, n_colon)
            hits = self.query_allele(truncated)
            if hits:
                label = _field_label(n_colon)
                # Scope describes the matched catalogue record's own
                # Resolution; broadening is recorded separately.
                scope = (
                    "allele"
                    if n_colon >= 3
                    else ("subtype" if n_colon == 2 else "allele_group")
                )
                broadened = n_colon < _colon_group_count(allele)
                # Index_siblings counts keys in THIS GWAS index under the
                # Truncated prefix — not IPD-IMGT/HLA alleles.
                expansion = (
                    1
                    if n_colon == _colon_group_count(allele)
                    else sum(
                        1
                        for key in self._by_allele
                        if key == truncated or key.startswith(truncated + ":")
                    )
                )
                # Copy rather than mutate: the hit objects live in the
                # Client's index and are shared by every query, so assigning
                # to them corrupts results already returned to a caller.
                hits = [
                    dataclasses.replace(
                        hit,
                        annotation_scope=scope,
                        index_siblings=max(1, expansion),
                        matched_allele=truncated,
                        match_was_broadened=broadened,
                    )
                    for hit in hits
                ]
                logger.debug(
                    "GWAS fallback: %r → %r (%s, scope=%s, exp=%d) %d hit(s)",
                    allele,
                    truncated,
                    label,
                    scope,
                    expansion,
                    len(hits),
                )
                return hits, label
        return [], RESOLUTION_LABEL_NONE

    # ---- Internal helpers ----

    def _locate_tsv(self) -> Optional[Path]:
        """
        Locate the HLA subset (preferred) or the full TSV.
        """
        subset = self.local_dir / GWAS_HLA_SUBSET_FILENAME
        if subset.is_file():
            return subset
        full = self.local_dir / GWAS_TSV_FILENAME_DEFAULT
        if full.is_file():
            return full
        for match in self.local_dir.rglob("*associations*full*.tsv"):
            return match
        for match in self.local_dir.rglob("*.tsv"):
            return match
        return None

    def _parse_tsv(self, path: Path) -> Dict[str, List[GWASHit]]:
        """
        Parse the GWAS TSV and group hits by allele.
        """
        by_allele: Dict[str, List[GWASHit]] = {}
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames is None:
                return by_allele
            cols = {name.strip().upper(): name for name in reader.fieldnames}

            def pick(*candidates: str) -> Optional[str]:
                for c in candidates:
                    real = cols.get(c.upper())
                    if real:
                        return real
                return None

            f_allele = pick(_COL_STRONGEST)
            f_disease = pick(_COL_DISEASE)
            f_mapped_trait = pick(_COL_MAPPED_TRAIT)
            f_pvalue = pick(_COL_PVALUE)
            f_or = pick(_COL_OR_BETA)
            f_pmid = pick(_COL_PUBMED)
            f_study = pick(_COL_STUDY_ACC)

            if f_allele is None or f_disease is None:
                raise GWASDatabaseError(
                    f"GWAS TSV is missing required columns: {reader.fieldnames}"
                )

            for row in reader:
                allele_name = _extract_allele_name(row.get(f_allele, ""))
                if allele_name is None:
                    continue
                raw_trait = (row.get(f_mapped_trait) or "" if f_mapped_trait else "").strip() or (
                    row.get(f_disease) or ""
                ).strip()
                if not raw_trait:
                    continue
                # Remap deprecated EFO terms at ingestion time.
                trait, trait_was_deprecated = _remap_trait(raw_trait)
                p_value = _coerce_float(row.get(f_pvalue) if f_pvalue else None)
                or_val = _coerce_float(row.get(f_or) if f_or else None)
                # Flag extreme / likely-quantitative effect sizes.
                effect_size_warning = _classify_effect_size(or_val, trait)
                pmid = ((row.get(f_pmid) or "").strip() if f_pmid else "") or None
                study = ((row.get(f_study) or "").strip() if f_study else "") or None
                url = f"https://www.ebi.ac.uk/gwas/studies/{study}" if study else None
                hit = GWASHit(
                    trait=trait,
                    p_value=p_value,
                    odds_ratio=or_val,
                    pmid=pmid,
                    study_accession=study,
                    allele=allele_name,
                    url=url,
                    trait_was_deprecated=trait_was_deprecated,
                    effect_size_warning=effect_size_warning,
                )
                by_allele.setdefault(allele_name, []).append(hit)
        return by_allele

    def _write_hla_subset(self) -> Optional[Path]:
        """
        After :meth:`update`, write a smaller subset containing only
        rows that mention an HLA allele. Later :meth:`load` calls will
        read this instead of the 59 MB full TSV.
        """
        full = self.local_dir / GWAS_TSV_FILENAME_DEFAULT
        # Fall back to whatever TSV name the archive actually produced.
        if not full.is_file():
            for match in self.local_dir.rglob("*associations*.tsv"):
                full = match
                break
        if not full.is_file():
            logger.warning("Full GWAS TSV not found; HLA subset was not written.")
            return None

        subset = self.local_dir / GWAS_HLA_SUBSET_FILENAME
        written = 0
        with full.open("r", encoding="utf-8") as src, subset.open("w", encoding="utf-8") as dst:
            header = src.readline()
            dst.write(header)
            for line in src:
                if "HLA-" in line or re.search(r"\t[A-Z]+\d*\*\d", line):
                    dst.write(line)
                    written += 1
        logger.info("GWAS HLA subset written: %s (%d row(s))", subset, written)
        return subset

    def _download_with_retry(self, url: str) -> bytes:
        """
        Download the bulk zip with exponential backoff on failure.
        """
        last_exc: Optional[BaseException] = None
        for attempt in range(self.max_retries):
            try:
                return self._fetcher(url)
            except (
                urllib.error.URLError,
                urllib.error.HTTPError,
                TimeoutError,
                OSError,
            ) as exc:
                last_exc = exc
                if attempt + 1 < self.max_retries:
                    delay = 2**attempt
                    logger.warning(
                        "GWAS attempt %d/%d failed (%s); retrying in %ds.",
                        attempt + 1,
                        self.max_retries,
                        exc,
                        delay,
                    )
                    self._sleep(float(delay))
        raise GWASDownloadError(
            f"GWAS zip could not be downloaded ({url}): {last_exc}"
        ) from last_exc

    def _default_fetcher(self, url: str) -> bytes:
        req = urllib.request.Request(url, headers={"User-Agent": "HLAnte/0.1"})
        with urllib.request.urlopen(req, timeout=self.timeout) as resp:  # noqa: S310
            return cast(bytes, resp.read())


__all__ = [
    "GWASClient",
    "GWASHit",
    "GWASDatabaseError",
    "GWASDownloadError",
    "GWAS_CATALOG_FULL_ZIP_URL",
    "GWAS_CATALOG_API_BASE",
    "DEFAULT_P_VALUE_THRESHOLD",
    "DEFAULT_CACHE_DIR",
    "DEFAULT_LOCAL_DIR",
    "DEFAULT_MAX_RETRIES",
    "DEFAULT_RATE_PER_SEC",
    "RESOLUTION_LABEL_NONE",
    "OBSOLETE_EFO_MAP",
]
