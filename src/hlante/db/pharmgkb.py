"""
hlante.db.pharmgkb
=================

PharmGKB interface.

Loads the PharmGKB ``clinicalAnnotations.zip`` bulk dump and serves
HLA-drug associations. By default only CPIC evidence levels 1A and 1B
are returned — the highest-confidence, action-grade annotations. Lower
levels (2A/2B/3) can be opted in via ``evidence_levels=``. In offline mode the bulk dump download is skipped;
a local copy is required.
"""

from __future__ import annotations

import csv
import io
import logging
import re
import urllib.error
import urllib.request
import zipfile
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Callable, Dict, FrozenSet, List, Optional, Set, Tuple, cast

logger: logging.Logger = logging.getLogger(__name__)


PHARMGKB_DOWNLOAD_BASE: str = "https://api.pharmgkb.org/v1/download/file/data"
PHARMGKB_CLINICAL_ANN_URL: str = f"{PHARMGKB_DOWNLOAD_BASE}/clinicalAnnotations.zip"

DEFAULT_LOCAL_DIR: Path = Path.home() / ".hlante" / "pharmgkb"
CLINICAL_ANN_FILENAME: str = "clinical_annotations.tsv"

DEFAULT_EVIDENCE_LEVELS: FrozenSet[str] = frozenset({"1A", "1B"})

PHARMGKB_ANNOTATION_URL_TEMPLATE: str = (
    "https://www.pharmgkb.org/clinicalAnnotation/{annotation_id}"
)


# ---------------------------------------------------------------------------
# Exceptions
# ---------------------------------------------------------------------------


class PharmGKBDatabaseError(Exception):
    """
    Base class for PharmGKB errors.
    """

    pass


class PharmGKBDownloadError(PharmGKBDatabaseError):
    """
    Raised when the PharmGKB bulk zip cannot be downloaded.
    """

    pass


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------


@dataclass
class PharmAnnotation:
    """
    A single PharmGKB clinical annotation.

    Attributes
    ----------
    drug : str
        Related drug (e.g., ``"abacavir"``).
    phenotype : str
        Phenotype / reaction text.
    evidence_level : str
        PharmGKB evidence level (``"1A"``, ``"1B"``, ``"2A"``, ...).
    pmid : list of str
        Supporting PubMed IDs (at least one required).
    allele : str
        HLA allele the annotation refers to.
    annotation_id : str, optional
        PharmGKB ``CA`` identifier.
    gene : str, optional
        Associated HLA gene (e.g., ``"HLA-B"``).
    cpic_url : str, optional
        Linked CPIC guideline URL, if any.
    pharmgkb_url : str, optional
        Source PharmGKB page.
    matched_form : str, optional
        Which candidate allele form produced the match (e.g., when the
        query is ``"A*02:01:01"`` but PharmGKB only stores ``"A*02:01"``,
        ``matched_form="A*02:01"``).
    """

    drug: str
    phenotype: str
    evidence_level: str
    pmid: List[str]
    allele: str
    annotation_id: Optional[str] = None
    gene: Optional[str] = None
    cpic_url: Optional[str] = None
    pharmgkb_url: Optional[str] = None
    matched_form: Optional[str] = None


# ---------------------------------------------------------------------------
# PharmGKB client
# ---------------------------------------------------------------------------


class PharmGKBClient:
    """
    PharmGKB bulk-data client.

    Parameters
    ----------
    local_dir : Path, optional
        Directory in which the bulk zip has been extracted.
        Defaults to :data:`DEFAULT_LOCAL_DIR`.
    evidence_levels : frozenset of str, optional
        Accepted evidence levels.
    offline : bool, optional
        If ``True``, remote download is not attempted; a local copy is
        required.
    fetcher : callable, optional
        URL → bytes callable (injectable for tests). Uses urllib by
        default.
    timeout : float, optional
        HTTP timeout.
    """

    def __init__(
        self,
        local_dir: Optional[Path] = None,
        *,
        evidence_levels: FrozenSet[str] = DEFAULT_EVIDENCE_LEVELS,
        offline: bool = False,
        fetcher: Optional[Callable[[str], bytes]] = None,
        timeout: float = 60.0,
    ) -> None:
        self.local_dir: Path = Path(local_dir) if local_dir else DEFAULT_LOCAL_DIR
        self.evidence_levels: FrozenSet[str] = frozenset(level.upper() for level in evidence_levels)
        self.offline: bool = offline
        self._fetcher: Callable[[str], bytes] = fetcher or self._default_fetcher
        self.timeout: float = timeout
        self._by_allele: Dict[str, List[PharmAnnotation]] = {}
        self._loaded: bool = False

    # ---- Public API ----

    def update(self) -> Path:
        """
        Download the PharmGKB ``clinicalAnnotations.zip`` and extract it.

        Returns
        -------
        Path
            Extraction directory.

        Raises
        ------
        PharmGKBDownloadError
            If downloading or extraction fails.
        """
        if self.offline:
            raise PharmGKBDownloadError("Cannot download in offline mode.")
        self.local_dir.mkdir(parents=True, exist_ok=True)
        try:
            raw = self._fetcher(PHARMGKB_CLINICAL_ANN_URL)
        except (urllib.error.URLError, urllib.error.HTTPError, OSError, TimeoutError) as exc:
            raise PharmGKBDownloadError(f"PharmGKB download failed: {exc}") from exc
        try:
            with zipfile.ZipFile(io.BytesIO(raw)) as zf:
                zf.extractall(self.local_dir)
        except (zipfile.BadZipFile, OSError) as exc:
            raise PharmGKBDownloadError(f"PharmGKB zip extraction failed: {exc}") from exc
        self._loaded = False
        logger.info("PharmGKB updated: %s", self.local_dir)
        return self.local_dir

    def load(self) -> None:
        """
        Load the local TSV into memory (idempotent).

        Raises
        ------
        PharmGKBDatabaseError
            If the expected TSV is missing.
        """
        if self._loaded:
            return
        tsv_path = self._locate_tsv()
        if tsv_path is None:
            raise PharmGKBDatabaseError(
                f"PharmGKB TSV not found in {self.local_dir}. "
                "Call update() first or pass the correct local_dir."
            )
        self._by_allele = self._parse_tsv(tsv_path)
        self._loaded = True
        logger.info(
            "PharmGKB loaded: %d allele(s), %d annotation(s)",
            len(self._by_allele),
            sum(len(v) for v in self._by_allele.values()),
        )

    def query_allele(self, allele: str) -> List[PharmAnnotation]:
        """
        Return all PharmGKB annotations for an HLA allele.

        The query is tried against several candidate forms (full,
        two- and one-field truncations, ``*``-prefixed forms). Results are
        de-duplicated by (annotation_id, drug, evidence_level, allele);
        each returned record stores the matched form in
        ``matched_form``.

        Parameters
        ----------
        allele : str
            Allele in IPD-IMGT/HLA format (``HLA-`` prefix optional).

        Returns
        -------
        list of PharmAnnotation
            Records passing the evidence-level filter.
        """
        self.load()
        seen: Set[Tuple[Optional[str], str, str, str]] = set()
        results: List[PharmAnnotation] = []
        for form in _candidate_allele_forms(allele):
            lookup_key = _strip_hla(form)
            if not lookup_key:
                continue
            for ann in self._by_allele.get(lookup_key, []):
                dedupe_key = (
                    ann.annotation_id,
                    ann.drug,
                    ann.evidence_level,
                    ann.allele,
                )
                if dedupe_key in seen:
                    continue
                seen.add(dedupe_key)
                results.append(replace(ann, matched_form=form))
        return results

    # ---- Internal helpers ----

    # File names present in the PharmGKB bulk dump that are *not* the
    # Main clinical_annotations table. They must not be picked up by
    # Glob fallbacks.
    _SIDE_FILENAMES: FrozenSet[str] = frozenset(
        {
            "clinical_ann_alleles.tsv",
            "clinical_ann_evidence.tsv",
            "clinical_ann_history.tsv",
        }
    )
    EVIDENCE_FILENAME: str = "clinical_ann_evidence.tsv"

    def _locate_tsv(self) -> Optional[Path]:
        """
        Locate the main clinical annotations TSV.

        The PharmGKB ``clinicalAnnotations.zip`` archive includes
        multiple TSV files; the primary one is
        ``clinical_annotations.tsv``. Exact name matches are tried
        first; glob fallbacks exclude side files.
        """
        candidates = [
            CLINICAL_ANN_FILENAME,
            "clinical_ann.tsv",
            "clinicalAnnotations.tsv",
        ]
        for name in candidates:
            p = self.local_dir / name
            if p.is_file():
                return p
        # Glob search in subdirectories (ignore side files).
        for match in self.local_dir.rglob("*.tsv"):
            if match.name in self._SIDE_FILENAMES:
                continue
            if match.name.startswith("clinical_annotations") or match.name in {
                "clinical_ann.tsv",
                "clinicalAnnotations.tsv",
            }:
                return match
        return None

    def _locate_evidence(self) -> Optional[Path]:
        """
        Locate the evidence file containing PMID ↔ annotation_id pairs.
        """
        primary = self.local_dir / self.EVIDENCE_FILENAME
        if primary.is_file():
            return primary
        for match in self.local_dir.rglob("clinical_ann_evidence*.tsv"):
            return match
        return None

    def _load_evidence_pmids(self, path: Path) -> Dict[str, List[str]]:
        """
        Extract the PMID list per Clinical Annotation ID from
        ``clinical_ann_evidence.tsv``.
        """
        ca_to_pmids: Dict[str, List[str]] = {}
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if not reader.fieldnames:
                return ca_to_pmids
            fields = {(name or "").strip().lower(): name for name in reader.fieldnames}
            f_ca = fields.get("clinical annotation id")
            f_pmid = fields.get("pmid")
            if not f_ca or not f_pmid:
                return ca_to_pmids
            for row in reader:
                ca = (row.get(f_ca) or "").strip()
                pmid = (row.get(f_pmid) or "").strip()
                if ca and pmid and pmid.lower() != "none":
                    bucket = ca_to_pmids.setdefault(ca, [])
                    if pmid not in bucket:
                        bucket.append(pmid)
        logger.info(
            "PharmGKB evidence loaded: %d annotation(s) with PMIDs (%s)",
            len(ca_to_pmids),
            path,
        )
        return ca_to_pmids

    def _parse_tsv(self, path: Path) -> Dict[str, List[PharmAnnotation]]:
        """
        Parse the PharmGKB TSV and group annotations by allele.

        PMIDs may come from two sources:
        1. A ``PMID`` / ``PMIDs`` column on the main TSV (HLAnte fixture
           layout).
        2. A separate ``clinical_ann_evidence.tsv`` (real PharmGKB dump).
        When the first source is absent, the second is JOIN-ed via
        Clinical Annotation ID.
        """
        evidence_path = self._locate_evidence()
        evidence_pmids: Dict[str, List[str]] = (
            self._load_evidence_pmids(evidence_path) if evidence_path else {}
        )

        by_allele: Dict[str, List[PharmAnnotation]] = {}
        dropped_no_pmid = 0
        accepted = 0
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            fieldnames = {(name or "").strip().lower(): name for name in (reader.fieldnames or [])}

            def col(*candidates: str) -> Optional[str]:
                for cand in candidates:
                    real = fieldnames.get(cand.lower())
                    if real:
                        return real
                return None

            f_allele = col("variant/haplotypes", "haplotype", "variant", "alleles", "allele")
            f_gene = col("gene", "genes")
            f_drug = col("drug(s)", "drug", "drugs", "related chemicals")
            f_phenotype = col("phenotype(s)", "phenotype category", "phenotype", "phenotypes")
            f_evidence = col("level of evidence", "evidence level")
            f_pmid = col("pmid", "pmids", "pubmed id", "pubmed ids")
            f_ann_id = col("clinical annotation id", "annotation id", "ca id")
            f_cpic = col("cpic url", "cpic guideline", "cpic")
            f_url = col("url", "pharmgkb url")

            if f_allele is None or f_drug is None or f_evidence is None:
                raise PharmGKBDatabaseError(
                    f"PharmGKB TSV is missing expected columns: {list(fieldnames)}"
                )

            for row in reader:
                evidence = (row.get(f_evidence) or "").strip().upper()
                if evidence not in self.evidence_levels:
                    continue

                gene = (row.get(f_gene) if f_gene else "") or ""
                variant_cell = (row.get(f_allele) or "").strip()
                if not _is_hla_related(gene, variant_cell):
                    continue

                annotation_id = (row.get(f_ann_id) or "").strip() if f_ann_id else ""

                # PMID: first from the main row, then from evidence.tsv
                pmids = _split_pmids(row.get(f_pmid) if f_pmid else "")
                if not pmids and annotation_id:
                    pmids = list(evidence_pmids.get(annotation_id, []))
                if not pmids:
                    dropped_no_pmid += 1
                    logger.debug(
                        "PharmGKB row has no PMID (neither main nor evidence): %s",
                        annotation_id or variant_cell,
                    )
                    continue

                drug = (row.get(f_drug) or "").strip()
                phenotype = (row.get(f_phenotype) or "").strip() if f_phenotype else ""
                cpic_url = (row.get(f_cpic) or "").strip() if f_cpic else ""
                custom_url = (row.get(f_url) or "").strip() if f_url else ""

                pharmgkb_url = custom_url or (
                    PHARMGKB_ANNOTATION_URL_TEMPLATE.format(annotation_id=annotation_id)
                    if annotation_id
                    else None
                )

                gene_clean = gene.strip() or None
                for allele_name in _expand_variant_field(variant_cell, gene_clean):
                    annotation = PharmAnnotation(
                        drug=drug,
                        phenotype=phenotype,
                        evidence_level=evidence,
                        pmid=pmids,
                        allele=allele_name,
                        annotation_id=annotation_id or None,
                        gene=gene_clean,
                        cpic_url=cpic_url or None,
                        pharmgkb_url=pharmgkb_url,
                    )
                    by_allele.setdefault(allele_name, []).append(annotation)
                    accepted += 1

        if dropped_no_pmid and not evidence_path:
            logger.warning(
                "PharmGKB: %d HLA row(s) were skipped because no PMID was "
                "available. Include clinical_ann_evidence.tsv in the local "
                "directory to recover those PMIDs.",
                dropped_no_pmid,
            )
        logger.info(
            "PharmGKB TSV parsed: %d HLA annotation(s), %d row(s) skipped for missing PMID",
            accepted,
            dropped_no_pmid,
        )
        return by_allele

    def _default_fetcher(self, url: str) -> bytes:
        req = urllib.request.Request(url, headers={"User-Agent": "HLAnte/0.1"})
        with urllib.request.urlopen(req, timeout=self.timeout) as resp:  # noqa: S310
            return cast(bytes, resp.read())


# ---------------------------------------------------------------------------
# Module helpers
# ---------------------------------------------------------------------------


_NOMENCLATURE_SUFFIX_RE: re.Pattern[str] = re.compile(r"[A-Z]$")


def _strip_hla(allele: str) -> str:
    if allele.upper().startswith("HLA-"):
        return allele[4:]
    return allele


def _candidate_allele_forms(allele: str) -> List[str]:
    """
    Produce the ordered set of candidate forms a query allele might
    match in the PharmGKB dump.

    PharmGKB typically stores two-field entries (``B*57:01``) while
    callers may pass three- or four-field identifiers or G/P-suffixed forms from
    arcasHLA (e.g. ``B*57:01G``). This helper generates the full
    truncation list with and without G/P suffix, plus gene-prefix-less
    and ``HLA-`` prefixed variants. The original form comes first.

    Example
    -------
    ``_candidate_allele_forms("B*57:01G")`` →
    ``["B*57:01G", "B*57G", "B*57:01", "B*57", "HLA-B*57:01G", ...]``.
    """
    bare = _strip_hla(allele).strip()
    if not bare or "*" not in bare:
        return [bare] if bare else []
    gene, rest = bare.split("*", 1)
    suffix_match = _NOMENCLATURE_SUFFIX_RE.search(rest)
    suffix = suffix_match.group() if suffix_match else ""
    body = rest[: -len(suffix)] if suffix else rest
    parts = body.split(":") if body else []

    # Build both suffixed and suffix-free truncation sets
    suffixes_to_try: List[str] = [suffix] if suffix else [""]
    if suffix:
        suffixes_to_try.append("")  # Also try without G/P suffix

    truncations: List[str] = []
    for sfx in suffixes_to_try:
        for n in range(len(parts), 0, -1):
            form = f"{gene}*{':'.join(parts[:n])}{sfx}"
            if form not in truncations:
                truncations.append(form)

    # HLA- prefixed forms
    hla_prefixed = [f"HLA-{t}" for t in truncations]

    # Gene-prefix-less forms (PharmGKB sometimes stores "*15:02")
    gene_prefixless: List[str] = []
    for sfx in suffixes_to_try:
        for n in range(len(parts), 0, -1):
            form = f"*{':'.join(parts[:n])}{sfx}"
            if form not in gene_prefixless:
                gene_prefixless.append(form)

    seen: List[str] = []
    for candidate in truncations + hla_prefixed + gene_prefixless:
        if candidate not in seen:
            seen.append(candidate)
    return seen


def _is_hla_related(gene: str, variant_cell: str) -> bool:
    """
    Return ``True`` when a row is HLA-related based on either the gene
    or the variant cell.
    """
    haystack = f"{gene} {variant_cell}".upper()
    return "HLA" in haystack or "HLA-" in haystack


def _expand_variant_field(cell: str, gene: Optional[str] = None) -> List[str]:
    """
    Extract allele names from a ``Variant/Haplotypes`` cell.

    A single cell may contain multiple alleles separated by ``,``,
    ``;``, or ``" + "``. Commonly observed forms:

    - ``"HLA-B*57:01"``   → ``"B*57:01"``
    - ``"B*57:01"``       → ``"B*57:01"``
    - ``"*15:02"`` + ``gene="HLA-B"`` → ``"B*15:02"``

    Parameters
    ----------
    cell : str
        Raw cell value.
    gene : str, optional
        ``Gene`` column from the same row, used to expand gene-less
        allele tokens (``"*15:02"``) back into a gene-prefixed form.
    """
    if not cell:
        return []
    for sep in (" + ", "+", ",", ";"):
        cell = cell.replace(sep, "|")
    tokens = [t.strip() for t in cell.split("|") if t.strip()]
    bare_gene = _strip_hla(gene) if gene else None  # "HLA-B" → "B"
    seen: List[str] = []
    for tok in tokens:
        cleaned = _strip_hla(tok)
        # Expand forms like "*15:02" using the Gene column.
        if cleaned.startswith("*") and bare_gene:
            cleaned = f"{bare_gene}{cleaned}"
        if cleaned and cleaned not in seen:
            seen.append(cleaned)
    return seen


def _split_pmids(cell: Optional[str]) -> List[str]:
    if not cell:
        return []
    out: List[str] = []
    for part in cell.replace(";", ",").split(","):
        pid = part.strip()
        if pid and pid.lower() != "none":
            out.append(pid)
    return out


__all__ = [
    "PharmGKBClient",
    "PharmAnnotation",
    "PharmGKBDatabaseError",
    "PharmGKBDownloadError",
    "PHARMGKB_CLINICAL_ANN_URL",
    "DEFAULT_EVIDENCE_LEVELS",
    "DEFAULT_LOCAL_DIR",
    "CLINICAL_ANN_FILENAME",
]
