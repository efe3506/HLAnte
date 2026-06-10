"""
hlante.db.curated
=================

Built-in curated HLA–disease association table.

Contains ~40 well-established allele–disease associations derived from
peer-reviewed literature (CPIC guidelines, PharmGKB Level 1A/1B, large
GWAS meta-analyses).

The table is bundled with the package as ``curated_hla_disease.tsv``
and requires no download or internet access.  Entries carry a
``review_status`` of ``"curated (HLAnte built-in)"`` so downstream
consumers can identify the provenance.

Coverage
--------
- CPIC Level 1A pharmacogenomics (B*57:01, B*58:01, B*15:02, A*31:01 …)
- Autoimmune Class I  (B*27:05 AS, B*51:01 Behçet, C*06:02 psoriasis …)
- Autoimmune Class II (DRB1*15:01 MS, DQB1*06:02 narcolepsy, DRB1*04
  RA/T1D, DQB1*02:01 celiac …)

Allele resolution fallback: 8 → 6 → 4 → 2 field.
"""

from __future__ import annotations

import hashlib
import logging
from pathlib import Path
from typing import List, Optional

import pandas as pd

from hlante.types import DiseaseEntry

logger: logging.Logger = logging.getLogger(__name__)


def _stable_curated_id(allele: str, condition: str) -> str:
    """
    Deterministic identifier for a curated (allele, condition) pair.

    Uses a stable hash (SHA-1) rather than the built-in ``hash()``, which
    is salted per process (``PYTHONHASHSEED``) and would otherwise yield a
    different ``variation_id`` on every run, breaking output reproducibility
    and any downstream join on the id.
    """
    digest = hashlib.sha1(f"{allele}|{condition}".encode("utf-8")).hexdigest()
    return f"CURATED_{int(digest, 16) % 10**6:06d}"

BUILTIN_CURATED_TSV: Path = Path(__file__).parent / "curated_hla_disease.tsv"

# Marker used in review_status so reporter can identify source
CURATED_REVIEW_STATUS: str = "curated (HLAnte built-in)"


class CuratedDiseaseClient:
    """
    Built-in curated HLA–disease association client.

    Reads :data:`BUILTIN_CURATED_TSV` and returns
    :class:`~hlante.types.DiseaseEntry` records.

    Parameters
    ----------
    tsv_path : Path, optional
        Override the default built-in TSV path.

    Notes
    -----
    ``load()`` is idempotent.  The client never raises; on any loading
    failure it silently returns empty results.
    """

    def __init__(self, tsv_path: Optional[Path] = None) -> None:
        self._tsv_path: Path = tsv_path or BUILTIN_CURATED_TSV
        self._df: Optional[pd.DataFrame] = None
        self._loaded: bool = False

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def load(self) -> None:
        """Load the curated TSV (idempotent)."""
        if self._loaded:
            return
        if not self._tsv_path.is_file():
            logger.warning("Curated HLA-disease TSV not found: %s", self._tsv_path)
            self._df = pd.DataFrame()
            self._loaded = True
            return
        try:
            self._df = pd.read_csv(self._tsv_path, sep="\t", comment="#", dtype=str).fillna("")
            logger.info(
                "Curated HLA-disease: %d entries loaded from %s",
                len(self._df),
                self._tsv_path,
            )
        except Exception as exc:  # noqa: BLE001
            logger.warning(
                "Curated HLA-disease: failed to load (%s): %s",
                self._tsv_path,
                exc,
            )
            self._df = pd.DataFrame()
        self._loaded = True

    def query_allele(self, allele: str) -> List[DiseaseEntry]:
        """
        Return curated entries for *allele*, with 8→6→4→2-field fallback.

        Parameters
        ----------
        allele : str
            HLA allele — ``HLA-`` prefix is optional.

        Returns
        -------
        list of DiseaseEntry
            Empty list when no match is found.
        """
        if not self._loaded:
            self.load()
        if self._df is None or self._df.empty:
            return []

        bare = allele[4:] if allele.upper().startswith("HLA-") else allele

        for candidate in _resolution_cascade(bare):
            hits = self._lookup(candidate)
            if hits:
                return hits
        return []

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _lookup(self, allele: str) -> List[DiseaseEntry]:
        """Exact case-insensitive match against the Allele column."""
        assert self._df is not None
        mask = self._df["Allele"].str.strip().str.upper() == allele.upper()
        results: List[DiseaseEntry] = []
        for _, row in self._df[mask].iterrows():
            pmids = [p.strip() for p in str(row.get("PMID", "")).split(";") if p.strip()]
            or_val = row.get("OR", "")
            condition = str(row.get("Disease", "")).strip()
            sig = str(row.get("ClinicalSignificance", "risk factor")).strip().lower()
            evidence = str(row.get("Evidence", "")).strip()
            pop = str(row.get("Population", "")).strip()
            citation = str(row.get("Citation", "")).strip()

            # Embed OR and population into condition text so reporter
            # can surface it without schema changes
            detail_parts: List[str] = []
            if or_val:
                detail_parts.append(f"OR={or_val}")
            if pop:
                detail_parts.append(pop)
            if evidence:
                detail_parts.append(f"{evidence} evidence")
            if citation:
                detail_parts.append(citation)
            detail = f" [{'; '.join(detail_parts)}]" if detail_parts else ""

            results.append(
                DiseaseEntry(
                    variation_id=_stable_curated_id(allele, condition),
                    significance=sig,
                    condition=f"{condition}{detail}",
                    review_status=CURATED_REVIEW_STATUS,
                    allele=allele,
                    pmid=pmids,
                )
            )
        return results


# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------


def _resolution_cascade(bare: str) -> List[str]:
    """
    Return allele forms from highest to lowest resolution.

    ``"B*57:01:01:01"`` → ``["B*57:01:01:01", "B*57:01:01", "B*57:01", "B*57"]``
    """
    if "*" not in bare:
        return [bare]
    gene, fields_str = bare.split("*", 1)
    parts = fields_str.split(":")
    return [f"{gene}*{':'.join(parts[:n])}" for n in range(len(parts), 0, -1)]


__all__ = [
    "CuratedDiseaseClient",
    "BUILTIN_CURATED_TSV",
    "CURATED_REVIEW_STATUS",
]
