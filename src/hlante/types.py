"""
hlante.types
============

Shared types used across the HLAnte pipeline.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Iterable, List, Optional, Tuple


class InputSource(str, Enum):
    """Provenance of HLA allele calls passed to HLAnte.

    Different input sources carry different intrinsic uncertainty.
    This enum drives input-quality-score adjustment in the annotator:

    - Resolution penalties always apply (a two-field allele is less specific
      than a four-field allele regardless of how it was obtained).
    - The ambiguity penalty (×0.75) reflects tool-level inability to
      discriminate sub-allele variants from read data; it does NOT apply
      to lab-validated calls that are exactly correct at their reported
      resolution.
    """

    TYPING_TOOL = "typing_tool"
    # Default. Output from arcasHLA / T1K / HLA-HD / OptiType.
    # Resolution and ambiguity penalties apply normally.

    VALIDATED = "validated"
    # Sanger-sequenced, PCR-SBT, or otherwise lab-validated alleles
    # (e.g. 1000 Genomes Project HLA typing, IHIW reference panels).
    # Resolution penalties still apply (a two-field validated allele is
    # Still less specific than a four-field validated allele), but the
    # Ambiguity penalty for "two-field = tool-ambiguous" is suppressed —
    # The call itself is not uncertain, only the resolution is coarse.

    SIMULATED = "simulated"
    # Synthetic alleles for testing/benchmarking.
    # Penalties as TYPING_TOOL but a warning is logged.

    UNKNOWN = "unknown"
    # User did not specify; treated as TYPING_TOOL with a warning.


@dataclass
class DiseaseEntry:
    """
    A curated or inferred disease association for an HLA allele.

    Attributes
    ----------
    variation_id : str
        Source record identifier (e.g. ``"CURATED_001234"``).
    significance : str
        Association label (e.g. ``"pathogenic"``, ``"risk factor"``).
    condition : str
        Disease or phenotype name.
    review_status : str
        Provenance / review level (e.g. ``"curated (HLAnte built-in)"``).
    allele : str
        HLA allele used during the lookup.
    pmid : list of str
        Supporting PMIDs (may be empty).
    accession : str, optional
        External record accession.
    url : str, optional
        Source URL.
    """

    variation_id: str
    significance: str
    condition: str
    review_status: str
    allele: str
    pmid: List[str] = field(default_factory=list)
    accession: Optional[str] = None
    url: Optional[str] = None


__all__ = ["InputSource", "DiseaseEntry"]


#: Provenance tokens for the evidence layer that produced a statement.
#: Reported in the ``significance_basis`` column so that a guideline-grade
#: pharmacogenomic recommendation, a database-derived association and the
#: authors' curated convenience table are never presented as equivalent.
LAYER_GUIDELINE_CPIC: str = "guideline:CPIC"
LAYER_DB_PHARMGKB: str = "database:PharmGKB"
LAYER_DB_GWAS: str = "database:GWAS"
LAYER_CURATED_BUILTIN: str = "curated:HLAnte"
LAYER_CURATED_USER: str = "curated:user"
LAYER_NONE: str = "none"

#: Emission order: strongest provenance first, curated transcription last.
LAYER_RANK: Tuple[str, ...] = (
    LAYER_GUIDELINE_CPIC,
    LAYER_DB_PHARMGKB,
    LAYER_DB_GWAS,
    LAYER_CURATED_USER,
    LAYER_CURATED_BUILTIN,
)


def order_layers(tokens: Iterable[str]) -> Tuple[str, ...]:
    """
    De-duplicate and order provenance tokens deterministically.
    """
    seen = {t for t in tokens if t}
    return tuple(t for t in LAYER_RANK if t in seen)

