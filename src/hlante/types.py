"""
hlante.types
============

Shared types used across the HLAnte pipeline.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import List, Optional


class InputSource(str, Enum):
    """Provenance of HLA allele calls passed to HLAnte.

    Different input sources carry different intrinsic uncertainty.
    This enum drives confidence-score adjustment in the annotator:

    - Resolution penalties always apply (a 2-field allele is less specific
      than a 4-field allele regardless of how it was obtained).
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
    # Resolution penalties still apply (a 2-field validated allele is
    # still less specific than a 4-field validated allele), but the
    # ambiguity penalty for "2-field = tool-ambiguous" is suppressed —
    # the call itself is not uncertain, only the resolution is coarse.

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
