"""
hlante.annotator
===============

Annotation engine that orchestrates GWAS / PharmGKB / AFND
queries for a list of :class:`NormalizedAllele` records.

For every input allele the engine queries the configured database
clients and produces a unified :class:`AnnotatedHLA` record that
includes a disease-risk summary, drug-response summary, overall
clinical significance, and a confidence score.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, FrozenSet, List, Optional, Sequence, Tuple, cast

from hlante.db.afnd import (
    AFNDClient,
    AllelFrequency,
    DEFAULT_POPULATION_GROUP as AFND_DEFAULT_POPULATION,
)
from hlante.db.nmdp import NMDPClient
from hlante.db.curated import CuratedDiseaseClient
from hlante.types import DiseaseEntry
from hlante.db.gwas import (
    DEFAULT_P_VALUE_THRESHOLD,
    GWASClient,
    GWASHit,
    RESOLUTION_LABEL_NONE,
)
from hlante.db.pharmgkb import (
    DEFAULT_EVIDENCE_LEVELS,
    PharmAnnotation,
    PharmGKBClient,
)
from hlante.normalizer import NormalizedAllele
from hlante.types import InputSource

# Association-strength prefixes and no-hit placeholders. These strings
# are also re-exported from ``hlante.reporter``; any change must be
# mirrored there so TSV / Markdown / JSON consumers stay in sync.
# P0-8: values were reworded from risk-magnitude language
# ("High risk") to association-strength language ("Strong association").
RISK_PREFIX_HIGH: str = "Strong association"
RISK_PREFIX_MODERATE: str = "Moderate association"
RISK_PREFIX_PROTECTIVE: str = "Inverse association"
RISK_PREFIX_ASSOCIATION: str = "Reported association"
NO_DISEASE_SUMMARY: str = "No disease association reported"
NO_DRUG_SUMMARY: str = "No drug response reported"

logger: logging.Logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Exceptions
# ---------------------------------------------------------------------------


class HLAAnnotationError(Exception):
    """
    Base class for errors raised during annotation.
    """

    pass


# ---------------------------------------------------------------------------
# Evidence-strength labels (P0-8)
# ---------------------------------------------------------------------------
#
# HLAnte does NOT implement ACMG/AMP criteria. The labels below are
# evidence-strength descriptors, not clinical classifications. The
# constant names (``SIGNIFICANCE_PATHOGENIC`` etc.) are preserved for
# backward compatibility with imports; the string values were revised
# to remove ACMG-adjacent language
# (``Pathogenic`` / ``Likely Pathogenic`` / ``VUS`` / ``Benign`` /
# ``Novel``) and replace them with neutral evidence-strength phrasing.


SIGNIFICANCE_NOVEL: str = "Not in IMGT"
SIGNIFICANCE_NULL_ALLELE: str = "Null allele (not expressed) — risk not assessed"
# NOTE: these are evidence-strength / actionability descriptors, deliberately
# phrased to AVOID ACMG/AMP variant-classification vocabulary ("pathogenic",
# "likely pathogenic", "VUS"), which does not apply to statistical HLA
# disease/PGx associations.
SIGNIFICANCE_PATHOGENIC: str = "Actionable pharmacogenomic risk (CPIC 1A — avoid)"
SIGNIFICANCE_LIKELY_PATHOGENIC: str = "Strong pharmacogenomic risk association"
SIGNIFICANCE_RISK_FACTOR: str = "Suggestive risk factor"
SIGNIFICANCE_BENIGN: str = "No reported risk"
SIGNIFICANCE_BENIGN_LIMITED: str = "Not assessed — insufficient coverage"
SIGNIFICANCE_VUS: str = "Inconclusive evidence"


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------


@dataclass
class AnnotatorConfig:
    """
    Runtime configuration for :func:`annotate_genotype`.

    Attributes
    ----------
    offline : bool
        If ``True``, no database HTTP calls are made; only local
        caches / bulk dumps are consulted.
    cache_root : Path, optional
        Root directory for per-database caches.
    pharmgkb_local_dir : Path, optional
        Directory of the extracted PharmGKB bulk dump.
    gwas_local_dir : Path, optional
        Directory of the extracted GWAS Catalog bulk dump.
    afnd_local_dir : Path, optional
        Directory of the AFND allele-frequency TSV.
    imgt_db_path : Path, optional
        Local IPD-IMGT/HLA directory.
    enable_gwas : bool
        Toggle GWAS Catalog queries.
    enable_pharmgkb : bool
        Toggle PharmGKB queries.
    enable_afnd : bool
        Toggle AFND lookups (used for confidence scoring).
    gwas_p_threshold : float
        p-value threshold for GWAS hits.
    pharmgkb_evidence_levels : frozenset of str
        Accepted PharmGKB evidence levels.
    population_group : str
        AFND population group code (``EUR``/``AFR``/``EAS``/``SAS``/
        ``MID``/``AMR``/``OCE``/``global``; ``ASN`` is an alias for
        ``EAS``) or custom substring.
    afnd_min_sample_size : int
        Minimum per-study sample size accepted by AFND aggregation.
    ncbi_api_key : str, optional
        NCBI Entrez API key.
    """

    offline: bool = False
    cache_root: Optional[Path] = None
    pharmgkb_local_dir: Optional[Path] = None
    gwas_local_dir: Optional[Path] = None
    afnd_local_dir: Optional[Path] = None
    imgt_db_path: Optional[Path] = None
    enable_gwas: bool = True
    enable_pharmgkb: bool = True
    enable_afnd: bool = True
    gwas_p_threshold: float = DEFAULT_P_VALUE_THRESHOLD
    pharmgkb_evidence_levels: FrozenSet[str] = field(
        default_factory=lambda: frozenset(DEFAULT_EVIDENCE_LEVELS)
    )
    population_group: str = AFND_DEFAULT_POPULATION
    afnd_min_sample_size: int = 50
    ncbi_api_key: Optional[str] = None
    input_source: InputSource = InputSource.TYPING_TOOL
    curated_tsv_path: Optional[Path] = None


@dataclass
class AnnotatedHLA:
    """
    Unified annotation record for a single :class:`NormalizedAllele`.

    Attributes
    ----------
    normalized_allele : NormalizedAllele
        Source normalized allele.
    gwas_hits : list of GWASHit
        GWAS Catalog hits.
    pharm_annotations : list of PharmAnnotation
        PharmGKB clinical annotations.
    disease_entries : list of DiseaseEntry
        Curated built-in disease association records.
    disease_risk_summary : str
        Human-readable disease-risk summary (e.g., ``"High risk: RA (OR=4.2)"``).
    drug_response_summary : str
        Human-readable drug-response summary.
    clinical_significance : str
        Overall evidence-strength label — one of the ``SIGNIFICANCE_*``
        constants. These are evidence-strength descriptors, NOT
        ACMG/AMP classifications. See the module-level constant block
        for definitions.
    gwas_resolution_used : str
        Resolution tier at which GWAS hits were found (``"2-field"``,
        ``"4-field"``, ``"6-field"``, ``"8-field"``) or ``"none"``.
    allele_frequency : float or None
        Population-specific frequency from AFND.
    frequency_population : str or None
        Label of the population/group used for the frequency lookup.
    frequency_sample_size : int or None
        Aggregate sample size behind the frequency estimate.
    frequency_is_estimated : bool
        ``True`` when fallback resolution was required to find a frequency.
    confidence_score : float
        Confidence score in ``[0.0, 1.0]``. Lowered for novel, rare,
        low-resolution, or ambiguous calls.
    confidence_rationale : str
        Pipe-delimited reason codes explaining the score
        (e.g., ``"rare_allele(0.0003)|novel_allele"``).
    """

    normalized_allele: NormalizedAllele
    gwas_hits: List[GWASHit]
    pharm_annotations: List[PharmAnnotation]
    disease_entries: List[DiseaseEntry]
    disease_risk_summary: str
    drug_response_summary: str
    clinical_significance: str
    gwas_resolution_used: str = RESOLUTION_LABEL_NONE
    allele_frequency: Optional[float] = None
    frequency_population: Optional[str] = None
    frequency_sample_size: Optional[int] = None
    frequency_is_estimated: bool = False
    confidence_score: float = 1.0
    confidence_rationale: str = "standard"
    confidence_tier: str = "HIGH"


# ---------------------------------------------------------------------------
# Client bundle
# ---------------------------------------------------------------------------


@dataclass
class AnnotatorClients:
    """
    Bundle of database clients used by :func:`annotate_genotype`.

    Tests may inject stub clients; :func:`annotate_genotype` accepts
    this bundle optionally.

    Attributes
    ----------
    gwas : GWASClient, optional
    pharmgkb : PharmGKBClient, optional
    afnd : AFNDClient, optional
    nmdp : NMDPClient, optional
    """

    gwas: Optional[GWASClient] = None
    pharmgkb: Optional[PharmGKBClient] = None
    afnd: Optional[AFNDClient] = None
    nmdp: Optional[NMDPClient] = None
    curated: Optional[CuratedDiseaseClient] = None


def build_clients(config: AnnotatorConfig) -> AnnotatorClients:
    """
    Build database clients from :class:`AnnotatorConfig`.

    Parameters
    ----------
    config : AnnotatorConfig
        Runtime configuration.

    Returns
    -------
    AnnotatorClients
        Instantiated client bundle.
    """
    gwas_client: Optional[GWASClient] = None
    pharm_client: Optional[PharmGKBClient] = None
    afnd_client: Optional[AFNDClient] = None
    nmdp_client: Optional[NMDPClient] = None

    if config.enable_gwas:
        gwas_client = GWASClient(
            local_dir=getattr(config, "gwas_local_dir", None),
            offline=config.offline,
            p_value_threshold=config.gwas_p_threshold,
        )

    if config.enable_pharmgkb:
        pharm_client = PharmGKBClient(
            local_dir=config.pharmgkb_local_dir,
            offline=config.offline,
            evidence_levels=config.pharmgkb_evidence_levels,
        )

    if config.enable_afnd:
        afnd_client = AFNDClient(
            local_dir=config.afnd_local_dir,
            population_group=config.population_group,
            min_sample_size=config.afnd_min_sample_size,
        )
        nmdp_client = NMDPClient(
            population_group=config.population_group,
            min_sample_size=config.afnd_min_sample_size,
        )

    curated_client = CuratedDiseaseClient(tsv_path=config.curated_tsv_path)
    curated_client.load()

    return AnnotatorClients(
        gwas=gwas_client,
        pharmgkb=pharm_client,
        afnd=afnd_client,
        nmdp=nmdp_client,
        curated=curated_client,
    )


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------


def annotate_genotype(
    normalized: Sequence[NormalizedAllele],
    config: AnnotatorConfig,
    *,
    clients: Optional[AnnotatorClients] = None,
) -> List[AnnotatedHLA]:
    """
    Produce clinical annotations for a list of normalized alleles.

    Parameters
    ----------
    normalized : sequence of NormalizedAllele
        Alleles to annotate.
    config : AnnotatorConfig
        Runtime configuration.
    clients : AnnotatorClients, optional
        Pre-built client bundle. If ``None``, :func:`build_clients`
        constructs one. Used primarily by tests.

    Returns
    -------
    list of AnnotatedHLA
        One annotation record per input allele.
    """
    if clients is None:
        clients = build_clients(config)

    if config.input_source == InputSource.SIMULATED:
        logger.warning(
            "Input source is SIMULATED — penalties as TYPING_TOOL. "
            "Scores do not reflect real biological data."
        )
    elif config.input_source == InputSource.UNKNOWN:
        logger.warning(
            "Input source is UNKNOWN — treating as TYPING_TOOL. "
            "Use --input-source to specify provenance explicitly."
        )

    results: List[AnnotatedHLA] = []
    for allele in normalized:
        gwas_hits: List[GWASHit]
        pharm_anns: List[PharmAnnotation]
        disease_entries: List[DiseaseEntry]
        if allele.is_null:
            # Null (non-expressed) allele: the antigen is absent from the
            # cell surface, so the surface-expression-dependent HLA disease
            # and drug-hypersensitivity associations do not apply. Suppress
            # those lookups to avoid emitting a false risk alert. The allele
            # frequency is still meaningful and is retained below.
            gwas_hits, gwas_resolution = [], RESOLUTION_LABEL_NONE
            pharm_anns = []
            disease_entries = []
        else:
            gwas_hits, gwas_resolution = _query_gwas_with_fallback(
                clients.gwas, allele.allele_name
            )
            pharm_anns = (
                _safe_query(
                    clients.pharmgkb,
                    lambda c: c.query_allele(allele.allele_name),
                    context=f"PharmGKB {allele.allele_name}",
                )
                or []
            )

            # Curated built-in disease associations
            disease_entries = (
                _safe_query(
                    clients.curated,
                    lambda c: c.query_allele(allele.allele_name),
                    context=f"Curated {allele.allele_name}",
                )
                or []
            )

        freq = _safe_query(
            clients.afnd,
            lambda c: c.get_frequency_with_fallback(allele.allele_name),
            context=f"AFND {allele.allele_name}",
        )
        if freq is None:
            freq = _safe_query(
                clients.nmdp,
                lambda c: c.get_frequency_with_fallback(allele.allele_name),
                context=f"NMDP {allele.allele_name}",
            )

        # Source / PMID is required — filter missing
        gwas_hits = [h for h in gwas_hits if h.pmid]
        pharm_anns = [p for p in pharm_anns if p.pmid]

        confidence, rationale = _compute_confidence_score(allele, freq, config.input_source)
        tier = _confidence_tier(confidence)

        if allele.is_null:
            clinical_significance = SIGNIFICANCE_NULL_ALLELE
            disease_risk_summary = (
                "Null allele (not expressed); surface-expression-dependent "
                "HLA disease associations not applicable."
            )
            drug_response_summary = (
                "Null allele (not expressed); HLA drug-hypersensitivity "
                "associations not applicable."
            )
            # A non-expressed allele must never be presented at HIGH
            # confidence, regardless of its frequency/resolution typing score.
            if tier == "HIGH":
                tier = "MODERATE"
        else:
            clinical_significance = _classify_significance(
                allele,
                gwas_hits,
                pharm_anns,
                disease_entries,
                gwas_resolution=gwas_resolution,
                p_threshold=config.gwas_p_threshold,
            )
            disease_risk_summary = _build_disease_risk_summary(gwas_hits, disease_entries)
            drug_response_summary = _build_drug_response_summary(pharm_anns)

        results.append(
            AnnotatedHLA(
                normalized_allele=allele,
                gwas_hits=gwas_hits,
                pharm_annotations=pharm_anns,
                disease_entries=disease_entries,
                disease_risk_summary=disease_risk_summary,
                drug_response_summary=drug_response_summary,
                clinical_significance=clinical_significance,
                gwas_resolution_used=(gwas_resolution if gwas_hits else RESOLUTION_LABEL_NONE),
                allele_frequency=freq.frequency if freq else None,
                frequency_population=(f"{freq.population} ({freq.source})" if freq else None),
                frequency_sample_size=freq.sample_size if freq else None,
                frequency_is_estimated=freq.is_estimated if freq else False,
                confidence_score=confidence,
                confidence_rationale=rationale,
                confidence_tier=tier,
            )
        )
    return results


def _compute_confidence_score(
    normalized: NormalizedAllele,
    freq: Optional[AllelFrequency],
    input_source: InputSource = InputSource.TYPING_TOOL,
) -> Tuple[float, str]:
    """
    Compute ``(confidence_score, rationale)`` for an allele.

    Deterministic, pure function. The score starts at 1.0 and is
    multiplied by stated penalties that reflect the annotation
    uncertainty introduced by each signal. Penalty values align with
    published HLA-typing concordance gaps (arcasHLA/HLA-HD ~99 % at
    2-field vs ~92 % at 4-field; Illing et al. 2022 Frontiers Immunol)
    and AFND rarity classifications. Because HLAnte is a research
    annotation aid — not a clinical decision tool — these criteria are
    explicitly declared rather than empirically calibrated against a
    gold-standard cohort.

    ================================== ========================= ============
    Signal                              Multiplier                Applies to
    ================================== ========================= ============
    Novel allele (not in IMGT)          × 0.30                   all sources
    Rare allele (freq < 0.001)          × 0.50                   all sources
    Uncommon (0.001 ≤ freq < 0.01)      × 0.80                   all sources
    Frequency unknown (no AFND/NMDP)    × 0.85                   all sources
    Low resolution (2-field)            × 0.70                   all sources
    Medium resolution (4-field)         × 0.90                   all sources
    Ambiguous (is_ambiguous)            × 0.75                   TYPING_TOOL,
                                                                 SIMULATED,
                                                                 UNKNOWN only
    ================================== ========================= ============

    Ambiguity penalty and source-awareness
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    For ``TYPING_TOOL`` / ``SIMULATED`` / ``UNKNOWN`` inputs the
    ambiguity penalty (×0.75) reflects the tool's inability to
    discriminate sub-allele variants from read data.

    For ``VALIDATED`` inputs (Sanger-sequenced, PCR-SBT, or IHIW
    reference panels such as the 1000 Genomes Project HLA types) the
    call is exactly correct at its reported resolution — only the
    resolution penalty applies; the ambiguity penalty is suppressed.

    Parameters
    ----------
    normalized : NormalizedAllele
        Allele to score.
    freq : AllelFrequency or None
        AFND frequency (``None`` if unavailable).
    input_source : InputSource
        Provenance of the allele call.  Defaults to
        :attr:`~hlante.types.InputSource.TYPING_TOOL`.

    Returns
    -------
    (score, rationale) : tuple of (float, str)
        ``score`` rounded to 4 decimals. ``rationale`` is a pipe-joined
        reason code list (``"standard"`` when no penalty applies).
    """
    score: float = 1.0
    reasons: List[str] = []

    # 1. IMGT recognition
    if normalized.is_novel:
        score *= 0.3
        reasons.append("novel_allele")

    # 2. Allele frequency
    if freq is not None:
        if freq.frequency < 0.001:
            score *= 0.5
            reasons.append(f"rare_allele(freq={freq.frequency:.4f})")
        elif freq.frequency < 0.01:
            score *= 0.8
            reasons.append(f"uncommon_allele(freq={freq.frequency:.4f})")
    else:
        score *= 0.85
        reasons.append("freq_unknown")

    # 3. Resolution
    if normalized.resolution_level == 2:
        score *= 0.7
        reasons.append("low_resolution(2-field)")
    elif normalized.resolution_level == 4:
        score *= 0.9
        reasons.append("medium_resolution(4-field)")

    # 4. Ambiguous allele — penalty depends on input source
    if normalized.is_ambiguous:
        if input_source == InputSource.VALIDATED:
            # Validated calls are exactly correct at their reported resolution.
            # The resolution penalty already captures the loss of specificity;
            # the ambiguity penalty (which reflects tool uncertainty) does not
            # apply.
            reasons.append("ambiguity_suppressed(validated_source)")
        else:
            score *= 0.75
            reasons.append("ambiguous")

    # 5. Expression status (IPD-IMGT/HLA suffix). A null allele is not
    # expressed; its disease/drug annotations are suppressed upstream and it
    # is never shown at HIGH tier, so the rationale records the status
    # without an additional numeric penalty (the call itself may be correct).
    # Reduced/aberrant-expression alleles (L/S/C/A/Q) are annotated but
    # down-weighted, as the clinical relevance of their altered expression is
    # uncertain.
    if normalized.is_null:
        reasons.append("null_allele(not_expressed)")
    elif normalized.is_low_or_aberrant_expression:
        score *= 0.85
        reasons.append(f"low_or_aberrant_expression({normalized.expression_suffix})")

    rationale = "|".join(reasons) if reasons else "standard"
    return round(score, 4), rationale


def _confidence_tier(score: Optional[float]) -> str:
    """
    Convert a numeric confidence score to an interpretable tier.

    Parameters
    ----------
    score : float or None

    Returns
    -------
    str
        ``"HIGH"`` (score ≥ 0.85), ``"MODERATE"`` (0.70 ≤ score < 0.85),
        ``"LOW"`` (score < 0.70), or ``"NA"`` (score is ``None``).

    Notes
    -----
    Tier thresholds are stated criteria for this research tool:
    HIGH (≥ 0.85) — well-characterised allele, ≥ 4-field, known frequency;
    MODERATE (0.70–0.85) — minor uncertainty (ambiguous or freq unknown);
    LOW (< 0.70) — substantial uncertainty (2-field, novel, or very rare).
    """
    if score is None:
        return "NA"
    if score >= 0.85:
        return "HIGH"
    if score >= 0.70:
        return "MODERATE"
    return "LOW"


def _query_gwas_with_fallback(
    client: Any,
    allele_name: str,
) -> Tuple[List[Any], str]:
    """
    Query GWAS with fallback support. Safely handles older clients
    (only ``query_allele``) and any exceptions from the client.
    """
    if client is None:
        return [], RESOLUTION_LABEL_NONE
    try:
        fallback = getattr(client, "query_allele_with_fallback", None)
        if callable(fallback):
            return cast(Tuple[List[Any], str], fallback(allele_name))
        hits = client.query_allele(allele_name)
        return hits, ("4-field" if hits else RESOLUTION_LABEL_NONE)
    except Exception as exc:  # noqa: BLE001 — DB exception swallowing
        logger.warning("GWAS query for %s failed: %s", allele_name, exc)
        return [], RESOLUTION_LABEL_NONE


def _safe_query(
    client: Any,
    runner: Any,
    *,
    context: str,
) -> Optional[Any]:
    """
    Wrapper that returns ``None`` when the client is absent or raises,
    logging the failure.
    """
    if client is None:
        return None
    try:
        return runner(client)
    except Exception as exc:  # noqa: BLE001 — DB exception swallowing
        logger.warning("%s query failed: %s", context, exc)
        return None


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _condition_key(condition: str) -> str:
    """Normalise a condition string for duplicate detection.

    Strips bracketed context and comma sub-qualifiers so that semantically
    identical disease names from different sources compare equal:
      'Multiple sclerosis, susceptibility to 1'  →  'multiple sclerosis'
      'Multiple sclerosis [OR=3.1; ...]'          →  'multiple sclerosis'
    """
    return condition.split("[")[0].split(",")[0].strip().lower()


# ---------------------------------------------------------------------------
# Summary generation
# ---------------------------------------------------------------------------


def _build_disease_risk_summary(
    gwas_hits: Sequence[GWASHit],
    disease_entries: Sequence[DiseaseEntry],
) -> str:
    """
    Turn the lowest-p-value GWAS hit and the leading curated disease entry
    into a one-line summary.
    """
    segments: List[str] = []

    if gwas_hits:
        top = min(
            gwas_hits,
            key=lambda h: h.p_value if h.p_value is not None else 1.0,
        )
        severity = _severity_label(top.odds_ratio)
        or_part = f" (OR={top.odds_ratio:.2f})" if top.odds_ratio else ""
        segments.append(f"{severity}: {top.trait}{or_part}")

    if disease_entries:
        cv = disease_entries[0]
        condition = cv.condition or "disease association"
        segments.append(f"Curated {cv.significance}: {condition}")

    if not segments:
        return NO_DISEASE_SUMMARY
    return "; ".join(segments)


def _build_drug_response_summary(
    pharm_annotations: Sequence[PharmAnnotation],
) -> str:
    """
    Summarize the strongest-evidence PharmGKB annotation.
    """
    if not pharm_annotations:
        return NO_DRUG_SUMMARY
    ordered = sorted(
        pharm_annotations,
        key=lambda p: _evidence_rank(p.evidence_level),
    )
    top = ordered[0]
    phenotype = top.phenotype or "reaction"
    level = top.evidence_level or ""
    suffix = " (low evidence)" if level.upper() in {"3", "4"} else ""
    return f"{top.drug} {phenotype}: {level} evidence{suffix}"


def _severity_label(odds_ratio: Optional[float]) -> str:
    # Direction of effect is decided first: any OR < 1.0 is an inverse
    # (protective) association, not a risk association. The
    # previous threshold (OR <= 0.5) mislabeled weakly protective alleles
    # (0.5 < OR < 1.0) as a plain "reported association".
    if odds_ratio is None:
        return RISK_PREFIX_ASSOCIATION
    if odds_ratio < 1.0:
        return RISK_PREFIX_PROTECTIVE
    if odds_ratio >= 3.0:
        return RISK_PREFIX_HIGH
    if odds_ratio >= 1.5:
        return RISK_PREFIX_MODERATE
    return RISK_PREFIX_ASSOCIATION


_EVIDENCE_ORDER: Dict[str, int] = {
    "1A": 0,
    "1B": 1,
    "2A": 2,
    "2B": 3,
    "3": 4,
    "4": 5,
}


def _evidence_rank(level: str) -> int:
    return _EVIDENCE_ORDER.get((level or "").upper(), 99)


# ---------------------------------------------------------------------------
# Clinical-significance classification
# ---------------------------------------------------------------------------


def _classify_significance(
    allele: NormalizedAllele,
    gwas_hits: Sequence[GWASHit],
    pharm_annotations: Sequence[PharmAnnotation],
    disease_entries: Sequence[DiseaseEntry],
    gwas_resolution: str = RESOLUTION_LABEL_NONE,
    p_threshold: float = DEFAULT_P_VALUE_THRESHOLD,
) -> str:
    """
    Assign an overall clinical-significance label for a single allele.

    Notes
    -----
    P0-4: ``Benign`` is now split. The plain :data:`SIGNIFICANCE_BENIGN`
    label is reserved for alleles where *at least one database query
    succeeded* (GWAS resolution not ``none``, or any pharm / curated
    disease record) and no signal of risk was found. When the allele is known
    in IPD-IMGT/HLA but no database query returned *any* result (offline
    mode, novel locus, DB sparse for this gene), the more cautious
    :data:`SIGNIFICANCE_BENIGN_LIMITED` label is used instead.
    """
    if allele.is_novel:
        return SIGNIFICANCE_NOVEL

    cv_labels = [e.significance.strip().lower() for e in disease_entries]
    if any(label == "pathogenic" for label in cv_labels):
        return SIGNIFICANCE_PATHOGENIC
    if any(label == "likely pathogenic" for label in cv_labels):
        return SIGNIFICANCE_LIKELY_PATHOGENIC

    has_strong_pharm = any(
        (p.evidence_level or "").upper() in {"1A", "1B"} for p in pharm_annotations
    )
    has_significant_gwas = any(
        (h.p_value is not None and h.p_value <= p_threshold) for h in gwas_hits
    )
    if has_strong_pharm or has_significant_gwas:
        return SIGNIFICANCE_RISK_FACTOR

    if gwas_hits or pharm_annotations:
        return SIGNIFICANCE_VUS

    # At this point: no hits from any source.
    # Determine whether at least one DB was actually queried + returned
    # something (even if curated entries with non-qualifying significance)
    # vs. every DB came back empty.
    any_db_queried = (
        gwas_resolution not in (RESOLUTION_LABEL_NONE, "none")
        or bool(pharm_annotations)
        or bool(disease_entries)
    )

    if allele.imgt_accession is not None and not allele.is_ambiguous:
        if any_db_queried:
            return SIGNIFICANCE_BENIGN
        return SIGNIFICANCE_BENIGN_LIMITED
    return SIGNIFICANCE_VUS


__all__ = [
    "AnnotatedHLA",
    "AnnotatorClients",
    "AnnotatorConfig",
    "HLAAnnotationError",
    "InputSource",
    "SIGNIFICANCE_NOVEL",
    "SIGNIFICANCE_NULL_ALLELE",
    "SIGNIFICANCE_PATHOGENIC",
    "SIGNIFICANCE_LIKELY_PATHOGENIC",
    "SIGNIFICANCE_RISK_FACTOR",
    "SIGNIFICANCE_BENIGN",
    "SIGNIFICANCE_BENIGN_LIMITED",
    "SIGNIFICANCE_VUS",
    "annotate_genotype",
    "build_clients",
    "_compute_confidence_score",
]
