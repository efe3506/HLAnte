"""
hlante.reporter
==============

Generate annotation outputs in three formats:

1. **TSV/CSV** — main annotation table with 22+ columns
   (machine- and human-readable).
2. **Markdown** — per-sample clinical summary report.
3. **JSON** — fully nested structure for API consumption.

Every output embeds the HLAnte version, run timestamp, reference
database versions, and a disclaimer indicating that the results are
intended for research use only.
"""

from __future__ import annotations

import csv
import json
import logging
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, FrozenSet, Iterable, List, Optional, Sequence, Tuple

from hlante import __version__
from hlante.annotator import AnnotatedHLA
from hlante.normalizer import NormalizedAllele
from hlante.parser import RESOLUTION_LABELS

logger: logging.Logger = logging.getLogger(__name__)


try:  # pragma: no cover - optional dependency
    from tqdm import tqdm as _tqdm
except ImportError:  # pragma: no cover

    def _tqdm(iterable: Any, **_kwargs: Any) -> Any:
        return iterable


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

NA: str = "NA"
PIPE: str = "|"
PROGRESS_THRESHOLD: int = 100

FLOAT_DECIMALS: int = 4

DISCLAIMER: str = (
    "RESEARCH USE ONLY. This report is an annotation aid. Nothing in "
    "this output constitutes a clinical diagnosis, medical advice, or "
    "pharmacogenomic recommendation. The significance labels used here "
    "are evidence-strength descriptors, NOT ACMG/AMP classifications. "
    "HLAnte does not implement ACMG/AMP criteria. The input-quality score / "
    "tier is an UNCALIBRATED heuristic describing how completely the allele "
    "CALL is supported by reference data (novelty, frequency, number of "
    "reported fields, ambiguity) — it is NOT a measure of genotype accuracy, "
    "NOT a posterior probability, and NOT the clinical certainty or "
    "correctness of any associated risk. A limited tier never down-weights "
    "an actionable association. Any clinical decision based on an HLA allele "
    "must rely on a certified laboratory result interpreted by a qualified "
    "clinician."
)

NO_DISEASE_SUMMARY: str = "No disease association reported"
NO_DRUG_SUMMARY: str = "No drug response reported"

#: Explicit limitation note for diplotype/zygosity biology that HLAnte does not
#: model. Surfaced in the Markdown note and JSON metadata so
#: the simplification is transparent rather than silent.
DIPLOTYPE_CAVEAT: str = (
    "Alleles are annotated independently. HLAnte does NOT apply "
    "compound-heterozygote rules (for example the "
    "HLA-DRB1*03:01+HLA-DRB1*04:01 genotype in type 1 diabetes), which "
    "require both alleles at a locus to be interpreted jointly but do not "
    "require chromosomal phase; nor does it model the HLA-DQ heterodimer "
    "(DQ2.5, DQ8) cis-versus-trans configuration across HLA-DQA1 and "
    "HLA-DQB1, which does require haplotypic phase that NGS typing tools do "
    "not routinely report. A locus reporting a single allele (allele2 = NA) "
    "is hemizygous or not fully reported, and is NOT assumed homozygous — "
    "homozygosity is asserted only when the same allele is reported twice. "
    "DRB3/4/5 are present-or-absent genes whose copy number is not inferred."
)

# Association-strength prefixes used in disease-risk summaries.
# Values were reworded from risk-magnitude language
# ("High risk", "Moderate risk") to association-strength language
# ("Strong association", "Moderate association") to avoid implying a
# Clinical magnitude of risk that the tool does not measure.
RISK_PREFIX_HIGH: str = "Strong association"
RISK_PREFIX_MODERATE: str = "Moderate association"
RISK_PREFIX_PROTECTIVE: str = "Inverse association"
RISK_PREFIX_ASSOCIATION: str = "Reported association"

# TSV column layout.
#
# The locus-level aggregate columns (``gwas_traits`` / ``pharm_drugs``)
# Pipe-join allele1 + allele2 hits. They are retained for backward
# Compatibility and machine-friendly bulk parsing. They cannot answer
# The question "which allele contributed which hit?".
#
# The explicit per-allele columns
# (``gwas_traits_allele1`` / ``gwas_traits_allele2``, and analogous
# ``pharm_drugs_*``) so the attribution is unambiguous without forcing
# The user to emit per-allele rows.
TSV_COLUMNS: Tuple[str, ...] = (
    "sample_id",
    "locus",
    "allele1",
    "allele2",
    "resolution",
    "gl_string",
    "tool",
    "imgt_accession",
    "hla_class",
    "hla_serotype",
    "protein_group",
    "imgt_match_category",
    "imgt_match_candidates",
    "gwas_traits",
    "gwas_traits_allele1",
    "gwas_traits_allele2",
    "gwas_p_values",
    "gwas_odds_ratios",
    "gwas_pmids",
    "gwas_annotation_resolution",
    "gwas_annotation_scope",
    "gwas_matched_allele",
    "gwas_match_broadened",
    "gwas_index_siblings",
    "pharm_drugs",
    "pharm_drugs_allele1",
    "pharm_drugs_allele2",
    "pharm_evidence",
    "pharm_cpic_action",
    "pharm_pmids",
    "disease_risk_summary",
    "drug_response_summary",
    "clinical_significance",
    "significance_basis",
    "allele_frequency",
    "allele_freq_population",
    "input_quality_score",
    "input_quality_tier",
    "input_quality_rationale",
    "caller_allele_quality",
)


# ---------------------------------------------------------------------------
# Exceptions
# ---------------------------------------------------------------------------


class HLAReportError(Exception):
    """
    Base class for errors raised during report generation.
    """

    pass


class OutputFileExistsError(HLAReportError):
    """
    Raised when an output path already exists and ``overwrite=False``.
    """

    pass


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------


@dataclass
class ReportContext:
    """
    Metadata shown in every report header.

    Attributes
    ----------
    hlante_version : str
        HLAnte package version.
    generated_at : str
        ISO-8601 UTC timestamp.
    db_versions : dict of str to str
        Database name → version string (e.g., ``{"imgt": "3.55.0"}``).
    disclaimer : str
        Disclaimer text appended to the report.
    cli_invocation : str
        The command line that produced the report (for provenance /
        reproducibility); empty when not invoked via the CLI.
    """

    hlante_version: str = __version__
    generated_at: str = field(default_factory=lambda: datetime.now(timezone.utc).isoformat())
    db_versions: Dict[str, str] = field(default_factory=dict)
    disclaimer: str = DISCLAIMER
    input_source: str = "typing_tool"
    cli_invocation: str = ""


@dataclass
class GenotypeRow:
    """
    Groups two :class:`AnnotatedHLA` records (``allele1`` / ``allele2``)
    from the same ``(sample_id, locus)`` into a single report row.
    """

    sample_id: str
    locus: str
    tool: str
    resolution: str
    allele1: AnnotatedHLA
    allele2: Optional[AnnotatedHLA]


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def generate_tsv(
    annotated: Sequence[AnnotatedHLA],
    output_path: Path,
    *,
    overwrite: bool = False,
    context: Optional[ReportContext] = None,
    delimiter: str = "\t",
) -> Path:
    """
    Write the main annotation table as TSV.

    Parameters
    ----------
    annotated : sequence of AnnotatedHLA
        Annotations to include in the report.
    output_path : Path
        Destination file path.
    overwrite : bool, optional
        If ``False``, raise :class:`OutputFileExistsError` when the
        destination already exists.
    context : ReportContext, optional
        Header metadata; defaults to a freshly created context.
    delimiter : str, optional
        Field delimiter (``"\\t"`` by default; pass ``","`` for CSV).

    Returns
    -------
    Path
        Path of the written file.
    """
    _check_overwrite(output_path, overwrite)
    ctx = context or ReportContext()
    rows = _group_genotypes(annotated)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    show_progress = len(rows) > PROGRESS_THRESHOLD
    iterator = _tqdm(rows, desc="TSV", disable=not show_progress)

    with output_path.open("w", encoding="utf-8", newline="") as handle:
        _write_text_metadata(handle, ctx)
        writer = csv.writer(
            handle,
            delimiter=delimiter,
            lineterminator="\n",
            quoting=csv.QUOTE_MINIMAL,
        )
        writer.writerow(TSV_COLUMNS)
        for row in iterator:
            writer.writerow([_sanitize_cell(c) for c in _row_cells(row)])

    logger.info("TSV written: %s (%d row(s))", output_path, len(rows))
    return output_path


def generate_csv(
    annotated: Sequence[AnnotatedHLA],
    output_path: Path,
    *,
    overwrite: bool = False,
    context: Optional[ReportContext] = None,
) -> Path:
    """
    Comma-delimited variant of :func:`generate_tsv`.
    """
    return generate_tsv(
        annotated,
        output_path,
        overwrite=overwrite,
        context=context,
        delimiter=",",
    )


def generate_markdown_report(
    annotated: Sequence[AnnotatedHLA],
    output_path: Path,
    *,
    overwrite: bool = False,
    context: Optional[ReportContext] = None,
) -> Path:
    """
    Write a per-sample Markdown research-annotation summary report.

    Each sample section contains an HLA genotype table, reported
    disease associations, reported pharmacogenomic associations, and
    an auto-generated interpretation note.

    Report framing uses "research-annotation summary" rather than
    "clinical summary" to keep the wording consistent with the
    research-use disclaimer and scope.

    Parameters
    ----------
    annotated : sequence of AnnotatedHLA
        Annotations.
    output_path : Path
        Destination ``.md`` file.
    overwrite : bool, optional
        Overwrite an existing file when ``True``.
    context : ReportContext, optional
        Header metadata.

    Returns
    -------
    Path
        Path of the written file.
    """
    _check_overwrite(output_path, overwrite)
    ctx = context or ReportContext()
    rows = _group_genotypes(annotated)
    samples = _group_by_sample(rows)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    sample_ids = list(samples.keys())
    show_progress = len(sample_ids) > PROGRESS_THRESHOLD
    iterator = _tqdm(sample_ids, desc="Markdown", disable=not show_progress)

    with output_path.open("w", encoding="utf-8") as handle:
        _write_markdown_header(handle, ctx)
        for sid in iterator:
            _write_markdown_sample(handle, sid, samples[sid])
        handle.write("\n---\n\n")
        handle.write(f"> **Disclaimer**: {ctx.disclaimer}\n")

    logger.info("Markdown written: %s (%d sample(s))", output_path, len(sample_ids))
    return output_path


def generate_json(
    annotated: Sequence[AnnotatedHLA],
    output_path: Path,
    *,
    overwrite: bool = False,
    context: Optional[ReportContext] = None,
    indent: Optional[int] = 2,
) -> Path:
    """
    Write the full annotation data as nested JSON (for API consumers).

    Parameters
    ----------
    annotated : sequence of AnnotatedHLA
        Annotations.
    output_path : Path
        Destination file.
    overwrite : bool, optional
        Overwrite flag.
    context : ReportContext, optional
        Header metadata.
    indent : int, optional
        JSON indent level; ``None`` for single-line output.

    Returns
    -------
    Path
        Path of the written file.
    """
    _check_overwrite(output_path, overwrite)
    ctx = context or ReportContext()
    rows = _group_genotypes(annotated)
    samples = _group_by_sample(rows)

    payload = {
        "metadata": {
            "hlante_version": ctx.hlante_version,
            "generated_at": ctx.generated_at,
            "db_versions": dict(ctx.db_versions),
            "cli_invocation": ctx.cli_invocation,
            "disclaimer": ctx.disclaimer,
            "research_use_only": True,
            "input_quality_score_definition": (
                "Heuristic summary in [0,1] of how completely the submitted "
                "allele call is supported by reference data and of the "
                "characteristics of the call itself (presence in IPD-IMGT/HLA, "
                "population frequency availability and rarity, number of "
                "reported fields, and whether the typing tool flagged the call "
                "as ambiguous). It is NOT a measure of genotype accuracy, which "
                "is determined upstream by the typing tool and the sequence "
                "data, and it is NOT a posterior probability: the penalty "
                "factors are declared values, not likelihoods. A limited tier "
                "reflects sparse supporting data or a low-detail call; it never "
                "down-weights an actionable association."
            ),
            "diplotype_caveat": DIPLOTYPE_CAVEAT,
        },
        "samples": [
            {
                "sample_id": sid,
                # Multilocus unphased genotype across the reported loci.
                "gl_string": _gl_string_for_sample(samples[sid]),
                "loci": [_row_to_json(row) for row in samples[sid]],
                # Make absence of typing explicit. Actionable loci with
                # no genotype row are indeterminate ("not typed"), not
                # Negative — a missing row must not read as "no risk".
                "actionable_loci_not_typed": _loci_not_typed(samples[sid]),
            }
            for sid in samples
        ],
    }

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=indent, default=_json_default)
        handle.write("\n")

    logger.info("JSON written: %s (%d sample(s))", output_path, len(samples))
    return output_path


def generate_all(
    annotated: Sequence[AnnotatedHLA],
    output_dir: Path,
    prefix: str = "hlante_report",
    *,
    overwrite: bool = False,
    context: Optional[ReportContext] = None,
) -> Dict[str, Path]:
    """
    Generate TSV, Markdown, and JSON outputs in a single call.

    Parameters
    ----------
    annotated : sequence of AnnotatedHLA
        Annotations.
    output_dir : Path
        Output directory (created if it does not exist).
    prefix : str, optional
        Output file name prefix (``"{prefix}.tsv"``, ``"{prefix}.md"``,
        ``"{prefix}.json"``).
    overwrite : bool, optional
        Overwrite existing files.
    context : ReportContext, optional
        Header metadata.

    Returns
    -------
    dict of str to Path
        Mapping ``{"tsv": ..., "markdown": ..., "json": ...}``.
    """
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    ctx = context or ReportContext()
    return {
        "tsv": generate_tsv(
            annotated,
            out_dir / f"{prefix}.tsv",
            overwrite=overwrite,
            context=ctx,
        ),
        "markdown": generate_markdown_report(
            annotated,
            out_dir / f"{prefix}.md",
            overwrite=overwrite,
            context=ctx,
        ),
        "json": generate_json(
            annotated,
            out_dir / f"{prefix}.json",
            overwrite=overwrite,
            context=ctx,
        ),
    }


# ---------------------------------------------------------------------------
# Grouping helpers
# ---------------------------------------------------------------------------


def _group_genotypes(annotated: Sequence[AnnotatedHLA]) -> List[GenotypeRow]:
    """
    Group :class:`AnnotatedHLA` items by ``(sample_id, locus)`` and emit
    ``GenotypeRow`` records. Input order is preserved.
    """
    groups: Dict[Tuple[str, str], List[AnnotatedHLA]] = {}
    order: List[Tuple[str, str]] = []
    for item in annotated:
        na = item.normalized_allele
        sid = na.sample_id or "unknown"
        locus = na.source_locus or na.gene or "HLA-?"
        key = (sid, locus)
        if key not in groups:
            groups[key] = []
            order.append(key)
        groups[key].append(item)

    rows: List[GenotypeRow] = []
    for key in order:
        items = sorted(
            groups[key],
            key=lambda x: x.normalized_allele.allele_index or 0,
        )
        first_na = items[0].normalized_allele
        second_na = items[1].normalized_allele if len(items) > 1 else None
        min_res = min(
            first_na.resolution_level,
            second_na.resolution_level if second_na else first_na.resolution_level,
        )
        rows.append(
            GenotypeRow(
                sample_id=key[0],
                locus=key[1],
                tool=first_na.source_tool or NA,
                resolution=RESOLUTION_LABELS.get(min_res, f"{min_res}-field"),
                allele1=items[0],
                allele2=items[1] if len(items) > 1 else None,
            )
        )
    return rows


def _group_by_sample(rows: Sequence[GenotypeRow]) -> Dict[str, List[GenotypeRow]]:
    """
    Bucket :class:`GenotypeRow` records by ``sample_id`` while preserving
    first-seen order.
    """
    samples: Dict[str, List[GenotypeRow]] = {}
    for row in rows:
        samples.setdefault(row.sample_id, []).append(row)
    return samples


#: Loci for which HLAnte carries actionable pharmacogenomic / disease
#: evidence. Used to make *absence of typing* explicit: a locus on this
#: panel that produced no genotype row for a sample is reported as
#: "not typed — not assessed", so that a missing row is never silently
#: read as "no risk".
ACTIONABLE_LOCI_CLASS_I: Tuple[str, ...] = ("HLA-A", "HLA-B", "HLA-C")
ACTIONABLE_LOCI_CLASS_II: Tuple[str, ...] = (
    "HLA-DRB1",
    "HLA-DQA1",
    "HLA-DQB1",
    "HLA-DPB1",
)

#: Typing tools that report Class I loci only; Class II "not typed" is
#: expected by design for these and is therefore not flagged.
CLASS_I_ONLY_TOOLS: FrozenSet[str] = frozenset({"optitype"})


def _expected_actionable_loci(tool: Optional[str]) -> Tuple[str, ...]:
    """
    Return the actionable-locus panel expected for a typing tool. Class
    I-only tools (e.g. OptiType) are not expected to report Class II.
    """
    if tool and tool.strip().lower() in CLASS_I_ONLY_TOOLS:
        return ACTIONABLE_LOCI_CLASS_I
    return ACTIONABLE_LOCI_CLASS_I + ACTIONABLE_LOCI_CLASS_II


def _loci_not_typed(rows: Sequence[GenotypeRow]) -> List[str]:
    """
    Actionable loci that were *not* typed for a sample.

    A clinically-relevant locus with no genotype row is indeterminate,
    not negative: the allele could not have been detected because the
    locus was not assayed. Surfacing these makes the difference between
    "allele absent" and "locus not typed" explicit in the report.
    """
    if not rows:
        return []
    typed = {row.locus.strip().upper() for row in rows}
    tool = rows[0].tool
    return [
        locus for locus in _expected_actionable_loci(tool) if locus.upper() not in typed
    ]


# ---------------------------------------------------------------------------
# TSV row construction
# ---------------------------------------------------------------------------


def _row_cells(row: GenotypeRow) -> List[str]:
    """
    Convert a :class:`GenotypeRow` into cells ordered by
    :data:`TSV_COLUMNS`.
    """
    a1 = row.allele1
    a2 = row.allele2
    na1 = a1.normalized_allele
    na2 = a2.normalized_allele if a2 else None

    gwas_all = list(a1.gwas_hits) + (list(a2.gwas_hits) if a2 else [])
    pharm_all = list(a1.pharm_annotations) + (list(a2.pharm_annotations) if a2 else [])

    return [
        row.sample_id,
        row.locus,
        na1.allele_name,
        na2.allele_name if na2 else NA,
        row.resolution,
        _gl_string_for_locus(row) or NA,
        row.tool,
        _pipe([na1.imgt_accession, na2.imgt_accession if na2 else None]),
        na1.hla_class or NA,
        # Per-allele serotype lookup.
        _pipe(
            [
                _serotype(na1.allele_name),
                _serotype(na2.allele_name) if na2 else None,
            ]
        ),
        _pipe([na1.protein_group, na2.protein_group if na2 else None]),
        _pipe_slots(
            [na1.imgt_match_category, na2.imgt_match_category if na2 else None]
        ),
        _pipe_slots(
            [
                str(na1.imgt_match_candidates),
                str(na2.imgt_match_candidates) if na2 else None,
            ]
        ),
        _pipe([h.trait for h in gwas_all]),
        # Per-allele attribution.
        _pipe([h.trait for h in a1.gwas_hits]),
        _pipe([h.trait for h in a2.gwas_hits]) if a2 else NA,
        _pipe([_fmt_float(h.p_value) for h in gwas_all]),
        _pipe([_fmt_float(h.odds_ratio) for h in gwas_all]),
        _pipe([h.pmid for h in gwas_all]),
        _pipe(
            [
                a1.gwas_resolution_used,
                a2.gwas_resolution_used if a2 else None,
            ]
        ),
        # Annotation scope per allele (worst case across hits).
        _pipe(
            [
                _worst_annotation_scope(a1.gwas_hits),
                _worst_annotation_scope(a2.gwas_hits) if a2 else None,
            ]
        ),
        _pipe_slots(
            [
                _matched_alleles(a1.gwas_hits),
                _matched_alleles(a2.gwas_hits) if a2 else None,
            ]
        ),
        _pipe_slots(
            [
                _broadened_flag(a1.gwas_hits),
                _broadened_flag(a2.gwas_hits) if a2 else None,
            ]
        ),
        _pipe(
            [
                _max_expansion(a1.gwas_hits),
                _max_expansion(a2.gwas_hits) if a2 else None,
            ]
        ),
        _pipe([p.drug for p in pharm_all]),
        # Per-allele attribution.
        _pipe([p.drug for p in a1.pharm_annotations]),
        _pipe([p.drug for p in a2.pharm_annotations]) if a2 else NA,
        _pipe([_evidence_label(p.evidence_level) for p in pharm_all]),
        # CPIC standardised action verb per pharm annotation, keyed on the
        # Carried allele and the drug — the recommendation is
        # Allele-dependent, not drug-only.
        _pipe([_cpic_action(p.allele, p.drug) for p in pharm_all]),
        # Aggregate long PMID lists (44 → "top3 (+41 more)").
        _pipe([_aggregate_pmids(p.pmid) for p in pharm_all]),
        _pipe(
            [
                a1.disease_risk_summary,
                a2.disease_risk_summary if a2 else None,
            ]
        ),
        _pipe(
            [
                a1.drug_response_summary,
                a2.drug_response_summary if a2 else None,
            ]
        ),
        _pipe(
            [
                a1.clinical_significance,
                a2.clinical_significance if a2 else None,
            ]
        ),
        _pipe_slots(
            [
                a1.significance_basis,
                a2.significance_basis if a2 else None,
            ]
        ),
        _pipe_slots(
            [
                _fmt_float(a1.allele_frequency, decimals=6),
                _fmt_float(a2.allele_frequency, decimals=6) if a2 else None,
            ]
        ),
        _pipe_slots(
            [
                a1.frequency_population,
                a2.frequency_population if a2 else None,
            ]
        ),
        _pipe_slots(
            [
                _fmt_float(a1.input_quality_score, decimals=4),
                _fmt_float(a2.input_quality_score, decimals=4) if a2 else None,
            ]
        ),
        _pipe(
            [
                getattr(a1, "input_quality_tier", "NA"),
                getattr(a2, "input_quality_tier", "NA") if a2 else None,
            ]
        ),
        # Use ";;" between allele1 and allele2 to avoid pipe
        # collision with intra-rationale reason codes.
        _rationale_pipe(
            a1.input_quality_rationale,
            a2.input_quality_rationale if a2 else None,
        ),
        # Quality reported by the typing tool itself, where it provides one.
        # Distinct from input_quality_score, which HLAnte computes.
        _pipe_slots(
            [
                _caller_quality(na1),
                _caller_quality(na2) if na2 else None,
            ]
        ),
    ]


# Static map from lowercased drug name → CPIC-standardised action verb.
#
# PharmGKB dumps large PMID lists
# But exposes no short "what do I do?" verb. CPIC does, via its
# Guideline tables. This map is intentionally small and hand-curated
# — the CPIC corpus has only a few Level-A HLA pairs and each has a
# Single standard recommendation.
#
# The phrasing uses CPIC's five-verb vocabulary:
#
# * ``Contraindicated``
# * ``Test required before prescribing``
# * ``Use with caution — monitor``
# * ``Alternative therapy recommended``
# * ``No specific action required``
#
# Missing drugs fall through to ``NA``. Drug names are matched in a
# Case-insensitive, whitespace-trimmed form.
# Static allele → WHO/IMGT serotype lookup.
#
# Most HLA-autoimmune literature
# is written at serotype level (DR2 / DR3 / DR4 / DQ2 / DQ8) rather
# Than at two-field resolution. The map is two-field-granular; more
# Specific keys (``"DQB1*03:02"``) take precedence over first-field
# Fallbacks (``"DRB1*03"``).
#
# The list below is intentionally hand-curated and limited to
# Commonly referenced serotypes; alleles not in the map yield ``NA``.
HLA_SEROTYPE_MAP: Dict[str, str] = {
    # Class-II DRB1 serotypes (DR nomenclature)
    "DRB1*01": "DR1",
    "DRB1*03": "DR3",
    "DRB1*04": "DR4",
    "DRB1*07": "DR7",
    "DRB1*08": "DR8",
    "DRB1*09": "DR9",
    "DRB1*10": "DR10",
    "DRB1*11": "DR11",
    "DRB1*12": "DR12",
    "DRB1*13": "DR13",
    "DRB1*14": "DR14",
    "DRB1*15": "DR15",
    "DRB1*16": "DR16",
    # DQB1 → DQ serotypes (subtype-specific for DQ7/DQ8/DQ9)
    "DQB1*02": "DQ2",
    "DQB1*03:01": "DQ7",
    "DQB1*03:02": "DQ8",
    "DQB1*03:03": "DQ9",
    "DQB1*03:04": "DQ7",
    "DQB1*04": "DQ4",
    "DQB1*05": "DQ5",
    "DQB1*06": "DQ6",
    # DPB1
    "DPB1*04": "DP4",
    # Class-I A serotypes
    "A*01": "A1",
    "A*02": "A2",
    "A*03": "A3",
    "A*11": "A11",
    "A*23": "A23",
    "A*24": "A24",
    "A*25": "A25",
    "A*26": "A26",
    "A*29": "A29",
    "A*30": "A30",
    "A*31": "A31",
    "A*32": "A32",
    "A*33": "A33",
    "A*34": "A34",
    "A*36": "A36",
    "A*66": "A66",
    "A*68": "A68",
    "A*69": "A69",
    "A*74": "A74",
    "A*80": "A80",
    # Class-I B serotypes
    "B*07": "B7",
    "B*08": "B8",
    "B*13": "B13",
    "B*14": "B14",
    "B*15": "B15",
    "B*18": "B18",
    "B*27": "B27",
    "B*35": "B35",
    "B*37": "B37",
    "B*38": "B38",
    "B*39": "B39",
    "B*40": "B40",
    "B*41": "B41",
    "B*42": "B42",
    "B*44": "B44",
    "B*45": "B45",
    "B*46": "B46",
    "B*47": "B47",
    "B*48": "B48",
    "B*49": "B49",
    "B*50": "B50",
    "B*51": "B51",
    "B*52": "B52",
    "B*53": "B53",
    "B*54": "B54",
    "B*55": "B55",
    "B*56": "B56",
    "B*57": "B57",
    "B*58": "B58",
    "B*59": "B59",
    "B*67": "B67",
    "B*73": "B73",
    "B*78": "B78",
    # Class-I C serotypes (Cw nomenclature)
    "C*01": "Cw1",
    "C*02": "Cw2",
    "C*03": "Cw3",
    "C*04": "Cw4",
    "C*05": "Cw5",
    "C*06": "Cw6",
    "C*07": "Cw7",
    "C*08": "Cw8",
}


def _serotype(allele_name: Optional[str]) -> Optional[str]:
    """
    Look up the serotype label for an HLA allele.

    Tries a two-field key first (``DQB1*03:02`` → ``DQ8``) and falls
    back to a first-field key (``DRB1*03:02`` → ``DR3``). Returns
    ``None`` for alleles not covered by :data:`HLA_SEROTYPE_MAP`.
    """
    if not allele_name:
        return None
    clean = allele_name.replace("HLA-", "").strip()
    if ":" in clean and "*" in clean:
        gene, rest = clean.split("*", 1)
        fields = rest.split(":")
        two_field = f"{gene}*{fields[0]}:{fields[1]}"
        first_field = f"{gene}*{fields[0]}"
        if two_field in HLA_SEROTYPE_MAP:
            return HLA_SEROTYPE_MAP[two_field]
        if first_field in HLA_SEROTYPE_MAP:
            return HLA_SEROTYPE_MAP[first_field]
    return HLA_SEROTYPE_MAP.get(clean)


#: Drug-level fallback action, used only when the (allele, drug) pair is not
#: present in :data:`ALLELE_DRUG_ACTION_MAP`. These are coarse, allele-blind
#: defaults; the recommendation for a given drug frequently depends on which
#: HLA allele is carried, so the allele-specific map below takes precedence.
CPIC_ACTION_MAP: Dict[str, str] = {
    "abacavir": "Avoid",
    "allopurinol": "Use alternative",
    "carbamazepine": "Use with caution — see allele-specific guidance",
    "oxcarbazepine": "Use with caution — see allele-specific guidance",
    "phenytoin": "Use with caution — see allele-specific guidance",
    "fosphenytoin": "Use with caution — see allele-specific guidance",
    "dapsone": "Test before prescribing",
    "flucloxacillin": "Use with caution — monitor (DILI risk)",
    "lamotrigine": "Use with caution — monitor",
    "lapatinib": "Use with caution — monitor",
    "nevirapine": "Use alternative",
    "ziagen": "Avoid",
}

#: CPIC-aligned action keyed on ``(two-field allele, drug)``. The carried allele
#: materially changes the recommendation (e.g. ``HLA-A*31:01`` carbamazepine is
#: a weaker association than ``HLA-B*15:02`` and is NOT a flat
#: contraindication), so the verb must not be derived from the drug alone.
#: Population caveats are embedded where CPIC stratifies. Falls back to
#: :data:`CPIC_ACTION_MAP` for pairs not listed here.
ALLELE_DRUG_ACTION_MAP: Dict[Tuple[str, str], str] = {
    ("B*57:01", "abacavir"): "Contraindicated (do not use)",
    ("B*57:01", "ziagen"): "Contraindicated (do not use)",
    ("B*58:01", "allopurinol"): "Contraindicated — use alternative",
    ("B*15:02", "carbamazepine"): "Avoid — contraindicated if carbamazepine-naïve",
    ("B*15:02", "oxcarbazepine"): "Avoid — contraindicated if oxcarbazepine-naïve",
    ("B*15:02", "phenytoin"): "Use alternative; consider CYP2C9 status",
    ("B*15:02", "fosphenytoin"): "Use alternative; consider CYP2C9 status",
    ("A*31:01", "carbamazepine"): "Use alternative if available; otherwise increased monitoring",
    ("B*13:01", "dapsone"): "Use alternative; test before prescribing",
}


def _evidence_label(level: Optional[str]) -> Optional[str]:
    """
    Return evidence level with ``(low evidence)`` suffix for Level-3/4,
    so TSV consumers can distinguish case-report associations at a glance.
    """
    if not level:
        return None
    return f"{level} (low evidence)" if level.upper() in {"3", "4"} else level


def _two_field_key(allele: Optional[str]) -> Optional[str]:
    """
    Reduce an allele string to its ``GENE*field1:field2`` key for CPIC
    lookup (``HLA-B*15:02:01`` → ``B*15:02``), dropping the ``HLA-``
    prefix and any expression/group suffix on the second field. Returns
    ``None`` when the input has fewer than two fields.
    """
    if not allele:
        return None
    bare = allele.strip()
    if bare.upper().startswith("HLA-"):
        bare = bare[4:]
    if "*" not in bare:
        return None
    gene, rest = bare.split("*", 1)
    fields = rest.split(":")
    if len(fields) < 2:
        return None
    second = fields[1].rstrip("NLSQCAGP")  # Drop any expression/group suffix
    if not second:
        return None
    return f"{gene}*{fields[0]}:{second}"


def _cpic_action(allele: Optional[str], drug: Optional[str]) -> Optional[str]:
    """
    Return the CPIC-aligned action verb for an ``(allele, drug)`` pair.

    The carried HLA allele is required because the recommendation for a
    drug is allele-dependent (e.g. ``HLA-B*15:02`` vs ``HLA-A*31:01`` for
    carbamazepine). Resolution order: exact ``(two-field allele, drug)``
    match in :data:`ALLELE_DRUG_ACTION_MAP`, then the drug-level fallback
    in :data:`CPIC_ACTION_MAP`, then ``None``.
    """
    if not drug:
        return None
    drug_key = drug.strip().lower()
    two_field = _two_field_key(allele)
    if two_field is not None:
        specific = ALLELE_DRUG_ACTION_MAP.get((two_field, drug_key))
        if specific is not None:
            return specific
    return CPIC_ACTION_MAP.get(drug_key)


# Ranking used to pick the "worst" scope for the per-allele summary
# Cell: an allele carrying even one locus-level fallback hit should
# Surface as ``locus``, so a clinician scanning the TSV sees the
# Broadest scope at a glance.
_SCOPE_RANK: Dict[str, int] = {"allele": 0, "subtype": 1, "locus": 2}


def _aggregate_pmids(pmids: Sequence[str], top: int = 3) -> str:
    """
    Compact PMID list for cell display.

    Examples
    --------
    ``["18", "19"]`` → ``"18,19"``.
    ``["18", "19", "20", "21", "22"]`` → ``"18,19,20 (+2 more)"`` when
    ``top=3``.

    Notes
    -----
    the B*58:01/allopurinol record carries 44 PMIDs. Dumping all
    44 into a single TSV cell is evidence-heavy but not actionable.
    The aggregate form keeps the top three for quick cross-reference
    and surfaces the remainder count so the user knows a lookup is
    appropriate.
    """
    clean = [p for p in pmids if p and p.strip() and p.strip() != NA]
    if not clean:
        return ""
    if len(clean) <= top:
        return ",".join(clean)
    extra = len(clean) - top
    return f"{','.join(clean[:top])} (+{extra} more)"


def _worst_annotation_scope(hits: Iterable[Any]) -> Optional[str]:
    """
    Return the broadest (worst-case) annotation scope across a set of
    GWAS hits, or ``None`` when the iterable is empty.
    """
    scopes = [getattr(h, "annotation_scope", "allele") or "allele" for h in hits]
    if not scopes:
        return None
    return max(scopes, key=lambda s: _SCOPE_RANK.get(s, 0))


def _matched_alleles(hits: Sequence[Any]) -> Optional[str]:
    """
    Catalogue key(s) that actually matched, de-duplicated in order.

    A reader needs this to tell whether an association was reported for the
    allele they submitted or for a less specific name derived from it.
    """
    seen: List[str] = []
    for hit in hits:
        key = getattr(hit, "matched_allele", "") or ""
        if key and key not in seen:
            seen.append(key)
    return ";".join(seen) if seen else None


def _broadened_flag(hits: Sequence[Any]) -> Optional[str]:
    """
    ``yes`` when any association for this allele required truncation.
    """
    if not hits:
        return None
    return "yes" if any(getattr(h, "match_was_broadened", False) for h in hits) else "no"


def _max_expansion(hits: Iterable[Any]) -> Optional[str]:
    """
    Return the maximum ``index_siblings`` across a set of
    GWAS hits as a string, or ``None`` when the iterable is empty.
    """
    sizes = [int(getattr(h, "index_siblings", 1) or 1) for h in hits]
    if not sizes:
        return None
    return str(max(sizes))


def _rationale_pipe(r1: Optional[str], r2: Optional[str]) -> str:
    """
    Join two allele-level rationale strings with ``";;"`` so that the
    pipes inside each rationale remain unambiguous.

    Returns :data:`NA` when both sides are empty.
    """
    left = (r1 or "").strip()
    right = (r2 or "").strip() if r2 is not None else ""
    parts = [p for p in (left, right) if p and p != NA]
    if not parts:
        return NA
    if r2 is None:
        return parts[0]
    return f"{left or NA};;{right or NA}"


def _pipe(values: Iterable[Any]) -> str:
    """
    Join non-empty values with a ``|`` separator.

    Returns :data:`NA` when every value is missing.
    """
    cleaned = [str(v).strip() for v in values if v is not None and str(v).strip() not in ("", NA)]
    if not cleaned:
        return NA
    return PIPE.join(cleaned)


def _pipe_slots(values: Iterable[Any]) -> str:
    """
    Like :func:`_pipe` but preserves positional slots: each value occupies
    its own ``|``-separated slot, using ``NA`` for missing entries.

    Use this for per-allele numeric columns (frequency, score) where the
    slot position carries meaning — i.e. the caller must know which slot
    belongs to allele1 vs allele2.
    """
    slots = [
        str(v).strip() if (v is not None and str(v).strip() not in ("", NA)) else NA for v in values
    ]
    if all(s == NA for s in slots):
        return NA
    return PIPE.join(slots)


def _fmt_float(value: Optional[float], decimals: int = FLOAT_DECIMALS) -> Optional[str]:
    """
    Format a float as fixed-point with ``decimals`` digits (no scientific
    notation), per the original output specification.

    Notes
    -----
    Very small p-values (< 10⁻⁴) are lossily rounded to ``0.0000``. Use
    the JSON output when exact precision matters.
    """
    if value is None:
        return None
    return f"{value:.{decimals}f}"


# ---------------------------------------------------------------------------
# Header writers
# ---------------------------------------------------------------------------


def _write_text_metadata(handle: Any, ctx: ReportContext) -> None:
    """
    Write a ``#``-commented metadata block at the top of TSV/CSV outputs.
    """
    handle.write(f"# HLAnte {ctx.hlante_version}\n")
    handle.write(f"# Generated at: {ctx.generated_at}\n")
    for db, ver in sorted(ctx.db_versions.items()):
        handle.write(f"# {db} version: {ver}\n")
    handle.write(f"# input_source: {ctx.input_source}\n")
    handle.write(f"# Disclaimer: {ctx.disclaimer}\n")


def _write_markdown_header(handle: Any, ctx: ReportContext) -> None:
    """
    Write the top-level Markdown report header.

    Notes
    -----
    The title deliberately avoids the word ``"Clinical"``, which would
    commit the output to clinical standards the tool does not meet (no
    ACMG/AMP implementation). The heading is
    ``"HLAnte Research Annotation Report"`` and a prominent disclaimer
    banner is emitted immediately after the metadata block.
    """
    handle.write("# HLAnte Research Annotation Report\n\n")
    handle.write(f"- **HLAnte version**: {ctx.hlante_version}\n")
    handle.write(f"- **Generated at**: {ctx.generated_at}\n")
    for db, ver in sorted(ctx.db_versions.items()):
        handle.write(f"- **{db} version**: {ver}\n")
    handle.write("\n")
    handle.write(
        "> ⚠️ **Research use only.** The labels in this report are "
        "evidence-strength descriptors, not ACMG/AMP classifications. "
        "Do not use this output for clinical decision-making.\n\n"
    )


# ---------------------------------------------------------------------------
# Markdown — per-sample section
# ---------------------------------------------------------------------------


def _write_markdown_sample(handle: Any, sample_id: str, rows: Sequence[GenotypeRow]) -> None:
    """
    Write the Markdown section for a single sample.
    """
    handle.write(f"## Sample: {sample_id}\n\n")

    # --- HLA genotype table ---
    handle.write("### HLA Genotype\n\n")
    handle.write("| Locus | Allele 1 | Allele 2 | Resolution | Tool |\n")
    handle.write("|-------|----------|----------|------------|------|\n")
    for row in rows:
        allele2_name = row.allele2.normalized_allele.allele_name if row.allele2 else NA
        handle.write(
            f"| {row.locus} | {row.allele1.normalized_allele.allele_name} "
            f"| {allele2_name} | {row.resolution} | {row.tool} |\n"
        )
    handle.write("\n")

    # Make absence of typing explicit so a missing locus is never read
    # as "no risk". Actionable loci with no genotype row are indeterminate.
    not_typed = _loci_not_typed(rows)
    if not_typed:
        handle.write(
            "> ⚠️ **Loci not typed (not assessed for HLA-linked risk):** "
            f"{', '.join(not_typed)}. A locus that was not typed is "
            "*indeterminate*, not negative — the absence of an alert does "
            "**not** mean the risk allele is absent.\n\n"
        )

    # --- Disease-association summary ---
    # Heading was "Disease Risk Summary"; renamed to avoid
    # Implying the tool quantifies individual clinical risk.
    handle.write("### Reported Disease Associations\n\n")
    any_risk = False
    for row in rows:
        for annot in _iter_alleles(row):
            allele_name = annot.normalized_allele.allele_name
            for hit in annot.gwas_hits:
                or_txt = f"OR={hit.odds_ratio:.2f}" if hit.odds_ratio is not None else "OR=NA"
                handle.write(f"- **HLA-{allele_name}**: {hit.trait} ({or_txt})\n")
                study = hit.study_accession or "GWAS Catalog"
                handle.write(f"  - Source: GWAS Catalog ({study}), PMID: {hit.pmid or NA}\n")
                # Annotate deprecated EFO provenance
                if getattr(hit, "trait_was_deprecated", False):
                    handle.write(
                        "  - ⚠️ Note: trait name was remapped from a "
                        "deprecated EFO term. Verify current "
                        "classification at https://www.ebi.ac.uk/efo/\n"
                    )
                # Extreme / quantitative-trait effect warning
                esw = getattr(hit, "effect_size_warning", "") or ""
                if esw:
                    if "quantitative_trait_effect" in esw:
                        handle.write(
                            f"  - ⚠️ Effect size note: OR={hit.odds_ratio:.2f} "
                            "exceeds 10 for a trait whose text suggests a "
                            "continuous / quantitative measure. The reported "
                            "value may be a β coefficient or per-SD effect "
                            "size, not a binary odds ratio. Verify the source "
                            f"study (PMID: {hit.pmid or NA}) before clinical "
                            "interpretation.\n"
                        )
                    else:
                        handle.write(
                            f"  - ⚠️ Effect size note: OR={hit.odds_ratio:.2f} "
                            "exceeds 10 — unusually large for a binary-disease "
                            "odds ratio. Verify the source study "
                            f"(PMID: {hit.pmid or NA}).\n"
                        )
                # Annotation scope when fallback was used
                scope = getattr(hit, "annotation_scope", "allele")
                if scope != "allele":
                    exp = getattr(hit, "index_siblings", 1)
                    handle.write(
                        f"  - ⚠️ Annotation scope: this GWAS hit was "
                        f"matched at {scope}-level resolution "
                        f"({exp} IMGT allele(s) share the matched prefix). "
                        f"The effect size may not apply specifically to "
                        f"HLA-{allele_name}.\n"
                    )
                any_risk = True
            for entry in annot.disease_entries:
                handle.write(
                    f"- **HLA-{allele_name}**: {entry.condition} (Curated: {entry.significance})\n"
                )
                handle.write(
                    f"  - Variation ID: {entry.variation_id}, Review: {entry.review_status}\n"
                )
                any_risk = True
    if not any_risk:
        handle.write("- No allele with a reported disease association.\n")
    handle.write("\n")

    # --- Pharmacogenomic associations ---
    # Heading was "Drug Response Warnings". Reworded to
    # "Reported Pharmacogenomic Associations" to avoid framing
    # Research-level PharmGKB / CPIC records as clinical warnings.
    handle.write("### Reported Pharmacogenomic Associations\n\n")
    any_drug = False
    for row in rows:
        for annot in _iter_alleles(row):
            allele_name = annot.normalized_allele.allele_name
            for ann in annot.pharm_annotations:
                is_strong = (ann.evidence_level or "").upper() in {"1A", "1B"}
                # Previously "⚠️ Critical" / "ℹ️ Info";
                # Reworded to evidence-strength markers.
                icon = "⚠️ Strong evidence" if is_strong else "ℹ️ Moderate evidence"
                phenotype = ann.phenotype or "interaction"
                handle.write(
                    f"- {icon}: **HLA-{allele_name}** — {ann.drug} {phenotype} "
                    f"({ann.evidence_level} evidence)\n"
                )
                if ann.cpic_url:
                    handle.write(f"  - CPIC guideline: {ann.cpic_url}\n")
                if ann.pharmgkb_url:
                    handle.write(f"  - PharmGKB: {ann.pharmgkb_url}\n")
                any_drug = True
    if not any_drug:
        handle.write("- No reported drug-allele association.\n")
    handle.write("\n")

    # --- Interpretation note (research use) ---
    # Heading was "Clinical Interpretation Note". Reworded to
    # Decouple the auto-generated sentence from a clinical frame.
    handle.write("### Interpretation Note (research use only)\n\n")
    handle.write(_build_clinical_note(rows))
    handle.write("\n\n")

    # Surface the diplotype/zygosity simplification explicitly.
    handle.write(f"> **Diplotype/zygosity note:** {DIPLOTYPE_CAVEAT}\n\n")


def _iter_alleles(row: GenotypeRow) -> Iterable[AnnotatedHLA]:
    yield row.allele1
    if row.allele2 is not None:
        yield row.allele2


def _build_clinical_note(rows: Sequence[GenotypeRow]) -> str:
    """
    Compose an auto-generated sentence summarizing alleles with
    evidence-backed associations.

    Notes
    -----
    the membership set was historically the ACMG-adjacent triad
    ``{"Pathogenic", "Likely Pathogenic", "Risk factor"}``. It is now
    imported from :mod:`hlante.annotator` so the set tracks the
    evidence-strength labels without drift.
    """
    from hlante.annotator import (
        SIGNIFICANCE_LIKELY_PATHOGENIC,
        SIGNIFICANCE_PATHOGENIC,
        SIGNIFICANCE_RISK_FACTOR,
    )

    relevant = {
        SIGNIFICANCE_PATHOGENIC,
        SIGNIFICANCE_LIKELY_PATHOGENIC,
        SIGNIFICANCE_RISK_FACTOR,
    }
    sentences: List[str] = []
    for row in rows:
        for annot in _iter_alleles(row):
            if annot.clinical_significance not in relevant:
                continue
            name = annot.normalized_allele.allele_name
            driver = _select_driver(annot)
            sentences.append(driver.format(name=name))
    if not sentences:
        return (
            "No HLA allele in this sample carries an evidence-backed "
            "association in the current annotation sources. This is a "
            "research-use annotation and is not a clinical finding."
        )
    return (
        " ".join(sentences) + " This summary is research-use only; any follow-up belongs "
        "with a qualified clinician."
    )


def _select_driver(annot: AnnotatedHLA) -> str:
    """
    Pick the strongest finding to describe in the auto-generated note.
    """
    if annot.pharm_annotations:
        ann = annot.pharm_annotations[0]
        return (
            f"This sample carries HLA-{{name}}. This allele is associated "
            f"with {ann.drug} {ann.phenotype or 'reactions'} "
            f"({ann.evidence_level} evidence, CPIC guideline)."
        )
    if annot.gwas_hits:
        hit = annot.gwas_hits[0]
        or_txt = f"OR={hit.odds_ratio:.2f}" if hit.odds_ratio is not None else "OR=NA"
        return (
            f"This sample carries HLA-{{name}}. The GWAS Catalog "
            f"reports an association with {hit.trait} ({or_txt})."
        )
    return (
        f"This sample carries HLA-{{name}} (evidence-strength label: "
        f"{annot.clinical_significance})."
    )


# ---------------------------------------------------------------------------
# JSON helpers
# ---------------------------------------------------------------------------


#: GL String operators HLAnte is able to emit. The tool holds at most two
#: allele designations per locus per sample and models neither chromosomal
#: phase nor alternative genotypes, so only the genotype delimiter (``+``)
#: and the locus delimiter (``^``) are ever produced. Emitting ``/``
#: (possible alleles), ``~`` (phased genes), ``|`` (possible genotypes) or
#: ``?`` (possible loci) would assert information HLAnte does not hold.
#: Grammar: Milius et al. 2013; GL String 1.1, Mack et al. 2023.
GL_OPERATORS_EMITTED: Tuple[str, ...] = ("+", "^")
GL_OPERATORS_NOT_EMITTED: Tuple[str, ...] = ("/", "~", "|", "?")


def _gl_token(allele_name: str, locus: str) -> Optional[str]:
    """
    Fully qualified GL String allele name, e.g. ``HLA-A*01:01``.

    GL Strings name alleles with the gene included, so a bare ``A*01:01``
    from a typing tool is qualified with the locus it was reported under.
    """
    name = (allele_name or "").strip()
    if not name or name == NA:
        return None
    if name.startswith("HLA-"):
        return name
    gene = (locus or "").strip()
    if gene and not gene.startswith("HLA-"):
        gene = f"HLA-{gene}"
    if gene and name.startswith(gene.split("-", 1)[1] + "*"):
        return f"HLA-{name}"
    return f"HLA-{name}" if "*" in name else None


def _gl_string_for_locus(row: "GenotypeRow") -> Optional[str]:
    """
    GL String for one locus: the gene copies joined by ``+``.

    A locus that reports a single allele yields that allele alone — the
    missing copy is not invented, consistent with :func:`_zygosity`.
    """
    first = _gl_token(row.allele1.normalized_allele.allele_name, row.locus)
    second = (
        _gl_token(row.allele2.normalized_allele.allele_name, row.locus)
        if row.allele2
        else None
    )
    tokens = [t for t in (first, second) if t]
    if not tokens:
        return None
    return "+".join(tokens)


def _gl_string_for_sample(rows: Sequence["GenotypeRow"]) -> Optional[str]:
    """
    Multilocus unphased genotype: per-locus GL Strings joined by ``^``.

    Loci are emitted in the order they appear in the report so the string is
    reproducible; the grammar imposes no ordering rule.
    """
    parts = [g for g in (_gl_string_for_locus(r) for r in rows) if g]
    return "^".join(parts) if parts else None


def _caller_quality(na: Optional[NormalizedAllele]) -> Optional[str]:
    """
    Per-allele quality as reported by the typing tool, formatted for the TSV.

    Only T1K's native layout supplies this. arcasHLA and HLA-HD report no
    per-allele quality, and OptiType reports a solution-level objective rather
    than an allele quality, so those tools yield ``None``.
    """
    if na is None:
        return None
    value = getattr(na, "caller_quality", None)
    return None if value is None else f"{value:g}"


def _zygosity(row: GenotypeRow) -> str:
    """
    Explicit, non-inferred zygosity state for a locus.

    ``"single_allele_reported"`` is deliberately distinct from
    ``"homozygous"``: a missing second allele is NOT treated as a second
    copy of the first.
    """
    a1 = row.allele1.normalized_allele.allele_name
    a2 = row.allele2.normalized_allele.allele_name if row.allele2 else None
    if a2 is None:
        return "single_allele_reported"
    return "homozygous" if a1 == a2 else "heterozygous"


#: Leading characters that a spreadsheet (Excel / LibreOffice / Sheets)
#: interprets as the start of a formula. A DB-sourced cell beginning with one
#: of these is neutralised so it cannot trigger CSV/TSV formula injection
#: (CWE-1236) when the report is opened in a spreadsheet.
_FORMULA_LEAD_CHARS: FrozenSet[str] = frozenset({"=", "+", "-", "@", "\t", "\r"})


def _sanitize_cell(value: Any) -> Any:
    """
    Prefix a single quote to any string cell that begins with a
    spreadsheet formula-trigger character, so it is rendered as literal
    text rather than evaluated.
    """
    if isinstance(value, str) and value and value[0] in _FORMULA_LEAD_CHARS:
        return "'" + value
    return value


def _row_to_json(row: GenotypeRow) -> Dict[str, Any]:
    return {
        "locus": row.locus,
        "resolution": row.resolution,
        # Standard GL String representation (Milius 2013; Mack 2023).
        "gl_string": _gl_string_for_locus(row),
        "tool": row.tool,
        # Machine-readable zygosity; single_allele_reported != homozygous.
        "zygosity": _zygosity(row),
        # Per-record research-use flag. A programmatic consumer that reads
        # Only the loci array (and never sees the metadata footer) still gets
        # an unambiguous, machine-readable "not for clinical use" marker on
        # Every annotated record.
        "research_use_only": True,
        "allele1": _annotated_to_json(row.allele1),
        "allele2": _annotated_to_json(row.allele2) if row.allele2 else None,
    }


def _annotated_to_json(annot: AnnotatedHLA) -> Dict[str, Any]:
    return {
        "normalized_allele": asdict(annot.normalized_allele),
        "gwas_hits": [asdict(h) for h in annot.gwas_hits],
        "pharm_annotations": [asdict(p) for p in annot.pharm_annotations],
        "disease_risk_summary": annot.disease_risk_summary,
        "drug_response_summary": annot.drug_response_summary,
        "clinical_significance": annot.clinical_significance,
    }


def _json_default(obj: Any) -> Any:
    if isinstance(obj, Path):
        return str(obj)
    if isinstance(obj, (set, frozenset)):
        return sorted(obj)
    raise TypeError(f"Non-serializable type for JSON output: {type(obj).__name__}")


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------


def _check_overwrite(path: Path, overwrite: bool) -> None:
    """
    Raise when ``path`` already exists and overwriting is disabled.
    """
    if Path(path).exists() and not overwrite:
        raise OutputFileExistsError(
            f"Output file already exists: {path}. "
            "Pass overwrite=True (or --overwrite on the CLI) to replace it."
        )


__all__ = [
    "NA",
    "PIPE",
    "PROGRESS_THRESHOLD",
    "DISCLAIMER",
    "DIPLOTYPE_CAVEAT",
    "NO_DISEASE_SUMMARY",
    "NO_DRUG_SUMMARY",
    "RISK_PREFIX_HIGH",
    "RISK_PREFIX_MODERATE",
    "RISK_PREFIX_PROTECTIVE",
    "RISK_PREFIX_ASSOCIATION",
    "TSV_COLUMNS",
    "HLAReportError",
    "OutputFileExistsError",
    "ReportContext",
    "GenotypeRow",
    "generate_tsv",
    "generate_csv",
    "generate_markdown_report",
    "generate_json",
    "generate_all",
]
