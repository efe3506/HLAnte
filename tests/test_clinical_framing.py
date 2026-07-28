"""
tests.test_clinical_framing
===========================

Tests for the clinical-safety framing fixes C4–C8:

- the input-quality score is documented as an uncalibrated heuristic
  heuristic, not clinical certainty (disclaimer + JSON metadata).
- direction of effect — any OR < 1.0 is an inverse (protective) association.
- significance labels avoid ACMG/AMP vocabulary ("pathogenic", "VUS").
- every JSON record carries a machine-readable research_use_only flag.
- a GWAS hit with a missing p-value does not pass the genome-wide
  significance filter, nor escalate the significance label.
"""

from __future__ import annotations

import json
from pathlib import Path

from hlante.annotator import (
    SIGNIFICANCE_LIKELY_PATHOGENIC,
    SIGNIFICANCE_PATHOGENIC,
    SIGNIFICANCE_RISK_FACTOR,
    SIGNIFICANCE_VUS,
    AnnotatedHLA,
    RISK_PREFIX_ASSOCIATION,
    RISK_PREFIX_HIGH,
    RISK_PREFIX_MODERATE,
    RISK_PREFIX_PROTECTIVE,
    _classify_significance,
    _severity_label,
)
from hlante.db.gwas import GWASClient, GWASHit
from hlante.normalizer import NormalizedAllele
from hlante.reporter import DISCLAIMER, GenotypeRow, _zygosity, generate_json


def _allele() -> NormalizedAllele:
    return NormalizedAllele(
        allele_name="A*02:01",
        imgt_accession="HLA00001",
        protein_group=None,
        hla_class="I",
        gene="HLA-A",
        resolution_level=2,
        is_ambiguous=False,
        is_novel=False,
    )


def _hit(p_value, *, trait="Some disease", oddsr=2.0) -> GWASHit:
    return GWASHit(
        trait=trait,
        p_value=p_value,
        odds_ratio=oddsr,
        pmid="12345678",
        study_accession="GCST000001",
        allele="A*02:01",
    )


# ---------------------------------------------------------------------------
# Direction of effect
# ---------------------------------------------------------------------------


class TestSeverityDirection:
    def test_weakly_protective_is_inverse(self) -> None:
        # 0.5 < OR < 1.0 used to be mislabeled "Reported association".
        assert _severity_label(0.6) == RISK_PREFIX_PROTECTIVE
        assert _severity_label(0.95) == RISK_PREFIX_PROTECTIVE

    def test_strongly_protective_is_inverse(self) -> None:
        assert _severity_label(0.3) == RISK_PREFIX_PROTECTIVE

    def test_null_effect_and_risk(self) -> None:
        assert _severity_label(1.0) == RISK_PREFIX_ASSOCIATION
        assert _severity_label(1.2) == RISK_PREFIX_ASSOCIATION
        assert _severity_label(2.0) == RISK_PREFIX_MODERATE
        assert _severity_label(4.0) == RISK_PREFIX_HIGH

    def test_none_or(self) -> None:
        assert _severity_label(None) == RISK_PREFIX_ASSOCIATION


# ---------------------------------------------------------------------------
# Non-ACMG vocabulary
# ---------------------------------------------------------------------------


class TestSignificanceVocabulary:
    def test_no_acmg_words(self) -> None:
        for label in (SIGNIFICANCE_PATHOGENIC, SIGNIFICANCE_LIKELY_PATHOGENIC):
            assert "pathogenic" not in label.lower()
        assert "vus" not in SIGNIFICANCE_VUS.lower()
        assert "uncertain significance" not in SIGNIFICANCE_VUS.lower()

    def test_labels_still_convey_actionability(self) -> None:
        # The highest tier must still read as actionable/avoid.
        assert "avoid" in SIGNIFICANCE_PATHOGENIC.lower() or "1a" in SIGNIFICANCE_PATHOGENIC.lower()


# ---------------------------------------------------------------------------
# Missing p-value does not pass the significance filter
# ---------------------------------------------------------------------------


class TestNullPValueFilter:
    def test_query_allele_drops_null_p(self) -> None:
        client = GWASClient()
        client._by_allele = {"A*02:01": [_hit(None), _hit(1e-12)]}
        client._loaded = True
        hits = client.query_allele("A*02:01")
        assert len(hits) == 1
        assert hits[0].p_value == 1e-12

    def test_classify_significance_ignores_null_p(self) -> None:
        # A present-but-not-significant (null-p) GWAS hit must not escalate to
        # The risk-factor label.
        label = _classify_significance(_allele(), [_hit(None)], [], [], gwas_resolution="one-field")
        assert label != SIGNIFICANCE_RISK_FACTOR
        assert label == SIGNIFICANCE_VUS

    def test_classify_significance_accepts_real_significant_hit(self) -> None:
        label = _classify_significance(_allele(), [_hit(1e-12)], [], [], gwas_resolution="one-field")
        assert label == SIGNIFICANCE_RISK_FACTOR


# ---------------------------------------------------------------------------
# C4 + C7 — input-quality framing & per-record research-use flag
# ---------------------------------------------------------------------------


def test_disclaimer_documents_uncalibrated_input_quality() -> None:
    """
    The disclaimer must say what the input-quality score is not: a measure
    of genotype accuracy, or a probability. The former "confidence" wording
    invited reading it as a claim about the annotation.
    """
    low = DISCLAIMER.lower()
    assert "uncalibrated" in low
    assert "input-quality" in low
    assert "not a measure of genotype accuracy" in low
    assert "not a posterior probability" in low
    assert "never down-weights an actionable association" in low


def _annotated() -> AnnotatedHLA:
    na = NormalizedAllele(
        allele_name="A*02:01", imgt_accession="HLA00001", protein_group=None,
        hla_class="I", gene="HLA-A", resolution_level=2, is_ambiguous=False,
        is_novel=False, sample_id="S1", source_tool="t1k", source_locus="HLA-A",
        source_resolution="two-field", allele_index=0,
    )
    return AnnotatedHLA(
        normalized_allele=na, gwas_hits=[], pharm_annotations=[], disease_entries=[],
        disease_risk_summary="x", drug_response_summary="y", clinical_significance="No reported risk",
    )


def test_json_carries_research_use_flags(tmp_path: Path) -> None:
    out = tmp_path / "o.json"
    generate_json([_annotated()], out, overwrite=True)
    payload = json.loads(out.read_text(encoding="utf-8"))
    assert payload["metadata"]["research_use_only"] is True
    assert "input_quality_score_definition" in payload["metadata"]
    definition = payload["metadata"]["input_quality_score_definition"].lower()
    assert "heuristic" in definition
    assert "not a measure of genotype accuracy" in definition
    assert "not a posterior probability" in definition
    # Every locus record is individually flagged
    for sample in payload["samples"]:
        for locus in sample["loci"]:
            assert locus["research_use_only"] is True


# ---------------------------------------------------------------------------
# Diplotype / zygosity transparency
# ---------------------------------------------------------------------------


def _annot_named(name: str, index: int) -> AnnotatedHLA:
    na = NormalizedAllele(
        allele_name=name, imgt_accession="HLA00001", protein_group=None,
        hla_class="II", gene="HLA-DQB1", resolution_level=2, is_ambiguous=False,
        is_novel=False, sample_id="S1", source_tool="t1k", source_locus="HLA-DQB1",
        source_resolution="two-field", allele_index=index,
    )
    return AnnotatedHLA(
        normalized_allele=na, gwas_hits=[], pharm_annotations=[], disease_entries=[],
        disease_risk_summary="x", drug_response_summary="y", clinical_significance="No reported risk",
    )


def _row(a1: str, a2) -> GenotypeRow:
    return GenotypeRow(
        sample_id="S1", locus="HLA-DQB1", tool="t1k", resolution="two-field",
        allele1=_annot_named(a1, 0),
        allele2=_annot_named(a2, 1) if a2 is not None else None,
    )


class TestZygosity:
    def test_homozygous(self) -> None:
        assert _zygosity(_row("DQB1*06:02", "DQB1*06:02")) == "homozygous"

    def test_heterozygous(self) -> None:
        assert _zygosity(_row("DQB1*06:02", "DQB1*03:01")) == "heterozygous"

    def test_single_allele_not_assumed_homozygous(self) -> None:
        # The key guarantee: a missing 2nd allele is NOT homozygous.
        assert _zygosity(_row("DQB1*06:02", None)) == "single_allele_reported"


def test_json_has_zygosity_and_diplotype_caveat(tmp_path: Path) -> None:
    out = tmp_path / "z.json"
    generate_json([_annotated()], out, overwrite=True)
    payload = json.loads(out.read_text(encoding="utf-8"))
    caveat = payload["metadata"]["diplotype_caveat"].lower()
    assert "heterodimer" in caveat and "homozygous" in caveat
    assert payload["samples"][0]["loci"][0]["zygosity"] == "single_allele_reported"
