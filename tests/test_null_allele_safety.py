"""
tests.test_null_allele_safety
=============================

Clinical-safety tests for IPD-IMGT/HLA *expression suffix* handling
.

A null allele (``N`` suffix, e.g. ``A*24:09N``) carries a change that
abolishes cell-surface expression of the antigen. Because HLA disease and
drug-hypersensitivity associations depend on surface peptide presentation,
those associations must NOT be attached to a null allele — doing so would
emit a false clinical alert (e.g. an abacavir contraindication for a
non-expressed ``B*57:01N``). These tests pin that behaviour:

1. The normalizer captures the expression suffix and exposes ``is_null`` /
   ``is_low_or_aberrant_expression``.
2. The HLA-HD "gene not present" placeholder (``DRB3*00:00`` /
   ``DRB3*00:00N``) normalizes to ``None`` (NA), as the Supplementary
   Methods S1.3 claim requires.
3. ``annotate_genotype`` never even *queries* the disease/drug databases
   for a null allele, and labels it ``SIGNIFICANCE_NULL_ALLELE`` at a
   non-HIGH tier — while an otherwise identical expressed allele is
   annotated normally.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional

import pytest

from hlante.annotator import (
    SIGNIFICANCE_NULL_ALLELE,
    SIGNIFICANCE_PATHOGENIC,
    AnnotatorClients,
    AnnotatorConfig,
    annotate_genotype,
)
from hlante.db.gwas import GWASHit
from hlante.db.pharmgkb import PharmAnnotation
from hlante.normalizer import NormalizedAllele, load_imgt_db, normalize_allele
from hlante.types import DiseaseEntry

IMGT_MINI_DIR: Path = Path(__file__).parent / "fixtures" / "imgt_mini"


# ---------------------------------------------------------------------------
# Normalizer-level: expression suffix capture + not-present placeholder
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def imgt_db() -> Dict[str, object]:
    return load_imgt_db(IMGT_MINI_DIR)


class TestExpressionSuffixCapture:
    def test_null_allele_flagged(self, imgt_db: Dict[str, object]) -> None:
        # A*01:01:01:02N is a real null allele present in the mini fixture.
        norm = normalize_allele("A*01:01:01:02N", imgt_db)
        assert norm is not None
        assert norm.expression_suffix == "N"
        assert norm.is_null is True
        assert norm.is_low_or_aberrant_expression is False
        # It is a recognised IPD-IMGT/HLA allele, not novel.
        assert norm.is_novel is False
        assert norm.imgt_accession == "HLA00002"

    def test_low_expression_allele_flagged(self, imgt_db: Dict[str, object]) -> None:
        norm = normalize_allele("A*02:01:01:01L", imgt_db)
        assert norm is not None
        assert norm.expression_suffix == "L"
        assert norm.is_null is False
        assert norm.is_low_or_aberrant_expression is True

    def test_expressed_allele_has_no_suffix(self, imgt_db: Dict[str, object]) -> None:
        norm = normalize_allele("A*24:02:01:01", imgt_db)
        assert norm is not None
        assert norm.expression_suffix is None
        assert norm.is_null is False
        assert norm.is_low_or_aberrant_expression is False

    def test_group_suffix_is_not_an_expression_suffix(
        self, imgt_db: Dict[str, object]
    ) -> None:
        # A trailing G/P denotes a sequence group, not an expression state.
        norm = normalize_allele("A*02:01:01G", imgt_db)
        assert norm is not None
        assert norm.expression_suffix is None
        assert norm.is_null is False

    @pytest.mark.parametrize(
        "placeholder",
        ["DRB3*00:00", "DRB3*00:00N", "DRB4*00:00:00", "DRB5*00:00N"],
    )
    def test_gene_not_present_placeholder_returns_none(
        self, placeholder: str, imgt_db: Dict[str, object]
    ) -> None:
        # HLA-HD emits an all-zero placeholder for absent DRB3/4/5; it must
        # be treated as NA (Supplementary Methods S1.3), not annotated.
        assert normalize_allele(placeholder, imgt_db) is None


# ---------------------------------------------------------------------------
# Annotator-level: null alleles are not annotated with disease/drug risk
# ---------------------------------------------------------------------------


@dataclass
class _Freq:
    frequency: float = 0.05
    population: str = "EUR"
    source: str = "AFND"
    sample_size: int = 1000
    is_estimated: bool = False


class _RecordingGWAS:
    """Returns a hit for every allele form and records each query."""

    def __init__(self) -> None:
        self.calls: List[str] = []

    def query_allele_with_fallback(self, allele: str, min_resolution: int = 2):
        self.calls.append(allele)
        hit = GWASHit(
            trait="Some disease",
            p_value=1e-12,
            odds_ratio=5.0,
            pmid="99999999",
            study_accession="GCST999999",
            allele=allele,
        )
        return [hit], "two-field"


class _RecordingPharmGKB:
    def __init__(self) -> None:
        self.calls: List[str] = []

    def query_allele(self, allele: str) -> List[PharmAnnotation]:
        self.calls.append(allele)
        return [
            PharmAnnotation(
                drug="abacavir",
                phenotype="HSR risk",
                evidence_level="1A",
                pmid=["18256392"],
                allele=allele,
                annotation_id="CA000001",
                gene="HLA-B",
                cpic_url="https://cpicpgx.org/guidelines/",
            )
        ]


class _RecordingCurated:
    def __init__(self) -> None:
        self.calls: List[str] = []

    def query_allele(self, allele: str) -> List[DiseaseEntry]:
        self.calls.append(allele)
        return [
            DiseaseEntry(
                variation_id="111",
                significance="Pathogenic",
                condition="Abacavir hypersensitivity",
                review_status="curated (HLAnte built-in)",
                allele=allele,
                pmid=["18322448"],
            )
        ]


class _StubAFND:
    def get_frequency_with_fallback(self, allele: str) -> Optional[_Freq]:
        return _Freq()


def _make_allele(allele_name: str, expression_suffix: Optional[str]) -> NormalizedAllele:
    return NormalizedAllele(
        allele_name=allele_name,
        imgt_accession="HLA00001",
        protein_group=None,
        hla_class="I",
        gene="HLA-B",
        resolution_level=2,
        is_ambiguous=False,
        is_novel=False,
        expression_suffix=expression_suffix,
    )


def _clients() -> tuple[AnnotatorClients, _RecordingGWAS, _RecordingPharmGKB, _RecordingCurated]:
    gwas, pharm, curated = _RecordingGWAS(), _RecordingPharmGKB(), _RecordingCurated()
    clients = AnnotatorClients(
        gwas=gwas, pharmgkb=pharm, curated=curated, afnd=_StubAFND()
    )
    return clients, gwas, pharm, curated


class TestNullAlleleSuppression:
    def test_null_allele_disease_drug_suppressed(self) -> None:
        clients, gwas, pharm, curated = _clients()
        allele = _make_allele("B*57:01N", expression_suffix="N")
        result = annotate_genotype([allele], AnnotatorConfig(offline=True), clients=clients)[0]

        # The databases must not even be queried for a null allele.
        assert gwas.calls == []
        assert pharm.calls == []
        assert curated.calls == []

        # No disease/drug evidence is attached.
        assert result.gwas_hits == []
        assert result.pharm_annotations == []
        assert result.disease_entries == []

        # It is labelled as a null allele, never as a pathogenic/risk hit.
        assert result.clinical_significance == SIGNIFICANCE_NULL_ALLELE
        assert "Null allele" in result.disease_risk_summary
        assert "Null allele" in result.drug_response_summary

        # A non-expressed allele is never presented in the detailed tier.
        assert result.input_quality_tier != "detailed"
        assert "null_allele" in result.input_quality_rationale

        # Allele frequency is still meaningful and is retained.
        assert result.allele_frequency == pytest.approx(0.05)

    def test_expressed_allele_still_annotated(self) -> None:
        """Control: the same allele without the N suffix is annotated."""
        clients, gwas, pharm, curated = _clients()
        allele = _make_allele("B*57:01", expression_suffix=None)
        result = annotate_genotype([allele], AnnotatorConfig(offline=True), clients=clients)[0]

        assert gwas.calls and pharm.calls and curated.calls
        assert len(result.gwas_hits) == 1
        assert len(result.pharm_annotations) == 1
        assert len(result.disease_entries) == 1
        assert result.clinical_significance == SIGNIFICANCE_PATHOGENIC
