"""
tests.test_resolution_fallback
==============================

Unit tests for the resolution-aware fallback mechanism.

Coverage
--------
- The ``_truncate_to_fields`` and ``_get_resolution_levels`` helpers.
- ``GWASClient.query_allele_with_fallback`` behaviour.
- ``PharmGKBClient`` candidate-form matching + ``matched_form``.
- The ``gwas_annotation_resolution`` column in the TSV report.
"""

from __future__ import annotations

from pathlib import Path
from typing import List

import pytest

from hlante.annotator import (
    AnnotatedHLA,
    AnnotatorClients,
    AnnotatorConfig,
    annotate_genotype,
)
from hlante.types import DiseaseEntry
from hlante.db.gwas import (
    GWASClient,
    GWASHit,
    _colon_group_count,
    _field_label,
    _get_resolution_levels,
    _truncate_to_fields,
)
from hlante.db.pharmgkb import (
    PharmAnnotation,
    PharmGKBClient,
    _candidate_allele_forms,
)
from hlante.normalizer import NormalizedAllele
from hlante.reporter import TSV_COLUMNS, generate_tsv


FIXTURES_DIR: Path = Path(__file__).parent / "fixtures"
GWAS_FIXTURE_DIR: Path = FIXTURES_DIR / "gwas"
PHARMGKB_FIXTURE_DIR: Path = FIXTURES_DIR / "pharmgkb"


# ---------------------------------------------------------------------------
# _truncate_to_fields + _get_resolution_levels + _colon_group_count
# ---------------------------------------------------------------------------


class TestTruncation:
    """
    Correctness of the pure helpers.
    """

    @pytest.mark.parametrize(
        "allele,n_fields,expected",
        [
            ("A*02:01:01:01", 2, "A*02:01"),
            ("A*02:01:01:01", 1, "A*02"),
            ("A*02:01:01:01", 4, "A*02:01:01:01"),  # already ≤4 groups
            ("DRB1*03:01:01", 4, "DRB1*03:01:01"),  # already ≤4, unchanged
            ("DRB1*03:01:01", 2, "DRB1*03:01"),
            ("B*57:01G", 2, "B*57:01G"),  # G-group suffix preserved
            ("B*57:01P", 2, "B*57:01P"),  # P-group suffix preserved
            ("A*02:01N", 1, "A*02N"),     # null-allele suffix preserved
            ("HLA-A*02:01:01:01", 2, "A*02:01"),  # HLA- prefix stripped
        ],
    )
    def test_truncate_cases(self, allele: str, n_fields: int, expected: str) -> None:
        assert _truncate_to_fields(allele, n_fields) == expected

    @pytest.mark.parametrize(
        "allele,expected",
        [
            ("A*02:01:01:01", [4, 3, 2, 1]),
            ("A*02:01:01", [3, 2, 1]),
            ("A*02:01", [2, 1]),
            ("A*02", [1]),
            ("B*57:01G", [2, 1]),  # suffix is not counted
            ("HLA-DRB1*03:01:01", [3, 2, 1]),
        ],
    )
    def test_resolution_levels(self, allele: str, expected: List[int]) -> None:
        assert _get_resolution_levels(allele) == expected

    @pytest.mark.parametrize(
        "allele,expected",
        [
            ("A*02:01:01:01", 4),
            ("A*02:01", 2),
            ("A*02", 1),
            ("B*57:01G", 2),
        ],
    )
    def test_colon_group_count(self, allele: str, expected: int) -> None:
        assert _colon_group_count(allele) == expected

    @pytest.mark.parametrize(
        "colon_groups,expected_label",
        [
            (1, "2-field"),
            (2, "4-field"),
            (3, "6-field"),
            (4, "8-field"),
        ],
    )
    def test_field_label_mapping(
        self, colon_groups: int, expected_label: str
    ) -> None:
        assert _field_label(colon_groups) == expected_label


# ---------------------------------------------------------------------------
# GWASClient.query_allele_with_fallback
# ---------------------------------------------------------------------------


class TestGWASFallback:
    """
    Behaviour of ``GWASClient.query_allele_with_fallback`` — the
    fixture TSV only contains 4-field (``A*02:01``, ``B*57:01``)
    records.
    """

    def test_6field_falls_back_to_4field(self) -> None:
        """
        The fixture has ``B*57:01`` but not ``B*57:01:01:01``. The
        fallback should descend from 8→6→4 and return the 4-field hit.
        """
        client = GWASClient(local_dir=GWAS_FIXTURE_DIR)
        hits, resolution = client.query_allele_with_fallback("B*57:01:01:01")
        assert hits, "Fallback should find a hit at B*57:01."
        assert resolution == "4-field"

    def test_exact_match_returns_same_resolution(self) -> None:
        """
        An exact-resolution query should return without falling back.
        """
        client = GWASClient(local_dir=GWAS_FIXTURE_DIR)
        hits, resolution = client.query_allele_with_fallback("B*57:01")
        assert hits
        assert resolution == "4-field"

    def test_no_match_at_any_level(self) -> None:
        """
        Without a match at any level the function must return ``"none"``.
        """
        client = GWASClient(local_dir=GWAS_FIXTURE_DIR)
        hits, resolution = client.query_allele_with_fallback("Z*99:99:99:99")
        assert hits == []
        assert resolution == "none"

    def test_min_resolution_limits_descent(self) -> None:
        """
        ``min_resolution=4`` must prevent descent to the 2-field
        (``A*02``) level.
        """
        client = GWASClient(local_dir=GWAS_FIXTURE_DIR)
        # The fixture has no entry for ``A*02`` or ``A*02:01`` → none.
        hits, resolution = client.query_allele_with_fallback(
            "A*02:01:01:01", min_resolution=4
        )
        assert resolution == "none"

    def test_hla_prefix_equivalent(self) -> None:
        """
        ``HLA-`` prefix must not affect the lookup path.
        """
        client = GWASClient(local_dir=GWAS_FIXTURE_DIR)
        a = client.query_allele_with_fallback("HLA-B*57:01")
        b = client.query_allele_with_fallback("B*57:01")
        assert a[1] == b[1]
        assert len(a[0]) == len(b[0])


# ---------------------------------------------------------------------------
# PharmGKB candidate-form matching
# ---------------------------------------------------------------------------


class TestPharmGKBCandidateForms:
    """
    Behaviour of :func:`_candidate_allele_forms` and the multi-form
    matching in :meth:`PharmGKBClient.query_allele`.
    """

    def test_candidate_forms_for_4field(self) -> None:
        """
        ``A*02:01`` must produce full, 2-field, HLA-prefixed, and
        gene-prefixless variants.
        """
        forms = _candidate_allele_forms("A*02:01")
        assert "A*02:01" in forms
        assert "A*02" in forms
        assert "HLA-A*02:01" in forms
        assert "HLA-A*02" in forms
        assert "*02:01" in forms
        assert "*02" in forms
        # The original form must come first (used for fallback order).
        assert forms[0] == "A*02:01"

    def test_candidate_forms_preserve_g_suffix(self) -> None:
        """
        The G-group suffix must be preserved in every variant.
        """
        forms = _candidate_allele_forms("B*57:01G")
        assert "B*57:01G" in forms
        assert "B*57G" in forms

    def test_query_with_6field_input_matches_4field_record(
        self, tmp_path: Path
    ) -> None:
        """
        PharmGKB stores 4-field alleles. A 6-field query must match a
        truncated variant, and ``matched_form`` must be populated.
        """
        client = PharmGKBClient(local_dir=PHARMGKB_FIXTURE_DIR)
        hits = client.query_allele("B*57:01:01")
        assert hits, "A 6-field query should match the 4-field record."
        drugs = {h.drug for h in hits}
        assert "abacavir" in drugs
        # matched_form must be the truncated 4-field form.
        assert any(h.matched_form == "B*57:01" for h in hits)

    def test_exact_match_records_exact_form(self) -> None:
        """
        ``matched_form`` must equal the query when an exact match
        exists.
        """
        client = PharmGKBClient(local_dir=PHARMGKB_FIXTURE_DIR)
        hits = client.query_allele("B*57:01")
        assert hits
        assert any(h.matched_form == "B*57:01" for h in hits)

    def test_dedupe_across_candidate_forms(self) -> None:
        """
        The same record should not be returned twice even when multiple
        candidate forms match it.
        """
        client = PharmGKBClient(local_dir=PHARMGKB_FIXTURE_DIR)
        hits = client.query_allele("B*57:01")
        ids = [(h.annotation_id, h.drug) for h in hits]
        assert len(ids) == len(set(ids)), "Duplicates detected after dedupe."


# ---------------------------------------------------------------------------
# New TSV column in the report
# ---------------------------------------------------------------------------


class TestResolutionColumnInReport:
    """
    The ``gwas_annotation_resolution`` column must sit in the right
    position and carry the expected values.
    """

    def _make_annotated(
        self,
        allele_name: str,
        gwas_hits,
        gwas_resolution: str,
    ) -> AnnotatedHLA:
        norm = NormalizedAllele(
            allele_name=allele_name,
            imgt_accession="HLA00001",
            protein_group=None,
            hla_class="I",
            gene="HLA-B",
            resolution_level=4,
            is_ambiguous=False,
            is_novel=False,
            sample_id="S1",
            source_tool="t1k",
            source_locus="HLA-B",
            source_resolution="4-field",
            allele_index=0,
        )
        return AnnotatedHLA(
            normalized_allele=norm,
            gwas_hits=gwas_hits,
            pharm_annotations=[],
            disease_entries=[],
            disease_risk_summary="",
            drug_response_summary="",
            clinical_significance="Benign",
            gwas_resolution_used=gwas_resolution,
        )

    def test_column_exists_in_schema(self) -> None:
        """
        The new column must appear immediately after ``gwas_pmids``.
        """
        idx_pmid = TSV_COLUMNS.index("gwas_pmids")
        idx_res = TSV_COLUMNS.index("gwas_annotation_resolution")
        assert idx_res == idx_pmid + 1

    def test_tsv_writes_resolution_label(self, tmp_path: Path) -> None:
        """
        The TSV row must carry the correct resolution label.
        """
        annot = self._make_annotated(
            "B*57:01",
            [GWASHit(
                trait="Abacavir HSR",
                p_value=1e-10,
                odds_ratio=4.2,
                pmid="18322448",
                study_accession="GCST000001",
                allele="B*57:01",
            )],
            gwas_resolution="4-field",
        )
        out = tmp_path / "r.tsv"
        generate_tsv([annot], out)
        text = out.read_text(encoding="utf-8")
        data_line = next(
            ln for ln in text.splitlines()
            if not ln.startswith("#") and "S1" in ln
        )
        cells = data_line.split("\t")
        idx = TSV_COLUMNS.index("gwas_annotation_resolution")
        assert cells[idx] == "4-field"

    def test_tsv_resolution_none_when_no_hits(self, tmp_path: Path) -> None:
        """
        When no hit is found the cell must be ``"none"`` — explicitly
        distinct from the ``NA`` placeholder used for absent data.
        """
        annot = self._make_annotated(
            "A*99:99",
            gwas_hits=[],
            gwas_resolution="none",
        )
        out = tmp_path / "r.tsv"
        generate_tsv([annot], out)
        text = out.read_text(encoding="utf-8")
        data_line = next(
            ln for ln in text.splitlines()
            if not ln.startswith("#") and "S1" in ln
        )
        cells = data_line.split("\t")
        idx = TSV_COLUMNS.index("gwas_annotation_resolution")
        assert cells[idx] == "none"


# ---------------------------------------------------------------------------
# Verify the annotator propagates the fallback resolution
# ---------------------------------------------------------------------------


class TestAnnotatorPropagatesResolution:
    """
    :func:`annotate_genotype` must carry the fallback result into
    ``AnnotatedHLA.gwas_resolution_used``.
    """

    def test_annotator_picks_up_fallback_resolution(self) -> None:
        """
        Using the real ``GWASClient`` + fixture, a 6-field query must
        fall back and produce ``gwas_resolution_used == "4-field"``.
        """
        norm = NormalizedAllele(
            allele_name="B*57:01:01:01",
            imgt_accession="HLA00009",
            protein_group=None,
            hla_class="I",
            gene="HLA-B",
            resolution_level=8,
            is_ambiguous=False,
            is_novel=False,
            sample_id="S1",
            source_tool="t1k",
            source_locus="HLA-B",
            source_resolution="8-field",
            allele_index=0,
        )
        gwas_client = GWASClient(local_dir=GWAS_FIXTURE_DIR)
        config = AnnotatorConfig(offline=True,
                                 enable_pharmgkb=False, enable_afnd=False)
        clients = AnnotatorClients(gwas=gwas_client, pharmgkb=None)
        result = annotate_genotype([norm], config, clients=clients)[0]
        assert result.gwas_resolution_used == "4-field"
        assert result.gwas_hits, "Fallback should have produced hits."
