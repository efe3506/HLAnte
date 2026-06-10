"""
tests.test_coverage_visibility
==============================

Tests that *absence of typing* is made explicit.

A clinically-actionable locus that produced no genotype row is
indeterminate ("not typed"), not negative. If a missing row is silently
read as "no risk", a clinician could conclude a patient is clear of, say,
``HLA-B*57:01`` when ``HLA-B`` was simply never typed. These tests pin:

1. The indeterminate label is unambiguous ("Not assessed — insufficient
   coverage").
2. ``_loci_not_typed`` reports the actionable loci a sample lacks, scoped
   to the typing tool's capability (Class I-only tools do not flag Class
   II).
3. The JSON and Markdown reports surface the not-typed loci.
"""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

from hlante.annotator import (
    SIGNIFICANCE_BENIGN,
    SIGNIFICANCE_BENIGN_LIMITED,
    AnnotatedHLA,
)
from hlante.normalizer import NormalizedAllele
from hlante.reporter import (
    _expected_actionable_loci,
    _loci_not_typed,
    generate_json,
    generate_markdown_report,
)


# ---------------------------------------------------------------------------
# Label clarity
# ---------------------------------------------------------------------------


def test_limited_data_label_is_unambiguous() -> None:
    # The indeterminate case must not be phrased as a kind of "no risk".
    assert SIGNIFICANCE_BENIGN_LIMITED == "Not assessed — insufficient coverage"
    assert SIGNIFICANCE_BENIGN != SIGNIFICANCE_BENIGN_LIMITED
    assert "no reported risk" not in SIGNIFICANCE_BENIGN_LIMITED.lower()


# ---------------------------------------------------------------------------
# _loci_not_typed / _expected_actionable_loci
# ---------------------------------------------------------------------------


def _row(locus: str, tool: str = "arcashla") -> SimpleNamespace:
    return SimpleNamespace(locus=locus, tool=tool)


class TestLociNotTyped:
    def test_missing_class_ii_flagged_for_full_tool(self) -> None:
        rows = [_row("HLA-A"), _row("HLA-B"), _row("HLA-C")]
        missing = _loci_not_typed(rows)
        assert "HLA-DRB1" in missing
        assert "HLA-DQB1" in missing
        assert "HLA-DQA1" in missing
        assert "HLA-DPB1" in missing
        # Class I that *was* typed is not flagged.
        assert "HLA-A" not in missing

    def test_optitype_class_i_only_does_not_flag_class_ii(self) -> None:
        rows = [_row("HLA-A", "optitype"), _row("HLA-B", "optitype"), _row("HLA-C", "optitype")]
        assert _loci_not_typed(rows) == []

    def test_single_locus_typed_flags_the_rest(self) -> None:
        rows = [_row("HLA-A")]
        missing = _loci_not_typed(rows)
        assert "HLA-B" in missing and "HLA-C" in missing
        assert "HLA-DRB1" in missing

    def test_case_insensitive_locus_match(self) -> None:
        rows = [_row("hla-a"), _row("HLA-B"), _row("HLA-C")]
        assert "HLA-A" not in _loci_not_typed(rows)

    def test_empty_rows(self) -> None:
        assert _loci_not_typed([]) == []

    def test_expected_panel_scopes_by_tool(self) -> None:
        assert _expected_actionable_loci("optitype") == ("HLA-A", "HLA-B", "HLA-C")
        assert "HLA-DQB1" in _expected_actionable_loci("t1k")


# ---------------------------------------------------------------------------
# End-to-end: JSON + Markdown surface the not-typed loci
# ---------------------------------------------------------------------------


def _make_annotated(locus: str, tool: str = "arcashla") -> AnnotatedHLA:
    na = NormalizedAllele(
        allele_name="A*02:01",
        imgt_accession="HLA00001",
        protein_group=None,
        hla_class="I",
        gene=locus,
        resolution_level=4,
        is_ambiguous=False,
        is_novel=False,
        sample_id="S1",
        source_tool=tool,
        source_locus=locus,
        source_resolution="4-field",
        allele_index=0,
    )
    return AnnotatedHLA(
        normalized_allele=na,
        gwas_hits=[],
        pharm_annotations=[],
        disease_entries=[],
        disease_risk_summary="No disease association reported",
        drug_response_summary="No drug response reported",
        clinical_significance=SIGNIFICANCE_BENIGN,
    )


def test_json_lists_not_typed_loci(tmp_path: Path) -> None:
    import json

    annotated = [_make_annotated("HLA-A")]  # only HLA-A typed
    out = tmp_path / "out.json"
    generate_json(annotated, out, overwrite=True)
    payload = json.loads(out.read_text(encoding="utf-8"))
    not_typed = payload["samples"][0]["actionable_loci_not_typed"]
    assert "HLA-B" in not_typed
    assert "HLA-DQB1" in not_typed


def test_markdown_warns_about_not_typed_loci(tmp_path: Path) -> None:
    annotated = [_make_annotated("HLA-A")]
    out = tmp_path / "out.md"
    generate_markdown_report(annotated, out, overwrite=True)
    text = out.read_text(encoding="utf-8")
    assert "Loci not typed" in text
    assert "HLA-DQB1" in text
    assert "indeterminate" in text.lower()
