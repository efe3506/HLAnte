"""
tests.test_reporter
===================

Unit tests for :mod:`hlante.reporter`.

Coverage
--------
- TSV / CSV: column order, pipe joining, float formatting, NA
  handling, metadata header.
- Markdown: per-sample sections, risk / drug blocks, auto-generated
  clinical note, disclaimer, emoji.
- JSON: valid and nested structure, metadata.
- ``generate_all``: all three formats produced together.
- Overwrite guard.
- The tqdm wrapper does not break for large cohorts.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import List, Optional

import pytest

from hlante import __version__
from hlante.annotator import (
    SIGNIFICANCE_BENIGN,
    SIGNIFICANCE_PATHOGENIC,
    SIGNIFICANCE_RISK_FACTOR,
    AnnotatedHLA,
)
from hlante.types import DiseaseEntry
from hlante.db.gwas import GWASHit
from hlante.db.pharmgkb import PharmAnnotation
from hlante.normalizer import NormalizedAllele
from hlante.reporter import (
    DISCLAIMER,
    NA,
    TSV_COLUMNS,
    OutputFileExistsError,
    ReportContext,
    generate_all,
    generate_csv,
    generate_json,
    generate_markdown_report,
    generate_tsv,
)


# ---------------------------------------------------------------------------
# Factory helpers
# ---------------------------------------------------------------------------


def _normalized(
    allele_name: str,
    *,
    sample_id: str = "S1",
    tool: str = "t1k",
    locus: str = "HLA-B",
    resolution: str = "two-field",
    allele_index: int = 0,
    imgt_accession: Optional[str] = "HLA00001",
    protein_group: Optional[str] = "B*57:01G",
    hla_class: str = "I",
    resolution_level: int = 2,
    is_ambiguous: bool = False,
    is_novel: bool = False,
) -> NormalizedAllele:
    return NormalizedAllele(
        allele_name=allele_name,
        imgt_accession=imgt_accession,
        protein_group=protein_group,
        hla_class=hla_class,
        gene=locus,
        resolution_level=resolution_level,
        is_ambiguous=is_ambiguous,
        is_novel=is_novel,
        sample_id=sample_id,
        source_tool=tool,
        source_locus=locus,
        source_resolution=resolution,
        allele_index=allele_index,
    )


def _annotated(
    normalized: NormalizedAllele,
    *,
    gwas: Optional[List[GWASHit]] = None,
    pharm: Optional[List[PharmAnnotation]] = None,
    entries: Optional[List[DiseaseEntry]] = None,
    disease: str = "No disease association reported",
    drug: str = "No drug response reported",
    significance: str = SIGNIFICANCE_BENIGN,
) -> AnnotatedHLA:
    return AnnotatedHLA(
        normalized_allele=normalized,
        gwas_hits=gwas or [],
        pharm_annotations=pharm or [],
        disease_entries=entries or [],
        disease_risk_summary=disease,
        drug_response_summary=drug,
        clinical_significance=significance,
    )


@pytest.fixture()
def sample_dataset() -> List[AnnotatedHLA]:
    """
    Two samples with several HLA loci and rich annotations.
    """
    # Sample S1 — HLA-B*57:01 + B*07:02
    b57 = _annotated(
        _normalized(
            "B*57:01",
            sample_id="S1",
            allele_index=0,
            locus="HLA-B",
            tool="t1k",
            resolution="two-field",
        ),
        gwas=[
            GWASHit(
                trait="Abacavir hypersensitivity",
                p_value=1.2e-10,
                odds_ratio=4.2,
                pmid="18322448",
                study_accession="GCST000001",
                allele="B*57:01",
            )
        ],
        pharm=[
            PharmAnnotation(
                drug="abacavir",
                phenotype="HSR risk",
                evidence_level="1A",
                pmid=["18256392"],
                allele="B*57:01",
                annotation_id="CA000001",
                gene="HLA-B",
                cpic_url="https://cpicpgx.org/guidelines/abacavir",
                pharmgkb_url="https://www.pharmgkb.org/clinicalAnnotation/CA000001",
            )
        ],
        entries=[
            DiseaseEntry(
                variation_id="111",
                significance="Pathogenic",
                condition="Stevens-Johnson syndrome",
                review_status="reviewed by expert panel",
                allele="B*57:01",
                pmid=["18322448"],
            )
        ],
        disease="Strong association: Abacavir hypersensitivity (OR=4.20)",
        drug="abacavir HSR risk: 1A evidence",
        significance=SIGNIFICANCE_PATHOGENIC,
    )
    b0702 = _annotated(
        _normalized(
            "B*07:02",
            sample_id="S1",
            allele_index=1,
            locus="HLA-B",
            tool="t1k",
            imgt_accession="HLA00009",
            protein_group="B*07:02:01G",
        ),
        significance=SIGNIFICANCE_BENIGN,
    )

    # Sample S2 — DRB1*15:01 single allele (homozygous)
    drb1 = _annotated(
        _normalized(
            "DRB1*15:01",
            sample_id="S2",
            allele_index=0,
            locus="HLA-DRB1",
            tool="hlahd",
            resolution="two-field",
            hla_class="II",
            imgt_accession="HLA00100",
            protein_group="DRB1*15:01:01G",
        ),
        gwas=[
            GWASHit(
                trait="Multiple sclerosis",
                p_value=1e-20,
                odds_ratio=3.1,
                pmid="17660530",
                study_accession="GCST000010",
                allele="DRB1*15:01",
            )
        ],
        disease="Strong association: Multiple sclerosis (OR=3.10)",
        drug="No drug response reported",
        significance=SIGNIFICANCE_RISK_FACTOR,
    )

    return [b57, b0702, drb1]


# ---------------------------------------------------------------------------
# TSV
# ---------------------------------------------------------------------------


class TestGenerateTSV:
    """
    Behaviour of :func:`generate_tsv`.
    """

    def _read_tsv(self, path: Path):
        """
        Skip ``#`` metadata lines and parse header + data rows.
        """
        lines = path.read_text(encoding="utf-8").splitlines()
        meta = [ln for ln in lines if ln.startswith("#")]
        data_lines = [ln for ln in lines if not ln.startswith("#")]
        reader = csv.reader(data_lines, delimiter="\t")
        rows = list(reader)
        return meta, rows[0], rows[1:]

    def test_column_order(
        self, sample_dataset: List[AnnotatedHLA], tmp_path: Path
    ) -> None:
        """
        The header row must match :data:`TSV_COLUMNS` exactly.
        """
        out = tmp_path / "report.tsv"
        generate_tsv(sample_dataset, out)
        _, header, _ = self._read_tsv(out)
        assert tuple(header) == TSV_COLUMNS

    def test_rows_grouped_by_locus(
        self, sample_dataset: List[AnnotatedHLA], tmp_path: Path
    ) -> None:
        """
        Same ``(sample, locus)`` records must merge into one row; an
        unpaired allele must render ``NA``.
        """
        out = tmp_path / "report.tsv"
        generate_tsv(sample_dataset, out)
        _, header, data = self._read_tsv(out)
        assert len(data) == 2  # S1/HLA-B and S2/HLA-DRB1

        idx = {col: i for i, col in enumerate(header)}
        s1_row = next(r for r in data if r[idx["sample_id"]] == "S1")
        s2_row = next(r for r in data if r[idx["sample_id"]] == "S2")

        assert s1_row[idx["locus"]] == "HLA-B"
        assert s1_row[idx["allele1"]] == "B*57:01"
        assert s1_row[idx["allele2"]] == "B*07:02"
        assert s1_row[idx["tool"]] == "t1k"
        assert s1_row[idx["hla_class"]] == "I"

        assert s2_row[idx["allele2"]] == NA
        assert s2_row[idx["hla_class"]] == "II"

    def test_pipe_joined_annotations(
        self, sample_dataset: List[AnnotatedHLA], tmp_path: Path
    ) -> None:
        """
        Multiple allele-specific values must join with a pipe.
        """
        out = tmp_path / "report.tsv"
        generate_tsv(sample_dataset, out)
        _, header, data = self._read_tsv(out)
        idx = {col: i for i, col in enumerate(header)}
        s1_row = next(r for r in data if r[idx["sample_id"]] == "S1")

        assert s1_row[idx["imgt_accession"]] == "HLA00001|HLA00009"
        assert s1_row[idx["protein_group"]] == "B*57:01G|B*07:02:01G"
        assert s1_row[idx["gwas_traits"]] == "Abacavir hypersensitivity"
        assert s1_row[idx["pharm_drugs"]] == "abacavir"
        assert s1_row[idx["pharm_evidence"]] == "1A"
        assert s1_row[idx["pharm_pmids"]] == "18256392"
        # Clinical_significance is now composed of the
        # Evidence-strength labels from hlante.annotator — imported as
        # Constants to avoid hard-coding display text here.
        assert s1_row[idx["clinical_significance"]] == (
            f"{SIGNIFICANCE_PATHOGENIC}|{SIGNIFICANCE_BENIGN}"
        )

    def test_float_fixed_decimals(
        self, sample_dataset: List[AnnotatedHLA], tmp_path: Path
    ) -> None:
        """
        Effect sizes must be written as fixed-point with 4 decimals, no
        scientific notation.
        """
        out = tmp_path / "report.tsv"
        generate_tsv(sample_dataset, out)
        _, header, data = self._read_tsv(out)
        idx = {col: i for i, col in enumerate(header)}
        s2_row = next(r for r in data if r[idx["sample_id"]] == "S2")
        # OR=3.1 → "3.1000"
        assert s2_row[idx["gwas_odds_ratios"]] == "3.1000"
        assert "e" not in s2_row[idx["gwas_odds_ratios"]].lower()

    def test_per_allele_columns_keep_their_slot(
        self, sample_dataset: List[AnnotatedHLA], tmp_path: Path
    ) -> None:
        """
        A per-allele value must stay in its own slot when the other allele
        has none.

        ``imgt_accession``, ``hla_serotype`` and ``protein_group`` used to
        drop missing entries rather than reserve them, so a row where only
        allele2 carried a value printed that value alone — and a reader
        following the positional convention attributed it to allele1. Real
        two-field calls hit this constantly: an exactly-matched allele has an
        accession but no G-group, and a prefix-matched one has the reverse.
        """
        dataset = [
            _annotated(
                _normalized(
                    "DPA1*01:03:01",
                    sample_id="SLOT",
                    locus="HLA-DPA1",
                    allele_index=0,
                    hla_class="II",
                    imgt_accession=None,
                    protein_group="DPA1*01:03:01G",
                )
            ),
            _annotated(
                _normalized(
                    "DPA1*01:03:04",
                    sample_id="SLOT",
                    locus="HLA-DPA1",
                    allele_index=1,
                    hla_class="II",
                    imgt_accession="HLA03224",
                    protein_group=None,
                )
            ),
        ]
        out = tmp_path / "report.tsv"
        generate_tsv(dataset, out)
        _, header, data = self._read_tsv(out)
        idx = {col: i for i, col in enumerate(header)}
        row = next(r for r in data if r[idx["sample_id"]] == "SLOT")

        # The accession belongs to allele2, the protein group to allele1.
        assert row[idx["imgt_accession"]] == "NA|HLA03224"
        assert row[idx["protein_group"]] == "DPA1*01:03:01G|NA"

    def test_pvalue_survives_the_tsv(
        self, sample_dataset: List[AnnotatedHLA], tmp_path: Path
    ) -> None:
        """
        p-values must reach the TSV intact.

        Only associations at or below 5e-8 are retained, so fixed-point
        formatting collapsed every p-value in every report to ``0.0000``.
        Scientific notation is therefore not a preference here — it is the
        only format that can carry the value at all.
        """
        out = tmp_path / "report.tsv"
        generate_tsv(sample_dataset, out)
        _, header, data = self._read_tsv(out)
        idx = {col: i for i, col in enumerate(header)}
        s2_row = next(r for r in data if r[idx["sample_id"]] == "S2")
        cell = s2_row[idx["gwas_p_values"]]
        assert cell == "1.00e-20", cell
        assert float(cell) == pytest.approx(1e-20)

    def test_na_for_missing_values(
        self, sample_dataset: List[AnnotatedHLA], tmp_path: Path
    ) -> None:
        """
        Cells with no source data must render as NA.
        """
        out = tmp_path / "report.tsv"
        generate_tsv(sample_dataset, out)
        _, header, data = self._read_tsv(out)
        idx = {col: i for i, col in enumerate(header)}
        s2_row = next(r for r in data if r[idx["sample_id"]] == "S2")
        assert s2_row[idx["pharm_drugs"]] == NA

    def test_metadata_header_written(
        self, sample_dataset: List[AnnotatedHLA], tmp_path: Path
    ) -> None:
        """
        The HLAnte version, timestamp, DB versions, and disclaimer must
        appear in the header.
        """
        out = tmp_path / "report.tsv"
        ctx = ReportContext(db_versions={"imgt": "3.55.0", "gwas": "2026-01"})
        generate_tsv(sample_dataset, out, context=ctx)
        text = out.read_text(encoding="utf-8")
        assert f"# HLAnte {__version__}" in text
        assert "# imgt version: 3.55.0" in text
        assert "# gwas version: 2026-01" in text
        assert "# Disclaimer:" in text
        assert DISCLAIMER in text

    def test_overwrite_guard(
        self, sample_dataset: List[AnnotatedHLA], tmp_path: Path
    ) -> None:
        """
        Writing to an existing file must raise unless overwrite=True.
        """
        out = tmp_path / "report.tsv"
        generate_tsv(sample_dataset, out)
        with pytest.raises(OutputFileExistsError, match="already exists"):
            generate_tsv(sample_dataset, out)
        # Overwrite=True succeeds
        generate_tsv(sample_dataset, out, overwrite=True)

    def test_csv_wrapper_uses_comma(
        self, sample_dataset: List[AnnotatedHLA], tmp_path: Path
    ) -> None:
        """
        :func:`generate_csv` must use a comma delimiter.
        """
        out = tmp_path / "report.csv"
        generate_csv(sample_dataset, out)
        text = out.read_text(encoding="utf-8")
        data_line = next(
            ln for ln in text.splitlines()
            if not ln.startswith("#") and "S1" in ln
        )
        assert "," in data_line
        assert "\t" not in data_line


# ---------------------------------------------------------------------------
# Markdown
# ---------------------------------------------------------------------------


class TestGenerateMarkdown:
    """
    Behaviour of :func:`generate_markdown_report`.
    """

    def test_sections_present(
        self, sample_dataset: List[AnnotatedHLA], tmp_path: Path
    ) -> None:
        """
        All expected section headings must be present.
        """
        out = tmp_path / "report.md"
        generate_markdown_report(sample_dataset, out)
        text = out.read_text(encoding="utf-8")
        # The title omits the word "Clinical" (which would imply
        # Clinical standards the tool does not meet). The three section
        # Headings use association-/evidence-strength framing.
        assert "# HLAnte Research Annotation Report" in text
        assert "## Sample: S1" in text
        assert "## Sample: S2" in text
        assert "### HLA Genotype" in text
        assert "### Reported Disease Associations" in text
        assert "### Reported Pharmacogenomic Associations" in text
        assert "### Interpretation Note (research use only)" in text

    def test_genotype_table(
        self, sample_dataset: List[AnnotatedHLA], tmp_path: Path
    ) -> None:
        """
        The genotype table must list the allele values.
        """
        out = tmp_path / "report.md"
        generate_markdown_report(sample_dataset, out)
        text = out.read_text(encoding="utf-8")
        assert "| HLA-B | B*57:01 | B*07:02 | two-field | t1k |" in text
        assert "| HLA-DRB1 | DRB1*15:01 | NA | two-field | hlahd |" in text

    def test_risk_summary_includes_gwas(
        self, sample_dataset: List[AnnotatedHLA], tmp_path: Path
    ) -> None:
        """
        GWAS hits must be listed with a source and PMID.
        """
        out = tmp_path / "report.md"
        generate_markdown_report(sample_dataset, out)
        text = out.read_text(encoding="utf-8")
        assert "Abacavir hypersensitivity" in text
        assert "OR=4.20" in text
        assert "GCST000001" in text
        assert "PMID: 18322448" in text

    def test_drug_warning_with_emoji(
        self, sample_dataset: List[AnnotatedHLA], tmp_path: Path
    ) -> None:
        """
        Strong (1A/1B) PharmGKB evidence must produce a ⚠️ marker with
        a CPIC link.
        """
        out = tmp_path / "report.md"
        generate_markdown_report(sample_dataset, out)
        text = out.read_text(encoding="utf-8")
        # "⚠️ Critical" was reworded to "⚠️ Strong evidence" to
        # Avoid implying a clinical severity judgement.
        assert "⚠️ Strong evidence" in text
        assert "abacavir HSR risk (1A evidence)" in text
        assert "cpicpgx.org" in text
        assert "pharmgkb.org" in text

    def test_clinical_note_autogenerated(
        self, sample_dataset: List[AnnotatedHLA], tmp_path: Path
    ) -> None:
        """
        The auto-generated interpretation note must follow the expected
        shape. The closing sentence reads "Clinician
        review is recommended" to a research-use reminder.
        """
        out = tmp_path / "report.md"
        generate_markdown_report(sample_dataset, out)
        text = out.read_text(encoding="utf-8")
        assert "HLA-B*57:01" in text
        assert "abacavir" in text
        assert "CPIC" in text
        assert "research-use only" in text
        assert "qualified clinician" in text

    def test_disclaimer_appended(
        self, sample_dataset: List[AnnotatedHLA], tmp_path: Path
    ) -> None:
        """
        The report must end with a blockquoted disclaimer.
        """
        out = tmp_path / "report.md"
        generate_markdown_report(sample_dataset, out)
        text = out.read_text(encoding="utf-8")
        assert f"> **Disclaimer**: {DISCLAIMER}" in text

    def test_empty_dataset(self, tmp_path: Path) -> None:
        """
        Even with no inputs the header and disclaimer must be written.
        """
        out = tmp_path / "empty.md"
        generate_markdown_report([], out)
        text = out.read_text(encoding="utf-8")
        assert "# HLAnte Research Annotation Report" in text
        assert DISCLAIMER in text
        assert "## Sample:" not in text

    def test_overwrite_guard(
        self, sample_dataset: List[AnnotatedHLA], tmp_path: Path
    ) -> None:
        out = tmp_path / "report.md"
        generate_markdown_report(sample_dataset, out)
        with pytest.raises(OutputFileExistsError):
            generate_markdown_report(sample_dataset, out)
        generate_markdown_report(sample_dataset, out, overwrite=True)


# ---------------------------------------------------------------------------
# JSON
# ---------------------------------------------------------------------------


class TestGenerateJSON:
    """
    Behaviour of :func:`generate_json`.
    """

    def test_structure_and_metadata(
        self, sample_dataset: List[AnnotatedHLA], tmp_path: Path
    ) -> None:
        """
        The generated JSON must be valid and contain metadata +
        ``samples[*]/loci[*]/allele*`` nesting.
        """
        out = tmp_path / "report.json"
        ctx = ReportContext(db_versions={"imgt": "3.55.0"})
        generate_json(sample_dataset, out, context=ctx)
        payload = json.loads(out.read_text(encoding="utf-8"))

        assert payload["metadata"]["hlante_version"] == __version__
        assert payload["metadata"]["db_versions"] == {"imgt": "3.55.0"}
        assert payload["metadata"]["disclaimer"] == DISCLAIMER

        samples = {s["sample_id"]: s for s in payload["samples"]}
        assert set(samples) == {"S1", "S2"}

        s1_loci = samples["S1"]["loci"]
        assert len(s1_loci) == 1
        locus = s1_loci[0]
        assert locus["locus"] == "HLA-B"
        assert locus["tool"] == "t1k"
        assert locus["allele1"]["normalized_allele"]["allele_name"] == "B*57:01"
        assert locus["allele2"]["normalized_allele"]["allele_name"] == "B*07:02"
        assert locus["allele1"]["clinical_significance"] == SIGNIFICANCE_PATHOGENIC

        s2_locus = samples["S2"]["loci"][0]
        assert s2_locus["allele2"] is None
        assert s2_locus["allele1"]["gwas_hits"][0]["trait"] == "Multiple sclerosis"

    def test_overwrite_guard(
        self, sample_dataset: List[AnnotatedHLA], tmp_path: Path
    ) -> None:
        out = tmp_path / "r.json"
        generate_json(sample_dataset, out)
        with pytest.raises(OutputFileExistsError):
            generate_json(sample_dataset, out)
        generate_json(sample_dataset, out, overwrite=True)


# ---------------------------------------------------------------------------
# Generate_all
# ---------------------------------------------------------------------------


class TestGenerateAll:
    """
    :func:`generate_all` must produce all three formats.
    """

    def test_all_three_outputs(
        self, sample_dataset: List[AnnotatedHLA], tmp_path: Path
    ) -> None:
        """
        TSV + Markdown + JSON files must be written with the same prefix.
        """
        out_dir = tmp_path / "out"
        paths = generate_all(sample_dataset, out_dir, prefix="run42")
        assert set(paths) == {"tsv", "markdown", "json"}
        for key, expected_suffix in [
            ("tsv", ".tsv"),
            ("markdown", ".md"),
            ("json", ".json"),
        ]:
            assert paths[key].name == f"run42{expected_suffix}"
            assert paths[key].is_file()
        # JSON must be valid.
        data = json.loads(paths["json"].read_text(encoding="utf-8"))
        assert "metadata" in data
        # TSV must contain the disclaimer.
        assert DISCLAIMER in paths["tsv"].read_text(encoding="utf-8")
        # Markdown must contain a Sample heading.
        assert "## Sample: S1" in paths["markdown"].read_text(encoding="utf-8")

    def test_overwrite_propagates(
        self, sample_dataset: List[AnnotatedHLA], tmp_path: Path
    ) -> None:
        """
        ``overwrite=True`` must apply to all three formats.
        """
        out_dir = tmp_path / "out"
        generate_all(sample_dataset, out_dir, prefix="run")
        # Second call with overwrite=False must raise.
        with pytest.raises(OutputFileExistsError):
            generate_all(sample_dataset, out_dir, prefix="run")
        # Overwrite=True must succeed.
        generate_all(sample_dataset, out_dir, prefix="run", overwrite=True)


# ---------------------------------------------------------------------------
# Scale / progress bar
# ---------------------------------------------------------------------------


class TestLargeCohort:
    """
    Verify the tqdm wrapper does not regress with 100+ samples.
    """

    def test_large_cohort_produces_all_rows(self, tmp_path: Path) -> None:
        """
        A 150-sample input must emit 150 data rows in the TSV.
        """
        many: List[AnnotatedHLA] = []
        for i in range(150):
            sample_id = f"S{i:03d}"
            many.append(
                _annotated(
                    _normalized(
                        "A*01:01:01:01",
                        sample_id=sample_id,
                        allele_index=0,
                        locus="HLA-A",
                        imgt_accession="HLA00001",
                        protein_group="A*01:01:01G",
                    ),
                    significance=SIGNIFICANCE_BENIGN,
                )
            )
        out = tmp_path / "big.tsv"
        generate_tsv(many, out)
        lines = [
            ln for ln in out.read_text(encoding="utf-8").splitlines()
            if not ln.startswith("#")
        ]
        # Header + 150 rows
        assert len(lines) == 151


# ---------------------------------------------------------------------------
# Input_quality_tier helper
# ---------------------------------------------------------------------------


class TestConfidenceTier:
    """
    ``_input_quality_tier`` thresholds: HIGH ≥ 0.85, MODERATE ≥ 0.70,
    else LOW; ``None`` → ``NA``.
    """

    def test_high(self) -> None:
        from hlante.annotator import _input_quality_tier
        assert _input_quality_tier(0.90) == "detailed"
        assert _input_quality_tier(0.85) == "detailed"

    def test_moderate(self) -> None:
        from hlante.annotator import _input_quality_tier
        assert _input_quality_tier(0.75) == "partial"
        assert _input_quality_tier(0.70) == "partial"

    def test_low(self) -> None:
        from hlante.annotator import _input_quality_tier
        assert _input_quality_tier(0.64) == "limited"
        assert _input_quality_tier(0.0) == "limited"

    def test_none(self) -> None:
        from hlante.annotator import _input_quality_tier
        assert _input_quality_tier(None) == "NA"


# ---------------------------------------------------------------------------
# Input_quality_rationale uses ";;" between alleles
# ---------------------------------------------------------------------------


class TestConfidenceRationaleSeparator:
    """
    The rationale column must use ``";;"`` between allele1 and allele2
    so intra-rationale pipes remain unambiguous.
    """

    def test_both_alleles_present(self, tmp_path: Path) -> None:
        norm1 = _normalized(
            "B*57:01", sample_id="S1", allele_index=0,
            imgt_accession="HLA00001", protein_group=None,
        )
        norm2 = _normalized(
            "B*07:02", sample_id="S1", allele_index=1,
            imgt_accession="HLA00009", protein_group=None,
        )
        a1 = _annotated(norm1)
        a1.input_quality_rationale = "freq_unknown|ambiguous"
        a2 = _annotated(norm2)
        a2.input_quality_rationale = "freq_unknown|ambiguous"

        out = tmp_path / "rationale.tsv"
        generate_tsv([a1, a2], out)
        text = out.read_text(encoding="utf-8")
        data_line = next(
            ln for ln in text.splitlines()
            if not ln.startswith("#") and "S1" in ln
        )
        cells = data_line.split("\t")
        from hlante.reporter import TSV_COLUMNS
        idx = TSV_COLUMNS.index("input_quality_rationale")
        cell = cells[idx]
        parts = cell.split(";;")
        assert len(parts) == 2, (
            f"Expected ;; separator between alleles, got cell {cell!r}"
        )
        assert "freq_unknown" in parts[0]
        assert "freq_unknown" in parts[1]


class TestGLString:
    """
    GL String output (Milius et al. 2013; GL String 1.1, Mack et al. 2023).

    HLAnte holds at most two allele designations per locus and models neither
    chromosomal phase nor alternative genotypes, so it must emit only the
    genotype delimiter (+) and the locus delimiter (^). Emitting /, ~, | or ?
    would assert information the tool does not have.
    """

    def _read(self, path: Path):
        lines = path.read_text(encoding="utf-8").splitlines()
        rows = list(csv.reader([ln for ln in lines if not ln.startswith("#")], delimiter="\t"))
        return rows[0], rows[1:]

    def test_locus_genotype_is_fully_qualified_and_plus_joined(
        self, sample_dataset: List[AnnotatedHLA], tmp_path: Path
    ) -> None:
        out = tmp_path / "report.tsv"
        generate_tsv(sample_dataset, out)
        header, data = self._read(out)
        idx = header.index("gl_string")
        values = [r[idx] for r in data if r[idx] != "NA"]
        assert values, "no GL String produced"
        for value in values:
            for allele in value.split("+"):
                assert allele.startswith("HLA-") and "*" in allele, value

    def test_never_emits_unsupported_operators(
        self, sample_dataset: List[AnnotatedHLA], tmp_path: Path
    ) -> None:
        out = tmp_path / "report.tsv"
        generate_tsv(sample_dataset, out)
        header, data = self._read(out)
        idx = header.index("gl_string")
        for row in data:
            for operator in ("/", "~", "|", "?"):
                assert operator not in row[idx], (operator, row[idx])

    def test_single_reported_allele_is_not_duplicated(
        self, sample_dataset: List[AnnotatedHLA], tmp_path: Path
    ) -> None:
        """
        The S2/HLA-DRB1 record in the fixture reports one allele; the missing
        gene copy must not be invented.
        """
        out = tmp_path / "report.tsv"
        generate_tsv(sample_dataset, out)
        header, data = self._read(out)
        idx = header.index("gl_string")
        singles = [r[idx] for r in data if r[header.index("allele2")] == "NA"]
        assert singles, "fixture must contain a single-allele locus"
        for value in singles:
            assert "+" not in value, value

    def test_sample_level_string_joins_loci(
        self, sample_dataset: List[AnnotatedHLA], tmp_path: Path
    ) -> None:
        out = tmp_path / "report.json"
        generate_json(sample_dataset, out)
        payload = json.loads(out.read_text(encoding="utf-8"))
        for sample in payload["samples"]:
            assert sample["gl_string"]
            for locus in sample["loci"]:
                assert locus["gl_string"]
