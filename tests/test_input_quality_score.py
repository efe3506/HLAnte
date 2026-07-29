"""
tests.test_input_quality_score
===========================

Unit tests for AFND integration and confidence score computation.

Coverage
--------
- ``_compute_input_quality_score`` deterministic across all combinations.
- ``AFNDClient`` population-group filtering (universal, no country
  aliasing).
- ``get_frequency_with_fallback`` sets ``is_estimated`` when needed.
- ``annotate_genotype`` end-to-end produces the correct confidence.
- The four new TSV columns are present and written correctly.
- ``POPULATION_GROUPS`` contains no country-level entries.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import pytest

from hlante.annotator import (
    AnnotatedHLA,
    AnnotatorClients,
    AnnotatorConfig,
    _compute_input_quality_score,
    annotate_genotype,
    build_clients,
)
from hlante.types import InputSource
from unittest.mock import patch

from hlante.db.afnd import (
    AFNDClient,
    AFNDDatabaseError,
    AllelFrequency,
    POPULATION_GROUPS,
)
from hlante.normalizer import NormalizedAllele
from hlante.reporter import TSV_COLUMNS, generate_tsv


FIXTURES_DIR: Path = Path(__file__).parent / "fixtures"
AFND_FIXTURE_DIR: Path = FIXTURES_DIR / "afnd"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_normalized(
    allele_name: str = "B*57:01",
    *,
    resolution_level: int = 4,
    is_ambiguous: bool = False,
    is_novel: bool = False,
    imgt_accession: Optional[str] = "HLA00009",
    gene: str = "HLA-B",
    hla_class: str = "I",
) -> NormalizedAllele:
    return NormalizedAllele(
        allele_name=allele_name,
        imgt_accession=imgt_accession,
        protein_group=None,
        hla_class=hla_class,
        gene=gene,
        resolution_level=resolution_level,
        is_ambiguous=is_ambiguous,
        is_novel=is_novel,
        sample_id="S1",
        source_tool="t1k",
        source_locus=gene,
        source_resolution=f"{resolution_level}-field",
        allele_index=0,
    )


def _make_freq(
    frequency: float,
    *,
    allele: str = "B*57:01",
    population_group: str = "EUR",
    sample_size: int = 500,
    resolution: int = 2,
    estimated: bool = False,
) -> AllelFrequency:
    return AllelFrequency(
        allele=allele,
        frequency=frequency,
        population_group=population_group,
        sample_size=sample_size,
        source_resolution=resolution,
        is_estimated=estimated,
    )


# ---------------------------------------------------------------------------
# _compute_input_quality_score
# ---------------------------------------------------------------------------


class TestComputeConfidenceScore:
    """
    Deterministic score computation — expected multipliers across all
    combinations.
    """

    def test_standard_allele_full_confidence(self) -> None:
        """
        four-field + known in IMGT + common allele → score close to 1.0.
        """
        norm = _make_normalized(resolution_level=4)
        freq = _make_freq(0.15)
        score, rationale = _compute_input_quality_score(norm, freq)
        assert score == 1.0
        assert rationale == "standard"

    def test_novel_allele_severe_penalty(self) -> None:
        """
        Novel → ×0.3 → score < 0.5 (spec assertion).
        """
        norm = _make_normalized(
            resolution_level=4, is_novel=True, imgt_accession=None
        )
        score, rationale = _compute_input_quality_score(norm, _make_freq(0.15))
        assert score < 0.5
        assert score == pytest.approx(0.3)
        assert "novel_allele" in rationale

    def test_rare_allele_penalty(self) -> None:
        """
        freq < 0.001 → ×0.5 → score < 0.6.
        """
        norm = _make_normalized(resolution_level=4)
        score, rationale = _compute_input_quality_score(norm, _make_freq(0.0005))
        assert score < 0.6
        assert score == pytest.approx(0.5)
        assert "rare_allele" in rationale

    def test_uncommon_allele_penalty(self) -> None:
        """
        0.001 ≤ freq < 0.01 → ×0.8.
        """
        norm = _make_normalized(resolution_level=4)
        score, _ = _compute_input_quality_score(norm, _make_freq(0.005))
        assert score == pytest.approx(0.8)

    def test_unknown_frequency_small_penalty(self) -> None:
        """
        Unknown frequency → ×0.85.
        """
        norm = _make_normalized(resolution_level=4)
        score, rationale = _compute_input_quality_score(norm, None)
        assert score == pytest.approx(0.85)
        assert rationale == "freq_unknown"

    def test_common_allele_4field_imgt_high_confidence(self) -> None:
        """
        freq > 0.1 + 4-field + known in IMGT → score > 0.8
        (spec assertion).
        """
        norm = _make_normalized(resolution_level=2)
        score, rationale = _compute_input_quality_score(norm, _make_freq(0.15))
        assert score > 0.8
        assert score == pytest.approx(0.9)
        assert "medium_resolution(two-field)" in rationale

    def test_low_resolution_penalty(self) -> None:
        """
        one-field → ×0.7.
        """
        norm = _make_normalized(resolution_level=1)
        score, rationale = _compute_input_quality_score(norm, _make_freq(0.15))
        assert score == pytest.approx(0.7)
        assert "low_resolution" in rationale

    def test_ambiguous_penalty(self) -> None:
        """
        is_ambiguous=True → ×0.75.
        """
        norm = _make_normalized(resolution_level=4, is_ambiguous=True)
        score, rationale = _compute_input_quality_score(norm, _make_freq(0.15))
        assert score == pytest.approx(0.75)
        assert "ambiguous" in rationale

    def test_combined_penalties_multiply(self) -> None:
        """
        Novel + rare + one-field + ambiguous penalties must compose:
        1.0 × 0.3 × 0.5 × 0.7 × 0.75 = 0.07875 → 0.0788.
        """
        norm = _make_normalized(
            resolution_level=1,
            is_ambiguous=True,
            is_novel=True,
            imgt_accession=None,
        )
        score, rationale = _compute_input_quality_score(
            norm, _make_freq(0.0005)
        )
        assert score == pytest.approx(0.0788, abs=1e-4)
        for code in (
            "novel_allele", "rare_allele", "low_resolution", "ambiguous"
        ):
            assert code in rationale

    def test_determinism(self) -> None:
        """
        The same inputs must produce the same score / rationale on
        every call.
        """
        norm = _make_normalized(resolution_level=2, is_ambiguous=True)
        freq = _make_freq(0.005)
        for _ in range(5):
            s, r = _compute_input_quality_score(norm, freq)
            assert s == pytest.approx(0.8 * 0.9 * 0.75)  # 0.54
            assert r == (
                "uncommon_allele(freq=0.0050)|medium_resolution(two-field)"
                "|ambiguous"
            )


# ---------------------------------------------------------------------------
# Source-aware confidence scoring
# ---------------------------------------------------------------------------


class TestSourceAwareConfidence:
    """Ambiguity penalty is source-dependent; resolution penalty is always applied."""

    def test_typing_tool_2field_keeps_ambiguity_penalty(self) -> None:
        """one-field with typing_tool source: ambiguity ×0.75 applies."""
        norm = _make_normalized(resolution_level=1, is_ambiguous=True)
        freq = _make_freq(0.05)
        score, rationale = _compute_input_quality_score(
            norm, freq, InputSource.TYPING_TOOL
        )
        assert "ambiguous" in rationale
        # 1.0 × 0.70 (one-field) × 0.75 (ambiguous) = 0.525
        assert score == pytest.approx(0.525, abs=1e-4)
        assert score < 0.60

    def test_validated_2field_suppresses_ambiguity_penalty(self) -> None:
        """one-field with validated source: ambiguity penalty NOT applied."""
        norm = _make_normalized(resolution_level=1, is_ambiguous=True)
        freq = _make_freq(0.05)
        score, rationale = _compute_input_quality_score(
            norm, freq, InputSource.VALIDATED
        )
        assert "ambiguity_suppressed(validated_source)" in rationale
        # "ambiguous" must not appear outside the suppressed message
        bare = rationale.replace("ambiguity_suppressed(validated_source)", "")
        assert "ambiguous" not in bare
        # 1.0 × 0.70 (one-field only) = 0.70
        assert score == pytest.approx(0.70, abs=1e-4)
        assert 0.65 <= score <= 0.75

    def test_validated_4field_no_ambiguity_high_tier(self) -> None:
        """two-field validated, not ambiguous: only medium-resolution penalty applies."""
        norm = _make_normalized(resolution_level=2, is_ambiguous=False)
        freq = _make_freq(0.05)
        score, _ = _compute_input_quality_score(norm, freq, InputSource.VALIDATED)
        # 1.0 × 0.90 (two-field) = 0.90
        assert score == pytest.approx(0.90, abs=1e-4)
        assert score >= 0.85  # HIGH tier

    def test_validated_does_not_suppress_resolution(self) -> None:
        """Resolution penalty must still apply for validated one-field."""
        freq = _make_freq(0.05)
        score_2f, _ = _compute_input_quality_score(
            _make_normalized(resolution_level=1, is_ambiguous=True),
            freq,
            InputSource.VALIDATED,
        )
        score_4f, _ = _compute_input_quality_score(
            _make_normalized(resolution_level=2, is_ambiguous=False),
            freq,
            InputSource.VALIDATED,
        )
        assert score_2f < score_4f  # 0.70 < 0.90

    def test_validated_does_not_suppress_novel(self) -> None:
        """Novel penalty must still apply regardless of source."""
        norm = _make_normalized(
            resolution_level=2, is_ambiguous=False,
            is_novel=True, imgt_accession=None,
        )
        score, rationale = _compute_input_quality_score(
            norm, None, InputSource.VALIDATED
        )
        assert "novel_allele" in rationale
        assert score <= 0.30

    def test_simulated_applies_ambiguity_penalty(self) -> None:
        """SIMULATED is treated identically to TYPING_TOOL for penalties."""
        norm = _make_normalized(resolution_level=1, is_ambiguous=True)
        freq = _make_freq(0.05)
        score_sim, rationale_sim = _compute_input_quality_score(
            norm, freq, InputSource.SIMULATED
        )
        score_tt, rationale_tt = _compute_input_quality_score(
            norm, freq, InputSource.TYPING_TOOL
        )
        assert score_sim == score_tt
        assert rationale_sim == rationale_tt

    def test_unknown_applies_ambiguity_penalty(self) -> None:
        """UNKNOWN is treated identically to TYPING_TOOL for penalties."""
        norm = _make_normalized(resolution_level=1, is_ambiguous=True)
        freq = _make_freq(0.05)
        score_unk, rationale_unk = _compute_input_quality_score(
            norm, freq, InputSource.UNKNOWN
        )
        score_tt, rationale_tt = _compute_input_quality_score(
            norm, freq, InputSource.TYPING_TOOL
        )
        assert score_unk == score_tt
        assert rationale_unk == rationale_tt

    def test_default_is_typing_tool(self) -> None:
        """Default input_source must be TYPING_TOOL (backward-compatible)."""
        norm = _make_normalized(resolution_level=1, is_ambiguous=True)
        freq = _make_freq(0.05)
        score_default, rationale_default = _compute_input_quality_score(norm, freq)
        score_tt, rationale_tt = _compute_input_quality_score(
            norm, freq, InputSource.TYPING_TOOL
        )
        assert score_default == score_tt
        assert rationale_default == rationale_tt


# ---------------------------------------------------------------------------
# POPULATION_GROUPS taxonomy
# ---------------------------------------------------------------------------


class TestPopulationGroupTaxonomy:
    """
    The :data:`POPULATION_GROUPS` mapping must be purely geographic.
    """

    def test_no_country_name_in_population_groups(self) -> None:
        """
        :data:`POPULATION_GROUPS` must not contain any country-level
        entries.
        """
        country_names = [
            "turkey", "turkish",
            "german", "germany",
            "french", "france",
            "british", "britain", "uk",
            "chinese", "china",
            "japanese", "japan",
            "indian", "india",
            "iranian", "iran",
            "saudi", "saudi arabia",
        ]
        all_keywords = [
            kw.lower()
            for keywords in POPULATION_GROUPS.values()
            for kw in keywords
        ]
        for country in country_names:
            assert country not in all_keywords, (
                f"Country name '{country}' must not appear in POPULATION_GROUPS"
            )

    def test_expected_group_codes_present(self) -> None:
        """
        Every expected universal group code must be a key.
        """
        expected = {"EUR", "AFR", "EAS", "SAS", "MID", "AMR", "OCE", "global"}
        assert expected == set(POPULATION_GROUPS.keys())

    def test_global_has_no_keywords(self) -> None:
        """
        ``"global"`` must be empty; it is interpreted as an all-pass
        filter.
        """
        assert POPULATION_GROUPS["global"] == []


# ---------------------------------------------------------------------------
# AFNDClient — population filter + fallback
# ---------------------------------------------------------------------------


class TestAFNDClient:
    """
    Local-TSV AFND client.
    """

    def test_loads_fixture(self) -> None:
        client = AFNDClient(local_dir=AFND_FIXTURE_DIR)
        client.load()
        # The fixture contains at least 20 rows.
        assert len(client._rows) >= 20

    def test_population_group_eur_aggregates_european(self) -> None:
        """
        ``population_group="EUR"`` must aggregate every row whose
        Population Group substring matches one of the EUR keywords.
        """
        client = AFNDClient(
            local_dir=AFND_FIXTURE_DIR,
            population_group="EUR",
        )
        freq = client.get_frequency("B*57:01")
        assert freq is not None
        # Fixture: Germany(400, 0.055) + UK(1000, 0.060) → weighted mean
        # (0.055*400 + 0.060*1000) / 1400 ≈ 0.0586
        assert freq.frequency == pytest.approx(0.0586, abs=1e-3)
        assert freq.sample_size == 1400
        assert freq.populations_aggregated == 2

    def test_population_group_mid(self) -> None:
        """
        ``MID`` group selects Middle Eastern / West Asian / Arab rows.
        """
        client = AFNDClient(
            local_dir=AFND_FIXTURE_DIR,
            population_group="MID",
        )
        freq = client.get_frequency("A*02:01")
        assert freq is not None
        assert freq.population_group == "MID"
        # Arab + Middle Eastern populations aggregated.
        assert freq.populations_aggregated >= 2

    def test_population_group_eas(self) -> None:
        """
        ``EAS`` selects East/Southeast Asian rows (Japanese mapped).
        """
        client = AFNDClient(
            local_dir=AFND_FIXTURE_DIR,
            population_group="EAS",
        )
        freq = client.get_frequency("A*02:01")
        assert freq is not None
        assert freq.population_group == "EAS"

    def test_population_group_global_aggregates_all(self) -> None:
        """
        ``global`` aggregates every population (weighted by size).
        """
        client = AFNDClient(
            local_dir=AFND_FIXTURE_DIR,
            population_group="global",
        )
        freq = client.get_frequency("A*02:01")
        assert freq is not None
        assert freq.population_group == "global"

    def test_min_sample_size_filter(self) -> None:
        """
        ``Small European Study`` (n=25) must be excluded when min=50.
        The A*02:01 EUR aggregate must therefore exclude that row.
        """
        client = AFNDClient(
            local_dir=AFND_FIXTURE_DIR,
            population_group="EUR",
            min_sample_size=50,
        )
        freq = client.get_frequency("A*02:01")
        assert freq is not None
        # Only Germany(500) + UK(1000) remain (Small study is dropped).
        assert freq.sample_size == 1500

    def test_no_match_returns_none(self) -> None:
        client = AFNDClient(
            local_dir=AFND_FIXTURE_DIR,
            population_group="EUR",
        )
        assert client.get_frequency("Z*99:99") is None

    def test_get_frequency_with_fallback_uses_2field(self) -> None:
        """
        Fixture holds only two-field alleles. A three-field query must fall
        back and set ``is_estimated=True``.
        """
        client = AFNDClient(
            local_dir=AFND_FIXTURE_DIR,
            population_group="EUR",
        )
        freq = client.get_frequency_with_fallback("B*57:01:01")
        assert freq is not None
        assert freq.is_estimated is True
        assert freq.allele == "B*57:01"
        assert freq.source_resolution == 2  # 2 colon-groups

    def test_fallback_exact_match_not_estimated(self) -> None:
        """
        An exact match must leave ``is_estimated=False``.
        """
        client = AFNDClient(
            local_dir=AFND_FIXTURE_DIR,
            population_group="EUR",
        )
        freq = client.get_frequency_with_fallback("B*57:01")
        assert freq is not None
        assert freq.is_estimated is False

    def test_hla_prefix_stripped(self) -> None:
        client = AFNDClient(local_dir=AFND_FIXTURE_DIR, population_group="EUR")
        a = client.get_frequency("HLA-B*57:01")
        b = client.get_frequency("B*57:01")
        assert a is not None and b is not None
        assert a.frequency == b.frequency

    def test_rare_allele_in_fixture(self) -> None:
        """
        ``B*99:99`` (freq=0.0005) — expected rare allele for
        confidence tests.
        """
        client = AFNDClient(local_dir=AFND_FIXTURE_DIR, population_group="EUR")
        freq = client.get_frequency("B*99:99")
        assert freq is not None
        assert freq.frequency < 0.001

    def test_missing_directory_falls_back_to_builtin(
        self, tmp_path: Path
    ) -> None:
        """When local_dir doesn't exist, load() uses the built-in TSV fallback."""
        client = AFNDClient(local_dir=tmp_path / "none")
        client.load()  # Must not raise
        assert client._loaded  # type: ignore[attr-defined]

    def test_missing_directory_and_builtin_raises(
        self, tmp_path: Path
    ) -> None:
        """Error is raised only when both local TSV and built-in fallback are absent."""
        client = AFNDClient(local_dir=tmp_path / "none")
        non_existent = tmp_path / "no_builtin.tsv"
        import hlante.db.afnd as _afnd_mod
        with patch.object(_afnd_mod, "BUILTIN_AFND_TSV", non_existent):
            with pytest.raises(AFNDDatabaseError, match="not found"):
                client.load()

    def test_update_raises_on_network_failure(self, tmp_path: Path) -> None:
        """AFNDDatabaseError is raised when the download fails."""
        import urllib.error

        client = AFNDClient(local_dir=tmp_path)
        with patch("hlante.db.afnd.urllib.request.urlopen", side_effect=urllib.error.URLError("simulated")):
            with pytest.raises(AFNDDatabaseError, match="download failed"):
                client.update()


# ---------------------------------------------------------------------------
# Annotate_genotype end to end
# ---------------------------------------------------------------------------


class TestEndToEnd:
    """
    ``annotate_genotype`` must invoke the AFND client and propagate the
    confidence score correctly.
    """

    def test_common_allele_in_european_population(self) -> None:
        """
        A*02:01 is common in European populations → high confidence.
        """
        norm = _make_normalized(
            "A*02:01", resolution_level=2, gene="HLA-A",
            imgt_accession="HLA00004",
        )
        config = AnnotatorConfig(
            offline=True,
            afnd_local_dir=AFND_FIXTURE_DIR,
            population_group="EUR",
            enable_gwas=False,
            enable_pharmgkb=False,

        )
        clients = AnnotatorClients(afnd=build_clients(config).afnd)
        result = annotate_genotype([norm], config, clients=clients)[0]
        # The query label is preserved; it is not hardcoded.
        assert result.frequency_population is not None
        assert result.allele_frequency is not None
        assert result.allele_frequency > 0.1
        assert result.input_quality_score > 0.8

    def test_rare_allele_lower_confidence(self) -> None:
        norm = _make_normalized(
            "B*99:99", resolution_level=2, gene="HLA-B",
            imgt_accession=None, is_novel=False,
        )
        config = AnnotatorConfig(
            offline=True,
            afnd_local_dir=AFND_FIXTURE_DIR,
            population_group="EUR",
            enable_gwas=False,
            enable_pharmgkb=False,

        )
        result = annotate_genotype(
            [norm], config, clients=build_clients(config)
        )[0]
        assert result.allele_frequency is not None
        assert result.allele_frequency < 0.001
        assert result.input_quality_score < 0.6
        assert "rare_allele" in result.input_quality_rationale

    def test_novel_allele_very_low_confidence(self) -> None:
        norm = _make_normalized(
            "A*99:99", resolution_level=2, gene="HLA-A",
            imgt_accession=None, is_novel=True,
        )
        config = AnnotatorConfig(
            offline=True,
            afnd_local_dir=AFND_FIXTURE_DIR,
            population_group="EUR",
            enable_gwas=False,
            enable_pharmgkb=False,

        )
        result = annotate_genotype(
            [norm], config, clients=build_clients(config)
        )[0]
        # Novel + no freq → 0.3 * 0.85 * 0.9 = 0.2295
        assert result.input_quality_score < 0.5
        assert "novel_allele" in result.input_quality_rationale

    def test_afnd_disabled_leaves_freq_none(self) -> None:
        norm = _make_normalized("A*02:01", resolution_level=2)
        config = AnnotatorConfig(
            offline=True,
            enable_gwas=False,
            enable_pharmgkb=False,

            enable_afnd=False,
        )
        clients = build_clients(config)
        assert clients.afnd is None
        result = annotate_genotype([norm], config, clients=clients)[0]
        assert result.allele_frequency is None
        assert "freq_unknown" in result.input_quality_rationale


# ---------------------------------------------------------------------------
# Reporter columns
# ---------------------------------------------------------------------------


class TestReporterColumns:
    """
    The four new TSV columns must be present with the right values.
    """

    def test_schema_contains_four_new_columns(self) -> None:
        for col in (
            "allele_frequency",
            "allele_freq_population",
            "input_quality_score",
            "input_quality_rationale",
        ):
            assert col in TSV_COLUMNS

    def test_columns_after_clinical_significance(self) -> None:
        idx_sig = TSV_COLUMNS.index("clinical_significance")
        for i, col in enumerate([
            "significance_basis",
            "allele_frequency",
            "allele_freq_population",
            "input_quality_score",
            "input_quality_tier",
            "input_quality_rationale",
        ], start=1):
            assert TSV_COLUMNS[idx_sig + i] == col

    def test_values_written_in_tsv(self, tmp_path: Path) -> None:
        norm = _make_normalized("B*57:01", resolution_level=2)
        annot = AnnotatedHLA(
            normalized_allele=norm,
            gwas_hits=[],
            pharm_annotations=[],
            disease_entries=[],
            disease_risk_summary="",
            drug_response_summary="",
            clinical_significance="Benign",
            allele_frequency=0.0509,
            frequency_population="EUR",
            frequency_sample_size=550,
            input_quality_score=0.9,
            input_quality_rationale="medium_resolution(two-field)",
        )
        out = tmp_path / "r.tsv"
        generate_tsv([annot], out)
        text = out.read_text(encoding="utf-8")
        data_line = next(
            ln for ln in text.splitlines()
            if not ln.startswith("#") and "S1" in ln
        )
        cells = data_line.split("\t")
        idx = {col: i for i, col in enumerate(TSV_COLUMNS)}
        assert cells[idx["allele_frequency"]] == "0.050900|NA"
        assert cells[idx["allele_freq_population"]] == "EUR|NA"
        assert cells[idx["input_quality_score"]] == "0.9000|NA"
        assert "medium_resolution" in cells[idx["input_quality_rationale"]]
