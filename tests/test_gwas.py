"""
tests.test_gwas
===============

Unit tests for :mod:`hlante.db.gwas` helper functions — obsolete
EFO remapping, extreme effect-size classification, and the
fallback-expansion / annotation-scope logic.

These tests are file-free (no DB dump required).
"""

from __future__ import annotations

import pytest

from hlante.db.gwas import (
    GWASHit,
    OBSOLETE_EFO_MAP,
    _classify_effect_size,
    _remap_trait,
)


# ---------------------------------------------------------------------------
# Obsolete EFO remapping
# ---------------------------------------------------------------------------


class TestRemapTrait:
    """
    :func:`_remap_trait` must tag every ``obsolete_*`` trait and pass
    non-obsolete traits through unchanged.
    """

    def test_obsolete_trait_remapped(self) -> None:
        """
        A mapped ``obsolete_*`` trait returns the clean current EFO
        term. No ``[remapped from deprecated EFO]`` bracket suffix is
        appended — provenance lives on the ``was_deprecated`` flag.
        """
        trait, deprecated = _remap_trait("obsolete_myositis")
        assert deprecated is True
        assert "inflammatory myopathy" in trait
        assert "EFO:" in trait
        # The trait cell must not carry a trailing bracket annotation.
        assert "[remapped" not in trait
        assert "deprecated" not in trait.lower()

    def test_non_obsolete_unchanged(self) -> None:
        trait, deprecated = _remap_trait("rheumatoid arthritis")
        assert deprecated is False
        assert trait == "rheumatoid arthritis"

    def test_unknown_obsolete_flagged(self) -> None:
        """
        An unmapped ``obsolete_*`` trait flags ``was_deprecated=True``
        but the trait cell itself is returned unchanged (no bracket
        suffix). The ``obsolete_`` prefix already conveys the status,
        and the flag is the canonical provenance channel.
        """
        trait, deprecated = _remap_trait("obsolete_unknown_condition")
        assert deprecated is True
        assert trait == "obsolete_unknown_condition"
        assert "[deprecated" not in trait

    def test_all_known_obsolete_terms_covered(self) -> None:
        """
        Every ``obsolete_*`` key observed in prior real output is
        listed in :data:`OBSOLETE_EFO_MAP`.
        """
        required = {
            "obsolete_myositis",
            "obsolete_uveitis",
            "obsolete_juvenile idiopathic arthritis",
            "obsolete_sclerosing cholangitis",
            "obsolete_neuromyelitis optica",
            "obsolete_Autoimmune Hepatitis",
            "obsolete_late-onset myasthenia gravis",
        }
        missing = required - set(OBSOLETE_EFO_MAP)
        assert not missing, f"Missing remappings: {missing}"

    def test_empty_string_unchanged(self) -> None:
        trait, deprecated = _remap_trait("")
        assert deprecated is False
        assert trait == ""

    def test_compound_two_obsolete_both_remapped(self) -> None:
        """
        The GWAS Catalog sometimes packs multiple traits into one
        cell separated by ``", "``. Every obsolete sub-term must be
        remapped independently and the flag set. Pre-fix the full
        concatenated string was not in the map so the whole cell
        passed through raw.
        """
        raw = "obsolete_uveitis, obsolete_juvenile idiopathic arthritis"
        trait, deprecated = _remap_trait(raw)
        assert deprecated is True
        assert "uveitis (EFO:0004284)" in trait
        assert "juvenile idiopathic arthritis (EFO:0000685)" in trait
        assert "obsolete_" not in trait

    def test_compound_mixed_obsolete_and_current(self) -> None:
        """
        A compound cell with one obsolete and one non-obsolete part
        must remap the obsolete part and pass the non-obsolete part
        through verbatim. The rejoined cell retains ``", "`` order.
        """
        raw = (
            "obsolete_toxic epidermal necrolysis, "
            "obsolete_Stevens-Johnson syndrome, "
            "response to methazolamide"
        )
        trait, deprecated = _remap_trait(raw)
        assert deprecated is True
        assert "toxic epidermal necrolysis (EFO:0004197)" in trait
        assert "Stevens-Johnson syndrome (EFO:0004190)" in trait
        assert "response to methazolamide" in trait
        assert "obsolete_" not in trait

    def test_compound_non_obsolete_unchanged(self) -> None:
        """
        A compound cell with zero obsolete tokens is untouched and
        ``was_deprecated`` is ``False``. Guards against false
        positives when a non-obsolete trait legitimately contains
        ``", "`` (e.g. GxE interaction traits).
        """
        raw = (
            "High density lipoprotein cholesterol measurement, "
            "interaction with smoking behaviour"
        )
        trait, deprecated = _remap_trait(raw)
        assert deprecated is False
        assert trait == raw

    def test_compound_unmapped_obsolete_preserved(self) -> None:
        """
        A compound cell where an obsolete sub-term is NOT in
        :data:`OBSOLETE_EFO_MAP` keeps that sub-term raw (prefix
        retained) and still sets ``was_deprecated=True``.
        """
        raw = "obsolete_uveitis, obsolete_never_heard_of_this"
        trait, deprecated = _remap_trait(raw)
        assert deprecated is True
        assert "uveitis (EFO:0004284)" in trait
        assert "obsolete_never_heard_of_this" in trait


# ---------------------------------------------------------------------------
# Extreme / quantitative-trait effect-size classification
# ---------------------------------------------------------------------------


class TestClassifyEffectSize:
    """
    :func:`_classify_effect_size` must warn when OR > 10 and add a
    ``quantitative_trait_effect`` tag when the trait text matches
    continuous-trait keywords.
    """

    def test_extreme_or_quantitative_flagged(self) -> None:
        w = _classify_effect_size(205.0, "blood immunoglobulin amount")
        assert "extreme_value" in w
        assert "quantitative_trait_effect" in w

    def test_extreme_or_binary_disease_flagged(self) -> None:
        w = _classify_effect_size(10.5, "ankylosing spondylitis")
        assert "extreme_value" in w
        assert "quantitative_trait_effect" not in w

    def test_normal_or_not_flagged(self) -> None:
        w = _classify_effect_size(1.4, "rheumatoid arthritis")
        assert w == ""

    def test_none_or_not_flagged(self) -> None:
        w = _classify_effect_size(None, "whatever")
        assert w == ""

    @pytest.mark.parametrize(
        "keyword_trait",
        [
            "protein level",
            "blood cell count",
            "cytokine concentration",
            "antibody seropositivity",
            "vaccine response",
        ],
    )
    def test_quantitative_keywords_detected(self, keyword_trait: str) -> None:
        w = _classify_effect_size(50.0, keyword_trait)
        assert "quantitative_trait_effect" in w


# ---------------------------------------------------------------------------
# GWASHit dataclass — new fields default correctly
# ---------------------------------------------------------------------------


class TestGWASHitDefaults:
    """
    New fields added in the v0.1.1 fix sprint must carry safe defaults
    so every existing constructor call keeps working.
    """

    def test_new_fields_default_values(self) -> None:
        hit = GWASHit(
            trait="example",
            p_value=1e-9,
            odds_ratio=2.0,
            pmid="12345",
            study_accession="GCST000001",
            allele="B*57:01",
        )
        assert hit.trait_was_deprecated is False
        assert hit.effect_size_warning == ""
        assert hit.annotation_scope == "allele"
        assert hit.index_siblings == 1
