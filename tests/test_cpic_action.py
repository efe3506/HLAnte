"""
tests.test_cpic_action
======================

Tests for allele-aware CPIC action verbs.

The CPIC therapeutic recommendation for a drug depends on *which* HLA
allele is carried, not on the drug alone. The clearest example is
carbamazepine: ``HLA-B*15:02`` is a strong contraindication in
carbamazepine-naïve patients, whereas ``HLA-A*31:01`` is a weaker
association with a milder recommendation. Keying the action verb on the
drug alone (the previous behaviour) collapsed both to an identical
"Contraindicated", over-restricting A*31:01.
"""

from __future__ import annotations

import pytest

from hlante.reporter import (
    ALLELE_DRUG_ACTION_MAP,
    CPIC_ACTION_MAP,
    _cpic_action,
    _two_field_key,
)


class TestTwoFieldKey:
    @pytest.mark.parametrize(
        "allele,expected",
        [
            ("B*15:02", "B*15:02"),
            ("HLA-B*15:02", "B*15:02"),
            ("HLA-B*15:02:01", "B*15:02"),
            ("B*15:02:01:02N", "B*15:02"),
            ("A*31:01:01G", "A*31:01"),
            ("B*57", None),  # only one field
            ("not-an-allele", None),
            (None, None),
            ("", None),
        ],
    )
    def test_two_field_key(self, allele, expected) -> None:
        assert _two_field_key(allele) == expected


class TestAlleleAwareCpicAction:
    def test_b1502_carbamazepine_avoid(self) -> None:
        action = _cpic_action("B*15:02", "carbamazepine")
        assert action is not None
        assert "Avoid" in action or "ontraindicated" in action

    def test_a3101_carbamazepine_is_milder(self) -> None:
        action = _cpic_action("A*31:01", "carbamazepine")
        assert action is not None
        # A*31:01 must NOT be a flat contraindication.
        assert "Contraindicated" not in action
        assert "alternative" in action.lower() or "monitoring" in action.lower()

    def test_carbamazepine_recommendation_is_allele_specific(self) -> None:
        """The core B2 fix: the two carbamazepine alleles differ."""
        b1502 = _cpic_action("B*15:02", "carbamazepine")
        a3101 = _cpic_action("A*31:01", "carbamazepine")
        assert b1502 != a3101

    def test_b5701_abacavir_contraindicated(self) -> None:
        assert _cpic_action("B*57:01", "abacavir") == "Contraindicated (do not use)"

    def test_b5801_allopurinol_contraindicated(self) -> None:
        action = _cpic_action("B*58:01", "allopurinol")
        assert action is not None and "ontraindicated" in action

    def test_two_field_reduction_on_high_resolution_input(self) -> None:
        # A 4-field, HLA-prefixed call still resolves to the 2-field rule.
        assert _cpic_action("HLA-B*15:02:01", "carbamazepine") == _cpic_action(
            "B*15:02", "carbamazepine"
        )

    def test_unlisted_pair_falls_back_to_drug_level(self) -> None:
        # C*07:01 has no specific lamotrigine rule → drug-level default.
        assert _cpic_action("C*07:01", "lamotrigine") == CPIC_ACTION_MAP["lamotrigine"]

    def test_missing_allele_uses_drug_fallback(self) -> None:
        assert _cpic_action(None, "abacavir") == CPIC_ACTION_MAP["abacavir"]

    def test_missing_drug_returns_none(self) -> None:
        assert _cpic_action("B*57:01", None) is None

    def test_unknown_drug_returns_none(self) -> None:
        assert _cpic_action("B*57:01", "aspirin") is None

    def test_every_specific_allele_drug_pair_has_a_nonempty_verb(self) -> None:
        for (allele, drug), verb in ALLELE_DRUG_ACTION_MAP.items():
            assert verb and _cpic_action(allele, drug) == verb
