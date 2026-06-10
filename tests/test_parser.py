"""
tests.test_parser
=================

Unit tests for :mod:`hlante.parser`.

Coverage
--------
- Allele regex validation.
- Resolution detection.
- Dispatcher behaviour (``parse_hla_output``).
- Happy-path for each tool-specific parser.
- Missing-file and malformed-input exceptions.
"""

from __future__ import annotations

import json
import textwrap
from pathlib import Path
from typing import List

import pytest

from hlante.parser import (
    ALLELE_REGEX,
    HLAGenotype,
    HLAnteParseError,
    SUPPORTED_TOOLS,
    TOOL_ARCASHLA,
    TOOL_HLAHD,
    TOOL_OPTITYPE,
    TOOL_T1K,
    UnsupportedToolError,
    _determine_resolution,
    _strip_trailing_asterisk,
    _validate_allele,
    parse_arcashla,
    parse_hla_output,
    parse_hlahd,
    parse_optitype,
    parse_t1k,
)


FIXTURES_DIR: Path = Path(__file__).parent / "fixtures"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _by_locus(records: List[HLAGenotype]) -> dict:
    """
    Build a ``{locus: record}`` map to simplify assertions.
    """
    return {rec.locus: rec for rec in records}


# ---------------------------------------------------------------------------
# Regex and resolution
# ---------------------------------------------------------------------------


class TestAlleleRegex:
    """
    Behaviour of :data:`ALLELE_REGEX` and :func:`_validate_allele`.
    """

    @pytest.mark.parametrize(
        "allele",
        [
            "A*02:01",
            "A*02:01:01",
            "A*02:01:01:01",
            "A*02:01:01G",
            "A*02:01P",
            "B*15:01:01",
            "DRB1*15:01",
            "DRB1*15:01:01G",
            "DQB1*06:02:01",
            "DPB1*04:01",
            "C*07:01:01:01",
            "HLA-A*02:01",
            "HLA-DRB1*04:01",
        ],
    )
    def test_valid_alleles(self, allele: str) -> None:
        """
        Accepted allele expressions must pass validation.
        """
        assert _validate_allele(allele), f"{allele!r} was rejected"

    @pytest.mark.parametrize(
        "allele",
        [
            "",
            "A",
            "A*",
            "*02:01",
            "A:02:01",
            "a*02:01",
            "A*2:1",
            "A**02:01",
            "A*02-01",
            "A*02:01:",
        ],
    )
    def test_invalid_alleles(self, allele: str) -> None:
        """
        Invalid expressions must be rejected.
        """
        assert not _validate_allele(allele), f"{allele!r} should be invalid"

    def test_pattern_is_compiled(self) -> None:
        """
        The module-level regex must be pre-compiled.
        """
        assert hasattr(ALLELE_REGEX, "match")


class TestResolution:
    """
    :func:`_determine_resolution` output.
    """

    @pytest.mark.parametrize(
        "allele,expected",
        [
            ("A*02", "2-field"),
            ("A*02:01", "4-field"),
            ("A*02:01:01", "6-field"),
            ("A*02:01:01:01", "8-field"),
            ("A*02:01:01G", "G-group"),
            ("A*02:01P", "P-group"),
            ("DRB1*15:01", "4-field"),
            ("HLA-A*02:01", "4-field"),
            # Three-digit first field: three colon-groups → 6-field.
            # The pre-fix digit-count rule labelled this 8-field
            # because 104+01+01 is 7 digits > 6.
            ("DPB1*104:01:01", "6-field"),
            ("DPB1*02:01:02", "6-field"),
            ("DPB1*104:01", "4-field"),
            ("DPB1*104", "2-field"),
        ],
    )
    def test_resolution(self, allele: str, expected: str) -> None:
        """
        The resolution mapping must match the expected table.
        """
        assert _determine_resolution(allele) == expected

    def test_parser_and_normalizer_agree_on_dpb1_three_digit(self) -> None:
        """
        The parser-side resolution label (string) and the normalizer's
        integer resolution must agree on colon-group semantics. This
        guards the DPB1\\*104:01:01 regression where the parser still
        returned ``"8-field"`` after the normalizer was fixed to use
        colon-group count.
        """
        from hlante.normalizer import _resolution_of
        cases = {
            "DPB1*104:01:01": ("6-field", 6),
            "A*02:01:01:01":  ("8-field", 8),
            "A*02:01":        ("4-field", 4),
            "A*02":           ("2-field", 2),
        }
        for allele, (label, level) in cases.items():
            assert _determine_resolution(allele) == label
            assert _resolution_of(allele) == level


# ---------------------------------------------------------------------------
# ARCAS-HLA
# ---------------------------------------------------------------------------


class TestParseArcasHLA:
    """
    Tests for :func:`parse_arcashla`.
    """

    def test_flat_schema(self) -> None:
        """
        A flat-dictionary JSON must parse successfully.
        """
        records = parse_arcashla(FIXTURES_DIR / "sample.genotype.json")
        assert len(records) == 5

        by_locus = _by_locus(records)
        assert set(by_locus) == {
            "HLA-A",
            "HLA-B",
            "HLA-C",
            "HLA-DRB1",
            "HLA-DQB1",
        }

        a_record = by_locus["HLA-A"]
        assert a_record.allele1 == "A*02:01:01G"
        assert a_record.allele2 == "A*24:02:01G"
        assert a_record.resolution == "G-group"
        assert a_record.tool == TOOL_ARCASHLA
        assert a_record.sample_id == "sample"
        assert a_record.quality_score is None
        assert "HLA-A" in a_record.raw_line

    def test_nested_alleles_schema(self) -> None:
        """
        Nested schema under ``"alleles"`` must also parse.
        """
        records = parse_arcashla(FIXTURES_DIR / "sample_nested.genotype.json")
        assert len(records) == 3
        loci = {r.locus for r in records}
        assert loci == {"HLA-A", "HLA-B", "HLA-C"}

    def test_missing_file(self, tmp_path: Path) -> None:
        """
        A missing path must raise :class:`FileNotFoundError`.
        """
        with pytest.raises(FileNotFoundError, match="ARCAS-HLA"):
            parse_arcashla(tmp_path / "nope.json")

    def test_malformed_json(self) -> None:
        """
        Malformed JSON must raise :class:`HLAnteParseError`.
        """
        with pytest.raises(HLAnteParseError, match="JSON"):
            parse_arcashla(FIXTURES_DIR / "malformed.json")

    def test_invalid_allele_raises(self, tmp_path: Path) -> None:
        """
        JSON with an allele that fails the regex must raise.
        """
        bad = tmp_path / "bad.genotype.json"
        bad.write_text(json.dumps({"HLA-A": ["A_02_01", "A*24:02"]}))
        with pytest.raises(HLAnteParseError, match="Invalid allele"):
            parse_arcashla(bad)

    def test_non_dict_root(self, tmp_path: Path) -> None:
        """
        A non-dict root element must raise.
        """
        bad = tmp_path / "list.json"
        bad.write_text(json.dumps(["A*02:01"]))
        with pytest.raises(HLAnteParseError, match="object"):
            parse_arcashla(bad)

    def test_trailing_asterisk_stripped(self, tmp_path: Path) -> None:
        """
        arcasHLA uncertain-call notation (trailing ``*``) must be stripped
        and the allele parsed normally.
        """
        f = tmp_path / "trailing.genotype.json"
        f.write_text(json.dumps({"DRB1": ["DRB1*04:92*", "DRB1*13:01"]}))
        records = parse_arcashla(f)
        assert len(records) == 1
        assert records[0].allele1 == "DRB1*04:92"
        assert records[0].allele2 == "DRB1*13:01"

    def test_space_separated_ambiguous_pair(self, tmp_path: Path) -> None:
        """
        arcasHLA space-separated ambiguous pair (``"B*49:01 50:01"``) must
        resolve to the first token.
        """
        f = tmp_path / "ambig.genotype.json"
        f.write_text(json.dumps({"B": ["B*49:01 50:01", "B*07:02"]}))
        records = parse_arcashla(f)
        assert len(records) == 1
        assert records[0].allele1 == "B*49:01"
        assert records[0].allele2 == "B*07:02"


# ---------------------------------------------------------------------------
# T1K
# ---------------------------------------------------------------------------


class TestParseT1K:
    """
    Tests for :func:`parse_t1k`.
    """

    def test_basic_parse(self) -> None:
        """
        A T1K TSV with a header row must parse successfully.
        """
        records = parse_t1k(FIXTURES_DIR / "sample_t1k_genotype.tsv")
        assert len(records) == 4

        by_locus = _by_locus(records)
        a = by_locus["HLA-A"]
        assert a.allele1 == "A*02:01"
        assert a.allele2 == "A*24:02"
        assert a.quality_score == pytest.approx(98.5)
        assert a.resolution == "4-field"
        assert a.tool == TOOL_T1K

    def test_missing_second_allele(self) -> None:
        """
        A second allele marked ``*`` must resolve to :data:`None`.
        """
        records = parse_t1k(FIXTURES_DIR / "sample_t1k_genotype.tsv")
        c = _by_locus(records)["HLA-C"]
        assert c.allele1 == "C*07:01"
        assert c.allele2 is None

    def test_missing_file(self, tmp_path: Path) -> None:
        """
        A missing path must raise :class:`FileNotFoundError`.
        """
        with pytest.raises(FileNotFoundError, match="T1K"):
            parse_t1k(tmp_path / "none.tsv")

    def test_native_headerless_layout(self, tmp_path: Path) -> None:
        """
        Real T1K releases emit a headerless 8-column TSV
        (``gene count a1 s1 q1 a2 s2 q2``). The parser must accept
        that format without requiring a ``gene`` header row. Earlier
        versions rejected native T1K output because they assumed the
        5-column headered layout; this regression fires on any real
        user file.
        """
        native = tmp_path / "native_genotype.tsv"
        native.write_text(
            "HLA-A\t2\tHLA-A*11:01:01\t14598.13\t60\tHLA-A*02:01:01\t14112.31\t60\t\n"
            "HLA-J\t1\tHLA-J*01:01:01\t5192.47\t60\t.\t0\t-1\t\n"
        )
        records = parse_t1k(native)
        by_locus = _by_locus(records)
        assert by_locus["HLA-A"].allele1 == "A*11:01:01"
        assert by_locus["HLA-A"].allele2 == "A*02:01:01"
        # count=1 with "." for the second allele must map to None.
        assert by_locus["HLA-J"].allele1 == "J*01:01:01"
        assert by_locus["HLA-J"].allele2 is None

    def test_native_ambiguity_list_uses_primary_call(
        self, tmp_path: Path
    ) -> None:
        """
        T1K can emit an equivalence-class ambiguity list as a single
        comma-separated allele token (``"HLA-DRA*01:01:01,HLA-DRA*01:01:02,
        …"``). The parser must take the first element as the primary
        call instead of rejecting the comma as malformed.
        """
        native = tmp_path / "ambig_genotype.tsv"
        native.write_text(
            "HLA-DRA\t1\t"
            "HLA-DRA*01:01:01,HLA-DRA*01:01:02,HLA-DRA*01:02:02\t"
            "1.23\t1\t.\t0\t-1\t\n"
        )
        records = parse_t1k(native)
        assert records[0].allele1 == "DRA*01:01:01"
        assert records[0].allele2 is None

    def test_short_row_raises(self, tmp_path: Path) -> None:
        """
        A data row with missing columns must raise.
        """
        bad = tmp_path / "short.tsv"
        bad.write_text(
            textwrap.dedent(
                """\
                gene\tallele1\tallele2\tscore1\tscore2
                HLA-A\tA*02:01
                """
            )
        )
        with pytest.raises(HLAnteParseError, match="missing columns"):
            parse_t1k(bad)

    def test_trailing_asterisk_stripped(self, tmp_path: Path) -> None:
        """
        T1K uncertain-call notation (trailing ``*``) must be stripped
        and the allele parsed normally, both in headered and native layout.
        """
        f = tmp_path / "trailing.tsv"
        f.write_text(
            "gene\tallele1\tallele2\tscore1\tscore2\n"
            "B\tB*40:02*\tB*40:04\t100\t100\n"
        )
        records = parse_t1k(f)
        assert records[0].allele1 == "B*40:02"
        assert records[0].allele2 == "B*40:04"

    def test_trailing_asterisk_stripped_native(self, tmp_path: Path) -> None:
        """
        Trailing ``*`` must also be stripped in the native headerless layout.
        """
        f = tmp_path / "trailing_native.tsv"
        f.write_text("HLA-B\t2\tHLA-B*40:02*\t100\t60\tHLA-B*40:04\t100\t60\n")
        records = parse_t1k(f)
        assert records[0].allele1 == "B*40:02"
        assert records[0].allele2 == "B*40:04"


# ---------------------------------------------------------------------------
# HLA-HD
# ---------------------------------------------------------------------------


class TestParseHLAHD:
    """
    Tests for :func:`parse_hlahd`.
    """

    def test_basic_parse(self) -> None:
        """
        An HLA-HD final.result.txt file must parse successfully.
        """
        records = parse_hlahd(FIXTURES_DIR / "sample_final.result.txt")
        # DPB1 is fully "-" → skipped.
        assert len(records) == 5

        by_locus = _by_locus(records)
        assert "HLA-DPB1" not in by_locus
        assert by_locus["HLA-A"].allele1 == "A*02:01:01"
        assert by_locus["HLA-A"].allele2 == "A*24:02:01"
        assert by_locus["HLA-A"].resolution == "6-field"
        assert by_locus["HLA-A"].tool == TOOL_HLAHD

    def test_ambiguous_second_allele(self) -> None:
        """
        A ``-`` second allele must become :data:`None`.
        """
        records = parse_hlahd(FIXTURES_DIR / "sample_final.result.txt")
        c = _by_locus(records)["HLA-C"]
        assert c.allele1 == "C*07:01:01"
        assert c.allele2 is None

    def test_sample_id_strips_suffix(self) -> None:
        """
        For ``sample_final.result.txt`` the derived sample_id must be
        ``"sample"``.
        """
        records = parse_hlahd(FIXTURES_DIR / "sample_final.result.txt")
        assert all(r.sample_id == "sample" for r in records)

    def test_missing_file(self, tmp_path: Path) -> None:
        """
        A missing path must raise :class:`FileNotFoundError`.
        """
        with pytest.raises(FileNotFoundError, match="HLA-HD"):
            parse_hlahd(tmp_path / "none.txt")

    def test_invalid_allele_raises(self, tmp_path: Path) -> None:
        """
        A row with an invalid allele must raise
        :class:`HLAnteParseError`.
        """
        bad = tmp_path / "bad_final.result.txt"
        bad.write_text("HLA-A\tBOGUS\t-\n")
        with pytest.raises(HLAnteParseError, match="Invalid allele"):
            parse_hlahd(bad)

    def test_empty_file_raises(self, tmp_path: Path) -> None:
        """
        A file where every locus is ``-`` must raise.
        """
        bad = tmp_path / "empty_final.result.txt"
        bad.write_text("HLA-A\t-\t-\nHLA-B\t-\t-\n")
        with pytest.raises(HLAnteParseError, match="No typed locus"):
            parse_hlahd(bad)

    def test_not_typed_null_token_with_embedded_space(
        self, tmp_path: Path
    ) -> None:
        """
        HLA-HD writes the null marker as the two-word token
        ``"Not typed"`` in tab-delimited files. The earlier regex
        ``r"[\\t ]+"`` split the interior space as a column boundary
        and produced ``['DRB5', 'Not', 'typed', 'Not', 'typed']``,
        then tripped on ``"Not"`` as an allele. The row must be
        recognised as an untyped locus and skipped cleanly.
        """
        f = tmp_path / "sample_final.result.txt"
        f.write_text(
            "A\tHLA-A*01:01:01\tHLA-A*11:01:01\n"
            "DRB5\tNot typed\tNot typed\n"
            "DRB6\tNot typed\tNot typed\n"
        )
        records = parse_hlahd(f)
        assert len(records) == 1
        assert records[0].locus == "HLA-A"
        assert records[0].allele1 == "A*01:01:01"

    def test_not_consistent_null_token(self, tmp_path: Path) -> None:
        """
        HLA-HD also emits ``"Not consistent"`` for ambiguous/unresolved
        loci; it must be recognised as a null token and skipped, not
        parsed as an allele.
        """
        f = tmp_path / "nc_final.result.txt"
        f.write_text(
            "A\tHLA-A*02:01:01\tHLA-A*24:02:01\n"
            "DRB1\tNot consistent\tNot consistent\n"
        )
        records = parse_hlahd(f)
        assert len(records) == 1
        assert records[0].locus == "HLA-A"

    def test_trailing_asterisk_stripped(self, tmp_path: Path) -> None:
        """
        HLA-HD uncertain-call notation (trailing ``*``) must be stripped
        and the allele parsed normally.
        """
        f = tmp_path / "trailing_final.result.txt"
        f.write_text("A\tHLA-A*02:01*\tHLA-A*24:02\n")
        records = parse_hlahd(f)
        assert records[0].allele1 == "A*02:01"
        assert records[0].allele2 == "A*24:02"


# ---------------------------------------------------------------------------
# _strip_trailing_asterisk
# ---------------------------------------------------------------------------


class TestStripTrailingAsterisk:
    """
    Unit tests for the shared :func:`_strip_trailing_asterisk` helper.
    """

    def test_strips_trailing_asterisk(self) -> None:
        assert _strip_trailing_asterisk("B*40:02*", "t1k") == "B*40:02"

    def test_no_change_when_absent(self) -> None:
        assert _strip_trailing_asterisk("B*40:02", "t1k") == "B*40:02"

    def test_strips_only_trailing(self) -> None:
        # The asterisk between gene and fields must NOT be removed.
        assert _strip_trailing_asterisk("DRB1*04:92*", "hlahd") == "DRB1*04:92"

    def test_empty_string(self) -> None:
        assert _strip_trailing_asterisk("", "t1k") == ""


# ---------------------------------------------------------------------------
# OptiType
# ---------------------------------------------------------------------------


class TestParseOptiType:
    """
    Tests for :func:`parse_optitype`.
    """

    def test_basic_parse(self) -> None:
        """
        An OptiType result file must parse successfully.
        """
        records = parse_optitype(FIXTURES_DIR / "sample_result.tsv")
        assert len(records) == 3

        by_locus = _by_locus(records)
        assert set(by_locus) == {"HLA-A", "HLA-B", "HLA-C"}

        a = by_locus["HLA-A"]
        assert a.allele1 == "A*02:01"
        assert a.allele2 == "A*24:02"
        assert a.resolution == "4-field"
        assert a.tool == TOOL_OPTITYPE
        assert a.quality_score == pytest.approx(1456.78)
        assert a.sample_id == "sample"

    def test_missing_file(self, tmp_path: Path) -> None:
        """
        A missing path must raise :class:`FileNotFoundError`.
        """
        with pytest.raises(FileNotFoundError, match="OptiType"):
            parse_optitype(tmp_path / "none.tsv")

    def test_missing_required_columns(self, tmp_path: Path) -> None:
        """
        A missing required column must raise.
        """
        bad = tmp_path / "opti_result.tsv"
        bad.write_text("\tA1\tA2\tReads\n0\tA*02:01\tA*24:02\t100\n")
        with pytest.raises(HLAnteParseError, match="missing required"):
            parse_optitype(bad)

    def test_insufficient_rows(self, tmp_path: Path) -> None:
        """
        A file with only a header row must raise.
        """
        bad = tmp_path / "opti_result.tsv"
        bad.write_text("\tA1\tA2\tB1\tB2\tC1\tC2\tReads\tObjective\n")
        with pytest.raises(HLAnteParseError, match="header"):
            parse_optitype(bad)


# ---------------------------------------------------------------------------
# Dispatcher
# ---------------------------------------------------------------------------


class TestDispatcher:
    """
    Behaviour of :func:`parse_hla_output`.
    """

    @pytest.mark.parametrize(
        "tool,fixture,count",
        [
            ("arcashla", "sample.genotype.json", 5),
            ("arcas-hla", "sample.genotype.json", 5),
            ("t1k", "sample_t1k_genotype.tsv", 4),
            ("hla-hd", "sample_final.result.txt", 5),
            ("HLAHD", "sample_final.result.txt", 5),
            ("optitype", "sample_result.tsv", 3),
        ],
    )
    def test_dispatch(self, tool: str, fixture: str, count: int) -> None:
        """
        Different aliases must dispatch to the correct parser.
        """
        records = parse_hla_output(FIXTURES_DIR / fixture, tool)
        assert len(records) == count

    def test_unsupported_tool_raises(self, tmp_path: Path) -> None:
        """
        An unsupported tool must raise :class:`UnsupportedToolError`.
        """
        dummy = tmp_path / "anything.txt"
        dummy.write_text("foo")
        with pytest.raises(UnsupportedToolError, match="Unsupported tool"):
            parse_hla_output(dummy, "madeuptool")

    def test_missing_file_raises(self, tmp_path: Path) -> None:
        """
        A valid tool but missing file must raise :class:`FileNotFoundError`.
        """
        with pytest.raises(FileNotFoundError):
            parse_hla_output(tmp_path / "ghost.tsv", "t1k")

    def test_supported_tools_set(self) -> None:
        """
        :data:`SUPPORTED_TOOLS` must contain the expected four tools.
        """
        assert SUPPORTED_TOOLS == {
            TOOL_ARCASHLA,
            TOOL_T1K,
            TOOL_HLAHD,
            TOOL_OPTITYPE,
        }
