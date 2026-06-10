"""
tests.test_normalizer
=====================

Unit tests for :mod:`hlante.normalizer`.

Coverage
--------
- ``load_imgt_db`` load behaviour and missing-file errors.
- ``normalize_allele``: exact match, prefix-ambiguous, G/P group,
  novel, and null alleles.
- ``resolve_ambiguity``: group expansion and prefix search.
- ``batch_normalize``: thread-pool parallelism, order preservation,
  and ``None`` skipping.
- Stale-copy (>6 months) warning.
"""

from __future__ import annotations

import json
import logging
import shutil
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Dict

import pytest

from hlante.db.imgt import VERSION_FILENAME
from hlante.normalizer import (
    IMGTDatabaseMissingError,
    InvalidAlleleError,
    NormalizedAllele,
    STALE_THRESHOLD_DAYS,
    batch_normalize,
    load_imgt_db,
    normalize_allele,
    resolve_ambiguity,
)
from hlante.normalizer import _resolution_of
from hlante.parser import HLAGenotype, parse_arcashla


FIXTURES_DIR: Path = Path(__file__).parent / "fixtures"
IMGT_MINI_DIR: Path = FIXTURES_DIR / "imgt_mini"


# ---------------------------------------------------------------------------
# P0-3 — resolution based on colon-group count
# ---------------------------------------------------------------------------


class TestResolutionColonGroupCount:
    """
    :func:`_resolution_of` must return field-level resolution based
    on colon-group count, not digit count. Before this fix,
    ``DPB1*104:01:01`` was wrongly labelled 8-field because its
    three colon-groups contain 7 digits (3+2+2).
    """

    def test_dpb1_three_digit_first_field_resolution(self) -> None:
        assert _resolution_of("DPB1*104:01:01") == 6

    @pytest.mark.parametrize(
        "allele,expected",
        [
            ("A*02", 2),
            ("A*02:01", 4),
            ("A*02:01:01", 6),
            ("A*02:01:01:01", 8),
            ("B*57:01G", 4),        # G-suffix stripped → still 2 fields
            ("A*02:01P", 4),        # P-suffix stripped
            ("A*02:01N", 4),        # nomenclature suffix stripped
            ("DRB1*08:04:01", 6),
        ],
    )
    def test_resolution_colon_group_count(
        self, allele: str, expected: int
    ) -> None:
        assert _resolution_of(allele) == expected


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def imgt_db() -> Dict[str, object]:
    """
    Load the mini IPD-IMGT/HLA fixture once per module.
    """
    return load_imgt_db(IMGT_MINI_DIR)


@pytest.fixture()
def imgt_db_dir(tmp_path: Path) -> Path:
    """
    A writable per-test copy of the fixture (used for version.json tests).
    """
    dest = tmp_path / "imgt"
    shutil.copytree(IMGT_MINI_DIR, dest)
    return dest


# ---------------------------------------------------------------------------
# load_imgt_db
# ---------------------------------------------------------------------------


class TestLoadIMGTDB:
    """
    Behaviour of :func:`load_imgt_db`.
    """

    def test_loads_alleles(self, imgt_db: Dict[str, object]) -> None:
        """
        The allele dictionary must contain the expected accessions.
        """
        alleles = imgt_db["alleles"]
        assert alleles["A*01:01:01:01"] == "HLA00001"
        assert alleles["DRB1*15:01:01:01"] == "HLA00100"
        assert len(alleles) == 20

    def test_loads_g_groups(self, imgt_db: Dict[str, object]) -> None:
        """
        The G-group file must be parsed and the reverse map populated.
        """
        g_groups = imgt_db["g_groups"]
        assert "A*02:01:01G" in g_groups
        assert "A*02:01:01:01" in g_groups["A*02:01:01G"]

        allele_to_g = imgt_db["allele_to_g"]
        assert allele_to_g["A*02:01:01:01"] == "A*02:01:01G"

    def test_loads_p_groups(self, imgt_db: Dict[str, object]) -> None:
        """
        The P-group file must be parsed.
        """
        p_groups = imgt_db["p_groups"]
        assert "A*02:01P" in p_groups

    def test_version_and_path(self, imgt_db: Dict[str, object]) -> None:
        """
        The version must be extracted from the Allelelist header.
        """
        assert imgt_db["version"] == "IPD-IMGT/HLA 3.55.0"
        assert Path(imgt_db["path"]) == IMGT_MINI_DIR

    def test_missing_dir_raises(self, tmp_path: Path) -> None:
        """
        Missing files must raise :class:`IMGTDatabaseMissingError`.
        """
        with pytest.raises(IMGTDatabaseMissingError, match="Allelelist"):
            load_imgt_db(tmp_path / "none")

    def test_accepts_file_path(self) -> None:
        """
        ``db_path`` may point directly at an Allelelist.txt file.
        """
        db = load_imgt_db(IMGT_MINI_DIR / "Allelelist.txt")
        assert db["alleles"]["A*01:01:01:01"] == "HLA00001"

    def test_stale_warning(
        self,
        imgt_db_dir: Path,
        caplog: pytest.LogCaptureFixture,
    ) -> None:
        """
        A 200-day-old ``downloaded_at`` must emit a warning log.
        """
        stale_time = datetime.now(timezone.utc) - timedelta(
            days=STALE_THRESHOLD_DAYS + 20
        )
        (imgt_db_dir / VERSION_FILENAME).write_text(
            json.dumps(
                {
                    "version": "3.54.0",
                    "downloaded_at": stale_time.isoformat(),
                }
            ),
            encoding="utf-8",
        )
        with caplog.at_level(logging.WARNING, logger="hlante.normalizer"):
            db = load_imgt_db(imgt_db_dir)
        assert db["is_stale"] is True
        assert any("old" in rec.message.lower() for rec in caplog.records)

    def test_fresh_no_warning(
        self,
        imgt_db_dir: Path,
        caplog: pytest.LogCaptureFixture,
    ) -> None:
        """
        A fresh download must not produce a staleness warning.
        """
        fresh = datetime.now(timezone.utc) - timedelta(days=1)
        (imgt_db_dir / VERSION_FILENAME).write_text(
            json.dumps(
                {
                    "version": "3.55.0",
                    "downloaded_at": fresh.isoformat(),
                }
            ),
            encoding="utf-8",
        )
        with caplog.at_level(logging.WARNING, logger="hlante.normalizer"):
            db = load_imgt_db(imgt_db_dir)
        assert db["is_stale"] is False
        assert not any("old" in rec.message.lower() for rec in caplog.records)


# ---------------------------------------------------------------------------
# normalize_allele
# ---------------------------------------------------------------------------


class TestNormalizeAllele:
    """
    Behaviour of :func:`normalize_allele`.
    """

    def test_exact_match(self, imgt_db: Dict[str, object]) -> None:
        """
        Exact IMGT match must populate ``imgt_accession``.
        """
        norm = normalize_allele("A*02:01:01:01", imgt_db)
        assert isinstance(norm, NormalizedAllele)
        assert norm.imgt_accession == "HLA00004"
        assert norm.gene == "HLA-A"
        assert norm.hla_class == "I"
        assert norm.resolution_level == 8
        assert norm.is_ambiguous is False
        assert norm.is_novel is False
        assert norm.protein_group == "A*02:01:01G"

    def test_class_ii_gene(self, imgt_db: Dict[str, object]) -> None:
        """
        DRB1-prefixed alleles must be classified as Class II.
        """
        norm = normalize_allele("DRB1*15:01:01:01", imgt_db)
        assert norm is not None
        assert norm.hla_class == "II"
        assert norm.gene == "HLA-DRB1"

    def test_strips_hla_prefix(self, imgt_db: Dict[str, object]) -> None:
        """
        The ``HLA-`` prefix must be stripped before lookup.
        """
        norm = normalize_allele("HLA-A*01:01:01:01", imgt_db)
        assert norm is not None
        assert norm.allele_name == "A*01:01:01:01"
        assert norm.imgt_accession == "HLA00001"

    def test_low_resolution_is_ambiguous(
        self, imgt_db: Dict[str, object]
    ) -> None:
        """
        An ``A*02`` 2-field allele must be flagged as ambiguous.
        """
        norm = normalize_allele("A*02", imgt_db)
        assert norm is not None
        assert norm.resolution_level == 2
        assert norm.is_ambiguous is True
        assert norm.imgt_accession is None
        assert norm.is_novel is False

    def test_prefix_match_without_exact(
        self, imgt_db: Dict[str, object]
    ) -> None:
        """
        ``A*02:01`` has no exact record but prefix matches exist →
        ambiguous but not novel.
        """
        norm = normalize_allele("A*02:01", imgt_db)
        assert norm is not None
        assert norm.resolution_level == 4
        assert norm.is_ambiguous is True
        assert norm.is_novel is False

    def test_g_group_notation(self, imgt_db: Dict[str, object]) -> None:
        """
        G-group inputs must populate ``protein_group``.
        """
        norm = normalize_allele("A*02:01:01G", imgt_db)
        assert norm is not None
        assert norm.protein_group == "A*02:01:01G"
        assert norm.is_ambiguous is True
        assert norm.imgt_accession == "HLA00004"  # first group member

    def test_p_group_notation(self, imgt_db: Dict[str, object]) -> None:
        """
        P-group inputs must be ambiguous and carry a member accession.
        """
        norm = normalize_allele("A*02:01P", imgt_db)
        assert norm is not None
        assert norm.is_ambiguous is True
        assert norm.imgt_accession == "HLA00004"

    def test_novel_allele(self, imgt_db: Dict[str, object]) -> None:
        """
        No match anywhere → ``is_novel=True``.
        """
        norm = normalize_allele("A*99:99:99:99", imgt_db)
        assert norm is not None
        assert norm.is_novel is True
        assert norm.is_ambiguous is True
        assert norm.imgt_accession is None

    @pytest.mark.parametrize(
        "token",
        ["*", "-", "Not typed", "NA", "", None, "  "],
    )
    def test_null_tokens_return_none(
        self,
        imgt_db: Dict[str, object],
        token,
        caplog: pytest.LogCaptureFixture,
    ) -> None:
        """
        Null-token inputs must return ``None`` with a debug log.
        """
        with caplog.at_level(logging.DEBUG, logger="hlante.normalizer"):
            result = normalize_allele(token, imgt_db)
        assert result is None

    def test_invalid_allele_raises(
        self, imgt_db: Dict[str, object]
    ) -> None:
        """
        Non-IMGT expressions must raise :class:`InvalidAlleleError`.
        """
        with pytest.raises(InvalidAlleleError):
            normalize_allele("A_02_01", imgt_db)


# ---------------------------------------------------------------------------
# resolve_ambiguity
# ---------------------------------------------------------------------------


class TestResolveAmbiguity:
    """
    Behaviour of :func:`resolve_ambiguity`.
    """

    def test_g_group_expands_to_members(
        self, imgt_db: Dict[str, object]
    ) -> None:
        """
        A G-group name must expand to its member list.
        """
        members = resolve_ambiguity("A*02:01:01G", imgt_db)
        assert set(members) == {
            "A*02:01:01:01",
            "A*02:01:01:02",
            "A*02:01:02:01",
        }

    def test_p_group_expands(self, imgt_db: Dict[str, object]) -> None:
        """
        A P-group name must expand to its members.
        """
        members = resolve_ambiguity("A*02:01P", imgt_db)
        assert "A*02:01:01:01" in members
        assert len(members) == 3

    def test_prefix_match(self, imgt_db: Dict[str, object]) -> None:
        """
        ``A*02`` must return every ``A*02:*`` full allele.
        """
        members = resolve_ambiguity("A*02", imgt_db)
        assert len(members) >= 3
        assert all(m.startswith("A*02:") for m in members)

    def test_exact_full_allele_returned_as_is(
        self, imgt_db: Dict[str, object]
    ) -> None:
        """
        A fully known allele must be returned as a single-item list.
        """
        assert resolve_ambiguity("A*01:01:01:01", imgt_db) == ["A*01:01:01:01"]

    def test_unknown_prefix_returns_empty(
        self, imgt_db: Dict[str, object]
    ) -> None:
        """
        An unknown prefix must return an empty list.
        """
        assert resolve_ambiguity("A*99", imgt_db) == []

    def test_null_token_returns_empty(
        self, imgt_db: Dict[str, object]
    ) -> None:
        """
        A null-token input must return an empty list.
        """
        assert resolve_ambiguity("-", imgt_db) == []


# ---------------------------------------------------------------------------
# batch_normalize
# ---------------------------------------------------------------------------


class TestBatchNormalize:
    """
    Behaviour of :func:`batch_normalize`.
    """

    def test_from_arcas_fixture(self, imgt_db: Dict[str, object]) -> None:
        """
        Genotypes parsed from the ARCAS fixture must batch-normalize.
        """
        genotypes = parse_arcashla(FIXTURES_DIR / "sample.genotype.json")
        normalized = batch_normalize(genotypes, imgt_db=imgt_db)
        # Five loci × two alleles = 10 normalized alleles.
        assert len(normalized) == 10
        assert all(isinstance(n, NormalizedAllele) for n in normalized)
        genes = {n.gene for n in normalized}
        assert {"HLA-A", "HLA-B", "HLA-C", "HLA-DRB1", "HLA-DQB1"} <= genes

    def test_skips_null_allele2(self, imgt_db: Dict[str, object]) -> None:
        """
        Inputs where ``allele2`` is :data:`None` must produce one result.
        """
        genotype = HLAGenotype(
            sample_id="s1",
            locus="HLA-A",
            allele1="A*02:01:01:01",
            allele2=None,
            resolution="8-field",
            quality_score=None,
            tool="t1k",
            raw_line="",
        )
        normalized = batch_normalize([genotype], imgt_db=imgt_db)
        assert len(normalized) == 1
        assert normalized[0].imgt_accession == "HLA00004"

    def test_empty_input(self, imgt_db: Dict[str, object]) -> None:
        """
        An empty input must yield an empty list.
        """
        assert batch_normalize([], imgt_db=imgt_db) == []

    def test_loads_db_when_not_provided(
        self, imgt_db_dir: Path
    ) -> None:
        """
        The DB must be loaded from ``db_path`` when ``imgt_db`` is absent.
        """
        genotype = HLAGenotype(
            sample_id="s1",
            locus="HLA-A",
            allele1="A*01:01:01:01",
            allele2=None,
            resolution="8-field",
            quality_score=None,
            tool="t1k",
            raw_line="",
        )
        normalized = batch_normalize([genotype], db_path=imgt_db_dir)
        assert len(normalized) == 1
        assert normalized[0].imgt_accession == "HLA00001"

    def test_parallel_deterministic_order(
        self, imgt_db: Dict[str, object]
    ) -> None:
        """
        Output order must follow the allele1 → allele2 sequence.
        """
        genotypes = [
            HLAGenotype("s1", "HLA-A", "A*01:01:01:01", "A*02:01:01:01",
                        "8-field", None, "t1k", ""),
            HLAGenotype("s1", "HLA-B", "B*07:02:01:01", "B*15:01:01:01",
                        "8-field", None, "t1k", ""),
        ]
        normalized = batch_normalize(genotypes, imgt_db=imgt_db, max_workers=4)
        assert [n.allele_name for n in normalized] == [
            "A*01:01:01:01",
            "A*02:01:01:01",
            "B*07:02:01:01",
            "B*15:01:01:01",
        ]

    def test_invalid_allele_propagates(
        self, imgt_db: Dict[str, object]
    ) -> None:
        """
        Any invalid allele must abort the batch with an exception.
        """
        genotypes = [
            HLAGenotype("s1", "HLA-A", "BOGUS", None,
                        "8-field", None, "t1k", ""),
        ]
        with pytest.raises(InvalidAlleleError):
            batch_normalize(genotypes, imgt_db=imgt_db)
