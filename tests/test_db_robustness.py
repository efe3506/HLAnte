"""
tests.test_db_robustness
========================

Software-phase DB-layer fixes:

- curated variation_id is deterministic across runs (stable hash).
- AFND/NMDP population matching uses precedence-based canonical
  classification, so "African American" maps to AFR only (never AMR).
"""

from __future__ import annotations

import hashlib
from pathlib import Path

import pytest

from hlante.cli import _parse_all
from hlante.db import DatabaseIntegrityError, atomic_install, sha256_file
from hlante.db.afnd import AFNDClient, _classify_population_group
from hlante.db.curated import _stable_curated_id
from hlante.parser import HLAnteParseError, _read_text, parse_hla_output
from hlante.reporter import _sanitize_cell

FIXTURES = Path(__file__).parent / "fixtures"


class TestStableCuratedId:
    def test_deterministic(self) -> None:
        a = _stable_curated_id("B*57:01", "Abacavir hypersensitivity")
        b = _stable_curated_id("B*57:01", "Abacavir hypersensitivity")
        assert a == b
        assert a.startswith("CURATED_")

    def test_distinct_inputs_distinct_ids(self) -> None:
        assert _stable_curated_id("B*57:01", "X") != _stable_curated_id("B*58:01", "X")


class TestAfndGroupMatching:
    def setup_method(self) -> None:
        self.client = AFNDClient.__new__(AFNDClient)

    def test_african_american_classifies_as_afr(self) -> None:
        assert _classify_population_group("African American") == "AFR"

    def test_african_american_does_not_leak_into_amr(self) -> None:
        # The C2 bug: "american" substring of "African American" matched AMR.
        assert self.client._matches_group("AFR African American", "AMR") is False
        assert self.client._matches_group("AFR African American", "AFR") is True

    def test_true_amr_still_matches(self) -> None:
        assert self.client._matches_group("AMR Mexico Mestizo", "AMR") is True
        assert self.client._matches_group("USA Hispanic", "AMR") is True

    def test_east_asian_matches_asn(self) -> None:
        assert self.client._matches_group("Japan pop 1 East Asian", "ASN") is True
        assert self.client._matches_group("Japan pop 1 East Asian", "AMR") is False

    def test_global_matches_everything(self) -> None:
        assert self.client._matches_group("anything at all", "global") is True


# ---------------------------------------------------------------------------
# Encoding robustness (BOM / UTF-16)
# ---------------------------------------------------------------------------


class TestEncodingRobustness:
    def test_utf16_decoded(self, tmp_path: Path) -> None:
        p = tmp_path / "x.json"
        p.write_bytes('{"A": ["A*02:01"]}'.encode("utf-16"))
        assert "A*02:01" in _read_text(p)

    def test_utf8_bom_stripped(self, tmp_path: Path) -> None:
        p = tmp_path / "y.txt"
        p.write_bytes(b"\xef\xbb\xbf" + "hello".encode("utf-8"))
        assert _read_text(p) == "hello"

    def test_undecodable_raises_parse_error(self, tmp_path: Path) -> None:
        p = tmp_path / "z.bin"
        p.write_bytes(b"\x80\x81\x82 not utf8")
        with pytest.raises(HLAnteParseError):
            _read_text(p)

    def test_arcashla_parses_utf16_without_crash(self, tmp_path: Path) -> None:
        p = tmp_path / "s.genotype.json"
        p.write_bytes('{"A": ["A*02:01:01:01", "A*01:01:01:01"]}'.encode("utf-16"))
        genos = parse_hla_output(p, "arcashla")
        assert genos


# ---------------------------------------------------------------------------
# CSV/TSV formula-injection neutralisation
# ---------------------------------------------------------------------------


class TestFormulaInjection:
    @pytest.mark.parametrize(
        "payload",
        ['=cmd|"/c calc"!A1', "+1+1", "-2+3", "@SUM(A1)", "\tTAB", "\rCR"],
    )
    def test_formula_cells_neutralised(self, payload: str) -> None:
        assert _sanitize_cell(payload).startswith("'")

    @pytest.mark.parametrize("benign", ["B*57:01", "Abacavir hypersensitivity", "0.05", "NA"])
    def test_benign_cells_unchanged(self, benign: str) -> None:
        assert _sanitize_cell(benign) == benign


# ---------------------------------------------------------------------------
# One bad file does not abort the cohort
# ---------------------------------------------------------------------------


class TestCohortResilience:
    def test_skips_bad_file_keeps_good(self) -> None:
        files = [FIXTURES / "sample.genotype.json", FIXTURES / "malformed.json"]
        genos, failures = _parse_all(files, "arcashla", strict=False)
        assert genos, "the valid file should still be parsed"
        assert len(failures) == 1
        assert failures[0][0].name == "malformed.json"

    def test_strict_mode_raises(self) -> None:
        files = [FIXTURES / "sample.genotype.json", FIXTURES / "malformed.json"]
        with pytest.raises(HLAnteParseError):
            _parse_all(files, "arcashla", strict=True)

    def test_all_bad_yields_no_genotypes(self) -> None:
        genos, failures = _parse_all([FIXTURES / "malformed.json"], "arcashla", strict=False)
        assert genos == []
        assert len(failures) == 1


# ---------------------------------------------------------------------------
# db-update integrity (checksum + .bak rollback + atomic verify)
# ---------------------------------------------------------------------------


class TestDbIntegrity:
    def test_sha256_file(self, tmp_path: Path) -> None:
        p = tmp_path / "f"
        p.write_bytes(b"hello")
        assert sha256_file(p) == hashlib.sha256(b"hello").hexdigest()

    def test_atomic_install_backs_up_and_swaps(self, tmp_path: Path) -> None:
        dest = tmp_path / "db.txt"
        dest.write_text("OLD")
        staged = tmp_path / "db.txt.part"
        staged.write_text("NEW")
        digest = atomic_install(staged, dest)
        assert dest.read_text() == "NEW"
        assert digest == hashlib.sha256(b"NEW").hexdigest()
        assert (tmp_path / "db.txt.bak").read_text() == "OLD"
        assert not staged.exists()

    def test_atomic_install_no_existing_dest(self, tmp_path: Path) -> None:
        dest = tmp_path / "new.txt"
        staged = tmp_path / "new.txt.part"
        staged.write_text("DATA")
        atomic_install(staged, dest)
        assert dest.read_text() == "DATA"
        assert not (tmp_path / "new.txt.bak").exists()

    def test_atomic_install_rolls_back_on_checksum_mismatch(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        dest = tmp_path / "db.txt"
        dest.write_text("OLD")
        staged = tmp_path / "db.txt.part"
        staged.write_text("NEW")
        # Force the staged vs installed checksums to differ.
        seq = iter(["staged-digest", "installed-digest"])
        monkeypatch.setattr("hlante.db.sha256_file", lambda _p: next(seq))
        with pytest.raises(DatabaseIntegrityError):
            atomic_install(staged, dest)
        # The previous file must be restored from the .bak.
        assert dest.read_text() == "OLD"


class TestIMGTRefNormalization:
    """
    The ANHIG mirror names release branches without separators, so the
    release number quoted in the manuscript and the documentation has to be
    translated before it is fetched.
    """

    def test_dotted_release_maps_to_branch_name(self) -> None:
        from hlante.db.imgt import normalize_imgt_ref

        assert normalize_imgt_ref("3.64.0") == "3640"
        assert normalize_imgt_ref("3.65.0") == "3650"

    def test_other_refs_pass_through_unchanged(self) -> None:
        from hlante.db.imgt import normalize_imgt_ref

        for ref in ("Latest", "3640", "v3.64.0-alpha", "a1b2c3d"):
            assert normalize_imgt_ref(ref) == ref


class TestIMGTRefPinning:
    """
    ``--imgt-ref`` has to replace the installed release, and ``version.json``
    has to describe the files that are actually on disk.

    The failure these cover was silent in both directions: asking for a
    different release left the cached copy in place, and the metadata was
    rewritten with the requested ref anyway, so a snapshot recorded as 3.65.0
    could hold 3.64.0 data together with 3.64.0 checksums.
    """

    @staticmethod
    def _fake_fetch(monkeypatch: pytest.MonkeyPatch, version: str) -> list:
        """Record every fetch and write an Allelelist carrying ``version``."""
        from hlante.db import imgt

        calls: list = []

        def _fetch(url: str, dest: Path, timeout: int = 60) -> Path:
            calls.append(url)
            if dest.name == "Allelelist.txt":
                dest.write_text(f"# version: {version}\nA*01:01:01:01,HLA00001\n")
            else:
                dest.write_text(f"# {version}\n")
            return dest

        monkeypatch.setattr(imgt, "_http_fetch", _fetch)
        return calls

    def _meta(self, root: Path) -> dict:
        import json

        return json.loads((root / "version.json").read_text())

    def test_changing_ref_refreshes_without_force(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        from hlante.db.imgt import download_imgt_db

        self._fake_fetch(monkeypatch, "IPD-IMGT/HLA 3.64.0")
        download_imgt_db(tmp_path, ref="3.64.0")
        first = self._meta(tmp_path)

        calls = self._fake_fetch(monkeypatch, "IPD-IMGT/HLA 3.65.0")
        download_imgt_db(tmp_path, ref="3.65.0")
        second = self._meta(tmp_path)

        assert calls, "a differing ref must re-fetch even without force"
        assert all("/3650/" in url for url in calls)
        assert second["ref"] == "3650"
        assert second["version"] == "IPD-IMGT/HLA 3.65.0"
        assert second["sha256"] != first["sha256"]

    def test_metadata_never_records_an_unfetched_ref(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        from hlante.db.imgt import download_imgt_db

        self._fake_fetch(monkeypatch, "IPD-IMGT/HLA 3.64.0")
        download_imgt_db(tmp_path, ref="3.64.0")

        # The requested release is never fetched, so it must not be recorded.
        self._fake_fetch(monkeypatch, "IPD-IMGT/HLA 3.65.0")
        download_imgt_db(tmp_path, ref="3.65.0")
        meta = self._meta(tmp_path)

        recorded_release = meta["ref"].replace("3650", "3.65.0").replace("3640", "3.64.0")
        assert recorded_release in str(meta["version"]), (
            "version.json records a ref that does not match the installed data"
        )
        assert meta["source_base"].endswith(meta["ref"])

    def test_same_ref_skips_and_keeps_the_acquisition_date(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        from hlante.db.imgt import download_imgt_db

        self._fake_fetch(monkeypatch, "IPD-IMGT/HLA 3.64.0")
        download_imgt_db(tmp_path, ref="3.64.0")
        first = self._meta(tmp_path)

        calls = self._fake_fetch(monkeypatch, "IPD-IMGT/HLA 3.64.0")
        download_imgt_db(tmp_path, ref="3.64.0")
        second = self._meta(tmp_path)

        assert calls == [], "an unchanged ref must not re-download"
        assert second["downloaded_at"] == first["downloaded_at"]
        assert second["sha256"] == first["sha256"]

    def test_force_redownloads_the_same_ref(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        from hlante.db.imgt import download_imgt_db

        self._fake_fetch(monkeypatch, "IPD-IMGT/HLA 3.64.0")
        download_imgt_db(tmp_path, ref="3.64.0")

        calls = self._fake_fetch(monkeypatch, "IPD-IMGT/HLA 3.64.0")
        download_imgt_db(tmp_path, ref="3.64.0", force=True)

        assert calls, "force must re-fetch even when the ref is unchanged"

    def test_files_without_metadata_are_replaced(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """A copy we cannot vouch for is not a copy we may describe."""
        from hlante.db.imgt import download_imgt_db

        (tmp_path / "Allelelist.txt").write_text("# version: unknown\n")
        calls = self._fake_fetch(monkeypatch, "IPD-IMGT/HLA 3.64.0")
        download_imgt_db(tmp_path, ref="3.64.0")

        assert calls, "an install with no version.json must be re-fetched"
        assert self._meta(tmp_path)["version"] == "IPD-IMGT/HLA 3.64.0"
