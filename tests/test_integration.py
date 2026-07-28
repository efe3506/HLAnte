"""
tests.test_integration
======================

End-to-end integration tests driven through :mod:`hlante.cli`.

Approach
--------
- Every command runs via :class:`click.testing.CliRunner`.
- Tests run in ``--offline`` mode; no network access.
- IPD-IMGT/HLA and PharmGKB use the mini fixtures shipped with the
  repository.
- Network-dependent aspects of ``db-update`` are covered by
  monkeypatching the download functions inside ``hlante.cli``.
"""

from __future__ import annotations

import csv
import json
import shutil
from pathlib import Path
from typing import List

import pytest
from click.testing import CliRunner

from hlante import __version__
from hlante.cli import cli
from hlante.db.imgt import ALLELE_LIST_FILENAME


FIXTURES_DIR: Path = Path(__file__).parent / "fixtures"
IMGT_MINI: Path = FIXTURES_DIR / "imgt_mini"
PHARMGKB_FIXTURE: Path = FIXTURES_DIR / "pharmgkb"
T1K_FIXTURE: Path = FIXTURES_DIR / "sample_t1k_genotype.tsv"
HLAHD_FIXTURE: Path = FIXTURES_DIR / "sample_final.result.txt"
ARCAS_FIXTURE: Path = FIXTURES_DIR / "sample.genotype.json"
OPTITYPE_FIXTURE: Path = FIXTURES_DIR / "sample_result.tsv"


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture()
def runner() -> CliRunner:
    """
    CliRunner that captures stderr separately on every supported Click
    release.

    Click <= 8.1 merges stderr into stdout unless ``mix_stderr=False`` is
    passed, and accessing ``result.stderr`` then raises ``ValueError``.
    Click 8.3 removed the argument and separates the streams by default.
    Requesting the flag and falling back keeps the suite green on both.
    """
    try:
        return CliRunner(mix_stderr=False)  # type: ignore[call-arg]
    except TypeError:
        return CliRunner()


@pytest.fixture()
def annotate_args(tmp_path: Path) -> List[str]:
    """
    CLI flags shared by the annotate tests.
    """
    return [
        "--offline",
        "--imgt-db-path", str(IMGT_MINI),
        "--pharmgkb-dir", str(PHARMGKB_FIXTURE),
        "--cache-dir", str(tmp_path / "cache"),
    ]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _read_tsv(path: Path):
    """
    Skip ``#`` metadata lines and return the header + data rows.
    """
    lines = path.read_text(encoding="utf-8").splitlines()
    data_lines = [ln for ln in lines if not ln.startswith("#")]
    reader = csv.reader(data_lines, delimiter="\t")
    rows = list(reader)
    return rows[0], rows[1:]


# ---------------------------------------------------------------------------
# Version + validate (no network required)
# ---------------------------------------------------------------------------


class TestVersionCommand:
    """
    Output of ``hlante version``.
    """

    def test_prints_version_and_dbs(self, runner: CliRunner, tmp_path: Path) -> None:
        """
        The HLAnte version and the fixture IPD-IMGT/HLA version must be
        listed.
        """
        result = runner.invoke(
            cli,
            [
                "version",
                "--imgt-db-path", str(IMGT_MINI),
                "--pharmgkb-dir", str(PHARMGKB_FIXTURE),
                "--cache-dir", str(tmp_path / "cache"),
            ],
        )
        assert result.exit_code == 0, result.output + result.stderr
        assert f"HLAnte v{__version__}" in result.output
        assert "IPD-IMGT/HLA 3.55.0" in result.output
        assert "PharmGKB" in result.output

    def test_no_imgt_installed(self, runner: CliRunner, tmp_path: Path) -> None:
        """
        When the IPD-IMGT/HLA directory is missing the status must read
        ``not installed``.
        """
        result = runner.invoke(
            cli,
            [
                "version",
                "--imgt-db-path", str(tmp_path / "missing"),
                "--pharmgkb-dir", str(tmp_path / "missing_pharm"),
                "--cache-dir", str(tmp_path / "cache"),
            ],
        )
        assert result.exit_code == 0
        assert "IPD-IMGT/HLA : not installed" in result.output
        assert "PharmGKB     : not installed" in result.output


class TestValidateCommand:
    """
    Basic behaviour of ``hlante validate``.
    """

    def test_valid_t1k(self, runner: CliRunner) -> None:
        """
        A valid T1K input must pass validation.
        """
        result = runner.invoke(
            cli, ["validate", "-i", str(T1K_FIXTURE), "-t", "t1k"]
        )
        assert result.exit_code == 0, result.output + result.stderr
        assert "Input is valid" in result.output
        assert "Samples: 1" in result.output
        assert "Locus calls: 4" in result.output

    def test_bad_file_exits_nonzero(
        self, runner: CliRunner, tmp_path: Path
    ) -> None:
        """
        A malformed T1K file must produce a non-zero exit code.
        """
        bad = tmp_path / "bad.tsv"
        bad.write_text("not a tsv header\njunk\tjunk\n")
        result = runner.invoke(cli, ["validate", "-i", str(bad), "-t", "t1k"])
        assert result.exit_code == 1
        assert "ERROR" in result.stderr

    def test_missing_input_exits(self, runner: CliRunner, tmp_path: Path) -> None:
        """
        Missing paths must be rejected by Click's own validation.
        """
        result = runner.invoke(
            cli,
            ["validate", "-i", str(tmp_path / "none.tsv"), "-t", "t1k"],
        )
        assert result.exit_code != 0

    def test_validate_directory(
        self, runner: CliRunner, tmp_path: Path
    ) -> None:
        """
        A directory input must discover per-tool matching files.
        """
        work = tmp_path / "many_t1k"
        work.mkdir()
        for idx in range(3):
            dest = work / f"sample{idx:02d}.tsv"
            shutil.copy(T1K_FIXTURE, dest)
        result = runner.invoke(
            cli, ["validate", "-i", str(work), "-t", "t1k"]
        )
        assert result.exit_code == 0, result.output + result.stderr
        assert "Files: 3" in result.output
        assert "Samples: 3" in result.output


# ---------------------------------------------------------------------------
# Annotate — end to end
# ---------------------------------------------------------------------------


class TestAnnotateEndToEnd:
    """
    End-to-end behaviour of ``hlante annotate`` (offline).
    """

    def test_t1k_produces_all_three_formats(
        self,
        runner: CliRunner,
        tmp_path: Path,
        annotate_args: List[str],
    ) -> None:
        """
        The T1K fixture must produce TSV + Markdown + JSON reports.
        """
        out_dir = tmp_path / "out"
        result = runner.invoke(
            cli,
            [
                "annotate",
                "-i", str(T1K_FIXTURE),
                "-t", "t1k",
                "-o", str(out_dir),
                "--format", "all",
                *annotate_args,
            ],
        )
        assert result.exit_code == 0, result.output + result.stderr

        tsv_path = out_dir / "hlante_report.tsv"
        md_path = out_dir / "hlante_report.md"
        json_path = out_dir / "hlante_report.json"
        for p in (tsv_path, md_path, json_path):
            assert p.is_file(), p

        # TSV — header and row count
        header, rows = _read_tsv(tsv_path)
        assert header[0] == "sample_id"
        assert "allele1" in header
        # The T1K fixture has 4 loci (A/B/C/DRB1).
        assert len(rows) == 4
        loci = {r[header.index("locus")] for r in rows}
        assert loci == {"HLA-A", "HLA-B", "HLA-C", "HLA-DRB1"}

        # JSON — metadata + sample structure
        payload = json.loads(json_path.read_text(encoding="utf-8"))
        assert payload["metadata"]["hlante_version"] == __version__
        assert payload["metadata"]["db_versions"].get("imgt") == "IPD-IMGT/HLA 3.55.0"
        assert len(payload["samples"]) == 1

        # Markdown — sample sections
        md_text = md_path.read_text(encoding="utf-8")
        assert "## Sample:" in md_text
        assert "### HLA Genotype" in md_text

    def test_format_tsv_only(
        self,
        runner: CliRunner,
        tmp_path: Path,
        annotate_args: List[str],
    ) -> None:
        """
        ``--format tsv`` must produce only the TSV output.
        """
        out_dir = tmp_path / "tsv_only"
        result = runner.invoke(
            cli,
            [
                "annotate",
                "-i", str(T1K_FIXTURE),
                "-t", "t1k",
                "-o", str(out_dir),
                "--format", "tsv",
                *annotate_args,
            ],
        )
        assert result.exit_code == 0, result.output + result.stderr
        assert (out_dir / "hlante_report.tsv").is_file()
        assert not (out_dir / "hlante_report.md").exists()
        assert not (out_dir / "hlante_report.json").exists()

    def test_resolution_filter_drops_low_res(
        self,
        runner: CliRunner,
        tmp_path: Path,
        annotate_args: List[str],
    ) -> None:
        """
        ``--resolution 6`` must drop two-field T1K records.
        """
        out_dir = tmp_path / "out"
        result = runner.invoke(
            cli,
            [
                "annotate",
                "-i", str(T1K_FIXTURE),
                "-t", "t1k",
                "-o", str(out_dir),
                "--format", "tsv",
                "--resolution", "6",
                *annotate_args,
            ],
        )
        assert result.exit_code == 0, result.output + result.stderr
        # T1K fixture has two-field calls → all dropped, TSV empty.
        header, rows = _read_tsv(out_dir / "hlante_report.tsv")
        assert rows == []

    def test_overwrite_guard(
        self,
        runner: CliRunner,
        tmp_path: Path,
        annotate_args: List[str],
    ) -> None:
        """
        The second run without ``--overwrite`` must fail; with it,
        succeed.
        """
        out_dir = tmp_path / "out"
        base_args = [
            "annotate", "-i", str(T1K_FIXTURE), "-t", "t1k",
            "-o", str(out_dir), "--format", "tsv",
            *annotate_args,
        ]
        r1 = runner.invoke(cli, base_args)
        assert r1.exit_code == 0, r1.output + r1.stderr

        r2 = runner.invoke(cli, base_args)
        assert r2.exit_code == 1
        assert "already exists" in r2.stderr

        r3 = runner.invoke(cli, base_args + ["--overwrite"])
        assert r3.exit_code == 0, r3.output + r3.stderr

    def test_hlahd_tool_works(
        self,
        runner: CliRunner,
        tmp_path: Path,
        annotate_args: List[str],
    ) -> None:
        """
        The HLA-HD format must work end-to-end.
        """
        out_dir = tmp_path / "hlahd_out"
        result = runner.invoke(
            cli,
            [
                "annotate",
                "-i", str(HLAHD_FIXTURE),
                "-t", "hla-hd",  # Alias
                "-o", str(out_dir),
                "--format", "tsv",
                *annotate_args,
            ],
        )
        assert result.exit_code == 0, result.output + result.stderr
        header, rows = _read_tsv(out_dir / "hlante_report.tsv")
        # HLA-HD fixture: 5 typed loci (DPB1 is skipped).
        assert len(rows) == 5

    def test_optitype_tool_works(
        self,
        runner: CliRunner,
        tmp_path: Path,
        annotate_args: List[str],
    ) -> None:
        """
        The OptiType (Class I) format must work end-to-end.
        """
        out_dir = tmp_path / "opt_out"
        result = runner.invoke(
            cli,
            [
                "annotate",
                "-i", str(OPTITYPE_FIXTURE),
                "-t", "optitype",
                "-o", str(out_dir),
                "--format", "tsv",
                *annotate_args,
            ],
        )
        assert result.exit_code == 0, result.output + result.stderr
        header, rows = _read_tsv(out_dir / "hlante_report.tsv")
        assert len(rows) == 3
        loci = {r[header.index("locus")] for r in rows}
        assert loci == {"HLA-A", "HLA-B", "HLA-C"}


# ---------------------------------------------------------------------------
# Offline behaviour
# ---------------------------------------------------------------------------


class TestOfflineBehaviour:
    """
    ``--offline`` must not make any HTTP requests.
    """

    def test_no_http_calls_in_offline(
        self,
        runner: CliRunner,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
        annotate_args: List[str],
    ) -> None:
        """
        When ``urllib.urlopen`` is replaced with a raiser and the
        pipeline still succeeds, no HTTP request was issued.
        """
        import urllib.request

        def _boom(*args, **kwargs):  # pragma: no cover — should not fire
            raise AssertionError("Offline mode must not issue HTTP calls.")

        monkeypatch.setattr(urllib.request, "urlopen", _boom)

        out_dir = tmp_path / "offline_out"
        result = runner.invoke(
            cli,
            [
                "annotate",
                "-i", str(T1K_FIXTURE),
                "-t", "t1k",
                "-o", str(out_dir),
                "--format", "tsv",
                *annotate_args,
            ],
        )
        assert result.exit_code == 0, result.output + result.stderr
        assert (out_dir / "hlante_report.tsv").is_file()


# ---------------------------------------------------------------------------
# 10-sample batch
# ---------------------------------------------------------------------------


class TestBatchCohort:
    """
    Multi-sample integration via a directory input.
    """

    def test_ten_sample_batch(
        self,
        runner: CliRunner,
        tmp_path: Path,
        annotate_args: List[str],
    ) -> None:
        """
        Ten T1K files must be combined into a single report.
        """
        cohort_dir = tmp_path / "cohort"
        cohort_dir.mkdir()
        for idx in range(10):
            dest = cohort_dir / f"sample_{idx:02d}.tsv"
            shutil.copy(T1K_FIXTURE, dest)

        out_dir = tmp_path / "cohort_out"
        result = runner.invoke(
            cli,
            [
                "annotate",
                "-i", str(cohort_dir),
                "-t", "t1k",
                "-o", str(out_dir),
                "--format", "all",
                *annotate_args,
            ],
        )
        assert result.exit_code == 0, result.output + result.stderr

        # TSV — 10 samples × 4 loci = 40 rows
        header, rows = _read_tsv(out_dir / "hlante_report.tsv")
        assert len(rows) == 40

        idx_sample = header.index("sample_id")
        samples = {r[idx_sample] for r in rows}
        assert len(samples) == 10
        assert all(sid.startswith("sample_") for sid in samples)

        # JSON — same sample count
        payload = json.loads(
            (out_dir / "hlante_report.json").read_text(encoding="utf-8")
        )
        assert len(payload["samples"]) == 10
        for sample_entry in payload["samples"]:
            assert len(sample_entry["loci"]) == 4

        # Markdown — ten sample headings
        md_text = (out_dir / "hlante_report.md").read_text(encoding="utf-8")
        assert md_text.count("## Sample:") == 10


# ---------------------------------------------------------------------------
# db-update
# ---------------------------------------------------------------------------


class TestDBUpdateCommand:
    """
    Basic behaviour of ``hlante db-update``.
    """

    def test_gwas_update_calls_download(
        self,
        runner: CliRunner,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        """
        ``--db gwas`` must invoke ``GWASClient.update``.
        """
        import io
        import zipfile

        # Minimal valid zip payload
        buf = io.BytesIO()
        with zipfile.ZipFile(buf, "w") as zf:
            zf.writestr(
                "gwas-catalog-download-associations-alt-full.tsv",
                "STRONGEST SNP-RISK ALLELE\tDISEASE/TRAIT\tMAPPED_TRAIT\n",
            )
        good_zip = buf.getvalue()

        captured: dict = {}

        def fake_fetcher(url: str) -> bytes:
            captured["url"] = url
            return good_zip

        # Patch the GWASClient constructor used inside the CLI.
        real_client_cls = __import__("hlante.cli", fromlist=["GWASClient"]).GWASClient

        def patched_ctor(local_dir):
            return real_client_cls(local_dir=local_dir, fetcher=fake_fetcher)

        monkeypatch.setattr("hlante.cli.GWASClient", patched_ctor)

        gwas_dir = tmp_path / "gwas"
        result = runner.invoke(
            cli,
            [
                "db-update",
                "--db", "gwas",
                "--gwas-cache-dir", str(gwas_dir),
            ],
        )
        assert result.exit_code == 0, result.output + result.stderr
        assert "GWAS updated" in result.output
        assert "ftp.ebi.ac.uk" in captured.get("url", "")

    def test_imgt_uses_download_function(
        self,
        runner: CliRunner,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        """
        ``--db imgt`` must call ``download_imgt_db`` inside the CLI.
        """
        called = {}

        def fake_download(target_dir, *, include_groups=True, force=False, ref="Latest"):
            called["target_dir"] = Path(target_dir)
            called["force"] = force
            called["ref"] = ref
            Path(target_dir).mkdir(parents=True, exist_ok=True)
            (Path(target_dir) / ALLELE_LIST_FILENAME).write_text(
                "# version: test\nHLA00001,A*01:01\n", encoding="utf-8"
            )
            return Path(target_dir)

        monkeypatch.setattr("hlante.cli.download_imgt_db", fake_download)

        target = tmp_path / "imgt"
        result = runner.invoke(
            cli,
            [
                "db-update",
                "--db", "imgt",
                "--force",
                "--imgt-dir", str(target),
            ],
        )
        assert result.exit_code == 0, result.output + result.stderr
        assert called["target_dir"] == target
        assert called["force"] is True
        assert "IPD-IMGT/HLA updated" in result.output

    def test_imgt_failure_propagates_exit_code(
        self,
        runner: CliRunner,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        """
        A download failure must produce exit_code=1 with an ERROR line.
        """

        def broken_download(*_args, **_kwargs):
            raise RuntimeError("no network")

        monkeypatch.setattr("hlante.cli.download_imgt_db", broken_download)

        result = runner.invoke(
            cli,
            [
                "db-update",
                "--db", "imgt",
                "--imgt-dir", str(tmp_path / "imgt"),
            ],
        )
        assert result.exit_code == 1
        assert "IPD-IMGT/HLA update failed" in result.stderr


# ---------------------------------------------------------------------------
# Main() entry point
# ---------------------------------------------------------------------------


class TestMainEntryPoint:
    """
    :func:`hlante.cli.main` must not raise :class:`SystemExit` — it
    returns an int.
    """

    def test_main_returns_zero_on_help(self) -> None:
        from hlante.cli import main

        # --help triggers click.ClickException → returns 0.
        code = main(["--help"])
        assert code == 0

    def test_main_returns_nonzero_on_missing_required(self) -> None:
        from hlante.cli import main

        # Missing required -i/-t → Click UsageError (exit 2).
        code = main(["annotate"])
        assert code != 0
