"""
hlante.cli
=========

Click-based command-line interface.

Commands
--------
- ``annotate`` — parse HLA typing tool output and produce a clinical
  annotation report.
- ``db-update`` — refresh local copies of IPD-IMGT/HLA, PharmGKB,
  GWAS Catalog, and AFND.
- ``validate`` — sanity-check the format of an input file or directory.
- ``version`` — print HLAnte and installed database versions.

Entry point: :func:`main` (registered in ``pyproject.toml`` as
``hlante = hlante.cli:main``).
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import click

from hlante import __version__
from hlante.annotator import (
    AnnotatedHLA,
    AnnotatorConfig,
    annotate_genotype,
)
from hlante.types import InputSource
from hlante.db.afnd import (
    AFNDClient,
    AFNDDatabaseError,
    DEFAULT_LOCAL_DIR as AFND_DEFAULT_DIR,
    DEFAULT_POPULATION_GROUP as AFND_DEFAULT_POPULATION,
)
from hlante.db.gwas import (
    DEFAULT_LOCAL_DIR as GWAS_DEFAULT_DIR,
    GWAS_HLA_SUBSET_FILENAME,
    GWAS_TSV_FILENAME_DEFAULT,
    GWASClient,
    GWASDatabaseError,
    GWASDownloadError,
)
from hlante.db.imgt import (
    DEFAULT_LOCAL_DIR as IMGT_DEFAULT_DIR,
    IMGT_DEFAULT_REF,
    download_imgt_db,
)
from hlante.db.pharmgkb import (
    DEFAULT_LOCAL_DIR as PHARMGKB_DEFAULT_DIR,
    PharmGKBClient,
    PharmGKBDatabaseError,
    PharmGKBDownloadError,
)
from hlante.normalizer import (
    IMGTDatabaseMissingError,
    batch_normalize,
    load_imgt_db,
)
from hlante.parser import (
    HLAGenotype,
    HLAnteParseError,
    SUPPORTED_TOOLS,
    UnsupportedToolError,
    parse_hla_output,
)
from hlante.reporter import (
    OutputFileExistsError,
    ReportContext,
    generate_all,
    generate_json,
    generate_markdown_report,
    generate_tsv,
)


logger: logging.Logger = logging.getLogger("hlante.cli")


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

DEFAULT_OUTPUT_DIR: Path = Path("./hlante_output")
DEFAULT_CACHE_ROOT: Path = Path.home() / ".hlante" / "cache"
DEFAULT_PREFIX: str = "hlante_report"

TOOL_CHOICES: List[str] = sorted(SUPPORTED_TOOLS) + ["arcas-hla", "hla-hd"]

# Tool name → glob patterns searched when --input is a directory
_TOOL_GLOBS: Dict[str, Tuple[str, ...]] = {
    "arcashla": ("*.json",),
    "t1k": ("*_genotype.tsv", "*.tsv"),
    "hlahd": ("*final.result.txt", "*.txt"),
    "optitype": ("*_result.tsv", "*.tsv"),
}

_FORMAT_CHOICES: List[str] = ["tsv", "json", "markdown", "all"]
#: Accepted ``--resolution`` values. HLAnte reports resolution as a field
#: count (1..4) from v0.2.0; releases up to v0.1.0 used a digit scale
#: (2/4/6/8). The legacy values stay accepted so existing scripts keep
#: running, but 6 and 8 are not valid field counts and are translated with a
#: warning.
_RESOLUTION_CHOICES: List[str] = ["1", "2", "3", "4", "6", "8"]

#: Digit-scale values that cannot be field counts, and their field equivalent.
_LEGACY_RESOLUTION: Dict[str, int] = {"6": 3, "8": 4}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _configure_logging(level_name: str) -> None:
    """
    Configure the root logger to the requested level.
    """
    level = getattr(logging, level_name.upper(), logging.INFO)
    logging.basicConfig(
        level=level,
        format="%(asctime)s %(levelname)-7s %(name)s: %(message)s",
        stream=sys.stderr,
        force=True,
    )


def _normalize_tool(tool: str) -> str:
    """
    Map user-visible tool aliases (``arcas-hla``, ``hla-hd``) to the
    canonical keys used by the parser module.
    """
    key = tool.strip().lower()
    return {
        "arcas": "arcashla",
        "arcas-hla": "arcashla",
        "arcas_hla": "arcashla",
        "arcashla": "arcashla",
        "hla-hd": "hlahd",
        "hla_hd": "hlahd",
        "hlahd": "hlahd",
        "t1k": "t1k",
        "optitype": "optitype",
        "opti-type": "optitype",
    }.get(key, key)


def _collect_input_files(input_path: Path, tool: str) -> List[Path]:
    """
    Return a single file list when ``input_path`` is a file, or the
    tool-appropriate glob matches when it is a directory.
    """
    if input_path.is_file():
        return [input_path]
    if not input_path.is_dir():
        raise click.BadParameter(f"Input path is neither a file nor a directory: {input_path}")
    patterns = _TOOL_GLOBS.get(tool, ("*",))
    seen: Dict[Path, None] = {}
    for pattern in patterns:
        for match in sorted(input_path.glob(pattern)):
            if match.is_file():
                seen.setdefault(match, None)
    files = list(seen)
    if not files:
        raise click.BadParameter(
            f"No files matching {tool} patterns found in {input_path} (patterns: {patterns})."
        )
    return files


def _parse_all(
    files: Sequence[Path], tool: str, *, strict: bool = False
) -> Tuple[List[HLAGenotype], List[Tuple[Path, Exception]]]:
    """
    Parse all files with the tool parser and aggregate the results.

    By default a file that fails to parse is skipped and recorded in the
    returned failure list, so one malformed sample does not discard an
    entire cohort. When ``strict`` is ``True`` the first
    parse error is re-raised to the caller.
    """
    all_genotypes: List[HLAGenotype] = []
    failures: List[Tuple[Path, Exception]] = []
    for path in files:
        try:
            all_genotypes.extend(parse_hla_output(path, tool))
        except (HLAnteParseError, UnsupportedToolError, FileNotFoundError, OSError) as exc:
            if strict:
                raise
            logger.warning("Skipping unparseable file %s: %s", path, exc)
            failures.append((path, exc))
    return all_genotypes, failures


def _echo_success(message: str) -> None:
    click.echo(click.style(f"✓ {message}", fg="green"))


def _echo_info(message: str) -> None:
    click.echo(click.style(f"  {message}", fg="cyan"))


def _echo_warn(message: str) -> None:
    click.echo(click.style(f"⚠ {message}", fg="yellow"), err=True)


def _echo_error(message: str) -> None:
    click.echo(click.style(f"ERROR: {message}", fg="red", bold=True), err=True)


def _db_versions(
    imgt_db_path: Optional[Path],
    pharmgkb_dir: Optional[Path],
    gwas_dir: Optional[Path] = None,
) -> Dict[str, str]:
    """
    Detect local database versions (when available).

    Notes
    -----
    exposes the cache mtime for the GWAS Catalog bulk TSV as
    ``gwas_cache_date`` so every report header records *when* the
    underlying data was downloaded. This is important because the
    GWAS Catalog updates continuously; without a snapshot date the
    report is not independently reproducible.
    """
    from datetime import datetime, timezone

    versions: Dict[str, str] = {}
    try:
        imgt = load_imgt_db(imgt_db_path)
        if imgt.get("version"):
            versions["imgt"] = str(imgt["version"])
    except IMGTDatabaseMissingError:
        pass
    pharm_root = pharmgkb_dir or PHARMGKB_DEFAULT_DIR
    if Path(pharm_root).exists():
        versions["pharmgkb"] = "local"
    gwas_root = gwas_dir or GWAS_DEFAULT_DIR
    gwas_root = Path(gwas_root)
    if gwas_root.exists():
        mtimes = [p.stat().st_mtime for p in gwas_root.rglob("*") if p.is_file()]
        if mtimes:
            ts = datetime.fromtimestamp(max(mtimes), tz=timezone.utc)
            versions["gwas_cache_date"] = ts.date().isoformat()
    return versions


# ---------------------------------------------------------------------------
# Root group
# ---------------------------------------------------------------------------


@click.group(
    name="hlante",
    help="HLA genotype → disease/drug research annotation aid (research use only).",
    context_settings={"help_option_names": ["-h", "--help"]},
)
@click.version_option(version=__version__, prog_name="hlante")
@click.option(
    "--log-level",
    default="INFO",
    type=click.Choice(["DEBUG", "INFO", "WARNING", "ERROR"], case_sensitive=False),
    show_default=True,
    help="Logging verbosity.",
)
@click.pass_context
def cli(ctx: click.Context, log_level: str) -> None:
    """
    HLAnte root command group.
    """
    _configure_logging(log_level)
    ctx.ensure_object(dict)
    ctx.obj["log_level"] = log_level


# ---------------------------------------------------------------------------
# Annotate
# ---------------------------------------------------------------------------


@cli.command("annotate")
@click.option(
    "-i",
    "--input",
    "input_path",
    required=True,
    type=click.Path(exists=True, path_type=Path),
    help="HLA typing tool output file or directory.",
)
@click.option(
    "-t",
    "--tool",
    required=True,
    type=click.Choice(TOOL_CHOICES, case_sensitive=False),
    help="HLA typing tool that produced the input.",
)
@click.option(
    "-o",
    "--output-dir",
    default=DEFAULT_OUTPUT_DIR,
    show_default=True,
    type=click.Path(file_okay=False, path_type=Path),
    help="Report output directory.",
)
@click.option(
    "--format",
    "output_format",
    default="all",
    show_default=True,
    type=click.Choice(_FORMAT_CHOICES, case_sensitive=False),
    help="Report format to produce.",
)
@click.option(
    "--offline",
    is_flag=True,
    help="Do not make HTTP requests; use only local caches / bulk dumps.",
)
@click.option(
    "--threads",
    default=4,
    show_default=True,
    type=click.IntRange(min=1),
    help="Thread count for normalization.",
)
@click.option(
    "--resolution",
    default=None,
    type=click.Choice(_RESOLUTION_CHOICES, case_sensitive=False),
    help="Minimum accepted resolution, in fields (1-4).",
)
@click.option(
    "--imgt-db-path",
    default=None,
    type=click.Path(path_type=Path),
    help="IPD-IMGT/HLA local directory (default: ~/.hlante/imgt_hla).",
)
@click.option(
    "--pharmgkb-dir",
    default=None,
    type=click.Path(path_type=Path),
    help="PharmGKB bulk dump directory (default: ~/.hlante/pharmgkb).",
)
@click.option(
    "--cache-dir",
    default=None,
    type=click.Path(path_type=Path),
    help="Cache root for GWAS (default: ~/.hlante/cache).",
)
@click.option(
    "--prefix",
    default=DEFAULT_PREFIX,
    show_default=True,
    help="Output file name prefix.",
)
@click.option(
    "--overwrite",
    is_flag=True,
    help="Overwrite existing output files.",
)
@click.option(
    "--no-gwas",
    is_flag=True,
    help="Skip GWAS Catalog queries.",
)
@click.option(
    "--no-pharmgkb",
    is_flag=True,
    help="Skip PharmGKB queries.",
)
@click.option(
    "--no-afnd",
    is_flag=True,
    help="Skip AFND allele-frequency lookups (affects the input-quality score).",
)
@click.option(
    "-p",
    "--population",
    "population_group",
    default=AFND_DEFAULT_POPULATION,
    show_default=True,
    help=(
        "Population group for allele frequency lookup from AFND.\n"
        "Available groups:\n"
        "  EUR    — European (Caucasian, Western/Northern/Southern/Eastern European)\n"
        "  AFR    — African (Sub-Saharan, West African, East African)\n"
        "  EAS    — East/Southeast Asian (Japanese, Chinese, Korean; alias: ASN)\n"
        "  SAS    — South Asian (Indian, Pakistani, Bangladeshi)\n"
        "  MID    — Middle Eastern / North African / West Asian\n"
        "  AMR    — American / Hispanic / Latino / Native American\n"
        "  OCE    — Oceanian / Pacific Islander\n"
        "  global — aggregate across all populations (weighted by sample size)"
    ),
)
@click.option(
    "--afnd-dir",
    default=None,
    type=click.Path(path_type=Path),
    help="AFND local data directory (default: ~/.hlante/afnd).",
)
@click.option(
    "--curated-tsv",
    "curated_tsv_path",
    default=None,
    type=click.Path(exists=True, path_type=Path),
    help=(
        "Path to a user-supplied curated HLA–disease TSV file. "
        "When provided, replaces the built-in curated table. "
        "The file must follow HLAnte's curated-table schema "
        "(columns: Allele, Disease, ClinicalSignificance, Evidence, "
        "Population, OR, PMID, Citation)."
    ),
)
@click.option(
    "--input-source",
    "input_source",
    default="typing_tool",
    show_default=True,
    type=click.Choice(
        ["typing_tool", "validated", "simulated", "unknown"],
        case_sensitive=False,
    ),
    help=(
        "Provenance of input HLA allele calls.\n\n"
        "typing_tool  (default) — output from arcasHLA / T1K / HLA-HD / OptiType. "
        "Resolution and ambiguity penalties apply normally.\n\n"
        "validated — lab-validated alleles (e.g. 1000 Genomes Sanger types, "
        "IHIW reference panels). Two-field calls are not penalised for ambiguity; "
        "only the resolution penalty applies.\n\n"
        "simulated — synthetic alleles for testing. Same penalties as typing_tool.\n\n"
        "unknown — treated as typing_tool with a warning."
    ),
)
@click.option(
    "--strict",
    is_flag=True,
    default=False,
    help=(
        "Abort the entire run if ANY input file fails to parse. By default, "
        "unparseable files are skipped with a warning and the remaining "
        "samples are still annotated (the run fails only if NO file parses)."
    ),
)
def annotate_cmd(
    input_path: Path,
    tool: str,
    output_dir: Path,
    output_format: str,
    offline: bool,
    threads: int,
    resolution: Optional[str],
    imgt_db_path: Optional[Path],
    pharmgkb_dir: Optional[Path],
    cache_dir: Optional[Path],
    prefix: str,
    overwrite: bool,
    no_gwas: bool,
    no_pharmgkb: bool,
    no_afnd: bool,
    population_group: str,
    afnd_dir: Optional[Path],
    curated_tsv_path: Optional[Path],
    input_source: str,
    strict: bool,
) -> None:
    """
    Run the full annotation pipeline (parse → normalize → annotate → report).
    """
    tool_key = _normalize_tool(tool)

    # 1. Collect input files
    try:
        files = _collect_input_files(input_path, tool_key)
    except click.BadParameter as exc:
        _echo_error(str(exc))
        sys.exit(1)

    _echo_info(f"Tool: {tool_key} | File count: {len(files)}")

    # 2. Parse
    try:
        genotypes, parse_failures = _parse_all(files, tool_key, strict=strict)
    except (HLAnteParseError, UnsupportedToolError) as exc:
        _echo_error(f"Parse error: {exc}")
        sys.exit(1)
    except FileNotFoundError as exc:
        _echo_error(str(exc))
        sys.exit(1)

    if parse_failures:
        _echo_warn(f"Skipped {len(parse_failures)} unparseable file(s):")
        for fpath, ferr in parse_failures:
            _echo_warn(f"  - {fpath}: {ferr}")
    if not genotypes:
        _echo_error("No parseable samples found in input; nothing to annotate.")
        sys.exit(1)

    _echo_success(
        f"Parse complete: {len(genotypes)} locus call(s) "
        f"({len({g.sample_id for g in genotypes})} sample(s))"
    )

    # 3. Normalize (load IPD-IMGT/HLA)
    try:
        imgt_db = load_imgt_db(imgt_db_path)
    except IMGTDatabaseMissingError as exc:
        _echo_error(f"IPD-IMGT/HLA could not be loaded: {exc}")
        sys.exit(1)

    try:
        normalized = batch_normalize(genotypes, imgt_db=imgt_db, max_workers=threads)
    except Exception as exc:  # noqa: BLE001
        _echo_error(f"Normalization failed: {exc}")
        sys.exit(1)

    # Resolution filter
    if resolution is not None:
        if resolution in _LEGACY_RESOLUTION:
            min_res = _LEGACY_RESOLUTION[resolution]
            _echo_warn(
                f"--resolution {resolution} is a digit-scale value from v0.1.0 and is no "
                f"longer a valid field count; interpreting it as {min_res} fields. "
                f"Use --resolution {min_res}."
            )
        else:
            min_res = int(resolution)
            if resolution in {"2", "4"}:
                _echo_warn(
                    f"--resolution now counts fields, not digits: {min_res} means "
                    f"{min_res} field(s). Before v0.2.0 this value meant "
                    f"{min_res // 2} field(s)."
                )
        filtered = [n for n in normalized if n.resolution_level >= min_res]
        dropped = len(normalized) - len(filtered)
        normalized = filtered
        if dropped:
            noun = "field" if min_res == 1 else "fields"
            _echo_warn(
                f"Resolution filter: dropped {dropped} allele(s) below {min_res} {noun}."
            )

    _echo_success(f"Normalization complete: {len(normalized)} allele(s)")

    if not normalized:
        _echo_warn("No annotatable alleles remain; output will be empty.")

    # 4. Annotate
    src = InputSource(input_source.lower())
    if src == InputSource.VALIDATED:
        logger.info(
            "Input source: validated. Ambiguity penalty suppressed for "
            "two-field alleles. Resolution penalties unchanged."
        )
    if curated_tsv_path is not None:
        _echo_info(f"Custom curated table: {curated_tsv_path}")
    config = AnnotatorConfig(
        offline=offline,
        cache_root=cache_dir or DEFAULT_CACHE_ROOT,
        pharmgkb_local_dir=pharmgkb_dir,
        afnd_local_dir=afnd_dir,
        imgt_db_path=imgt_db_path,
        enable_gwas=not no_gwas,
        enable_pharmgkb=not no_pharmgkb,
        enable_afnd=not no_afnd,
        population_group=population_group,
        input_source=src,
        curated_tsv_path=curated_tsv_path,
    )
    annotated = annotate_genotype(normalized, config)
    _echo_success(f"Annotation complete: {len(annotated)} record(s) processed.")

    # 5. Report
    ctx = ReportContext(
        db_versions=_db_versions(
            imgt_db_path,
            pharmgkb_dir,
            gwas_dir=getattr(config, "gwas_local_dir", None),
        ),
        input_source=src.value,
        cli_invocation=" ".join(["hlante", *sys.argv[1:]]),
    )

    try:
        output_paths = _write_reports(
            annotated,
            output_dir=output_dir,
            output_format=output_format.lower(),
            prefix=prefix,
            overwrite=overwrite,
            context=ctx,
        )
    except OutputFileExistsError as exc:
        _echo_error(str(exc))
        sys.exit(1)

    _echo_success("Report(s) written:")
    for fmt, path in output_paths.items():
        _echo_info(f"{fmt}: {path}")


def _write_reports(
    annotated: List[AnnotatedHLA],
    *,
    output_dir: Path,
    output_format: str,
    prefix: str,
    overwrite: bool,
    context: ReportContext,
) -> Dict[str, Path]:
    """
    Dispatch to the appropriate reporter function(s).
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if output_format == "all":
        return generate_all(
            annotated,
            output_dir,
            prefix=prefix,
            overwrite=overwrite,
            context=context,
        )
    if output_format == "tsv":
        return {
            "tsv": generate_tsv(
                annotated,
                output_dir / f"{prefix}.tsv",
                overwrite=overwrite,
                context=context,
            )
        }
    if output_format == "markdown":
        return {
            "markdown": generate_markdown_report(
                annotated,
                output_dir / f"{prefix}.md",
                overwrite=overwrite,
                context=context,
            )
        }
    if output_format == "json":
        return {
            "json": generate_json(
                annotated,
                output_dir / f"{prefix}.json",
                overwrite=overwrite,
                context=context,
            )
        }
    raise click.BadParameter(f"Unknown format: {output_format}")


# ---------------------------------------------------------------------------
# Validate
# ---------------------------------------------------------------------------


@cli.command("validate")
@click.option(
    "-i",
    "--input",
    "input_path",
    required=True,
    type=click.Path(exists=True, path_type=Path),
    help="File or directory to validate.",
)
@click.option(
    "-t",
    "--tool",
    required=True,
    type=click.Choice(TOOL_CHOICES, case_sensitive=False),
    help="Input tool.",
)
def validate_cmd(input_path: Path, tool: str) -> None:
    """
    Verify that the input file or directory is parseable.
    """
    tool_key = _normalize_tool(tool)
    try:
        files = _collect_input_files(input_path, tool_key)
    except click.BadParameter as exc:
        _echo_error(str(exc))
        sys.exit(1)

    try:
        # Validate is a strict pre-flight check: any unparseable file fails.
        genotypes, _ = _parse_all(files, tool_key, strict=True)
    except (HLAnteParseError, UnsupportedToolError, FileNotFoundError) as exc:
        _echo_error(f"Invalid input: {exc}")
        sys.exit(1)

    sample_ids = {g.sample_id for g in genotypes}
    loci = {g.locus for g in genotypes}

    _echo_success("Input is valid")
    _echo_info(f"Tool: {tool_key}")
    _echo_info(f"Files: {len(files)}")
    _echo_info(f"Samples: {len(sample_ids)}")
    _echo_info(f"Locus calls: {len(genotypes)}")
    _echo_info(f"Unique loci: {', '.join(sorted(loci))}")


# ---------------------------------------------------------------------------
# db-update
# ---------------------------------------------------------------------------


@cli.command("db-update")
@click.option(
    "--db",
    "which",
    default="all",
    show_default=True,
    type=click.Choice(["imgt", "gwas", "pharmgkb", "afnd", "all"], case_sensitive=False),
    help="Target database.",
)
@click.option(
    "--force",
    is_flag=True,
    help="Ignore TTL; re-download / clear cache.",
)
@click.option(
    "--imgt-dir",
    default=None,
    type=click.Path(path_type=Path),
    help="IPD-IMGT/HLA directory (default: ~/.hlante/imgt_hla).",
)
@click.option(
    "--imgt-ref",
    default=IMGT_DEFAULT_REF,
    show_default=True,
    help=(
        "Git ref (release branch or tag) of the ANHIG/IMGTHLA mirror to "
        "download. 'Latest' tracks the moving release branch; pass an explicit "
        "release tag to pin a reproducible, checksum-recorded snapshot."
    ),
)
@click.option(
    "--pharmgkb-dir",
    default=None,
    type=click.Path(path_type=Path),
    help="PharmGKB directory (default: ~/.hlante/pharmgkb).",
)
@click.option(
    "--gwas-cache-dir",
    default=None,
    type=click.Path(path_type=Path),
    help="GWAS bulk dump directory (default: ~/.hlante/gwas).",
)
@click.option(
    "--afnd-dir",
    default=None,
    type=click.Path(path_type=Path),
    help="AFND directory (default: ~/.hlante/afnd).",
)
@click.option(
    "--afnd-url",
    default=None,
    help=(
        "Override the default Slowikowski AFND mirror URL "
        "(default: https://raw.githubusercontent.com/slowkow/allelefrequencies/main/afnd.tsv). "
        "Must point to either the Slowikowski 7-column schema or HLAnte's 5-column TSV."
    ),
)
def db_update_cmd(
    which: str,
    force: bool,
    imgt_dir: Optional[Path],
    imgt_ref: str,
    pharmgkb_dir: Optional[Path],
    gwas_cache_dir: Optional[Path],
    afnd_dir: Optional[Path],
    afnd_url: Optional[str],
) -> None:
    """
    Refresh local database copies.

    \b
    AFND (--db afnd):
      Downloads the machine-readable AFND mirror maintained by Slowikowski
      (2024) from GitHub (~6 MB; >3,000 population studies across 8 HLA loci).
      The Slowikowski 7-column schema is automatically transformed to HLAnte's
      internal 5-column format and saved to ~/.hlante/afnd/afnd_frequencies.tsv.
      Until this command is run, HLAnte uses a compact built-in fallback table
      (7 loci × 5 populations) derived from Gonzalez-Galarza et al. (2020).
    """
    which = which.lower()
    had_error = False

    if which in ("imgt", "all"):
        try:
            path = download_imgt_db(imgt_dir or IMGT_DEFAULT_DIR, force=force, ref=imgt_ref)
            _echo_success(f"IPD-IMGT/HLA updated → {path}")
        except Exception as exc:  # noqa: BLE001 — report network errors
            had_error = True
            _echo_error(f"IPD-IMGT/HLA update failed: {exc}")

    if which in ("pharmgkb", "all"):
        client = PharmGKBClient(local_dir=pharmgkb_dir or PHARMGKB_DEFAULT_DIR)
        try:
            path = client.update()
            _echo_success(f"PharmGKB updated → {path}")
        except PharmGKBDownloadError as exc:
            had_error = True
            _echo_error(f"PharmGKB update failed: {exc}")
        except PharmGKBDatabaseError as exc:
            had_error = True
            _echo_error(f"PharmGKB error: {exc}")

    if which in ("gwas", "all"):
        target = Path(gwas_cache_dir or GWAS_DEFAULT_DIR)
        gwas_client = GWASClient(local_dir=target)
        try:
            _echo_info(f"Downloading GWAS Catalog bulk dump → {target}")
            path = gwas_client.update()
            _echo_success(f"GWAS updated → {path}")
        except GWASDownloadError as exc:
            had_error = True
            _echo_error(f"GWAS download failed: {exc}")
        except GWASDatabaseError as exc:
            had_error = True
            _echo_error(f"GWAS error: {exc}")

    if which in ("afnd", "all"):
        target = Path(afnd_dir or AFND_DEFAULT_DIR)
        afnd_client = AFNDClient(local_dir=target)
        try:
            path = afnd_client.update(source_url=afnd_url)
            _echo_success(f"AFND updated → {path}")
        except AFNDDatabaseError as exc:
            # Network or transform failure; built-in fallback still covers basic use.
            if which == "afnd":
                had_error = True
                _echo_error(str(exc))
            else:
                _echo_warn(str(exc))

    if had_error:
        sys.exit(1)


# ---------------------------------------------------------------------------
# Version
# ---------------------------------------------------------------------------


@cli.command("version")
@click.option(
    "--imgt-db-path",
    default=None,
    type=click.Path(path_type=Path),
    help="IPD-IMGT/HLA directory (default: ~/.hlante/imgt_hla).",
)
@click.option(
    "--pharmgkb-dir",
    default=None,
    type=click.Path(path_type=Path),
    help="PharmGKB directory (default: ~/.hlante/pharmgkb).",
)
@click.option(
    "--gwas-cache-dir",
    default=None,
    type=click.Path(path_type=Path),
    help="GWAS bulk dump directory (default: ~/.hlante/gwas).",
)
def version_cmd(
    imgt_db_path: Optional[Path],
    pharmgkb_dir: Optional[Path],
    gwas_cache_dir: Optional[Path],
) -> None:
    """
    Print HLAnte and local database versions.
    """
    click.echo(f"HLAnte v{__version__}")

    # IPD-IMGT/HLA
    try:
        imgt = load_imgt_db(imgt_db_path)
        version = imgt.get("version") or "unknown"
        click.echo(f"  IPD-IMGT/HLA : {version}  ({imgt.get('path')})")
        if imgt.get("is_stale"):
            _echo_warn("  IPD-IMGT/HLA copy is >6 months old; db-update is recommended.")
    except IMGTDatabaseMissingError:
        click.echo("  IPD-IMGT/HLA : not installed")

    # PharmGKB
    pharm_dir = Path(pharmgkb_dir or PHARMGKB_DEFAULT_DIR)
    if pharm_dir.exists():
        tsv_candidates = list(pharm_dir.glob("clinical_ann*.tsv")) + list(
            pharm_dir.glob("clinical_annotations*.tsv")
        )
        status = "installed" if tsv_candidates else "partial"
        click.echo(f"  PharmGKB     : {status}  ({pharm_dir})")
    else:
        click.echo("  PharmGKB     : not installed")

    # GWAS bulk dump
    gwas_dir = Path(gwas_cache_dir or GWAS_DEFAULT_DIR)
    if (gwas_dir / GWAS_HLA_SUBSET_FILENAME).is_file() or (
        gwas_dir / GWAS_TSV_FILENAME_DEFAULT
    ).is_file():
        click.echo(f"  GWAS Catalog : installed  ({gwas_dir})")
    else:
        click.echo("  GWAS Catalog : not installed")

    # AFND
    afnd_dir_local = Path(AFND_DEFAULT_DIR)
    if afnd_dir_local.exists() and any(afnd_dir_local.glob("*.tsv")):
        click.echo(f"  AFND         : installed  ({afnd_dir_local})")
    else:
        click.echo("  AFND         : not installed")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def main(argv: Optional[Sequence[str]] = None) -> int:
    """
    Package entry point.

    Parameters
    ----------
    argv : sequence of str, optional
        Command-line arguments; ``sys.argv`` is used when ``None``.

    Returns
    -------
    int
        Process exit code.
    """
    try:
        cli.main(args=list(argv) if argv is not None else None, standalone_mode=False)
    except click.ClickException as exc:
        exc.show()
        return exc.exit_code
    except click.Abort:
        click.echo("Aborted.", err=True)
        return 130
    except SystemExit as exc:
        code = exc.code
        if code is None:
            return 0
        return int(code) if isinstance(code, int) else 1
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
