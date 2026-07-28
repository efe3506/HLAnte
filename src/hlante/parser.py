"""
hlante.parser
============

Parsers for HLA typing tool output files.

Supported tools
---------------
- ARCAS-HLA (``*.genotype.json``)
- T1K (``*_genotype.tsv`` / ``result_hla_genotype.tsv``)
- HLA-HD (``*_final.result.txt``)
- OptiType (``*_result.tsv``)

Each tool has a dedicated ``parse_*`` function; all return a list of
:class:`HLAGenotype` dataclasses so downstream modules can work against
a single schema.
"""

from __future__ import annotations

import io
import json
import logging
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Union

logger: logging.Logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

TOOL_ARCASHLA: str = "arcashla"
TOOL_T1K: str = "t1k"
TOOL_HLAHD: str = "hlahd"
TOOL_OPTITYPE: str = "optitype"

SUPPORTED_TOOLS: frozenset[str] = frozenset({TOOL_ARCASHLA, TOOL_T1K, TOOL_HLAHD, TOOL_OPTITYPE})

_TOOL_ALIASES: Dict[str, str] = {
    "arcas": TOOL_ARCASHLA,
    "arcashla": TOOL_ARCASHLA,
    "arcas-hla": TOOL_ARCASHLA,
    "arcas_hla": TOOL_ARCASHLA,
    "t1k": TOOL_T1K,
    "hlahd": TOOL_HLAHD,
    "hla-hd": TOOL_HLAHD,
    "hla_hd": TOOL_HLAHD,
    "optitype": TOOL_OPTITYPE,
    "opti-type": TOOL_OPTITYPE,
}

#: Field-count labels for the TSV ``resolution`` column, keyed by the
#: number of colon-separated fields.
RESOLUTION_LABELS: Dict[int, str] = {
    1: "one-field",
    2: "two-field",
    3: "three-field",
    4: "four-field",
}

RESOLUTION_G_GROUP: str = "G-group"
RESOLUTION_P_GROUP: str = "P-group"

# Spec regex: ``^[A-Z]+\*\d{2,3}(:\d{2,3})*[A-Z]?$``.
# Real HLA gene names end with digits (DRB1, DQB1, DPB1, DRB3, ...),
# so ``\d*`` is allowed between the letters and ``*``.
ALLELE_REGEX: re.Pattern[str] = re.compile(r"^[A-Z]+\d*\*\d{2,3}(:\d{2,3})*[A-Z]?$")

_HLA_PREFIX: str = "HLA-"


# ---------------------------------------------------------------------------
# Exceptions
# ---------------------------------------------------------------------------


class HLAnteParseError(Exception):
    """
    Raised for unrecognised or malformed tool outputs.

    Parameters
    ----------
    message : str
        Error message.
    """

    def __init__(self, message: str) -> None:
        super().__init__(message)


class UnsupportedToolError(HLAnteParseError):
    """
    Raised when an unsupported tool name is passed to the dispatcher.
    """

    pass


def _read_text(path: "Path") -> str:
    """
    Read a typing-tool output file as text, tolerating byte-order marks and
    UTF-16 (commonly emitted by Windows tools). A file that cannot be decoded
    is reported as a :class:`HLAnteParseError` rather than crashing with a raw
    ``UnicodeDecodeError`` that escapes the CLI's error handling.
    """
    raw = path.read_bytes()
    if raw.startswith((b"\xff\xfe", b"\xfe\xff")):
        encoding = "utf-16"
    elif raw.startswith(b"\xef\xbb\xbf"):
        encoding = "utf-8-sig"
    else:
        encoding = "utf-8"
    try:
        return raw.decode(encoding)
    except UnicodeDecodeError as exc:
        raise HLAnteParseError(
            f"Could not decode {path} as text (tried {encoding}); "
            f"re-save the file as UTF-8: {exc}"
        ) from exc


# ---------------------------------------------------------------------------
# Shared data structure
# ---------------------------------------------------------------------------


@dataclass
class HLAGenotype:
    """
    Shared dataclass representing a single HLA locus call.

    Attributes
    ----------
    sample_id : str
        Sample identifier (derived from file name).
    locus : str
        HLA locus with ``HLA-`` prefix (e.g., ``"HLA-A"``, ``"HLA-DRB1"``).
    allele1 : str
        First allele (``HLA-`` prefix stripped, e.g., ``"A*02:01"``).
    allele2 : str or None
        Second allele; ``None`` for homozygous or missing calls.
    resolution : str
        Resolution label: ``"one-field"``, ``"two-field"``,
        ``"three-field"``, ``"four-field"`` or
        ``"G-group"`` / ``"P-group"``.
    caller_quality1 : float or None
        Per-allele quality reported by the typing tool for allele 1, when the
        tool provides one. Populated only from the T1K native layout; arcasHLA
        and HLA-HD report no per-allele quality, and OptiType reports a
        solution-level objective rather than an allele quality.
    caller_quality2 : float or None
        Quality score reported by the tool, if any.
    tool : str
        Name of the tool that produced the call (``"arcashla"``,
        ``"t1k"``, ``"hlahd"``, ``"optitype"``).
    raw_line : str
        Original line / text for debugging.
    """

    sample_id: str
    locus: str
    allele1: str
    allele2: Optional[str]
    resolution: str
    tool: str
    raw_line: str
    caller_quality1: Optional[float] = None
    caller_quality2: Optional[float] = None


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _strip_hla_prefix(allele: str) -> str:
    """
    Strip the leading ``HLA-`` prefix from an allele name.

    Parameters
    ----------
    allele : str
        Input allele.

    Returns
    -------
    str
        Allele without the prefix.
    """
    if allele.upper().startswith(_HLA_PREFIX):
        return allele[len(_HLA_PREFIX) :]
    return allele


def _normalize_locus(locus: str) -> str:
    """
    Normalize a locus label to the ``HLA-<GENE>`` form.

    Parameters
    ----------
    locus : str
        Raw locus label (``"A"``, ``"HLA-A"``, ``"hla-drb1"``, ...).

    Returns
    -------
    str
        Normalized locus label.
    """
    stripped = _strip_hla_prefix(locus).upper()
    return f"{_HLA_PREFIX}{stripped}"


def _strip_trailing_asterisk(raw: str, tool: str) -> str:
    """
    Strip a trailing ``*`` from an allele string.

    Several HLA typing tools append a trailing ``*`` to uncertain calls
    (e.g. ``"B*40:02*"``).  This is non-standard but observed in real
    arcasHLA, T1K, and HLA-HD outputs.  The asterisk is stripped before
    regex validation; a DEBUG message records each occurrence.

    Parameters
    ----------
    raw : str
        Allele string (``HLA-`` prefix already stripped).
    tool : str
        Tool name used in the DEBUG log message.

    Returns
    -------
    str
        Allele with trailing ``*`` removed, or the original string
        unchanged when no trailing ``*`` is present.
    """
    if raw.endswith("*"):
        cleaned = raw.rstrip("*")
        logger.debug(
            "%s trailing-asterisk notation stripped: %r → %r",
            tool,
            raw,
            cleaned,
        )
        return cleaned
    return raw


def _sanitize_arcashla_allele(raw: str) -> str:
    """
    Normalise an arcasHLA allele string before regex validation.

    arcasHLA uses two non-standard notations:

    * **Trailing asterisk** — ``"DRB1*04:92*"`` — indicates an uncertain
      call.  The trailing ``*`` is stripped; the allele is otherwise valid.
    * **Space-separated pair** — ``"B*49:01 50:01"`` — arcasHLA
      ambiguity notation for two equally likely alleles.  The first token
      is taken as the primary call.

    In both cases a DEBUG message is logged.

    Parameters
    ----------
    raw : str
        Allele string from an arcasHLA genotype JSON, already stripped of
        the ``HLA-`` prefix.

    Returns
    -------
    str
        Cleaned allele string.
    """
    cleaned = raw

    # Space-separated pair: "B*49:01 50:01" → "B*49:01"
    if " " in cleaned:
        primary = cleaned.split()[0]
        logger.debug(
            "arcasHLA ambiguous pair %r → taking primary call %r",
            cleaned,
            primary,
        )
        cleaned = primary

    return _strip_trailing_asterisk(cleaned, "arcasHLA")


def _validate_allele(allele: str) -> bool:
    """
    Return whether an allele passes the accepted regex.

    Parameters
    ----------
    allele : str
        Allele to validate (``HLA-`` prefix is stripped automatically).

    Returns
    -------
    bool
        ``True`` when the regex matches.
    """
    if not allele:
        return False
    return bool(ALLELE_REGEX.match(_strip_hla_prefix(allele)))


def _ensure_valid_allele(allele: Optional[str], source: Path) -> None:
    """
    Raise :class:`HLAnteParseError` when ``allele`` is invalid.

    Parameters
    ----------
    allele : str or None
        Allele to check; no-op when ``None``.
    source : Path
        Source file path included in the error message.

    Raises
    ------
    HLAnteParseError
        When the allele does not pass the regex.
    """
    if allele is None:
        return
    if not _validate_allele(allele):
        raise HLAnteParseError(f"Invalid allele format: {allele!r} (source: {source})")


def _determine_resolution(allele: str) -> str:
    """
    Determine the resolution label of an allele.

    Returns ``"one-field"``, ``"two-field"``, ``"three-field"`` or
    ``"four-field"`` for the number of colon-separated fields,
    matching :func:`hlante.normalizer._resolution_of`. Returns
    ``"G-group"`` / ``"P-group"`` for group suffixes.

    Notes
    -----
    Releases up to v0.1.0 wrote these labels on a digit scale
    (``"2-field"`` for a one-field call, ``"4-field"`` for a two-field
    call), which conflated fields with digits — current HLA
    nomenclature admits at most four fields. The word forms used here
    cannot be confused with the old numeric labels, so a downstream
    filter written against v0.1.0 fails loudly rather than silently
    selecting the wrong resolution.

    Parameters
    ----------
    allele : str
        IMGT-like allele.

    Returns
    -------
    str
        Resolution label.
    """
    core = _strip_hla_prefix(allele)
    if core.endswith("G"):
        return RESOLUTION_G_GROUP
    if core.endswith("P"):
        return RESOLUTION_P_GROUP

    if "*" in core:
        _, rest = core.split("*", 1)
    else:
        rest = core
    # Strip a single trailing nomenclature letter (N/L/S/Q/C/A).
    rest = re.sub(r"[A-Z]$", "", rest)
    if not rest:
        return RESOLUTION_LABELS[1]
    n = len(rest.split(":"))
    return RESOLUTION_LABELS[min(n, 4)]


def _resolve_tool(tool: str) -> str:
    """
    Map a user-supplied tool name to its canonical form.

    Parameters
    ----------
    tool : str
        Raw tool name.

    Returns
    -------
    str
        One of the ``TOOL_*`` constants.

    Raises
    ------
    UnsupportedToolError
        When the name cannot be resolved.
    """
    key = tool.strip().lower()
    resolved = _TOOL_ALIASES.get(key)
    if resolved is None:
        raise UnsupportedToolError(
            f"Unsupported tool: {tool!r}. Accepted values: {sorted(SUPPORTED_TOOLS)}"
        )
    return resolved


def _ensure_file(filepath: Union[str, Path], tool_label: str) -> Path:
    """
    Ensure the input path is an existing file.

    Parameters
    ----------
    filepath : str or Path
        Path to check.
    tool_label : str
        Tool label used in the error message.

    Returns
    -------
    Path
        A verified :class:`Path` instance.

    Raises
    ------
    FileNotFoundError
        When the path does not exist or is not a regular file.
    """
    path = Path(filepath)
    if not path.exists():
        raise FileNotFoundError(f"{tool_label} output file not found: {filepath}")
    if not path.is_file():
        raise FileNotFoundError(f"{tool_label} input path is not a file: {filepath}")
    return path


# ---------------------------------------------------------------------------
# Tool-specific parsers
# ---------------------------------------------------------------------------


def parse_arcashla(filepath: Union[str, Path]) -> List[HLAGenotype]:
    """
    Parse an ARCAS-HLA ``.genotype.json`` output.

    Supported schemas:

    - Flat dictionary: ``{"HLA-A": ["A*02:01:01G", "A*24:02:01G"], ...}``
    - Nested under ``"alleles"``:
      ``{"alleles": {"HLA-A": [...], ...}}``

    Parameters
    ----------
    filepath : str or Path
        ARCAS-HLA JSON file.

    Returns
    -------
    list of HLAGenotype
        One record per locus.

    Raises
    ------
    FileNotFoundError
        When the file is missing.
    HLAnteParseError
        When the JSON is malformed or the schema is not recognised.
    """
    path = _ensure_file(filepath, "ARCAS-HLA")

    try:
        data: Any = json.loads(_read_text(path))
    except json.JSONDecodeError as exc:
        raise HLAnteParseError(f"ARCAS-HLA JSON parse error ({path}): {exc}") from exc

    if not isinstance(data, dict):
        raise HLAnteParseError(f"ARCAS-HLA JSON root must be an object (source: {path})")

    if "alleles" in data and isinstance(data["alleles"], dict):
        alleles_map: Dict[str, Any] = data["alleles"]
    else:
        alleles_map = data

    sample_id = path.stem
    if sample_id.endswith(".genotype"):
        sample_id = sample_id[: -len(".genotype")]

    results: List[HLAGenotype] = []
    for locus_key, calls in alleles_map.items():
        if not isinstance(calls, list):
            raise HLAnteParseError(
                f"ARCAS-HLA: expected a list for {locus_key!r}, got "
                f"{type(calls).__name__} (source: {path})"
            )
        if not calls:
            continue

        allele1_raw = str(calls[0]).strip()
        allele2_raw = str(calls[1]).strip() if len(calls) > 1 else None

        allele1 = _sanitize_arcashla_allele(_strip_hla_prefix(allele1_raw))
        allele2 = _sanitize_arcashla_allele(_strip_hla_prefix(allele2_raw)) if allele2_raw else None

        _ensure_valid_allele(allele1, path)
        _ensure_valid_allele(allele2, path)

        results.append(
            HLAGenotype(
                sample_id=sample_id,
                locus=_normalize_locus(locus_key),
                allele1=allele1,
                allele2=allele2,
                resolution=_determine_resolution(allele1),
                tool=TOOL_ARCASHLA,
                raw_line=json.dumps({locus_key: calls}, ensure_ascii=False),
            )
        )

    if not results:
        raise HLAnteParseError(f"No locus calls found in ARCAS-HLA JSON: {path}")

    return results


def parse_t1k(filepath: Union[str, Path]) -> List[HLAGenotype]:
    """
    Parse a T1K ``*_genotype.tsv`` output.

    Two layouts are supported.

    **Native (headerless) T1K output.** Current T1K releases emit 8
    tab-separated columns with no header row::

        gene  allele_count  allele1  score1  qual1  allele2  score2  qual2

    When ``allele_count`` is ``1`` (only one allele called), ``allele2``
    is written as ``"."`` and ``qual2`` as ``-1``. Missing alleles are
    marked with ``.`` or ``*``.

    **Legacy headered layout.** Earlier integrations (and the existing
    HLAnte test fixture) used a named 5-column layout whose first row
    reads ``gene  allele1  allele2  score1  score2``. This layout is
    still accepted for back-compatibility.

    Layout detection: if the first non-empty row starts with the
    literal ``gene`` token, the headered path is taken; otherwise the
    native 8-column path is used.

    Parameters
    ----------
    filepath : str or Path
        T1K TSV file.

    Returns
    -------
    list of HLAGenotype
        One record per locus.

    Raises
    ------
    FileNotFoundError
        When the file is missing.
    HLAnteParseError
        When the row layout cannot be interpreted.
    """
    path = _ensure_file(filepath, "T1K")

    raw_lines: List[str] = _read_text(path).splitlines(keepends=True)

    data_rows = [ln.rstrip("\n") for ln in raw_lines if ln.strip()]
    if not data_rows:
        raise HLAnteParseError(f"T1K output file is empty: {path}")

    first_tokens = [col.strip().lower() for col in data_rows[0].split("\t")]
    headered = first_tokens[:1] == ["gene"]

    if headered:
        required = ["gene", "allele1", "allele2", "score1", "score2"]
        if first_tokens[: len(required)] != required:
            raise HLAnteParseError(
                f"T1K header is missing required columns. "
                f"Expected {required}, got: {first_tokens} (source: {path})"
            )
        rows_iter = enumerate(data_rows[1:], start=2)
    else:
        # Headerless native layout; treat every row as data.
        rows_iter = enumerate(data_rows, start=1)

    def _is_t1k_null(token: str) -> bool:
        return token in ("", "*", ".", "-", "NA")

    def _primary_call(token: str) -> str:
        # T1K emits equivalence-class ambiguity lists as one
        # Comma-separated string, e.g.
        # "HLA-DRA*01:01:01,HLA-DRA*01:01:02,...". Take the first
        # Candidate as the primary call; the ambiguity itself is
        # Surfaced later by the normalizer's ambiguity handling.
        return token.split(",", 1)[0].strip() if "," in token else token

    sample_id = path.stem
    results: List[HLAGenotype] = []
    for line_no, row in rows_iter:
        if row.startswith("#"):
            continue
        parts = [p.strip() for p in row.split("\t")]

        if headered:
            if len(parts) < 5:
                raise HLAnteParseError(
                    f"T1K row has missing columns (line {line_no}): {row!r} (source: {path})"
                )
            gene, a1, a2, s1, s2 = parts[:5]
            q1 = q2 = None
        else:
            # Native: gene count a1 s1 q1 a2 s2 q2
            if len(parts) < 5:
                raise HLAnteParseError(
                    f"T1K row has missing columns (line {line_no}): {row!r} (source: {path})"
                )
            gene = parts[0]
            a1 = parts[2] if len(parts) > 2 else ""
            s1 = parts[3] if len(parts) > 3 else ""
            a2 = parts[5] if len(parts) > 5 else ""
            s2 = parts[6] if len(parts) > 6 else ""
            # Columns 4 and 7 are T1K's per-allele quality values; columns 3
            # And 6 (read above) are abundance. Only the quality values are
            # Reported, and only for the native layout — the legacy headered
            # Layout labels its two numbers "score" without saying which.
            q1 = _parse_float(parts[4]) if len(parts) > 4 else None
            q2 = _parse_float(parts[7]) if len(parts) > 7 else None

        allele1 = (
            None
            if _is_t1k_null(a1)
            else _strip_trailing_asterisk(_strip_hla_prefix(_primary_call(a1)), TOOL_T1K)
        )
        allele2 = (
            None
            if _is_t1k_null(a2)
            else _strip_trailing_asterisk(_strip_hla_prefix(_primary_call(a2)), TOOL_T1K)
        )

        if allele1 is None:
            logger.debug("T1K: locus %s has no allele1, skipping.", gene)
            continue

        _ensure_valid_allele(allele1, path)
        _ensure_valid_allele(allele2, path)

        score1 = _parse_float(s1)
        score2 = _parse_float(s2)

        results.append(
            HLAGenotype(
                sample_id=sample_id,
                locus=_normalize_locus(gene),
                allele1=allele1,
                allele2=allele2,
                resolution=_determine_resolution(allele1),
                # A quality is only meaningful for a reported allele.
                caller_quality1=q1 if allele1 else None,
                caller_quality2=q2 if allele2 else None,
                tool=TOOL_T1K,
                raw_line=row,
            )
        )

    if not results:
        raise HLAnteParseError(f"No valid locus calls found in T1K output: {path}")

    return results


def parse_hlahd(filepath: Union[str, Path]) -> List[HLAGenotype]:
    """
    Parse an HLA-HD ``*_final.result.txt`` output.

    Format: the first column is the locus (``HLA-A`` or ``A``), followed
    by two allele columns. Unresolved / missing alleles are indicated
    with ``-`` or ``Not typed``. Fields may be tab- or multi-space-
    delimited.

    Parameters
    ----------
    filepath : str or Path
        HLA-HD final result file.

    Returns
    -------
    list of HLAGenotype
        One record per typed locus (untyped loci are skipped).

    Raises
    ------
    FileNotFoundError
        When the file is missing.
    HLAnteParseError
        When rows do not match the expected schema.
    """
    path = _ensure_file(filepath, "HLA-HD")

    sample_id = path.stem
    for suffix in ("_final.result", ".final.result", "_final", ".final"):
        if sample_id.endswith(suffix):
            sample_id = sample_id[: -len(suffix)]
            break

    results: List[HLAGenotype] = []
    with io.StringIO(_read_text(path)) as handle:
        for line_no, raw in enumerate(handle, start=1):
            line = raw.rstrip("\n")
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            # Real HLA-HD final.result files are tab-delimited and contain
            # The null marker ``"Not typed"`` (with a literal space). The
            # Former ``r"[\t ]+"`` split also treated that interior space
            # as a column boundary, producing ``['DRB5', 'Not', 'typed',
            # 'Not', 'typed']`` and then tripping on ``"Not"`` as an
            # Allele. Prefer tab-only splitting when any tab is present;
            # Fall back to whitespace splitting for rare space-only
            # Variants (which never use the two-word null token).
            if "\t" in line:
                parts = [p.strip() for p in line.split("\t") if p.strip()]
            else:
                parts = re.split(r"\s+", line.strip())
            if len(parts) < 2:
                raise HLAnteParseError(
                    f"HLA-HD row has missing columns (line {line_no}): {line!r} (source: {path})"
                )
            locus = parts[0]
            a1_raw = parts[1] if len(parts) > 1 else ""
            a2_raw = parts[2] if len(parts) > 2 else ""

            def _clean(value: str) -> Optional[str]:
                v = value.strip()
                if not v or v == "-" or v.lower() in {"not typed", "nottyped", "not consistent"}:
                    return None
                return _strip_trailing_asterisk(_strip_hla_prefix(v), TOOL_HLAHD)

            allele1 = _clean(a1_raw)
            allele2 = _clean(a2_raw)

            if allele1 is None:
                logger.debug("HLA-HD: locus %s is not typed, skipping.", locus)
                continue

            _ensure_valid_allele(allele1, path)
            _ensure_valid_allele(allele2, path)

            results.append(
                HLAGenotype(
                    sample_id=sample_id,
                    locus=_normalize_locus(locus),
                    allele1=allele1,
                    allele2=allele2,
                    resolution=_determine_resolution(allele1),
                    tool=TOOL_HLAHD,
                    raw_line=line,
                )
            )

    if not results:
        raise HLAnteParseError(f"No typed locus found in HLA-HD output: {path}")

    return results


def parse_optitype(filepath: Union[str, Path]) -> List[HLAGenotype]:
    """
    Parse an OptiType ``*_result.tsv`` output.

    OptiType emits only HLA class I (A, B, C) calls at two-field
    resolution (e.g., ``A*02:01``). Expected header columns:
    ``A1``, ``A2``, ``B1``, ``B2``, ``C1``, ``C2``, ``Reads``,
    ``Objective``. A leading unnamed index column is tolerated.

    Parameters
    ----------
    filepath : str or Path
        OptiType TSV file.

    Returns
    -------
    list of HLAGenotype
        One record each for A / B / C.

    Raises
    ------
    FileNotFoundError
        When the file is missing.
    HLAnteParseError
        When required columns are missing or the file is malformed.
    """
    path = _ensure_file(filepath, "OptiType")

    sample_id = path.stem
    if sample_id.endswith("_result"):
        sample_id = sample_id[: -len("_result")]

    lines = [ln.rstrip("\n") for ln in _read_text(path).splitlines() if ln.strip()]

    if len(lines) < 2:
        raise HLAnteParseError(
            f"OptiType output must contain at least a header and one data row: {path}"
        )

    header = [col.strip().lower() for col in lines[0].split("\t")]
    data_row = lines[1].split("\t")

    # OptiType can enumerate several optimal solutions. Only the top-ranked
    # One is annotated; the alternatives are not retained, and downstream
    # Annotations could differ between them, so the collapse is announced
    # Rather than performed silently.
    solution_count = len(lines) - 1
    if solution_count > 1:
        logger.warning(
            "OptiType output %s contains %d enumerated solutions; only the "
            "top-ranked solution is annotated and the alternatives are not "
            "retained.",
            path,
            solution_count,
        )

    required = {"a1", "a2", "b1", "b2", "c1", "c2"}
    if not required.issubset(set(header)):
        raise HLAnteParseError(
            f"OptiType header is missing required columns. "
            f"Expected superset of: {sorted(required)}, got: {header} "
            f"(source: {path})"
        )

    idx: Dict[str, int] = {name: i for i, name in enumerate(header)}

    def _cell(name: str) -> str:
        pos = idx[name]
        return data_row[pos].strip() if pos < len(data_row) else ""

    objective: Optional[float] = None
    if "objective" in idx:
        objective = _parse_float(_cell("objective"))

    raw_line = lines[1]
    results: List[HLAGenotype] = []
    for locus_letter in ("a", "b", "c"):
        a1_col, a2_col = f"{locus_letter}1", f"{locus_letter}2"
        a1_val = _cell(a1_col)
        a2_val = _cell(a2_col)

        allele1 = None if not a1_val or a1_val == "*" else _strip_hla_prefix(a1_val)
        allele2 = None if not a2_val or a2_val == "*" else _strip_hla_prefix(a2_val)

        if allele1 is None:
            logger.debug("OptiType: locus %s has no allele1, skipping.", locus_letter)
            continue

        _ensure_valid_allele(allele1, path)
        _ensure_valid_allele(allele2, path)

        results.append(
            HLAGenotype(
                sample_id=sample_id,
                locus=f"{_HLA_PREFIX}{locus_letter.upper()}",
                allele1=allele1,
                allele2=allele2,
                resolution=_determine_resolution(allele1),
                tool=TOOL_OPTITYPE,
                raw_line=raw_line,
            )
        )

    if not results:
        raise HLAnteParseError(f"No valid A/B/C calls found in OptiType output: {path}")

    return results


# ---------------------------------------------------------------------------
# Dispatcher
# ---------------------------------------------------------------------------


_DISPATCH: Dict[str, Callable[[Union[str, Path]], List[HLAGenotype]]] = {
    TOOL_ARCASHLA: parse_arcashla,
    TOOL_T1K: parse_t1k,
    TOOL_HLAHD: parse_hlahd,
    TOOL_OPTITYPE: parse_optitype,
}


def parse_hla_output(
    filepath: Union[str, Path],
    tool: str,
) -> List[HLAGenotype]:
    """
    Dispatch to the appropriate parser based on ``tool``.

    Parameters
    ----------
    filepath : str or Path
        Input file to parse.
    tool : str
        Tool name. Accepted aliases: ``arcashla``/``arcas-hla``,
        ``t1k``, ``hlahd``/``hla-hd``, ``optitype``. Case-insensitive.

    Returns
    -------
    list of HLAGenotype
        One record per locus.

    Raises
    ------
    FileNotFoundError
        When the file is missing.
    UnsupportedToolError
        When ``tool`` is not recognised.
    HLAnteParseError
        When the file is malformed.
    """
    resolved_tool = _resolve_tool(tool)
    logger.debug(
        "parse_hla_output: tool=%r (resolved=%r) path=%s",
        tool,
        resolved_tool,
        filepath,
    )
    return _DISPATCH[resolved_tool](filepath)


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------


def _parse_float(value: str) -> Optional[float]:
    """
    Safely convert a string to float; return ``None`` on failure.

    Parameters
    ----------
    value : str
        Value to convert.

    Returns
    -------
    float or None
        Converted value, or ``None`` when the input is empty / unusable.
    """
    if value is None:
        return None
    stripped = value.strip()
    if not stripped or stripped in {"*", "-", "NA", "nan", "NaN"}:
        return None
    try:
        return float(stripped)
    except ValueError:
        return None


__all__ = [
    "HLAGenotype",
    "HLAnteParseError",
    "UnsupportedToolError",
    "ALLELE_REGEX",
    "SUPPORTED_TOOLS",
    "TOOL_ARCASHLA",
    "TOOL_T1K",
    "TOOL_HLAHD",
    "TOOL_OPTITYPE",
    "parse_arcashla",
    "parse_t1k",
    "parse_hlahd",
    "parse_optitype",
    "parse_hla_output",
]
