#!/usr/bin/env python3
"""
Convert 1000 Genomes Project HLA types to the 4 tool input formats expected by HLAnte.

Reads a TSV with columns:
  Region  Population  Sample ID  HLA-A 1  HLA-A 2  HLA-B 1  HLA-B 2  ...

and produces per-sample fixture files:
  <output-dir>/arcashla/{sample}.genotype.json
  <output-dir>/t1k/{sample}.t1k.tsv
  <output-dir>/hlahd/{sample}.hlahd.txt
  <output-dir>/optitype/{sample}.optitype.tsv  (Class I only)

Usage:
    python scripts/benchmark/convert_1000g_to_tool_formats.py \\
        --input benchmark/1000g_HLA_types.tsv \\
        --output-dir benchmark/fixtures_1000g/ \\
        --tools arcashla,t1k,hlahd,optitype
"""

import argparse
import csv
import json
import logging
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

_ALLELE_CELL_RE: re.Pattern = re.compile(r"^\d{2,3}(:\d{2,3})*$")

# Ordered list of loci in the 1000G file
LOCI: List[str] = ["A", "B", "C", "DQB1", "DRB1"]
CLASS_I_LOCI: List[str] = ["A", "B", "C"]

NULL_TOKENS = frozenset({"", "*", "-", "na", "not typed", "nottyped", "none", "."})
ALL_TOOLS = frozenset({"arcashla", "t1k", "hlahd", "optitype"})


def _is_null(value: str) -> bool:
    return value.strip().lower() in NULL_TOKENS


def _resolve_cell(
    raw: str,
    gene: str,
    sample_id: str,
    ambig_log: List[str],
    corrupt_log: List[str],
) -> Optional[str]:
    """
    Parse one allele cell from the 1000G TSV.

    '02:01/02' → 'GENE*02:01'  (keeps first option, logs ambiguity)
    '23:01'    → 'GENE*23:01'
    '' / '*'   → None
    """
    raw = raw.strip()
    if _is_null(raw):
        return None
    if "/" in raw:
        first = raw.split("/")[0].strip()
        ambig_log.append(f"{sample_id}\t{gene}\t{raw}\t{first}")
        raw = first
    if not _ALLELE_CELL_RE.match(raw):
        corrupt_log.append(f"{sample_id}\t{gene}\t{raw}")
        logger.warning("Corrupt allele value %r for %s/%s — skipped", raw, sample_id, gene)
        return None
    return f"{gene}*{raw}"


def _parse_row(
    row: Dict[str, str],
    ambig_log: List[str],
    missing_log: List[str],
    corrupt_log: List[str],
) -> Dict[str, Tuple[Optional[str], Optional[str]]]:
    """
    Extract per-locus allele pairs from a CSV row.

    Returns {gene: (allele1_or_None, allele2_or_None)}.
    """
    sample_id = row["Sample ID"].strip()
    alleles: Dict[str, Tuple[Optional[str], Optional[str]]] = {}
    for gene in LOCI:
        a1 = _resolve_cell(row.get(f"HLA-{gene} 1", ""), gene, sample_id, ambig_log, corrupt_log)
        a2 = _resolve_cell(row.get(f"HLA-{gene} 2", ""), gene, sample_id, ambig_log, corrupt_log)
        if a1 is None:
            missing_log.append(
                f"{sample_id}\t{gene}\t{row.get(f'HLA-{gene} 1', 'empty')}"
            )
        alleles[gene] = (a1, a2)
    return alleles


# ── Format writers ────────────────────────────────────────────────────────────

def _write_arcashla(
    sample_id: str,
    alleles: Dict[str, Tuple[Optional[str], Optional[str]]],
    out_dir: Path,
) -> bool:
    """
    Write ARCAS-HLA JSON (flat locus-keyed dict, no HLA- prefix on keys).

    {"A": ["A*23:01", "A*68:02"], "B": [...], ...}
    """
    data: Dict[str, List[str]] = {}
    for gene, (a1, a2) in alleles.items():
        calls = [a for a in (a1, a2) if a is not None]
        if calls:
            data[gene] = calls
    if not data:
        return False
    path = out_dir / f"{sample_id}.genotype.json"
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return True


def _write_t1k(
    sample_id: str,
    alleles: Dict[str, Tuple[Optional[str], Optional[str]]],
    out_dir: Path,
) -> bool:
    """
    Write T1K headered TSV.

    gene  allele1   allele2   score1  score2
    A     A*23:01   A*68:02   100     100
    """
    rows: List[str] = []
    for gene, (a1, a2) in alleles.items():
        if a1 is None:
            continue
        a2_val = a2 if a2 is not None else "."
        rows.append(f"{gene}\t{a1}\t{a2_val}\t100\t100")
    if not rows:
        return False
    path = out_dir / f"{sample_id}.t1k.tsv"
    path.write_text(
        "gene\tallele1\tallele2\tscore1\tscore2\n" + "\n".join(rows) + "\n",
        encoding="utf-8",
    )
    return True


def _write_hlahd(
    sample_id: str,
    alleles: Dict[str, Tuple[Optional[str], Optional[str]]],
    out_dir: Path,
) -> bool:
    """
    Write HLA-HD tab-delimited TXT (no header, HLA- prefix on locus column).

    HLA-A     A*23:01   A*68:02
    HLA-B     B*13:02   B*42:01
    """
    rows: List[str] = []
    for gene, (a1, a2) in alleles.items():
        if a1 is None:
            continue
        a2_val = a2 if a2 is not None else "-"
        rows.append(f"HLA-{gene}\t{a1}\t{a2_val}")
    if not rows:
        return False
    path = out_dir / f"{sample_id}.hlahd.txt"
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return True


def _write_optitype(
    sample_id: str,
    alleles: Dict[str, Tuple[Optional[str], Optional[str]]],
    out_dir: Path,
) -> bool:
    """
    Write OptiType TSV (Class I only: A, B, C).

    A1       A2       B1       B2       C1       C2       Reads  Objective
    A*23:01  A*68:02  B*13:02  B*42:01  C*08:04  C*17:01  1500   1485.0

    Skipped when no Class I allele is available.
    """
    has_class_i = any(alleles.get(g, (None, None))[0] is not None for g in CLASS_I_LOCI)
    if not has_class_i:
        return False

    def _slot(gene: str, idx: int) -> str:
        pair = alleles.get(gene, (None, None))
        val = pair[idx]
        return val if val is not None else "*"

    a1, a2 = _slot("A", 0), _slot("A", 1)
    b1, b2 = _slot("B", 0), _slot("B", 1)
    c1, c2 = _slot("C", 0), _slot("C", 1)

    path = out_dir / f"{sample_id}.optitype.tsv"
    path.write_text(
        "A1\tA2\tB1\tB2\tC1\tC2\tReads\tObjective\n"
        f"{a1}\t{a2}\t{b1}\t{b2}\t{c1}\t{c2}\t1500\t1485.0\n",
        encoding="utf-8",
    )
    return True


_WRITERS = {
    "arcashla": _write_arcashla,
    "t1k": _write_t1k,
    "hlahd": _write_hlahd,
    "optitype": _write_optitype,
}


# ── Main conversion ───────────────────────────────────────────────────────────

def convert(input_path: Path, output_dir: Path, tools: frozenset) -> None:
    for tool in tools:
        (output_dir / tool).mkdir(parents=True, exist_ok=True)

    ambig_log: List[str] = []
    missing_log: List[str] = []
    corrupt_log: List[str] = []
    total = samples_ambig = samples_null = files_written = 0

    with input_path.open(encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            sample_id = row.get("Sample ID", "").strip()
            if not sample_id:
                continue

            prev_ambig = len(ambig_log)
            prev_missing = len(missing_log)

            alleles = _parse_row(row, ambig_log, missing_log, corrupt_log)

            if len(ambig_log) > prev_ambig:
                samples_ambig += 1
            if len(missing_log) > prev_missing:
                samples_null += 1

            for tool in tools:
                written = _WRITERS[tool](sample_id, alleles, output_dir / tool)
                if written:
                    files_written += 1

            total += 1
            if total % 500 == 0:
                logger.info("  Processed %d samples …", total)

    # Write logs
    ambig_path = output_dir / "ambiguity_resolved.log"
    missing_path = output_dir / "missing_calls.log"
    corrupt_path = output_dir / "corrupt_alleles.log"

    with ambig_path.open("w", encoding="utf-8") as fh:
        fh.write("# sample_id\tgene\toriginal\tkept\n")
        if ambig_log:
            fh.write("\n".join(ambig_log) + "\n")

    with missing_path.open("w", encoding="utf-8") as fh:
        fh.write("# sample_id\tgene\toriginal_value\n")
        if missing_log:
            fh.write("\n".join(missing_log) + "\n")

    with corrupt_path.open("w", encoding="utf-8") as fh:
        fh.write("# sample_id\tgene\traw_value\n")
        if corrupt_log:
            fh.write("\n".join(corrupt_log) + "\n")

    print(f"\nTotal samples processed:   {total:>6}")
    print(f"Samples with ambiguity:    {samples_ambig:>6}")
    print(f"Samples with null calls:   {samples_null:>6}")
    print(f"Corrupt allele entries:    {len(corrupt_log):>6}")
    print(f"Files written:             {files_written:>6}  ({len(tools)} tool(s) × ~{total})")
    print(f"Ambiguity log:  {ambig_path}")
    print(f"Missing log:    {missing_path}")
    print(f"Corrupt log:    {corrupt_path}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("benchmark/1000g_HLA_types.tsv"),
        help="1000G HLA TSV file.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("benchmark/fixtures_1000g/"),
        help="Root directory for fixture output.",
    )
    parser.add_argument(
        "--tools",
        default="arcashla,t1k,hlahd,optitype",
        help="Comma-separated tools to convert (default: all four).",
    )
    args = parser.parse_args()

    tools = frozenset(t.strip().lower() for t in args.tools.split(",") if t.strip())
    unknown = tools - ALL_TOOLS
    if unknown:
        logger.error("Unknown tools: %s. Valid choices: %s", unknown, sorted(ALL_TOOLS))
        sys.exit(1)

    if not args.input.is_file():
        logger.error("Input file not found: %s", args.input)
        sys.exit(1)

    logger.info(
        "Converting %s  →  %s  (tools: %s)",
        args.input,
        args.output_dir,
        sorted(tools),
    )
    convert(args.input, args.output_dir, tools)


if __name__ == "__main__":
    main()
