#!/usr/bin/env python3
"""
How each benchmark allele call resolves against IPD-IMGT/HLA.

A name-match rate answers only whether a call could be matched at all, which
for a reference set of real IPD-IMGT/HLA names is close to tautological. What
matters downstream is *how* it matched: a two-field name matches unambiguously
while still denoting a group of alleles, and every annotation attached to it is
a property of that group rather than of one allele.

Reported over the benchmark cohort:

* **Category** — ``exact`` (the call is a listed allele), ``prefix_unique``
  (a prefix of exactly one listed allele), ``prefix_multiple`` (a prefix of
  several), ``g_group`` / ``p_group`` (group notation), or ``unmatched``.
* **Candidates** — how many release alleles a call maps to. The median over
  ``prefix_multiple`` calls is the honest measure of how much specificity a
  two-field reference call actually carries.
* **Rejected cells** — ground-truth cells that never reach the pipeline
  because ``convert_1000g_to_tool_formats.py`` cannot render them as a
  typing-tool call. Counted and broken down by cause so the denominator used
  here is reconcilable with the reference table it came from.

The denominator is the set of calls that entered the pipeline, which is the
reference table's non-empty cells minus the rejected ones. This script applies
the conversion rules itself rather than re-reading converted fixtures, so the
count cannot drift from the table it is derived from.

Usage
-----
    python scripts/benchmark/imgt_match_categories.py \\
        --ground-truth benchmarks/1000g_HLA_types.tsv \\
        --output-dir   benchmarks/1000g_imgt_match_categories/
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from collections import Counter
from pathlib import Path
from statistics import median
from typing import Dict, List, Optional, Tuple

from hlante.normalizer import MATCH_CATEGORIES, load_imgt_db, normalize_allele

LOCI = ("A", "B", "C", "DQB1", "DRB1")
NO_CALL = {"", "none", "na", "-", "*", "nan"}

#: Cell grammar accepted by ``convert_1000g_to_tool_formats.py``. Kept
#: identical here on purpose: this script must reject exactly what the
#: conversion rejects, or the denominator it reports would not be the one the
#: benchmark actually ran on.
CELL_RE = re.compile(r"^\d{2,3}(:\d{2,3})*$")


def rejection_cause(value: str) -> str:
    """Why the conversion cannot render *value* as a typing-tool call."""
    if value.endswith("*"):
        return "trailing_asterisk"
    if re.search(r"N$", value):
        return "null_allele_designation"
    if "," in value:
        return "malformed_numeric"
    if " " in value:
        return "two_alleles_in_one_cell"
    return "other"


def resolve_cell(raw: str) -> Tuple[Optional[str], Optional[str]]:
    """
    Apply the conversion rules to one reference-table cell.

    Returns ``(call, None)`` when the cell yields a call, or
    ``(None, cause)`` when the conversion would reject it. Empty cells yield
    ``(None, None)`` — absent, not rejected.
    """
    raw = raw.strip()
    if raw.lower() in NO_CALL:
        return None, None
    candidate = raw.split("/")[0].strip() if "/" in raw else raw
    if not CELL_RE.match(candidate):
        return None, rejection_cause(candidate)
    return candidate, None


def cohort_calls(path: Path) -> Tuple[List[Tuple[str, str, str]], Counter, int]:
    """
    Read the reference table into (sample, locus, call) triples.

    Also returns the rejected-cell tally and the number of non-empty cells,
    so the reported denominator can be reconciled with the source table.
    """
    calls: List[Tuple[str, str, str]] = []
    rejected: Counter = Counter()
    non_empty = 0
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            sample = (row.get("Sample ID") or "").strip()
            for gene in LOCI:
                for slot in ("1", "2"):
                    raw = (row.get(f"HLA-{gene} {slot}") or "").strip()
                    if raw.lower() in NO_CALL:
                        continue
                    non_empty += 1
                    call, cause = resolve_cell(raw)
                    if call is None:
                        rejected[cause] += 1
                        continue
                    calls.append((sample, f"HLA-{gene}", f"{gene}*{call}"))
    return calls, rejected, non_empty


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ground-truth", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--imgt-db-path",
        type=Path,
        default=None,
        help="Installed IPD-IMGT/HLA release; defaults to the local snapshot.",
    )
    args = parser.parse_args()

    imgt_db = load_imgt_db(args.imgt_db_path)
    calls, rejected, non_empty = cohort_calls(args.ground_truth)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    rows: List[Dict[str, object]] = []
    categories: Counter = Counter()
    candidates_by_category: Dict[str, List[int]] = {c: [] for c in MATCH_CATEGORIES}

    for sample, locus, call in calls:
        normalized = normalize_allele(call, imgt_db)
        if normalized is None:
            categories["unmatched"] += 1
            rows.append(
                {
                    "sample": sample,
                    "locus": locus,
                    "call": call,
                    "imgt_match_category": "unmatched",
                    "imgt_match_candidates": 0,
                }
            )
            continue
        category = normalized.imgt_match_category or "unmatched"
        candidates = normalized.imgt_match_candidates
        categories[category] += 1
        candidates_by_category.setdefault(category, []).append(candidates)
        rows.append(
            {
                "sample": sample,
                "locus": locus,
                "call": call,
                "imgt_match_category": category,
                "imgt_match_candidates": candidates,
            }
        )

    detail = args.output_dir / "imgt_match_detail.tsv"
    with detail.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["sample", "locus", "call", "imgt_match_category", "imgt_match_candidates"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)

    total = len(calls)
    single_allele = categories["exact"] + categories["prefix_unique"]
    multiple = candidates_by_category.get("prefix_multiple", [])
    summary = {
        "reference_table_non_empty_cells": non_empty,
        "rejected_by_conversion": sum(rejected.values()),
        "rejected_by_cause": dict(sorted(rejected.items())),
        "calls_entering_pipeline": total,
        "categories": {c: categories.get(c, 0) for c in MATCH_CATEGORIES},
        "categories_pct": {
            c: round(100 * categories.get(c, 0) / total, 2) if total else None
            for c in MATCH_CATEGORIES
        },
        "resolving_to_a_single_allele": single_allele,
        "resolving_to_a_single_allele_pct": round(100 * single_allele / total, 2) if total else None,
        "median_candidates_prefix_multiple": median(multiple) if multiple else None,
    }
    (args.output_dir / "imgt_match_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )

    print(f"{total} calls entering the pipeline from {non_empty} non-empty cells")
    for cause, count in sorted(rejected.items()):
        print(f"  rejected by conversion: {count:>5}  {cause}")
    for category in MATCH_CATEGORIES:
        count = categories.get(category, 0)
        pct = 100 * count / total if total else 0.0
        print(f"  {category:<16} {count:>6}  {pct:5.2f}%")
    if multiple:
        print(f"  median candidates for prefix_multiple: {median(multiple)}")
    print(f"Written to {args.output_dir}")


if __name__ == "__main__":
    main()
