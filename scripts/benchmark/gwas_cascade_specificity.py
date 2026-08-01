#!/usr/bin/env python3
"""
Specificity of the GWAS fallback cascade.

The cascade answers a query by truncating the submitted allele one field at a
time until the GWAS Catalog index yields a record, so an association can be
returned for a less specific name than the one submitted. This script
quantifies how often that happens over the benchmark cohort, and what is lost
when it does.

Reported per locus and overall:

* **Broadening rate** — the share of annotated alleles whose association was
  found only after truncation. An association returned without truncation is
  reported for exactly the allele submitted.
* **Scope** — the resolution of the catalogue record that matched, which is a
  property of the record and not of the query.
* **Trait dilution** — the number of distinct traits attached by a broadened
  match relative to an exact match. Broadening trades specificity for recall;
  this is the size of that trade.
* **Sibling ambiguity** — how many other alleles in the cohort share the
  truncated prefix that produced a broadened match, i.e. how many distinct
  alleles would receive the same annotation.
* **Negative controls** — alleles absent from the catalogue at every
  resolution. The cascade must return nothing for these; anything else would
  mean it invents associations.

Usage
-----
    python scripts/benchmark/gwas_cascade_specificity.py \\
        --ground-truth benchmarks/1000g_HLA_types.tsv \\
        --output-dir   benchmarks/1000g_cascade_specificity/
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

from hlante.db.gwas import GWASClient

LOCI = ("HLA-A", "HLA-B", "HLA-C", "HLA-DQB1", "HLA-DRB1")
NO_CALL = {"", "None", "NA", "-", "*", "."}

#: Cell grammar the conversion to typing-tool input accepts. Applied here too:
#: this analysis has to describe the same cohort the pipeline was run on, not
#: the raw reference table.
_CELL_RE = re.compile(r"^\d{2,3}(:\d{2,3})*$")


def _resolve_cell(raw: str) -> Optional[str]:
    """Normalise one reference-table cell the way the conversion does."""
    raw = raw.strip()
    if raw in NO_CALL:
        return None
    if "/" in raw:
        raw = raw.split("/")[0].strip()
    if raw.endswith("*"):
        raw = raw[:-1].strip()
    return raw if _CELL_RE.match(raw) else None


def cohort_alleles(path: Path) -> Dict[str, Counter]:
    """Distinct allele calls per locus, with the number of carriers."""
    per_locus: Dict[str, Counter] = {locus: Counter() for locus in LOCI}
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            for locus in LOCI:
                gene = locus.split("-", 1)[1]
                for slot in ("1", "2"):
                    value = _resolve_cell(row.get(f"{locus} {slot}") or "")
                    if value is None:
                        continue
                    per_locus[locus][f"{gene}*{value}"] += 1
    return per_locus


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ground-truth", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    client = GWASClient()
    client.load()

    per_locus = cohort_alleles(args.ground_truth)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    rows: List[Dict[str, object]] = []
    # Which cohort alleles would collapse onto the same truncated key.
    prefix_members: Dict[str, Set[str]] = defaultdict(set)

    for locus, counter in per_locus.items():
        for allele in counter:
            for depth in range(1, allele.count(":") + 2):
                key = ":".join(allele.split(":")[:depth])
                prefix_members[key].add(allele)

    for locus, counter in per_locus.items():
        for allele, carriers in sorted(counter.items()):
            hits, _label = client.query_allele_with_fallback(allele)
            if not hits:
                rows.append(
                    {
                        "locus": locus,
                        "allele": allele,
                        "carriers": carriers,
                        "annotated": "no",
                        "broadened": "NA",
                        "scope": "NA",
                        "matched_allele": "NA",
                        "n_traits": 0,
                        "cohort_alleles_sharing_match": 0,
                    }
                )
                continue
            first = hits[0]
            matched = getattr(first, "matched_allele", "") or allele
            rows.append(
                {
                    "locus": locus,
                    "allele": allele,
                    "carriers": carriers,
                    "annotated": "yes",
                    "broadened": "yes" if getattr(first, "match_was_broadened", False) else "no",
                    "scope": getattr(first, "annotation_scope", "NA"),
                    "matched_allele": matched,
                    "n_traits": len({h.trait for h in hits}),
                    "cohort_alleles_sharing_match": len(prefix_members.get(matched, {allele})),
                }
            )

    detail = args.output_dir / "cascade_detail.tsv"
    with detail.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    annotated = [r for r in rows if r["annotated"] == "yes"]
    broadened = [r for r in annotated if r["broadened"] == "yes"]
    exact = [r for r in annotated if r["broadened"] == "no"]

    def mean(values: List[int]) -> Optional[float]:
        return round(sum(values) / len(values), 2) if values else None

    # Two views, both reported: per distinct allele (how the catalogue covers
    # The allele space) and per carrier (what a cohort of samples actually
    # Receives, which is dominated by common alleles).
    def carriers(subset: List[Dict[str, object]]) -> int:
        return sum(int(r["carriers"]) for r in subset)

    carriers_annotated = carriers(annotated)

    summary = {
        "cohort_distinct_alleles": len(rows),
        "annotated": len(annotated),
        "not_annotated_negative_controls": len(rows) - len(annotated),
        "broadened": len(broadened),
        "broadened_pct": round(100 * len(broadened) / len(annotated), 1) if annotated else None,
        "exact": len(exact),
        "exact_pct": round(100 * len(exact) / len(annotated), 1) if annotated else None,
        "scope_distribution": dict(Counter(r["scope"] for r in annotated)),
        "mean_traits_exact_match": mean([int(r["n_traits"]) for r in exact]),
        "mean_traits_broadened_match": mean([int(r["n_traits"]) for r in broadened]),
        "mean_cohort_alleles_sharing_a_broadened_match": mean(
            [int(r["cohort_alleles_sharing_match"]) for r in broadened]
        ),
        "carrier_weighted": {
            "annotated_calls": carriers_annotated,
            "not_annotated_calls": carriers(rows) - carriers_annotated,
            "exact_calls": carriers(exact),
            "exact_pct": round(100 * carriers(exact) / carriers_annotated, 1)
            if carriers_annotated
            else None,
            "broadened_calls": carriers(broadened),
            "broadened_pct": round(100 * carriers(broadened) / carriers_annotated, 1)
            if carriers_annotated
            else None,
        },
        "per_locus": {
            locus: {
                "distinct_alleles": sum(1 for r in rows if r["locus"] == locus),
                "annotated": sum(
                    1 for r in rows if r["locus"] == locus and r["annotated"] == "yes"
                ),
                "broadened": sum(
                    1 for r in rows if r["locus"] == locus and r["broadened"] == "yes"
                ),
            }
            for locus in LOCI
        },
    }
    (args.output_dir / "cascade_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )

    print(f"distinct alleles in cohort : {summary['cohort_distinct_alleles']}")
    print(f"  annotated                : {summary['annotated']}")
    print(f"  returned nothing         : {summary['not_annotated_negative_controls']}")
    print(f"  matched as submitted     : {summary['exact']} ({summary['exact_pct']}%)")
    print(f"  broadened                : {summary['broadened']} ({summary['broadened_pct']}%)")
    print(f"  scope                    : {summary['scope_distribution']}")
    print(f"  mean traits, exact       : {summary['mean_traits_exact_match']}")
    print(f"  mean traits, broadened   : {summary['mean_traits_broadened_match']}")
    print(
        "  cohort alleles sharing a broadened match (mean): "
        f"{summary['mean_cohort_alleles_sharing_a_broadened_match']}"
    )
    cw = summary["carrier_weighted"]
    print("\ncarrier-weighted (what a cohort actually receives):")
    print(f"  annotated calls          : {cw['annotated_calls']}")
    print(f"  matched as submitted     : {cw['exact_calls']} ({cw['exact_pct']}%)")
    print(f"  broadened                : {cw['broadened_calls']} ({cw['broadened_pct']}%)")
    print(f"\nwritten: {detail}")


if __name__ == "__main__":
    main()
