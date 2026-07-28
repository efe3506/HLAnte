#!/usr/bin/env python3
"""
Build the 14-sample example batch used in the worked examples.

The batch is drawn from the Sanger sequence-based typing reference calls in
``benchmarks/1000g_HLA_types.tsv`` and written in the arcasHLA genotype JSON
layout, one file per sample. Selection is deterministic, so the batch can be
reproduced exactly:

  * the four carriers of a CPIC Level 1A sentinel allele that are quoted in the
    worked examples (HLA-B*57:01, HLA-B*58:01, HLA-B*15:02, HLA-A*31:01), and
  * the first two samples of each 1000 Genomes super-population, in sample-ID
    order, among those with a call at all five loci and not already selected.

Usage
-----
    python scripts/examples/build_example_batch.py \\
        --ground-truth benchmarks/1000g_HLA_types.tsv \\
        --output-dir   example_batch/
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path
from typing import Dict, List

LOCI = {
    "HLA-A": ("HLA-A 1", "HLA-A 2"),
    "HLA-B": ("HLA-B 1", "HLA-B 2"),
    "HLA-C": ("HLA-C 1", "HLA-C 2"),
    "HLA-DQB1": ("HLA-DQB1 1", "HLA-DQB1 2"),
    "HLA-DRB1": ("HLA-DRB1 1", "HLA-DRB1 2"),
}

# Carriers quoted in the worked examples, one per CPIC Level 1A sentinel.
SENTINEL_CARRIERS = ["HG02144", "HG01889", "HG01796", "HG02009"]

PER_POPULATION = 2
NO_CALL = {"", "None", "NA", "-"}


def is_complete(row: Dict[str, str]) -> bool:
    return all(
        row.get(col, "").strip() not in NO_CALL
        for cols in LOCI.values()
        for col in cols
    )


def to_arcashla(row: Dict[str, str]) -> Dict[str, List[str]]:
    """Render one reference row as an arcasHLA genotype object."""
    genotype: Dict[str, List[str]] = {}
    for locus, (first, second) in LOCI.items():
        gene = locus.split("-", 1)[1]
        genotype[locus] = [
            f"{gene}*{row[first].strip()}",
            f"{gene}*{row[second].strip()}",
        ]
    return genotype


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ground-truth", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    with args.ground_truth.open(newline="", encoding="utf-8") as handle:
        rows = [r for r in csv.DictReader(handle, delimiter="\t") if is_complete(r)]

    by_id = {r["Sample ID"]: r for r in rows}
    selected: List[str] = [s for s in SENTINEL_CARRIERS if s in by_id]

    by_population: Dict[str, List[str]] = defaultdict(list)
    for row in sorted(rows, key=lambda r: r["Sample ID"]):
        by_population[row["Region"]].append(row["Sample ID"])

    for population in sorted(by_population):
        added = 0
        for sample in by_population[population]:
            if added == PER_POPULATION:
                break
            if sample not in selected:
                selected.append(sample)
                added += 1

    args.output_dir.mkdir(parents=True, exist_ok=True)
    for sample in selected:
        destination = args.output_dir / f"{sample}.genotype.json"
        destination.write_text(
            json.dumps(to_arcashla(by_id[sample]), indent=2) + "\n",
            encoding="utf-8",
        )

    print(f"{len(selected)} samples written to {args.output_dir}")
    for sample in selected:
        marker = " (CPIC sentinel carrier)" if sample in SENTINEL_CARRIERS else ""
        print(f"  {sample} [{by_id[sample]['Region']}]{marker}")


if __name__ == "__main__":
    main()
