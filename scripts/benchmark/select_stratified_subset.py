#!/usr/bin/env python3
"""
Select a stratified N-sample subset from the 1000G HLA types file.

Selection rules:
  - N // 5 samples per super-population (AFR, AMR, EAS, EUR, SAS)
  - Round-robin across sub-populations within each super-pop
  - Priority for samples carrying clinically important alleles:
      B*57:01  (abacavir CPIC 1A),  B*58:01  (allopurinol CPIC 1A)
      B*15:02  (carbamazepine 1A),  A*31:01  (carbamazepine 1A)
      B*27:05  (ankylosing spondylitis), B*51:01 (Behcet disease)
      DRB1*03:01 (SLE / T1D),  DRB1*04:01 (RA)
      DRB1*15:01 (multiple sclerosis), DQB1*06:02 (narcolepsy)
      DQB1*02:01 (celiac disease)

Outputs:
  benchmark/subset_30.tsv           (subset in 1000G format)
  benchmark/subset_30_carriers.tsv  (clinical allele carriers)

Usage:
    python scripts/benchmark/select_stratified_subset.py \\
        --input benchmark/1000g_HLA_types.tsv \\
        --output benchmark/subset_30.tsv \\
        --n 30
"""

import argparse
import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Set, Tuple

SUPERPOPS = ["AFR", "AMR", "EAS", "EUR", "SAS"]

# Target allele (gene*fields) → minimum carriers in subset
CLINICAL_TARGETS: Dict[str, int] = {
    "B*57:01":    3,   # Abacavir CPIC 1A
    "B*58:01":    3,   # Allopurinol CPIC 1A
    "B*15:02":    2,   # Carbamazepine CPIC 1A
    "A*31:01":    2,   # Carbamazepine CPIC 1A
    "B*27:05":    3,   # Ankylosing spondylitis
    "B*51:01":    2,   # Behcet disease
    "DRB1*03:01": 3,   # SLE / type 1 diabetes
    "DRB1*04:01": 3,   # Rheumatoid arthritis
    "DRB1*15:01": 3,   # Multiple sclerosis
    "DQB1*06:02": 3,   # Narcolepsy
    "DQB1*02:01": 2,   # Celiac disease
}

LOCI = ["A", "B", "C", "DQB1", "DRB1"]


def _carried(row: Dict[str, str], allele_key: str) -> bool:
    """Return True if the sample carries allele_key (e.g. 'B*57:01')."""
    gene, fields = allele_key.split("*", 1)
    for slot in (1, 2):
        if row.get(f"HLA-{gene} {slot}", "").strip() == fields:
            return True
    return False


def _clinical_score(row: Dict[str, str]) -> int:
    return sum(1 for a in CLINICAL_TARGETS if _carried(row, a))


def _carried_list(row: Dict[str, str]) -> List[str]:
    return [a for a in CLINICAL_TARGETS if _carried(row, a)]


def _round_robin(candidates: List[Dict[str, str]], n: int) -> List[Dict[str, str]]:
    """
    Select up to n samples with round-robin across sub-populations.

    Within each sub-pop candidates are assumed pre-sorted by clinical score
    descending; the round-robin picks one per sub-pop per cycle.
    """
    by_subpop: Dict[str, List[Dict[str, str]]] = defaultdict(list)
    for s in candidates:
        by_subpop[s["Population"]].append(s)

    subpops = sorted(by_subpop)
    indices = {p: 0 for p in subpops}
    selected: List[Dict[str, str]] = []

    while len(selected) < n:
        advanced = False
        for pop in subpops:
            if indices[pop] < len(by_subpop[pop]):
                selected.append(by_subpop[pop][indices[pop]])
                indices[pop] += 1
                advanced = True
                if len(selected) >= n:
                    break
        if not advanced:
            break
    return selected[:n]


def select_subset(
    input_path: Path,
    n: int = 30,
) -> Tuple[List[Dict[str, str]], List[str]]:
    """
    Select a stratified subset of n samples.

    Strategy:
      Phase 1 — balanced core selection.
        For each super-pop, sort candidates by clinical score (descending),
        then pick n_per_pop samples via round-robin across sub-populations.
        This guarantees n_per_pop slots per super-pop while maximising the
        chance that clinical-allele carriers are included.

      Phase 2 — top up clinical minimums.
        If any CLINICAL_TARGETS minimum is still unmet after Phase 1, add
        the smallest number of additional carriers needed, choosing from the
        super-pop with the most remaining budget (or any, if all are full).

    Returns (selected_rows, fieldnames).
    """
    n_per_pop = max(1, n // len(SUPERPOPS))

    all_rows: List[Dict[str, str]] = []
    fieldnames: List[str] = []
    with input_path.open(encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fieldnames = list(reader.fieldnames or [])
        all_rows = list(reader)

    by_superpop: Dict[str, List[Dict[str, str]]] = defaultdict(list)
    for row in all_rows:
        by_superpop[row["Region"]].append(row)

    carriers: Dict[str, List[Dict[str, str]]] = {
        a: [r for r in all_rows if _carried(r, a)] for a in CLINICAL_TARGETS
    }

    # ── Phase 1: balanced core — n_per_pop per super-pop ─────────────────────
    selected: List[Dict[str, str]] = []
    selected_ids: Set[str] = set()

    for pop in SUPERPOPS:
        candidates = list(by_superpop[pop])
        # Sort: clinical score descending so high-value carriers rise to the top
        candidates.sort(key=_clinical_score, reverse=True)
        fill = _round_robin(candidates, n_per_pop)
        for s in fill:
            selected.append(s)
            selected_ids.add(s["Sample ID"])

    # ── Phase 2: top up unmet clinical minimums ───────────────────────────────
    # Process alleles by rarity (fewest carriers first) so that hard-to-find
    # Alleles are prioritised when the budget is tight.
    for allele, min_count in sorted(CLINICAL_TARGETS.items(), key=lambda x: len(carriers[x[0]])):
        have = sum(1 for s in selected if _carried(s, allele))
        if have >= min_count:
            continue
        for s in carriers[allele]:
            if have >= min_count:
                break
            if s["Sample ID"] not in selected_ids:
                selected.append(s)
                selected_ids.add(s["Sample ID"])
                have += 1

    return selected[:n], fieldnames


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("benchmark/1000g_HLA_types.tsv"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("benchmark/subset_30.tsv"),
    )
    parser.add_argument("--n", type=int, default=30, help="Subset size (default: 30).")
    args = parser.parse_args()

    if not args.input.is_file():
        print(f"ERROR: Input not found: {args.input}", file=sys.stderr)
        sys.exit(1)

    print(f"Selecting {args.n} stratified samples from {args.input} …")
    selected, fieldnames = select_subset(args.input, args.n)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(selected)

    carriers_path = args.output.with_name(args.output.stem + "_carriers.tsv")
    with carriers_path.open("w", encoding="utf-8") as fh:
        fh.write("Sample ID\tRegion\tPopulation\tCarried clinical alleles\n")
        for s in selected:
            carried = _carried_list(s)
            if carried:
                fh.write(
                    f"{s['Sample ID']}\t{s['Region']}\t{s['Population']}\t"
                    f"{', '.join(carried)}\n"
                )

    # ── Print summary ─────────────────────────────────────────────────────────
    pop_counts: Counter = Counter(s["Region"] for s in selected)
    allele_counts: Counter = Counter(
        a for s in selected for a in CLINICAL_TARGETS if _carried(s, a)
    )

    print(f"\nSelected {len(selected)} samples:")
    for pop in SUPERPOPS:
        print(f"  {pop}: {pop_counts.get(pop, 0)}")

    print("\nClinical allele coverage:")
    for allele, min_count in CLINICAL_TARGETS.items():
        count = allele_counts[allele]
        status = "✓" if count >= min_count else "✗"
        print(f"  {status} {allele:<16}  {count}/{min_count}")

    print(f"\nSubset TSV:    {args.output}  ({len(selected)} rows)")
    print(f"Carriers TSV:  {carriers_path}")


if __name__ == "__main__":
    main()
