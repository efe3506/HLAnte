#!/usr/bin/env python3
"""
Measure format-independence of HLAnte's annotation output.

CLAIM: the same biological alleles fed via different tool formats
(ARCAS-HLA, T1K, HLA-HD, OptiType) produce identical HLAnte annotations.

For every sample that has fixtures for at least two tools, this script:
  1. Runs HLAnte annotations for each available tool.
  2. Pairs annotations by locus (allele name matched).
  3. Computes concordance metrics per pair:
       - Allele identity (same normalised allele name?)
       - imgt_accession identity
       - GWAS Jaccard: |A∩B| / |A∪B| on trait sets
       - PharmGKB Jaccard: on drug sets
       - Disease Jaccard: on condition sets
       - Confidence tier match (identical tier?)
       - Clinical significance match
  4. Produces a concordance matrix and per-pair detail TSV.

Expected result: ≥ 99.5 % concordance for shared Class I+II loci
across arcashla / t1k / hlahd; OptiType compared on Class I only.

Usage:
    python scripts/benchmark/inter_format_concordance.py \\
        --fixtures-dir benchmark/fixtures_1000g/ \\
        --output docs/INTER_FORMAT_CONCORDANCE.md \\
        --max-samples 500 \\
        --threads 4
"""

import argparse
import csv
import logging
import sys
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Set, Tuple

logger = logging.getLogger(__name__)

ALL_TOOLS = ["arcashla", "t1k", "hlahd", "optitype"]
CLASS_I_ONLY_TOOLS = {"optitype"}
CLASS_I_LOCI = {"HLA-A", "HLA-B", "HLA-C"}

FIXTURE_EXT = {
    "arcashla": ".genotype.json",
    "t1k":      ".t1k.tsv",
    "hlahd":    ".hlahd.txt",
    "optitype": ".optitype.tsv",
}

TOOL_PAIRS = [
    ("arcashla", "t1k"),
    ("arcashla", "hlahd"),
    ("t1k",      "hlahd"),
    ("arcashla", "optitype"),
    ("t1k",      "optitype"),
    ("hlahd",    "optitype"),
]


# ── Jaccard helpers ───────────────────────────────────────────────────────────

def jaccard(a: Set[str], b: Set[str]) -> float:
    if not a and not b:
        return 1.0
    union = len(a | b)
    return len(a & b) / union if union else 1.0


def _trait_set(annotated: List[Any]) -> Set[str]:
    return {h.trait.strip().lower() for ann in annotated for h in ann.gwas_hits if h.trait}


def _drug_set(annotated: List[Any]) -> Set[str]:
    drugs: Set[str] = set()
    for ann in annotated:
        for p in ann.pharm_annotations:
            if p.drug:
                drugs.add(p.drug.strip().lower())
        for e in ann.disease_entries:
            if e.condition:
                # include pharmacogenomics conditions (contain drug name)
                drugs.add(e.condition.strip().lower().split("[")[0].strip())
    return drugs


def _cv_set(annotated: List[Any]) -> Set[str]:
    return {
        e.condition.strip().lower().split("[")[0].strip()
        for ann in annotated for e in ann.disease_entries
        if e.condition
    }


# ── Per-sample annotation builder ────────────────────────────────────────────

def _annotate_tool(
    sample_id: str,
    tool: str,
    fixture_path: Path,
    imgt_db: Dict[str, Any],
    clients: Any,
    config: Any,
) -> Optional[Dict[str, Any]]:
    """
    Run HLAnte on one fixture file.

    Returns {locus → AnnotatedHLA} keyed by normalised locus name (HLA-X),
    or None on failure.
    """
    from hlante.parser import parse_hla_output, HLAnteParseError
    from hlante.normalizer import batch_normalize
    from hlante.annotator import annotate_genotype

    try:
        genotypes = parse_hla_output(fixture_path, tool)
        normalized = batch_normalize(genotypes, imgt_db=imgt_db)
        if not normalized:
            return None
        annotated_list = annotate_genotype(normalized, config, clients=clients)
    except Exception as exc:  # noqa: BLE001
        logger.debug("%s/%s failed: %s", tool, sample_id, exc)
        return None

    # Index by (locus, allele_index) for allele-level comparison
    by_allele: Dict[Tuple[str, int], Any] = {}
    for ann in annotated_list:
        na = ann.normalized_allele
        locus = na.source_locus or na.gene or "?"
        idx = na.allele_index or 0
        by_allele[(locus, idx)] = ann
    return by_allele  # type: ignore[return-value]


# ── Concordance accumulator ───────────────────────────────────────────────────

class PairMetrics:
    def __init__(self, tool_a: str, tool_b: str) -> None:
        self.tool_a = tool_a
        self.tool_b = tool_b
        self.n_allele_pairs = 0
        self.n_allele_identical = 0
        self.n_accession_identical = 0
        self.gwas_jaccard_sum = 0.0
        self.pharm_jaccard_sum = 0.0
        self.cv_jaccard_sum = 0.0
        self.n_tier_match = 0
        self.n_sig_match = 0

    @property
    def allele_identity(self) -> float:
        return self.n_allele_identical / self.n_allele_pairs if self.n_allele_pairs else 1.0

    @property
    def accession_identity(self) -> float:
        return self.n_accession_identical / self.n_allele_pairs if self.n_allele_pairs else 1.0

    @property
    def mean_gwas_jaccard(self) -> float:
        return self.gwas_jaccard_sum / self.n_allele_pairs if self.n_allele_pairs else 1.0

    @property
    def mean_pharm_jaccard(self) -> float:
        return self.pharm_jaccard_sum / self.n_allele_pairs if self.n_allele_pairs else 1.0

    @property
    def mean_cv_jaccard(self) -> float:
        return self.cv_jaccard_sum / self.n_allele_pairs if self.n_allele_pairs else 1.0

    @property
    def tier_match_rate(self) -> float:
        return self.n_tier_match / self.n_allele_pairs if self.n_allele_pairs else 1.0

    def add_pair(self, ann_a: Any, ann_b: Any, shared_loci_only: bool) -> None:
        """Compare two AnnotatedHLA records for the same (locus, allele_index)."""
        na_a = ann_a.normalized_allele
        na_b = ann_b.normalized_allele

        # For optitype pairs, only compare Class I loci
        locus = na_a.source_locus or na_a.gene or ""
        if shared_loci_only and locus not in CLASS_I_LOCI:
            return

        self.n_allele_pairs += 1

        if na_a.allele_name == na_b.allele_name:
            self.n_allele_identical += 1
        if na_a.imgt_accession and na_a.imgt_accession == na_b.imgt_accession:
            self.n_accession_identical += 1

        self.gwas_jaccard_sum += jaccard(
            {h.trait for h in ann_a.gwas_hits if h.trait},
            {h.trait for h in ann_b.gwas_hits if h.trait},
        )
        self.pharm_jaccard_sum += jaccard(
            {p.drug for p in ann_a.pharm_annotations if p.drug},
            {p.drug for p in ann_b.pharm_annotations if p.drug},
        )
        self.cv_jaccard_sum += jaccard(
            {e.condition for e in ann_a.disease_entries if e.condition},
            {e.condition for e in ann_b.disease_entries if e.condition},
        )
        tier_a = getattr(ann_a, "input_quality_tier", "NA")
        tier_b = getattr(ann_b, "input_quality_tier", "NA")
        if tier_a == tier_b:
            self.n_tier_match += 1
        if ann_a.clinical_significance == ann_b.clinical_significance:
            self.n_sig_match += 1


# ── Markdown report ───────────────────────────────────────────────────────────

def _fmt(val: float) -> str:
    return f"{val:.4f}"


def write_concordance_report(
    pair_metrics: Dict[Tuple[str, str], PairMetrics],
    n_samples: int,
    output_path: Path,
) -> None:
    lines: List[str] = []
    lines.append("# Format-Independence: HLAnte Concordance Matrix\n")
    lines.append(f"Generated: {datetime.now(timezone.utc).date().isoformat()}  ")
    lines.append(f"N samples compared: {n_samples}\n")

    lines.append(
        "> **Claim**: HLAnte's annotation output is independent of the "
        "input tool format. The same biological alleles, fed via different "
        "tool output formats, produce identical HLAnte annotations.\n"
    )

    # Mean GWAS Jaccard matrix
    lines.append("\n## Mean GWAS Jaccard Similarity\n")
    lines.append("(Mean over all matched allele pairs; * = Class I loci only)\n")
    header = "|          | " + " | ".join(f"{t:<9}" for t in ALL_TOOLS) + " |"
    sep    = "|----------|" + "-----------|" * len(ALL_TOOLS)
    lines.append(header)
    lines.append(sep)
    for ta in ALL_TOOLS:
        row = f"| {ta:<8} |"
        for tb in ALL_TOOLS:
            if ta == tb:
                row += " 1.0000    |"
            else:
                key = (ta, tb) if (ta, tb) in pair_metrics else (tb, ta)
                pm = pair_metrics.get(key)
                suffix = "*" if "optitype" in (ta, tb) else " "
                val = _fmt(pm.mean_gwas_jaccard) if pm else "  N/A  "
                row += f" {val}{suffix}    |"
        lines.append(row)

    # Per-pair detail table
    lines.append("\n## Per-Pair Concordance Detail\n")
    lines.append(
        "| Pair | Allele pairs | Identity | Accession | GWAS J | PharmGKB J | Tier match |"
    )
    lines.append(
        "|------|-------------|---------|---------|--------|-----------|-----------|"
    )
    for (ta, tb), pm in pair_metrics.items():
        suffix = " *" if "optitype" in (ta, tb) else ""
        lines.append(
            f"| {ta} vs {tb}{suffix} "
            f"| {pm.n_allele_pairs:>11} "
            f"| {_fmt(pm.allele_identity):>7} "
            f"| {_fmt(pm.accession_identity):>7} "
            f"| {_fmt(pm.mean_gwas_jaccard):>6} "
            f"| {_fmt(pm.mean_pharm_jaccard):>9} "
            f"| {_fmt(pm.tier_match_rate):>9} |"
        )

    lines.append("\n## Notes\n")
    lines.append(
        "- Concordance < 1.0 between tool formats is expected for a small "
        "fraction of alleles where ambiguity resolution differs (e.g., "
        "ambiguity-class '02:01/02' is resolved to the first option by the "
        "converter, but the original tool may have reported a different allele).\n"
    )
    lines.append(
        "- OptiType comparisons (*) are restricted to Class I loci "
        "(HLA-A, HLA-B, HLA-C).\n"
    )
    lines.append(
        "- GWAS Jaccard = 1.0 when both alleles have no GWAS hits (the "
        "empty-intersection / empty-union convention).\n"
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--fixtures-dir",
        type=Path,
        default=Path("benchmark/fixtures_1000g/"),
        help="Root dir with per-tool fixture subdirs (same as convert step).",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("docs/INTER_FORMAT_CONCORDANCE.md"),
        help="Output Markdown report path.",
    )
    parser.add_argument(
        "--max-samples", type=int, default=500,
        help="Max samples per tool pair (0 = all; default: 500).",
    )
    parser.add_argument(
        "--threads", type=int, default=4,
        help="Parallel worker threads.",
    )
    parser.add_argument(
        "--imgt-db-path", type=Path, default=None,
    )
    parser.add_argument(
        "--offline", action="store_true",
    )
    parser.add_argument(
        "--log-level", default="WARNING",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
    )
    args = parser.parse_args()

    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s %(levelname)-7s %(name)s: %(message)s",
        stream=sys.stderr,
    )

    try:
        from hlante.normalizer import load_imgt_db, IMGTDatabaseMissingError
        from hlante.annotator import AnnotatorConfig, build_clients
    except ImportError as exc:
        logger.error("Cannot import hlante: %s", exc)
        sys.exit(1)

    print("Loading IMGT-HLA database …", end=" ", flush=True)
    try:
        imgt_db = load_imgt_db(args.imgt_db_path)
        print(f"OK (version={imgt_db.get('version', 'unknown')})")
    except IMGTDatabaseMissingError as exc:
        print(f"FAILED: {exc}")
        sys.exit(1)

    config = AnnotatorConfig(offline=args.offline, population_group="global")
    clients = build_clients(config)

    # Collect sample IDs that exist for ALL four tools (or at least two)
    sample_sets: Dict[str, Set[str]] = {}
    for tool in ALL_TOOLS:
        tool_dir = args.fixtures_dir / tool
        if not tool_dir.is_dir():
            continue
        ext = FIXTURE_EXT[tool]
        sids: Set[str] = set()
        for f in tool_dir.glob(f"*{ext}"):
            stem = f.stem
            for suffix in (".genotype", ".t1k", ".hlahd", ".optitype"):
                if stem.endswith(suffix):
                    stem = stem[: -len(suffix)]
                    break
            sids.add(stem)
        sample_sets[tool] = sids
        print(f"[{tool}] {len(sids)} fixture samples found.")

    if len(sample_sets) < 2:
        print("ERROR: Need fixtures for at least 2 tools.", file=sys.stderr)
        sys.exit(1)

    # Samples common to at least two tools (use union of all pairs)
    all_sids: Set[str] = set().union(*sample_sets.values())
    # For efficiency, limit to max_samples
    if args.max_samples and args.max_samples > 0:
        all_sids = set(sorted(all_sids)[: args.max_samples])

    print(f"Comparing {len(all_sids)} samples across {len(sample_sets)} tools …")

    pair_metrics: Dict[Tuple[str, str], PairMetrics] = {}
    for ta, tb in TOOL_PAIRS:
        if ta in sample_sets and tb in sample_sets:
            pair_metrics[(ta, tb)] = PairMetrics(ta, tb)

    # Annotate samples per tool
    annotations: Dict[str, Dict[str, Optional[Dict]]] = defaultdict(dict)

    for tool, sids in sample_sets.items():
        ext = FIXTURE_EXT[tool]
        tool_dir = args.fixtures_dir / tool
        to_run = sorted(sids & all_sids)

        def _run(sid: str, tl: str = tool) -> None:
            fixture_path = tool_dir / f"{sid}{ext}"
            if not fixture_path.exists():
                return
            result = _annotate_tool(sid, tl, fixture_path, imgt_db, clients, config)
            annotations[tl][sid] = result

        if args.threads > 1:
            with ThreadPoolExecutor(max_workers=args.threads) as pool:
                futures = [pool.submit(_run, sid) for sid in to_run]
                for fut in as_completed(futures):
                    fut.result()
        else:
            for sid in to_run:
                _run(sid)

        n_ok = sum(1 for v in annotations[tool].values() if v is not None)
        print(f"  [{tool}] annotated {n_ok}/{len(to_run)}")

    # Compute concordance per pair
    n_samples_compared = 0
    for sample_id in sorted(all_sids):
        for (ta, tb), pm in pair_metrics.items():
            ann_a = annotations.get(ta, {}).get(sample_id)
            ann_b = annotations.get(tb, {}).get(sample_id)
            if ann_a is None or ann_b is None:
                continue

            # Match alleles by (locus, allele_index) key
            shared_loci_only = "optitype" in (ta, tb)
            common_keys = set(ann_a.keys()) & set(ann_b.keys())
            for key in common_keys:
                locus, _idx = key
                if shared_loci_only and locus not in CLASS_I_LOCI:
                    continue
                pm.add_pair(ann_a[key], ann_b[key], shared_loci_only=shared_loci_only)

        n_samples_compared += 1

    print(f"\nConcordance computed over {n_samples_compared} samples.")

    # Write detail TSV
    detail_rows = []
    for (ta, tb), pm in pair_metrics.items():
        detail_rows.append({
            "tool_a": ta, "tool_b": tb,
            "n_allele_pairs": pm.n_allele_pairs,
            "allele_identity": f"{pm.allele_identity:.4f}",
            "accession_identity": f"{pm.accession_identity:.4f}",
            "gwas_jaccard": f"{pm.mean_gwas_jaccard:.4f}",
            "pharm_jaccard": f"{pm.mean_pharm_jaccard:.4f}",
            "cv_jaccard": f"{pm.mean_cv_jaccard:.4f}",
            "tier_match": f"{pm.tier_match_rate:.4f}",
        })
    tsv_path = args.output.with_suffix(".tsv")
    if detail_rows:
        tsv_path.parent.mkdir(parents=True, exist_ok=True)
        with tsv_path.open("w", encoding="utf-8", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(detail_rows[0].keys()), delimiter="\t")
            writer.writeheader()
            writer.writerows(detail_rows)

    write_concordance_report(pair_metrics, n_samples_compared, args.output)
    print(f"Report: {args.output}")
    print(f"Detail: {tsv_path}")


if __name__ == "__main__":
    main()
