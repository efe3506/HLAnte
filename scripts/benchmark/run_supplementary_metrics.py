#!/usr/bin/env python3
"""
Supplementary benchmark metrics for HLAnte.

Computes three metrics not captured by run_annotation_benchmark.py, reusing the
exact same parse -> normalize -> annotate path:

  R5  GWAS annotation-scope (cascade depth) distribution on the 1000G cohort
        -> gwas_scope_dist_1000g.tsv
  S2  Per-sentinel-allele confidence score by super-population
        -> sentinel_confidence_by_pop.tsv

S1 (typing_tool-mode metrics) is produced by re-running the main
run_annotation_benchmark.py with --input-source typing_tool; it is not
duplicated here.

This script imports helpers from run_annotation_benchmark.py so the annotation
logic is identical to the validated benchmark.

Usage:
    python scripts/benchmark/run_supplementary_metrics.py \\
        --fixtures-dir benchmarks/fixtures_1000g/ \\
        --ground-truth benchmarks/1000g_HLA_types.tsv \\
        --output-dir benchmarks/1000g_supplementary/ \\
        --tool arcashla \\
        --input-source validated \\
        --threads 8
"""

import argparse
import csv
import sys
import threading
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from statistics import mean
from typing import Dict, List

# Import the validated benchmark's helpers so annotation is byte-for-byte identical.
sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_annotation_benchmark import (  # noqa: E402
    CPIC_TARGETS,
    FIXTURE_EXT,
    GWAS_TARGETS,
    REGION_TO_AFND,
    _carried,
    annotate_sample,
    load_ground_truth,
)

# Sentinel alleles to profile for S2 (4 CPIC + 7 GWAS = 11)
SENTINELS: List[str] = list(CPIC_TARGETS.keys()) + list(GWAS_TARGETS.keys())

SUPERPOPS = ["EUR", "EAS", "AMR", "SAS", "AFR"]


def _strip_prefix(name: str) -> str:
    return name[4:] if name.startswith("HLA-") else name


def _matches_sentinel(allele_name: str, sentinel: str) -> bool:
    """True if a normalized allele_name corresponds to the sentinel key.

    Matches exact two-field (B*57:01) and higher-resolution extensions
    (B*57:01:01) but not sibling alleles (B*57:02)."""
    a = _strip_prefix(allele_name)
    return a == sentinel or a.startswith(sentinel + ":")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fixtures-dir", type=Path, default=Path("benchmarks/fixtures_1000g/"))
    ap.add_argument("--ground-truth", type=Path, default=Path("benchmarks/1000g_HLA_types.tsv"))
    ap.add_argument("--output-dir", type=Path, default=Path("benchmarks/1000g_supplementary/"))
    ap.add_argument("--tool", default="arcashla",
                    help="Representative tool to profile (default: arcashla).")
    ap.add_argument("--input-source", default="validated",
                    choices=["typing_tool", "validated", "simulated", "unknown"])
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--max-samples", type=int, default=0)
    args = ap.parse_args()

    from hlante.normalizer import load_imgt_db
    from hlante.annotator import AnnotatorConfig, build_clients
    from hlante.types import InputSource

    print("Loading IMGT-HLA database …", end=" ", flush=True)
    imgt_db = load_imgt_db(None)
    print(f"OK (version={imgt_db.get('version')})")

    gt = load_ground_truth(args.ground_truth)
    print(f"Ground truth: {len(gt)} samples; input-source={args.input_source}")

    input_src = InputSource(args.input_source)
    base_config = AnnotatorConfig(offline=False, population_group="global",
                                  input_source=input_src)
    base_clients = build_clients(base_config)

    ext = FIXTURE_EXT[args.tool]
    fixture_files = sorted((args.fixtures_dir / args.tool).glob(f"*{ext}"))
    if args.max_samples:
        fixture_files = fixture_files[: args.max_samples]
    print(f"[{args.tool}] {len(fixture_files)} fixtures")

    # R5: scope counter across all GWAS hits (deduplicated per allele record)
    scope_counter: Counter = Counter()
    # S2: sentinel -> region -> list of confidence scores; and tier counts
    sent_conf: Dict[str, Dict[str, List[float]]] = {
        s: defaultdict(list) for s in SENTINELS
    }
    sent_tier: Dict[str, Dict[str, Counter]] = {
        s: defaultdict(Counter) for s in SENTINELS
    }
    lock = threading.Lock()

    def _process(fpath: Path) -> None:
        stem = fpath.stem
        for suffix in (".genotype", ".t1k", ".hlahd", ".optitype"):
            if stem.endswith(suffix):
                stem = stem[: -len(suffix)]
                break
        gt_entry = gt.get(stem)
        if gt_entry is None:
            return
        region = gt_entry["region"]
        afnd_pop = REGION_TO_AFND.get(region, "global")
        cfg = AnnotatorConfig(offline=False, population_group=afnd_pop,
                              input_source=input_src)
        annotated = annotate_sample(fpath, args.tool, afnd_pop, imgt_db, base_clients, cfg)
        if not annotated:
            return

        local_scope: Counter = Counter()
        for ann in annotated:
            # R5: record the cascade scope of each GWAS hit
            for h in ann.gwas_hits:
                local_scope[getattr(h, "annotation_scope", "allele") or "allele"] += 1

        with lock:
            scope_counter.update(local_scope)
            # S2: for each sentinel the sample carries, record this allele's confidence
            for sentinel in SENTINELS:
                if not _carried(gt_entry["row"], sentinel):
                    continue
                for ann in annotated:
                    if _matches_sentinel(ann.normalized_allele.allele_name, sentinel):
                        sent_conf[sentinel][region].append(ann.input_quality_score or 0.0)
                        tier = getattr(ann, "input_quality_tier", "NA") or "NA"
                        sent_tier[sentinel][region][tier] += 1
                        break

    with ThreadPoolExecutor(max_workers=args.threads) as pool:
        futs = {pool.submit(_process, f): f for f in fixture_files}
        done = 0
        for fut in as_completed(futs):
            fut.result()
            done += 1
            if done % 500 == 0:
                print(f"  {done}/{len(fixture_files)} …", flush=True)

    args.output_dir.mkdir(parents=True, exist_ok=True)

    # ── R5: scope distribution ────────────────────────────────────────────────
    total_scope = sum(scope_counter.values())
    scope_path = args.output_dir / "gwas_scope_dist_1000g.tsv"
    # Report the known scopes in cascade order, then anything else the
    # annotator produced. Hard-coding the vocabulary once silently dropped a
    # whole tier: the one-field scope was renamed from "locus" to
    # "allele_group" and this loop kept writing a zero row for the old name
    # while the real counts vanished from the table.
    known = ("allele", "subtype", "allele_group")
    scopes = list(known) + sorted(set(scope_counter) - set(known))
    with scope_path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["scope", "count", "pct"])
        for scope in scopes:
            c = scope_counter.get(scope, 0)
            pct = f"{100 * c / total_scope:.1f}" if total_scope else "0.0"
            w.writerow([scope, c, pct])
    print(f"\nR5 → {scope_path}  (total GWAS hits: {total_scope})")
    for scope in scopes:
        c = scope_counter.get(scope, 0)
        pct = 100 * c / total_scope if total_scope else 0
        print(f"   {scope:<12} {c:>7}  {pct:5.1f}%")

    # ── S2: per-sentinel confidence by population ─────────────────────────────
    s2_path = args.output_dir / "sentinel_confidence_by_pop.tsv"
    with s2_path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sentinel", "population", "n_carriers", "mean_input_quality",
                    "min_input_quality", "tier_detailed", "tier_partial", "tier_limited"])
        for sentinel in SENTINELS:
            for pop in SUPERPOPS:
                scores = sent_conf[sentinel].get(pop, [])
                if not scores:
                    continue
                tiers = sent_tier[sentinel][pop]
                w.writerow([
                    sentinel, pop, len(scores),
                    f"{mean(scores):.3f}", f"{min(scores):.3f}",
                    tiers.get("detailed", 0), tiers.get("partial", 0),
                    tiers.get("limited", 0),
                ])
    print(f"S2 → {s2_path}")


if __name__ == "__main__":
    main()
