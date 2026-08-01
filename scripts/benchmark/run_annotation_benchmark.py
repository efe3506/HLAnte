#!/usr/bin/env python3
"""
Run the HLAnte annotation benchmark on 1000G fixture files.

Uses the HLAnte Python API directly (no subprocess) for speed:
  - IMGT DB is loaded once and shared across all samples.
  - DB clients are built once per population group configuration.

Metrics computed
----------------
  3a  Parser fidelity          — alleles parsed vs expected from ground truth
  3b  Normalisation success    — % alleles with imgt_accession != NA
  3c  CPIC Level 1A recall     — for carriers of B*57:01, B*58:01, B*15:02, A*31:01
  3d  GWAS / curated recall    — for carriers of key disease alleles
  3e  Population-stratified frequency coverage
  3f  Confidence tier distribution

Usage:
    python scripts/benchmark/run_annotation_benchmark.py \\
        --fixtures-dir benchmark/fixtures_1000g/ \\
        --ground-truth benchmark/1000g_HLA_types.tsv \\
        --output-dir benchmark/results/ \\
        --tools all \\
        --threads 8 \\
        --max-samples 0
"""

import argparse
import csv
import json
import logging
import sys
import threading
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional

logger = logging.getLogger(__name__)

LOCI = ["A", "B", "C", "DQB1", "DRB1"]
CLASS_I_LOCI = ["A", "B", "C"]

ALL_TOOLS = ["arcashla", "t1k", "hlahd", "optitype"]

# Fixture file extension per tool
FIXTURE_EXT = {
    "arcashla": ".genotype.json",
    "t1k":      ".t1k.tsv",
    "hlahd":    ".hlahd.txt",
    "optitype": ".optitype.tsv",
}

# Tools that provide Class I only (limits expected allele count)
CLASS_I_ONLY_TOOLS = {"optitype"}

# Super-population code mapping: 1000G → HLAnte AFND code
REGION_TO_AFND: Dict[str, str] = {
    "AFR": "AFR",
    "AMR": "AMR",
    "EAS": "ASN",
    "EUR": "EUR",
    "SAS": "SAS",
}

# ── CPIC Level 1A recall targets ─────────────────────────────────────────────
# Allele → expected drug name (case-insensitive substring in the pharm_drugs /
# drug_response_summary columns) that HLAnte should return for a carrier.
CPIC_TARGETS: Dict[str, str] = {
    "B*57:01": "abacavir",
    "B*58:01": "allopurinol",
    "B*15:02": "carbamazepine",
    "A*31:01": "carbamazepine",
}

# ── GWAS / curated disease recall targets ────────────────────────────────────
# Allele → list of keyword substrings (ANY match counts as a hit)
GWAS_TARGETS: Dict[str, List[str]] = {
    "DRB1*03:01": ["lupus", "diabetes", "type 1", "t1d", "sle"],
    "DRB1*04:01": ["arthritis", "rheumatoid"],
    "DRB1*15:01": ["multiple sclerosis", "sclerosis"],
    "B*27:05":    ["ankylosing", "spondylitis"],
    "B*51:01":    ["behcet", "behçet"],
    "DQB1*02:01": ["celiac", "coeliac"],
    "DQB1*06:02": ["narcolepsy"],
}


def _null(v: str) -> bool:
    return v.strip().lower() in {"", "*", "-", "na", "not typed", "nottyped", "."}


def _carried(row: Dict[str, str], allele_key: str) -> bool:
    gene, fields = allele_key.split("*", 1)
    for slot in (1, 2):
        if row.get(f"HLA-{gene} {slot}", "").strip() == fields:
            return True
    return False


def _allele_count_gt(row: Dict[str, str], loci: List[str]) -> int:
    """Count non-null allele slots in the ground-truth row for the given loci."""
    count = 0
    for gene in loci:
        for slot in (1, 2):
            val = row.get(f"HLA-{gene} {slot}", "")
            if not _null(val):
                count += 1
    return count


# ── Ground-truth loader ───────────────────────────────────────────────────────

def load_ground_truth(path: Path) -> Dict[str, Dict[str, Any]]:
    """
    Returns {sample_id: {"region": str, "population": str, "row": dict}}.

    Samples where every allele slot across all loci is null are excluded —
    they carry no typing data and no fixture files exist for them.
    """
    gt: Dict[str, Dict[str, Any]] = {}
    with path.open(encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            sid = row["Sample ID"].strip()
            if all(
                _null(row.get(f"HLA-{gene} {slot}", ""))
                for gene in LOCI
                for slot in (1, 2)
            ):
                logger.debug("Skipping fully-untyped sample: %s", sid)
                continue
            gt[sid] = {
                "region":     row["Region"].strip(),
                "population": row["Population"].strip(),
                "row":        dict(row),
            }
    return gt


# ── Annotation helpers ────────────────────────────────────────────────────────

def _has_drug_hit(annotated_hla_list: List[Any], drug: str) -> bool:
    """
    Return True if any AnnotatedHLA in the list has the drug in its
    pharm_annotations OR mentions it in disease_entries conditions.
    """
    drug_lower = drug.lower()
    for ann in annotated_hla_list:
        for p in ann.pharm_annotations:
            if p.drug and drug_lower in p.drug.lower():
                return True
        for e in ann.disease_entries:
            if e.condition and drug_lower in e.condition.lower():
                return True
    return False


def _has_trait_hit(annotated_hla_list: List[Any], keywords: List[str]) -> bool:
    """
    Return True if any AnnotatedHLA has any keyword in gwas_hits traits
    OR disease_entries conditions.
    """
    for ann in annotated_hla_list:
        for h in ann.gwas_hits:
            if any(kw in (h.trait or "").lower() for kw in keywords):
                return True
        for e in ann.disease_entries:
            if any(kw in (e.condition or "").lower() for kw in keywords):
                return True
    return False


# ── Per-sample annotation ─────────────────────────────────────────────────────

def annotate_sample(
    fixture_path: Path,
    tool: str,
    population: str,
    imgt_db: Dict[str, Any],
    clients: Any,
    config: Any,
) -> Optional[List[Any]]:
    """
    Parse → normalize → annotate one fixture file.
    Returns list of AnnotatedHLA or None on failure.
    """
    from hlante.parser import parse_hla_output, HLAnteParseError
    from hlante.normalizer import batch_normalize, InvalidAlleleError
    from hlante.annotator import annotate_genotype

    try:
        genotypes = parse_hla_output(fixture_path, tool)
    except (HLAnteParseError, FileNotFoundError) as exc:
        logger.debug("Parse failed %s: %s", fixture_path.name, exc)
        return None

    try:
        normalized = batch_normalize(genotypes, imgt_db=imgt_db)
    except InvalidAlleleError as exc:
        logger.debug("Normalize failed %s: %s", fixture_path.name, exc)
        return None

    if not normalized:
        return None

    try:
        return annotate_genotype(normalized, config, clients=clients)
    except Exception as exc:  # noqa: BLE001
        logger.debug("Annotate failed %s: %s", fixture_path.name, exc)
        return None


# ── Metrics aggregation ───────────────────────────────────────────────────────

class ToolMetrics:
    """Accumulates per-tool benchmark metrics."""

    def __init__(self, tool: str) -> None:
        self.tool = tool
        self.n_samples = 0
        self.alleles_expected = 0
        self.alleles_parsed = 0
        self.alleles_normalized = 0
        self.alleles_with_freq = 0
        self.tier_counts: Counter = Counter()
        self.cpic_carriers: Counter = Counter()
        self.cpic_hits: Counter = Counter()
        self.gwas_carriers: Counter = Counter()
        self.gwas_hits: Counter = Counter()
        self.failed_samples: List[str] = []
        self._lock = threading.Lock()

    def add(
        self,
        sample_id: str,
        gt_row: Dict[str, str],
        annotated: Optional[List[Any]],
        n_expected: int,
        n_parsed: int,
    ) -> None:
        with self._lock:
            self.n_samples += 1
            self.alleles_expected += n_expected
            self.alleles_parsed += n_parsed

            if annotated is None:
                self.failed_samples.append(sample_id)
                return

            for ann in annotated:
                # Count as normalised when IMGT recognises the allele (exact or
                # prefix match).  Two-field inputs always have imgt_accession=None
                # (ambiguous) but is_novel=False — both are "normalised" in the
                # sense that IMGT knows the allele.
                if not ann.normalized_allele.is_novel:
                    self.alleles_normalized += 1
                if ann.allele_frequency is not None:
                    self.alleles_with_freq += 1
                tier = getattr(ann, "input_quality_tier", "NA") or "NA"
                self.tier_counts[tier] += 1

            # CPIC recall
            for allele, drug in CPIC_TARGETS.items():
                if _carried(gt_row, allele):
                    self.cpic_carriers[allele] += 1
                    if _has_drug_hit(annotated, drug):
                        self.cpic_hits[allele] += 1

            # GWAS / curated recall
            for allele, keywords in GWAS_TARGETS.items():
                if _carried(gt_row, allele):
                    self.gwas_carriers[allele] += 1
                    if _has_trait_hit(annotated, keywords):
                        self.gwas_hits[allele] += 1


# ── Report writers ────────────────────────────────────────────────────────────

def _pct(num: int, denom: int) -> str:
    if denom == 0:
        return "N/A"
    return f"{100 * num / denom:.1f}%"


def write_tsv(rows: List[Dict[str, str]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("# no data\n", encoding="utf-8")
        return
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_markdown_summary(
    metrics_by_tool: Dict[str, ToolMetrics],
    pop_metrics: Dict[str, Dict[str, Any]],
    hlante_version: str,
    imgt_version: str,
    n_samples_total: int,
    output_path: Path,
    input_source: str = "validated",
) -> None:
    lines: List[str] = []
    lines.append("# HLAnte 1000 Genomes Annotation Benchmark\n")
    lines.append(f"Generated: {datetime.now(timezone.utc).date().isoformat()}  ")
    lines.append(f"HLAnte version: {hlante_version}  IMGT-HLA: {imgt_version}  ")
    lines.append(f"N samples: {n_samples_total}\n")

    lines.append("## Scope\n")
    lines.append(
        "This benchmark evaluates HLAnte's annotation pipeline using "
        "Sanger-validated HLA types from the 1000 Genomes Project (Abi-Rached et al. 2018) "
        "as input. HLAnte is an annotation tool, not a typing tool; "
        "typing accuracy is out of scope.\n"
    )

    # Parser fidelity
    lines.append("## Parser Fidelity\n")
    lines.append("| Tool | Alleles parsed | Expected | Success rate |")
    lines.append("|------|----------------|----------|--------------|")
    for tool in ALL_TOOLS:
        m = metrics_by_tool.get(tool)
        if m is None:
            continue
        lines.append(
            f"| {tool:<9} | {m.alleles_parsed:>14} | {m.alleles_expected:>8} "
            f"| {_pct(m.alleles_parsed, m.alleles_expected):>12} |"
        )

    # Normalisation
    lines.append(
        "\n## Normalisation Success (IMGT-HLA recognised; exact or prefix match)\n"
    )
    lines.append(
        "Note: two-field input alleles have no exact IMGT accession (ambiguous "
        "by design) but are counted as recognised when a prefix match exists "
        "(`is_novel=False`).  Novel alleles are those with no IMGT prefix.\n"
    )
    lines.append("| Tool | Allele count | IMGT recognised | Rate |")
    lines.append("|------|--------------|-----------------|------|")
    for tool in ALL_TOOLS:
        m = metrics_by_tool.get(tool)
        if m is None:
            continue
        lines.append(
            f"| {tool:<9} | {m.alleles_parsed:>12} | {m.alleles_normalized:>15} "
            f"| {_pct(m.alleles_normalized, m.alleles_parsed):>4} |"
        )

    # CPIC recall
    lines.append("\n## CPIC Level 1A Pharmacogenomic Recall\n")
    lines.append("Results shown for arcashla (representative; other tools similar).")
    lines.append("\n| Drug-allele | Carriers | Hits | Recall |")
    lines.append("|-------------|----------|------|--------|")
    ref = metrics_by_tool.get("arcashla") or next(iter(metrics_by_tool.values()), None)
    if ref:
        for allele, drug in CPIC_TARGETS.items():
            carriers = ref.cpic_carriers[allele]
            hits = ref.cpic_hits[allele]
            lines.append(
                f"| {allele} → {drug:<14} | {carriers:>8} | {hits:>4} "
                f"| {_pct(hits, carriers):>6} |"
            )

    # GWAS recall
    lines.append("\n## GWAS / Curated Disease Recall\n")
    lines.append("| Allele → Trait | Carriers | Hits | Recall |")
    lines.append("|----------------|----------|------|--------|")
    if ref:
        LABELS = {
            "DRB1*03:01": "SLE / T1D",
            "DRB1*04:01": "rheumatoid arthritis",
            "DRB1*15:01": "multiple sclerosis",
            "B*27:05":    "ankylosing spondylitis",
            "B*51:01":    "Behcet disease",
            "DQB1*02:01": "celiac disease",
            "DQB1*06:02": "narcolepsy",
        }
        for allele, label in LABELS.items():
            carriers = ref.gwas_carriers[allele]
            hits = ref.gwas_hits[allele]
            lines.append(
                f"| {allele} → {label:<22} | {carriers:>8} | {hits:>4} "
                f"| {_pct(hits, carriers):>6} |"
            )

    # Population-stratified performance (arcashla)
    lines.append("\n## Population-Stratified Performance (arcashla)\n")
    lines.append(
        "| Super-pop | N | AFND coverage | Mean input quality | limited tier % |"
    )
    lines.append("|-----------|---|---------------|-----------------|------------|")
    for pop in ["EUR", "AFR", "EAS", "SAS", "AMR"]:
        pm = pop_metrics.get(pop, {})
        lines.append(
            f"| {pop:<9} | {pm.get('n', 0):>3} "
            f"| {pm.get('freq_pct', 'N/A'):>13} "
            f"| {pm.get('mean_input_quality', 'N/A'):>15} "
            f"| {pm.get('low_pct', 'N/A'):>10} |"
        )

    # Input-quality tier distribution. The note has to follow the mode this run
    # actually used: printing the "validated" explanation above a typing_tool
    # table made the committed summary contradict its own numbers.
    lines.append("\n## Input-Quality Tier Distribution (arcashla)\n")
    if input_source == "validated":
        lines.append(
            "Note: this run used --input-source validated, because the 1000 "
            "Genomes types are Sanger-validated reference data. The ambiguity "
            "penalty (×0.75) is suppressed under that mode, leaving the "
            "two-field resolution penalty (×0.90), so most calls land in the "
            "detailed tier.\n"
        )
    else:
        lines.append(
            "Note: this run used --input-source typing_tool (the default). "
            "Two-field calls take the ambiguity penalty (×0.75) on top of the "
            "two-field resolution penalty (×0.90), so nearly all of them land "
            "in the limited tier. That is the intended conservative behaviour "
            "for unvalidated input, not an annotation failure.\n"
        )
    lines.append("| Tier | Count | Percentage |")
    lines.append("|------|-------|------------|")
    if ref:
        total_alleles = sum(ref.tier_counts.values())
        for tier in ("detailed", "partial", "limited", "NA"):
            count = ref.tier_counts.get(tier, 0)
            lines.append(
                f"| {tier:<8} | {count:>5} | {_pct(count, total_alleles):>10} |"
            )

    lines.append("\n---\n")
    lines.append(
        "> **Scope note**: This benchmark validates annotation pipeline "
        "behaviour, not typing accuracy. HLAnte does not perform HLA typing."
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
        "--fixtures-dir", type=Path, default=Path("benchmark/fixtures_1000g/"),
        help="Root directory containing per-tool fixture subdirectories.",
    )
    parser.add_argument(
        "--ground-truth", type=Path, default=Path("benchmark/1000g_HLA_types.tsv"),
        help="1000G HLA TSV file used as ground truth.",
    )
    parser.add_argument(
        "--output-dir", type=Path, default=Path("benchmark/results/"),
        help="Directory for benchmark result files.",
    )
    parser.add_argument(
        "--tools", default="all",
        help="Comma-separated tools or 'all' (default: all).",
    )
    parser.add_argument(
        "--threads", type=int, default=4,
        help="Parallel worker threads per tool.",
    )
    parser.add_argument(
        "--max-samples", type=int, default=0,
        help="Maximum samples to process per tool (0 = all).",
    )
    parser.add_argument(
        "--imgt-db-path", type=Path, default=None,
        help="Local IMGT-HLA directory (default: ~/.hlante/imgt_hla).",
    )
    parser.add_argument(
        "--offline", action="store_true",
        help="Disable live HTTP requests; use local dumps only.",
    )
    parser.add_argument(
        "--input-source",
        default="validated",
        choices=["typing_tool", "validated", "simulated", "unknown"],
        help=(
            "Provenance of the HLA allele calls. "
            "Default: 'validated' (1000G data is Sanger-typed). "
            "Use 'typing_tool' for arcasHLA/T1K/HLA-HD/OptiType output."
        ),
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

    if not args.ground_truth.is_file():
        logger.error("Ground truth not found: %s", args.ground_truth)
        sys.exit(1)

    tools_to_run = (
        ALL_TOOLS
        if args.tools.strip().lower() == "all"
        else [t.strip().lower() for t in args.tools.split(",")]
    )

    # ── Lazy HLAnte imports (after path check) ────────────────────────────────
    try:
        from hlante import __version__ as hlante_version
        from hlante.normalizer import load_imgt_db, IMGTDatabaseMissingError
        from hlante.annotator import AnnotatorConfig, build_clients
        from hlante.types import InputSource
    except ImportError as exc:
        logger.error("Cannot import hlante: %s", exc)
        sys.exit(1)

    # Load IMGT DB once
    print("Loading IMGT-HLA database …", end=" ", flush=True)
    try:
        imgt_db = load_imgt_db(args.imgt_db_path)
        imgt_version = imgt_db.get("version") or "unknown"
        print(f"OK (version={imgt_version})")
    except IMGTDatabaseMissingError as exc:
        print(f"FAILED: {exc}")
        print(
            "Run 'hlante db-update --db imgt' first, or pass --imgt-db-path.",
            file=sys.stderr,
        )
        sys.exit(1)

    gt = load_ground_truth(args.ground_truth)
    print(f"Ground truth: {len(gt)} samples")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    failed_log = (args.output_dir / "failed_samples.log").open("w", encoding="utf-8")

    all_metrics: Dict[str, ToolMetrics] = {}

    # Per-population tracking for arcashla (representative tool)
    pop_allele_counts: Dict[str, Counter] = defaultdict(Counter)  # Region → {stat: count}

    for tool in tools_to_run:
        tool_dir = args.fixtures_dir / tool
        if not tool_dir.is_dir():
            print(f"[{tool}] fixture directory not found: {tool_dir} — SKIPPED")
            continue

        ext = FIXTURE_EXT[tool]
        fixture_files = sorted(tool_dir.glob(f"*{ext}"))
        if args.max_samples and args.max_samples > 0:
            fixture_files = fixture_files[: args.max_samples]

        print(f"[{tool}] {len(fixture_files)} fixture files …", end=" ", flush=True)

        # Build clients for the default population; per-sample population is
        # handled by building a lightweight config overriding population_group.
        # For efficiency, we build one global config and override population per
        # Sample via a separate config instance (cheap to construct).
        input_src = InputSource(args.input_source)
        base_config = AnnotatorConfig(
            offline=args.offline,
            population_group="global",
            input_source=input_src,
        )
        base_clients = build_clients(base_config)

        m = ToolMetrics(tool)
        loci_for_tool = CLASS_I_LOCI if tool in CLASS_I_ONLY_TOOLS else LOCI

        def _process(fpath: Path) -> None:
            # Derive sample_id from filename by stripping known extensions
            stem = fpath.stem
            for suffix in (".genotype", ".t1k", ".hlahd", ".optitype"):
                if stem.endswith(suffix):
                    stem = stem[: -len(suffix)]
                    break
            sample_id = stem

            gt_entry = gt.get(sample_id)
            if gt_entry is None:
                logger.debug("Sample %s not in ground truth — skipped", sample_id)
                return

            region = gt_entry["region"]
            afnd_pop = REGION_TO_AFND.get(region, "global")

            # Build a sample-specific config only if population differs from base
            cfg = AnnotatorConfig(
                offline=args.offline,
                population_group=afnd_pop,
                input_source=input_src,
            )

            n_expected = _allele_count_gt(gt_entry["row"], loci_for_tool)
            annotated = annotate_sample(fpath, tool, afnd_pop, imgt_db, base_clients, cfg)

            n_parsed = 0
            if annotated is not None:
                n_parsed = len(annotated)

            m.add(sample_id, gt_entry["row"], annotated, n_expected, n_parsed)

            if annotated is None:
                failed_log.write(f"{tool}\t{sample_id}\tparse_or_annotate_failed\n")

            # Population-level tracking (arcashla only to avoid duplication)
            if tool == "arcashla" and annotated:
                c = pop_allele_counts[region]
                c["n_alleles"] += len(annotated)
                c["n_freq"] += sum(
                    1 for a in annotated if a.allele_frequency is not None
                )
                conf_sum = sum(a.input_quality_score or 0 for a in annotated)
                c["conf_sum"] += conf_sum
                c["n_low"] += sum(
                    1 for a in annotated
                    if getattr(a, "input_quality_tier", "NA") == "limited"
                )
                c["n_samples"] = c.get("n_samples", 0) + 1

        if args.threads > 1:
            with ThreadPoolExecutor(max_workers=args.threads) as pool:
                futures = {pool.submit(_process, f): f for f in fixture_files}
                done = 0
                for fut in as_completed(futures):
                    fut.result()
                    done += 1
                    if done % 200 == 0:
                        print(f"\n  [{tool}] {done}/{len(fixture_files)} …", end="", flush=True)
        else:
            for i, f in enumerate(fixture_files):
                _process(f)
                if (i + 1) % 200 == 0:
                    print(f"\n  [{tool}] {i+1}/{len(fixture_files)} …", end="", flush=True)

        all_metrics[tool] = m
        rate = _pct(m.alleles_parsed, m.alleles_expected)
        print(f" done  (parsed={m.alleles_parsed}/{m.alleles_expected} = {rate})")

    failed_log.close()

    # ── Build population metrics dict ─────────────────────────────────────────
    pop_metrics: Dict[str, Dict[str, Any]] = {}
    for region, c in pop_allele_counts.items():
        n = c.get("n_alleles", 0)
        pop_metrics[region] = {
            "n":        c.get("n_samples", 0),
            "freq_pct": _pct(c.get("n_freq", 0), n),
            "mean_input_quality": (
                f"{c.get('conf_sum', 0) / n:.3f}" if n else "N/A"
            ),
            "low_pct": _pct(c.get("n_low", 0), n),
        }

    # ── Write TSV outputs ─────────────────────────────────────────────────────
    # Per-tool metrics TSV
    tool_rows = []
    for tool, m in all_metrics.items():
        tool_rows.append({
            "tool":               tool,
            "n_samples":          m.n_samples,
            "alleles_expected":   m.alleles_expected,
            "alleles_parsed":     m.alleles_parsed,
            "parse_rate":         _pct(m.alleles_parsed, m.alleles_expected),
            "alleles_normalized": m.alleles_normalized,
            "norm_rate":          _pct(m.alleles_normalized, m.alleles_parsed),
            "alleles_with_freq":  m.alleles_with_freq,
            "freq_rate":          _pct(m.alleles_with_freq, m.alleles_parsed),
            "tier_detailed":          m.tier_counts.get("detailed", 0),
            "tier_partial":      m.tier_counts.get("partial", 0),
            "tier_limited":           m.tier_counts.get("limited", 0),
            "tier_NA":            m.tier_counts.get("NA", 0),
        })
    write_tsv(tool_rows, args.output_dir / "per_tool_metrics.tsv")

    # CPIC recall detail TSV
    cpic_rows = []
    for tool, m in all_metrics.items():
        for allele, drug in CPIC_TARGETS.items():
            carriers = m.cpic_carriers[allele]
            hits = m.cpic_hits[allele]
            cpic_rows.append({
                "tool": tool, "allele": allele, "expected_drug": drug,
                "carriers": carriers, "hits": hits,
                "recall": _pct(hits, carriers),
            })
    write_tsv(cpic_rows, args.output_dir / "pharm_recall_detail.tsv")

    # GWAS recall detail TSV
    gwas_rows = []
    GWAS_LABELS = {
        "DRB1*03:01": "SLE/T1D", "DRB1*04:01": "RA",
        "DRB1*15:01": "MS",      "B*27:05":    "AS",
        "B*51:01":    "Behcet",  "DQB1*02:01": "celiac",
        "DQB1*06:02": "narcolepsy",
    }
    for tool, m in all_metrics.items():
        for allele, keywords in GWAS_TARGETS.items():
            carriers = m.gwas_carriers[allele]
            hits = m.gwas_hits[allele]
            gwas_rows.append({
                "tool": tool, "allele": allele,
                "trait": GWAS_LABELS.get(allele, ""),
                "keywords": "|".join(keywords),
                "carriers": carriers, "hits": hits,
                "recall": _pct(hits, carriers),
            })
    write_tsv(gwas_rows, args.output_dir / "gwas_recall_detail.tsv")

    # Per-population metrics TSV (arcashla)
    pop_rows = [
        {
            "region": region,
            **{k: v for k, v in pm.items()},
        }
        for region, pm in pop_metrics.items()
    ]
    write_tsv(pop_rows, args.output_dir / "per_population_metrics.tsv")

    # Confidence distribution TSV (arcashla)
    ref_m = all_metrics.get("arcashla") or (next(iter(all_metrics.values()), None))
    if ref_m:
        total_a = sum(ref_m.tier_counts.values())
        conf_rows = [
            {
                "tier": tier,
                "count": ref_m.tier_counts.get(tier, 0),
                "pct": _pct(ref_m.tier_counts.get(tier, 0), total_a),
            }
            for tier in ("detailed", "partial", "limited", "NA")
        ]
        write_tsv(conf_rows, args.output_dir / "input_quality_distribution.tsv")

    # JSON summary
    summary = {
        "hlante_version": hlante_version,
        "imgt_version":   imgt_version,
        "generated_at":   datetime.now(timezone.utc).isoformat(),
        "n_samples_total": len(gt),
        "tools": {
            tool: {
                "n_samples":          m.n_samples,
                "parse_rate":         _pct(m.alleles_parsed, m.alleles_expected),
                "norm_rate":          _pct(m.alleles_normalized, m.alleles_parsed),
                "cpic_recall":        {
                    a: _pct(m.cpic_hits[a], m.cpic_carriers[a])
                    for a in CPIC_TARGETS
                },
                "gwas_recall":        {
                    a: _pct(m.gwas_hits[a], m.gwas_carriers[a])
                    for a in GWAS_TARGETS
                },
            }
            for tool, m in all_metrics.items()
        },
        "population_metrics": pop_metrics,
    }
    (args.output_dir / "benchmark_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )

    # Markdown summary
    write_markdown_summary(
        all_metrics,
        pop_metrics,
        hlante_version,
        imgt_version,
        len(gt),
        args.output_dir / "benchmark_summary.md",
        input_source=args.input_source,
    )

    print(f"\nResults written to {args.output_dir}/")
    print("  benchmark_summary.md")
    print("  per_tool_metrics.tsv")
    print("  pharm_recall_detail.tsv")
    print("  gwas_recall_detail.tsv")
    print("  input_quality_distribution.tsv")
    print("  failed_samples.log")


if __name__ == "__main__":
    main()
