#!/usr/bin/env python3
"""
Per-super-population AFND coverage under each of the two frequency snapshots.

Coverage is the proportion of annotated allele calls that received a population
allele-frequency record for the sample's super-population — the quantity in the
``AFND coverage`` column of the main benchmark's ``per_population_metrics.tsv``.
It is a property of the installed snapshot rather than of the tool, and this
script measures that directly by annotating the same cohort twice:

``installed``
    the table under ``~/.hlante/afnd`` (the full release fetched by
    ``hlante db-update --db afnd``);
``builtin``
    the compiled-in fallback that lets HLAnte run before any download, forced by
    pointing the annotator at a directory holding no AFND table.

The accumulation reproduces ``run_annotation_benchmark.py``: arcasHLA input,
``--input-source validated``, one AFND population per sample derived from its
1000 Genomes super-population. Running the ``installed`` arm therefore
reproduces the committed per-population coverage figures, which is what makes
the ``builtin`` arm comparable to them.

Usage
-----
    python scripts/benchmark/afnd_snapshot_coverage.py \\
        --ground-truth benchmarks/1000g_HLA_types.tsv \\
        --fixtures-dir <fixtures>/arcashla \\
        --output-dir benchmarks/1000g_afnd_snapshots
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
import tempfile
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Dict, List, Optional

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from run_annotation_benchmark import (  # noqa: E402
    REGION_TO_AFND,
    annotate_sample,
    load_ground_truth,
)

from hlante.annotator import AnnotatorConfig, build_clients  # noqa: E402
from hlante.db.afnd import BUILTIN_AFND_TSV, DEFAULT_LOCAL_DIR  # noqa: E402
from hlante.normalizer import load_imgt_db  # noqa: E402
from hlante.types import InputSource  # noqa: E402

#: Report the super-populations in the order used by main-text Table 4.
REGION_ORDER = ["EUR", "EAS", "AMR", "SAS", "AFR"]


def _pct(num: int, den: int) -> str:
    return f"{100.0 * num / den:.1f}%" if den else "N/A"


def _table_rows(path: Optional[Path]) -> int:
    """Data rows in a frequency table, header excluded."""
    if path is None or not path.is_file():
        return 0
    with path.open(encoding="utf-8") as handle:
        return max(0, sum(1 for _ in handle) - 1)


def measure(
    label: str,
    afnd_local_dir: Optional[Path],
    fixtures: List[Path],
    ground_truth: Dict[str, Dict[str, Any]],
    imgt_db: Dict[str, Any],
    threads: int,
) -> Dict[str, Dict[str, Any]]:
    """Annotate the cohort once and tally coverage per super-population."""
    base_config = AnnotatorConfig(
        offline=True,
        input_source=InputSource("validated"),
        afnd_local_dir=afnd_local_dir,
    )
    clients = build_clients(base_config)

    counts: Dict[str, Dict[str, int]] = defaultdict(
        lambda: {"n_samples": 0, "n_alleles": 0, "n_freq": 0}
    )

    def _process(path: Path) -> None:
        sample_id = path.stem
        if sample_id.endswith(".genotype"):
            sample_id = sample_id[: -len(".genotype")]
        entry = ground_truth.get(sample_id)
        if entry is None:
            return
        region = entry["region"]
        config = AnnotatorConfig(
            offline=True,
            population_group=REGION_TO_AFND.get(region, "global"),
            input_source=InputSource("validated"),
            afnd_local_dir=afnd_local_dir,
        )
        annotated = annotate_sample(
            path, "arcashla", config.population_group, imgt_db, clients, config
        )
        if not annotated:
            return
        bucket = counts[region]
        bucket["n_samples"] += 1
        bucket["n_alleles"] += len(annotated)
        bucket["n_freq"] += sum(1 for a in annotated if a.allele_frequency is not None)

    if threads > 1:
        with ThreadPoolExecutor(max_workers=threads) as pool:
            for future in as_completed({pool.submit(_process, f): f for f in fixtures}):
                future.result()
    else:
        for path in fixtures:
            _process(path)

    result: Dict[str, Dict[str, Any]] = {}
    for region in REGION_ORDER:
        bucket = counts.get(region, {"n_samples": 0, "n_alleles": 0, "n_freq": 0})
        result[region] = {
            "n_samples": bucket["n_samples"],
            "n_alleles": bucket["n_alleles"],
            "n_with_frequency": bucket["n_freq"],
            "coverage": _pct(bucket["n_freq"], bucket["n_alleles"]),
        }
    print(f"  {label}: " + "  ".join(
        f"{r}={result[r]['coverage']}" for r in REGION_ORDER
    ))
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ground-truth", type=Path, required=True)
    parser.add_argument(
        "--fixtures-dir",
        type=Path,
        required=True,
        help="Directory of arcasHLA fixtures (…/arcashla).",
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=8)
    args = parser.parse_args()

    ground_truth = load_ground_truth(args.ground_truth)
    fixtures = sorted(args.fixtures_dir.glob("*.json"))
    if not fixtures:
        raise SystemExit(f"no arcasHLA fixtures under {args.fixtures_dir}")
    print(f"{len(fixtures)} fixture(s); {len(ground_truth)} ground-truth sample(s)")

    imgt_db = load_imgt_db()

    installed_tsv = DEFAULT_LOCAL_DIR / "afnd_frequencies.tsv"
    snapshots: Dict[str, Dict[str, Any]] = {}

    # An empty directory holds no AFND table, so the client falls back to the
    # compiled-in one — the state a laboratory is in before its first db-update.
    with tempfile.TemporaryDirectory(prefix="hlante-afnd-empty-") as empty:
        for label, local_dir, source in (
            ("installed", None, installed_tsv),
            ("builtin", Path(empty), BUILTIN_AFND_TSV),
        ):
            snapshots[label] = {
                "source": str(source),
                "rows": _table_rows(source),
                "per_population": measure(
                    label, local_dir, fixtures, ground_truth, imgt_db, args.threads
                ),
            }

    args.output_dir.mkdir(parents=True, exist_ok=True)
    rows = [
        {
            "snapshot": label,
            "table_rows": snapshots[label]["rows"],
            "super_population": region,
            "n_samples": values["n_samples"],
            "annotated_alleles": values["n_alleles"],
            "alleles_with_frequency": values["n_with_frequency"],
            "afnd_coverage": values["coverage"],
        }
        for label in ("installed", "builtin")
        for region, values in snapshots[label]["per_population"].items()
    ]
    tsv_path = args.output_dir / "afnd_snapshot_coverage.tsv"
    with tsv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    (args.output_dir / "afnd_snapshot_coverage.json").write_text(
        json.dumps(snapshots, indent=2) + "\n", encoding="utf-8"
    )
    print(f"written → {tsv_path}")


if __name__ == "__main__":
    main()
