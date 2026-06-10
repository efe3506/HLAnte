"""
tests.test_qa_full_panel
========================

Drive every allele in :data:`tests.fixtures.qa_allele_panel.QA_PANEL`
through the full HLAnte pipeline
(normalize → annotate) and produce a correctness report.

Run
---
Explicit, does not run in default pytest pass:

.. code-block:: bash

    pytest -m qa tests/test_qa_full_panel.py -s --no-header -q
    # or directly:
    python -m tests.test_qa_full_panel

Output artifacts
----------------
On completion, two files are written in ``docs/``:

- ``docs/QA_PANEL_RESULTS.md`` — Markdown table for human review.
- ``docs/QA_PANEL_RESULTS.json`` — machine-readable structured log.

Default local databases are used
(``~/.hlante/imgt_hla``, ``~/.hlante/pharmgkb``, ``~/.hlante/gwas``,
``~/.hlante/afnd``). Missing databases are reported as pipeline
warnings, not fatal, so the panel can be run incrementally.
"""

from __future__ import annotations

import json
import logging
import sys
import traceback
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional

import pytest

from hlante.annotator import (
    AnnotatedHLA,
    AnnotatorClients,
    AnnotatorConfig,
    annotate_genotype,
    build_clients,
)
from hlante.normalizer import (
    IMGTDatabaseMissingError,
    InvalidAlleleError,
    NormalizedAllele,
    load_imgt_db,
    normalize_allele,
)

from tests.fixtures.qa_allele_panel import QA_PANEL


pytestmark = pytest.mark.qa


DOCS_DIR: Path = Path(__file__).resolve().parent.parent / "docs"
MD_OUT: Path = DOCS_DIR / "QA_PANEL_RESULTS.md"
JSON_OUT: Path = DOCS_DIR / "QA_PANEL_RESULTS.json"


# ---------------------------------------------------------------------------
# Result row model
# ---------------------------------------------------------------------------


@dataclass
class PanelRow:
    allele: str
    status: str = "PENDING"
    result: str = ""
    imgt_accession: Optional[str] = None
    is_novel: Optional[bool] = None
    is_ambiguous: Optional[bool] = None
    resolution_level: Optional[int] = None
    gwas_hits_count: int = 0
    pharm_hits_count: int = 0
    disease_hits_count: int = 0
    confidence_score: Optional[float] = None
    gwas_annotation_resolution: Optional[str] = None
    clinical_significance: Optional[str] = None
    gwas_traits_snippet: List[str] = field(default_factory=list)
    pharm_drugs_snippet: List[str] = field(default_factory=list)
    issue: str = ""
    note: str = ""
    expected: Dict[str, Any] = field(default_factory=dict)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _build_norm(allele_name: str, imgt_db: Dict[str, Any]) -> Optional[NormalizedAllele]:
    """
    Attempt to normalize; return :data:`None` for null tokens,
    :class:`NormalizedAllele` otherwise. Raises :class:`InvalidAlleleError`
    for malformed input.
    """
    return normalize_allele(allele_name, imgt_db)


def _synthetic_normalized(allele_name: str) -> NormalizedAllele:
    """
    Fallback normalization used when IMGT DB is unavailable or the
    allele could not be looked up — produces a coarse NormalizedAllele
    so the annotator path can still be exercised.
    """
    gene_letter = allele_name.split("*", 1)[0] if "*" in allele_name else allele_name
    return NormalizedAllele(
        allele_name=allele_name,
        imgt_accession=None,
        protein_group=None,
        hla_class="I" if not gene_letter.startswith(("DR", "DQ", "DP", "DM", "DO")) else "II",
        gene=f"HLA-{gene_letter}",
        resolution_level=4,
        is_ambiguous=True,
        is_novel=True,
        sample_id="qa_panel",
        source_tool="qa",
        source_locus=f"HLA-{gene_letter}",
        source_resolution="4-field",
        allele_index=0,
    )


def _evaluate_expectations(
    row: PanelRow, result: Optional[AnnotatedHLA]
) -> None:
    """
    Compare actual annotation against expected fields on the panel entry.
    """
    issues: List[str] = []

    exp_drug = row.expected.get("expect_drug")
    exp_ev = row.expected.get("expect_ev")
    if exp_drug is not None:
        pharm_drugs = [p.drug.lower() for p in (result.pharm_annotations if result else [])]
        if not any(exp_drug.lower() in d for d in pharm_drugs):
            issues.append(f"Expected pharm drug '{exp_drug}' not found")
        elif exp_ev:
            matching = [
                p for p in result.pharm_annotations
                if exp_drug.lower() in p.drug.lower()
            ]
            evs = {p.evidence_level for p in matching}
            if exp_ev not in evs:
                issues.append(f"Expected '{exp_drug}' evidence {exp_ev}, got {evs}")
    elif exp_drug is None and "expect_drug" in row.expected:
        # Explicitly expect NO pharm drug — any strong hit is a flag
        strong = [
            p for p in (result.pharm_annotations if result else [])
            if (p.evidence_level or "").upper() in {"1A", "1B"}
        ]
        if strong:
            issues.append(
                f"Unexpected strong PharmGKB hit: {[p.drug for p in strong]}"
            )

    exp_trait = row.expected.get("expect_trait")
    if exp_trait is not None:
        traits = [h.trait.lower() for h in (result.gwas_hits if result else [])]
        if not any(exp_trait.lower() in t for t in traits):
            issues.append(f"Expected GWAS trait '{exp_trait}' not found")
    elif exp_trait is None and "expect_trait" in row.expected:
        # Optional — if explicitly None, flag any GWAS hits
        if result and result.gwas_hits:
            # Informational only; not an automatic issue
            pass

    if issues:
        row.status = "FAIL"
        row.issue = "; ".join(issues)
    else:
        row.status = "PASS"


def _snippet(items: List[str], max_items: int = 3) -> List[str]:
    return items[:max_items]


# ---------------------------------------------------------------------------
# Main runner
# ---------------------------------------------------------------------------


def run_panel() -> List[PanelRow]:
    """
    Process every panel entry; populate :class:`PanelRow` for each.
    """
    logging.basicConfig(level=logging.WARNING)
    rows: List[PanelRow] = []

    try:
        imgt_db = load_imgt_db()
    except IMGTDatabaseMissingError as exc:
        print(f"[WARN] IMGT DB missing: {exc}", file=sys.stderr)
        imgt_db = None

    config = AnnotatorConfig(offline=True)
    try:
        clients = build_clients(config)
    except Exception as exc:  # pragma: no cover - defensive
        print(f"[WARN] Could not build clients: {exc}", file=sys.stderr)
        clients = AnnotatorClients()

    for entry in QA_PANEL:
        row = PanelRow(
            allele=entry["allele"],
            note=entry.get("note", ""),
            expected={k: v for k, v in entry.items() if k.startswith("expect_")},
        )

        # Step 1 — normalize
        norm: Optional[NormalizedAllele]
        try:
            if imgt_db is not None:
                norm = _build_norm(entry["allele"], imgt_db)
            else:
                norm = _synthetic_normalized(entry["allele"])
        except InvalidAlleleError as exc:
            row.status = "MALFORMED"
            row.issue = f"InvalidAlleleError: {exc}"
            rows.append(row)
            continue
        except Exception as exc:
            row.status = "CRASH"
            row.issue = f"{type(exc).__name__}: {exc}\n{traceback.format_exc(limit=2)}"
            rows.append(row)
            continue

        if norm is None:
            # Null token — do not annotate
            row.status = "NULL_TOKEN"
            row.result = "normalize_allele returned None (null token handled)"
            rows.append(row)
            continue

        row.imgt_accession = norm.imgt_accession
        row.is_novel = norm.is_novel
        row.is_ambiguous = norm.is_ambiguous
        row.resolution_level = norm.resolution_level

        # Attach minimal provenance for annotator's fallback
        norm.sample_id = "qa_panel"
        norm.source_tool = "qa"
        norm.source_locus = norm.gene
        norm.source_resolution = f"{norm.resolution_level}-field"
        norm.allele_index = 0

        # Step 2 — annotate
        try:
            results = annotate_genotype([norm], config, clients=clients)
        except Exception as exc:
            row.status = "CRASH"
            row.issue = f"annotate_genotype raised: {exc}"
            rows.append(row)
            continue

        if not results:
            row.status = "NO_RESULT"
            row.issue = "annotate_genotype returned empty list"
            rows.append(row)
            continue

        ann = results[0]
        row.gwas_hits_count = len(ann.gwas_hits)
        row.pharm_hits_count = len(ann.pharm_annotations)
        row.disease_hits_count = len(ann.disease_entries)
        row.confidence_score = ann.confidence_score
        row.gwas_annotation_resolution = ann.gwas_resolution_used
        row.clinical_significance = ann.clinical_significance
        row.gwas_traits_snippet = _snippet([h.trait for h in ann.gwas_hits])
        row.pharm_drugs_snippet = _snippet([
            f"{p.drug}({p.evidence_level})" for p in ann.pharm_annotations
        ])

        _evaluate_expectations(row, ann)
        rows.append(row)

    return rows


def write_markdown(rows: List[PanelRow], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)

    lines: List[str] = []
    lines.append("# HLAnte — QA Panel Results\n")
    lines.append(f"Total panel size: **{len(rows)}**\n")

    status_counts: Dict[str, int] = {}
    for r in rows:
        status_counts[r.status] = status_counts.get(r.status, 0) + 1
    lines.append("## Status summary\n")
    for status, n in sorted(status_counts.items(), key=lambda x: -x[1]):
        lines.append(f"- **{status}**: {n}")
    lines.append("")

    lines.append("## Per-allele results\n")
    lines.append(
        "| Allele | Status | IMGT | Novel | Ambig | GWAS | Pharm | Disease | Conf | GWAS-res | Sig | Issue / Note |"
    )
    lines.append(
        "|--------|--------|------|-------|-------|------|-------|---------|------|----------|-----|--------------|"
    )
    for r in rows:
        imgt = (r.imgt_accession or "")[:10]
        novel = "-" if r.is_novel is None else ("T" if r.is_novel else "F")
        ambig = "-" if r.is_ambiguous is None else ("T" if r.is_ambiguous else "F")
        conf = f"{r.confidence_score:.3f}" if r.confidence_score is not None else "-"
        gres = r.gwas_annotation_resolution or "-"
        sig = (r.clinical_significance or "-").replace("|", "/")
        msg = r.issue or r.note or ""
        msg = msg.replace("|", "/")[:80]
        lines.append(
            f"| {r.allele} | {r.status} | {imgt} | {novel} | {ambig} | "
            f"{r.gwas_hits_count} | {r.pharm_hits_count} | {r.disease_hits_count} | "
            f"{conf} | {gres} | {sig} | {msg} |"
        )

    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_json(rows: List[PanelRow], out_path: Path) -> None:
    payload = [asdict(r) for r in rows]
    out_path.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8"
    )


# ---------------------------------------------------------------------------
# pytest integration — does not run by default (needs -m qa)
# ---------------------------------------------------------------------------


def test_qa_full_panel() -> None:
    rows = run_panel()
    write_markdown(rows, MD_OUT)
    write_json(rows, JSON_OUT)

    # Verbose print — captured only with ``-s``
    print(f"\nQA panel complete. Results:\n  {MD_OUT}\n  {JSON_OUT}")
    for row in rows:
        if row.status in {"FAIL", "CRASH"}:
            print(f"  [{row.status}] {row.allele:25s} {row.issue}")

    # Do not fail the test on panel issues — this is a reporting run.
    # The findings are written to docs/QA_PANEL_RESULTS.{md,json} for review.


if __name__ == "__main__":
    rows = run_panel()
    write_markdown(rows, MD_OUT)
    write_json(rows, JSON_OUT)
    print(f"Wrote {MD_OUT}")
    print(f"Wrote {JSON_OUT}")
