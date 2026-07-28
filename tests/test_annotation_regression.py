"""
tests.test_annotation_regression
================================

Regression tests against the real local database dumps.

Execution
---------
Skipped during default pytest runs. Run explicitly with:

.. code-block:: bash

    pytest -m integration tests/test_annotation_regression.py

Prerequisites
-------------
- ``hlante db-update --db imgt``
- ``hlante db-update --db pharmgkb``
- ``hlante db-update --db gwas``

Tests depending on a missing dump are skipped via :func:`pytest.skip`.

Expected findings
-----------------
- ``HLA-B*57:01`` → PharmGKB 1A abacavir (and secondarily
  flucloxacillin)
- ``HLA-B*57:01`` → at least one GWAS hit (HIV progression, abacavir
  hypersensitivity, or psoriasis)
- ``HLA-DRB1*04:01`` → at least one GWAS hit at coarse resolution
"""

from __future__ import annotations

from pathlib import Path
from typing import List

import pytest

from hlante.db.gwas import GWASClient, GWASDatabaseError
from hlante.db.pharmgkb import PharmGKBClient, PharmGKBDatabaseError

pytestmark = pytest.mark.integration


IMGT_LOCAL = Path.home() / ".hlante" / "imgt_hla"
PHARMGKB_LOCAL = Path.home() / ".hlante" / "pharmgkb"
GWAS_LOCAL = Path.home() / ".hlante" / "gwas"


def _require_pharmgkb() -> PharmGKBClient:
    if not (PHARMGKB_LOCAL / "clinical_annotations.tsv").is_file():
        pytest.skip(
            f"PharmGKB dump is not present ({PHARMGKB_LOCAL}). "
            "Run `hlante db-update --db pharmgkb`."
        )
    return PharmGKBClient(local_dir=PHARMGKB_LOCAL)


def _require_gwas() -> GWASClient:
    # Either the subset or the full TSV is acceptable.
    subset = GWAS_LOCAL / "gwas-hla-subset.tsv"
    full = GWAS_LOCAL / "gwas-catalog-download-associations-alt-full.tsv"
    if not (subset.is_file() or full.is_file()):
        pytest.skip(
            f"GWAS dump is not present ({GWAS_LOCAL}). "
            "Run `hlante db-update --db gwas` (~59 MB)."
        )
    return GWASClient(local_dir=GWAS_LOCAL)


# ---------------------------------------------------------------------------
# PharmGKB regression
# ---------------------------------------------------------------------------


class TestPharmGKBRegression:
    """
    Minimum expected hits against the real PharmGKB dump.
    """

    def test_b5701_abacavir_1a(self) -> None:
        """
        HLA-B*57:01 must have a 1A abacavir annotation (PharmGKB
        CA000001, Tier 1 VIP).
        """
        client = _require_pharmgkb()
        hits = client.query_allele("B*57:01")
        assert hits, "No PharmGKB hit returned for B*57:01 — parser regression."

        abacavir = [h for h in hits if "abacavir" in h.drug.lower()]
        assert abacavir, (
            f"No abacavir hit. Drugs returned: "
            f"{sorted({h.drug for h in hits})}"
        )
        evidence_levels = {h.evidence_level for h in abacavir}
        assert "1A" in evidence_levels, (
            f"Missing 1A evidence for abacavir: {evidence_levels}"
        )

        # Every record must carry at least one PMID.
        assert all(h.pmid for h in abacavir), \
            "PharmGKB annotations are missing PMIDs — evidence join failed."

    def test_b5701_flucloxacillin_1a(self) -> None:
        """
        HLA-B*57:01 must also carry a 1A flucloxacillin annotation
        (CA 1184996860).
        """
        client = _require_pharmgkb()
        hits = client.query_allele("B*57:01")
        flucloxacillin = [
            h for h in hits if "flucloxacillin" in h.drug.lower()
        ]
        assert flucloxacillin, (
            f"No flucloxacillin hit. Drugs returned: "
            f"{sorted({h.drug for h in hits})}"
        )
        assert any(h.evidence_level == "1A" for h in flucloxacillin)

    def test_b1502_carbamazepine(self) -> None:
        """
        HLA-B*15:02 must have a carbamazepine annotation (the classic
        SJS risk, CPIC guideline).
        """
        client = _require_pharmgkb()
        hits = client.query_allele("B*15:02")
        carbamazepine = [
            h for h in hits if "carbamazepine" in h.drug.lower()
        ]
        assert carbamazepine, (
            f"No carbamazepine hit. Drugs returned: "
            f"{sorted({h.drug for h in hits})}"
        )


# ---------------------------------------------------------------------------
# GWAS regression
# ---------------------------------------------------------------------------


class TestGWASRegression:
    """
    Minimum expected hits against the real GWAS Catalog bulk dump.
    """

    def test_b5701_has_hits(self) -> None:
        """
        HLA-B*57:01 must have at least one GWAS hit (abacavir HSR,
        HIV progression, psoriasis, β-2 microglobulin, etc.).
        """
        client = _require_gwas()
        hits = client.query_allele("B*57:01")
        assert hits, (
            "No GWAS hit for B*57:01. If access to "
            "'STRONGEST SNP-RISK ALLELE' is broken the parser regressed."
        )
        # Traits must be populated.
        assert all(h.trait for h in hits)

    def test_drb1_0401_via_fallback(self) -> None:
        """
        DRB1*04:01 may not have a direct exact-resolution hit (the
        STRONGEST SNP-RISK ALLELE column is often rs-ID based). The
        fallback must pick it up at ``4-field`` or ``2-field``.
        """
        client = _require_gwas()
        hits, resolution = client.query_allele_with_fallback("DRB1*04:01")
        assert hits, (
            f"No hit for DRB1*04:01 at any fallback level. "
            f"resolution={resolution}"
        )
        assert resolution in ("two-field", "one-field"), (
            f"Unexpected fallback resolution: {resolution}"
        )

    def test_b2705_via_fallback(self) -> None:
        """
        B*27:05 must have at least one GWAS hit at some fallback level.
        """
        client = _require_gwas()
        hits, resolution = client.query_allele_with_fallback("B*27:05")
        assert hits, (
            f"No hit for B*27:05 at any fallback level. resolution={resolution}"
        )
        assert resolution in ("two-field", "one-field"), (
            f"Unexpected fallback resolution: {resolution}"
        )
        # Traits must be populated.
        assert all(h.trait for h in hits)


# ---------------------------------------------------------------------------
# End-to-end (annotate_genotype)
# ---------------------------------------------------------------------------


class TestFullAnnotationRegression:
    """
    Real-world clinical classification through annotate_genotype.
    """

    def test_b5701_significance(self) -> None:
        """
        Against the real databases, B*57:01 should be classified as
        ``Risk factor`` (or ``Pathogenic``).
        """
        _require_pharmgkb()  # Without PharmGKB this test is meaningless

        from hlante.annotator import (
            AnnotatorConfig,
            AnnotatorClients,
            annotate_genotype,
            build_clients,
            SIGNIFICANCE_PATHOGENIC,
            SIGNIFICANCE_RISK_FACTOR,
        )
        from hlante.normalizer import NormalizedAllele

        # Enable GWAS only if it is available.
        gwas_available = (
            (GWAS_LOCAL / "gwas-hla-subset.tsv").is_file()
            or (
                GWAS_LOCAL
                / "gwas-catalog-download-associations-alt-full.tsv"
            ).is_file()
        )

        norm = NormalizedAllele(
            allele_name="B*57:01",
            imgt_accession=None,
            protein_group=None,
            hla_class="I",
            gene="HLA-B",
            resolution_level=2,
            is_ambiguous=False,
            is_novel=False,
            sample_id="regression",
            source_tool="debug",
            source_locus="HLA-B",
            source_resolution="two-field",
            allele_index=0,
        )
        config = AnnotatorConfig(
            offline=True,
            pharmgkb_local_dir=PHARMGKB_LOCAL,
            gwas_local_dir=GWAS_LOCAL if gwas_available else None,
            enable_gwas=gwas_available,
            enable_afnd=False,
        )
        clients = build_clients(config)
        result = annotate_genotype([norm], config, clients=clients)[0]

        assert result.pharm_annotations, "Missing PharmGKB annotation."
        assert result.clinical_significance in (
            SIGNIFICANCE_RISK_FACTOR,
            SIGNIFICANCE_PATHOGENIC,
        ), (
            f"Expected Risk factor / Pathogenic, got: "
            f"{result.clinical_significance}"
        )
        assert "abacavir" in result.drug_response_summary.lower()
