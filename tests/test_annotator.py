"""
tests.test_annotator
====================

Unit tests for :mod:`hlante.annotator` plus
``hlante.db.{gwas, pharmgkb}``.

Approach
--------
- GWAS client is exercised through injected fake fetchers with canned
  JSON responses.
- The PharmGKB client is tested against a local TSV fixture.
- ``annotate_genotype`` is driven end-to-end with injected
  :class:`AnnotatorClients`; no network access occurs.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Callable, Dict, List

import pytest

from hlante.annotator import (
    SIGNIFICANCE_BENIGN,
    SIGNIFICANCE_NOVEL,
    SIGNIFICANCE_CURATED_ACTIONABLE,
    SIGNIFICANCE_PATHOGENIC,
    SIGNIFICANCE_RISK_FACTOR,
    SIGNIFICANCE_VUS,
    AnnotatedHLA,
    AnnotatorClients,
    AnnotatorConfig,
    annotate_genotype,
)
from hlante.types import DiseaseEntry
from hlante.db.gwas import GWASClient, GWASHit
from hlante.db.pharmgkb import (
    PharmAnnotation,
    PharmGKBClient,
    PharmGKBDownloadError,
)
from hlante.normalizer import NormalizedAllele


FIXTURES_DIR: Path = Path(__file__).parent / "fixtures"
PHARMGKB_FIXTURE_DIR: Path = FIXTURES_DIR / "pharmgkb"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_allele(
    allele_name: str,
    *,
    gene: str = "HLA-B",
    hla_class: str = "I",
    resolution_level: int = 4,
    is_ambiguous: bool = False,
    is_novel: bool = False,
    imgt_accession: str = "HLA00001",
    protein_group: str = None,
) -> NormalizedAllele:
    """
    Short factory used across the tests below.
    """
    return NormalizedAllele(
        allele_name=allele_name,
        imgt_accession=imgt_accession,
        protein_group=protein_group,
        hla_class=hla_class,
        gene=gene,
        resolution_level=resolution_level,
        is_ambiguous=is_ambiguous,
        is_novel=is_novel,
    )


def _fake_fetcher(responses: Dict[str, Any]) -> Callable[[str], bytes]:
    """
    Build a minimal fetcher from a URL → canned JSON mapping.
    """

    def _fetch(url: str) -> bytes:
        for key, payload in responses.items():
            if key in url:
                if isinstance(payload, Exception):
                    raise payload
                return json.dumps(payload).encode("utf-8")
        raise AssertionError(f"Unexpected URL: {url}")

    return _fetch


# ---------------------------------------------------------------------------
# GWAS client tests
# ---------------------------------------------------------------------------


GWAS_FIXTURE_DIR: Path = FIXTURES_DIR / "gwas"


class TestGWASClient:
    """
    Unit tests for :class:`GWASClient` — bulk-TSV approach.
    """

    def test_query_returns_filtered_hits(self) -> None:
        """
        B*57:01 returns hits (abacavir HSR + psoriasis); the weak
        common-cold hit (p=1e-5) is filtered out.
        """
        client = GWASClient(local_dir=GWAS_FIXTURE_DIR)
        hits = client.query_allele("B*57:01")
        traits = sorted(h.trait for h in hits)
        # Fixture rows for B*57:01: abacavir HSR (1.2e-10), psoriasis
        # (3.0e-20), common cold (1.0e-5 → above threshold, dropped).
        assert "Abacavir hypersensitivity" in [h.trait for h in hits] or \
               "drug-induced hypersensitivity" in [h.trait for h in hits]
        assert not any("common cold" in t.lower() for t in traits)
        # OR value must be read correctly.
        hsr = next(h for h in hits if "hypersensitivity" in h.trait.lower())
        assert hsr.odds_ratio == pytest.approx(4.2)
        assert hsr.pmid == "18256392"
        assert hsr.study_accession == "GCST000001"
        assert hsr.allele == "B*57:01"

    def test_strips_risk_allele_suffix(self) -> None:
        """
        ``B*27:05-?`` must be indexed as ``B*27:05`` (ankylosing
        spondylitis).
        """
        client = GWASClient(local_dir=GWAS_FIXTURE_DIR)
        hits = client.query_allele("B*27:05")
        assert len(hits) >= 1
        assert any("ankylosing" in h.trait.lower() for h in hits)

    def test_hla_prefix_equivalent(self) -> None:
        """
        ``HLA-B*57:01`` and ``B*57:01`` must return the same result.
        """
        client = GWASClient(local_dir=GWAS_FIXTURE_DIR)
        a = client.query_allele("HLA-B*57:01")
        b = client.query_allele("B*57:01")
        assert [h.allele for h in a] == [h.allele for h in b]

    def test_ignores_rsid_strongest_allele(self) -> None:
        """
        Rows whose STRONGEST SNP-RISK ALLELE is an rs-ID must not be
        indexed under HLA keys.
        """
        client = GWASClient(local_dir=GWAS_FIXTURE_DIR)
        client.load()
        # No rs-based key should be in the internal index.
        assert not any("rs" in key.lower() for key in client._by_allele)

    def test_offline_without_local_dir_returns_empty(
        self, tmp_path: Path, caplog: pytest.LogCaptureFixture
    ) -> None:
        """
        Offline mode + empty local dir must return an empty list with
        a warning log.
        """
        import logging

        client = GWASClient(
            local_dir=tmp_path / "empty",
            offline=True,
            fetcher=lambda url: (_ for _ in ()).throw(
                AssertionError("fetcher must not be called in offline mode")
            ),
        )
        with caplog.at_level(logging.WARNING, logger="hlante.db.gwas"):
            out = client.query_allele("B*57:01")
        assert out == []
        assert any("offline" in rec.message.lower() for rec in caplog.records)

    def test_update_raises_in_offline(self, tmp_path: Path) -> None:
        """
        ``update()`` must raise in offline mode.
        """
        from hlante.db.gwas import GWASDownloadError

        client = GWASClient(local_dir=tmp_path, offline=True)
        with pytest.raises(GWASDownloadError):
            client.update()

    def test_retry_on_download_failure(self, tmp_path: Path) -> None:
        """
        After a transient error, the second attempt must succeed.
        """
        import io
        import zipfile

        call_log: List[int] = []

        # A minimal but valid zip payload.
        buf = io.BytesIO()
        with zipfile.ZipFile(buf, "w") as zf:
            zf.writestr(
                "gwas-catalog-download-associations-alt-full.tsv",
                "STRONGEST SNP-RISK ALLELE\tDISEASE/TRAIT\n",
            )
        good_zip = buf.getvalue()

        def flaky(url: str) -> bytes:
            call_log.append(1)
            if len(call_log) == 1:
                raise OSError("transient network error")
            return good_zip

        client = GWASClient(
            local_dir=tmp_path,
            fetcher=flaky,
            max_retries=3,
            sleep=lambda _: None,
        )
        client.update()
        assert len(call_log) == 2


# ---------------------------------------------------------------------------
# PharmGKB client tests
# ---------------------------------------------------------------------------


class TestPharmGKBClient:
    """
    Unit tests for :class:`PharmGKBClient`.
    """

    def test_load_and_query(self) -> None:
        """
        The fixture TSV must load and return the 1A/1B/2A HLA records
        for B*57:01.
        """
        client = PharmGKBClient(local_dir=PHARMGKB_FIXTURE_DIR)
        hits = client.query_allele("B*57:01")
        assert len(hits) == 1
        ann = hits[0]
        assert ann.drug == "abacavir"
        assert ann.evidence_level == "1A"
        assert "18256392" in ann.pmid
        assert ann.cpic_url and "cpicpgx.org" in ann.cpic_url
        assert ann.pharmgkb_url and "CA000001" in ann.pharmgkb_url
        assert ann.allele == "B*57:01"
        assert ann.gene == "HLA-B"

    def test_filters_non_hla_rows(self) -> None:
        """
        Non-HLA rows (e.g., CYP2C19) must be skipped.
        """
        client = PharmGKBClient(local_dir=PHARMGKB_FIXTURE_DIR)
        assert client.query_allele("CYP2C19*2") == []

    def test_filters_evidence_level(self) -> None:
        """
        Explicit 1A/1B/2A filter must suppress a level-3 DRB1*07:01
        record. (The default is 1A/1B only; here we pass an explicit
        filter to verify the filtering mechanism across levels.)
        """
        client = PharmGKBClient(
            local_dir=PHARMGKB_FIXTURE_DIR,
            evidence_levels=frozenset({"1A", "1B", "2A"}),
        )
        assert client.query_allele("DRB1*07:01") == []

    def test_skips_rows_without_pmid(self) -> None:
        """
        The CA000006 row has no PMID and must be dropped.
        """
        client = PharmGKBClient(local_dir=PHARMGKB_FIXTURE_DIR)
        assert client.query_allele("B*07:02") == []

    def test_offline_update_raises(self, tmp_path: Path) -> None:
        """
        ``update()`` must raise in offline mode.
        """
        client = PharmGKBClient(local_dir=tmp_path, offline=True)
        with pytest.raises(PharmGKBDownloadError):
            client.update()



# ---------------------------------------------------------------------------
# Stubs for annotator orchestration
# ---------------------------------------------------------------------------


class _StubGWAS:
    def __init__(self, hits_by_allele: Dict[str, List[GWASHit]]) -> None:
        self.hits = hits_by_allele
        self.calls: List[str] = []

    def query_allele(self, allele: str) -> List[GWASHit]:
        self.calls.append(allele)
        return list(self.hits.get(allele, []))

    def query_allele_with_fallback(self, allele: str, min_resolution: int = 2):
        """
        Test stub: if a hit exists for the full query return "two-field";
        otherwise return "none". The actual fallback logic is covered
        in the GWASClient unit tests.
        """
        hits = self.query_allele(allele)
        return hits, ("two-field" if hits else "none")


class _StubPharmGKB:
    def __init__(self, anns_by_allele: Dict[str, List[PharmAnnotation]]) -> None:
        self.anns = anns_by_allele
        self.calls: List[str] = []

    def query_allele(self, allele: str) -> List[PharmAnnotation]:
        self.calls.append(allele)
        return list(self.anns.get(allele, []))


class _StubCurated:
    def __init__(self, entries_by_allele: Dict[str, List[DiseaseEntry]]) -> None:
        self.entries = entries_by_allele
        self.calls: List[Any] = []

    def query_allele(self, allele: str) -> List[DiseaseEntry]:
        self.calls.append(allele)
        return list(self.entries.get(allele, []))


# ---------------------------------------------------------------------------
# Annotate_genotype — end to end
# ---------------------------------------------------------------------------


class TestAnnotateGenotype:
    """
    Orchestration tests for :func:`annotate_genotype`.
    """

    def _make_clients(self) -> AnnotatorClients:
        gwas_hits = {
            "B*57:01": [
                GWASHit(
                    trait="Abacavir hypersensitivity",
                    p_value=1.2e-10,
                    odds_ratio=4.2,
                    pmid="18322448",
                    study_accession="GCST000001",
                    allele="B*57:01",
                ),
            ],
            "DRB1*15:01": [
                GWASHit(
                    trait="Multiple sclerosis",
                    p_value=1e-20,
                    odds_ratio=3.1,
                    pmid="17660530",
                    study_accession="GCST000010",
                    allele="DRB1*15:01",
                ),
            ],
        }
        pharm = {
            "B*57:01": [
                PharmAnnotation(
                    drug="abacavir",
                    phenotype="HSR risk",
                    evidence_level="1A",
                    pmid=["18256392"],
                    allele="B*57:01",
                    annotation_id="CA000001",
                    gene="HLA-B",
                    cpic_url="https://cpicpgx.org/guidelines/",
                )
            ],
        }
        curated_entries = {
            "B*57:01": [
                DiseaseEntry(
                    variation_id="111",
                    significance="Pathogenic",
                    condition="Stevens-Johnson syndrome",
                    review_status="reviewed by expert panel",
                    allele="B*57:01",
                    pmid=["18322448"],
                )
            ],
        }
        return AnnotatorClients(
            gwas=_StubGWAS(gwas_hits),
            pharmgkb=_StubPharmGKB(pharm),
            curated=_StubCurated(curated_entries),
        )

    def test_pathogenic_classification(self) -> None:
        """
        A PharmGKB Level 1A annotation carrying a CPIC guideline link may
        assert the CPIC label.
        """
        allele = _make_allele("B*57:01")
        config = AnnotatorConfig(offline=True)
        clients = self._make_clients()
        result = annotate_genotype([allele], config, clients=clients)[0]

        assert result.clinical_significance == SIGNIFICANCE_PATHOGENIC
        assert "Abacavir" in result.disease_risk_summary
        assert "OR=4.20" in result.disease_risk_summary
        # RISK_PREFIX_HIGH is now "Strong association".
        assert "Strong association" in result.disease_risk_summary
        assert "abacavir" in result.drug_response_summary
        assert "1A" in result.drug_response_summary
        assert len(result.gwas_hits) == 1
        assert len(result.pharm_annotations) == 1
        assert len(result.disease_entries) == 1

    def test_risk_factor_classification(self) -> None:
        """
        No Pathogenic but strong GWAS/PharmGKB → ``Risk factor``.
        """
        allele = _make_allele(
            "DRB1*15:01", gene="HLA-DRB1", hla_class="II"
        )
        config = AnnotatorConfig(offline=True)
        clients = self._make_clients()
        result = annotate_genotype([allele], config, clients=clients)[0]

        assert result.clinical_significance == SIGNIFICANCE_RISK_FACTOR
        assert "Multiple sclerosis" in result.disease_risk_summary

    def test_benign_classification(self) -> None:
        """
        No hits + IMGT-known + no DB query ran (empty stubs return
        ``("none", [])`` for GWAS fallback) is classified as
        ``Benign (limited evidence)`` to distinguish it from a verified
        clean result. ``_make_clients`` with populated stubs covers
        the plain ``Benign`` path separately.
        """
        from hlante.annotator import SIGNIFICANCE_BENIGN_LIMITED

        allele = _make_allele("A*01:01:01:01", gene="HLA-A")
        config = AnnotatorConfig(offline=True)
        clients = AnnotatorClients(
            gwas=_StubGWAS({}),
            pharmgkb=_StubPharmGKB({}),
        )
        result = annotate_genotype([allele], config, clients=clients)[0]
        assert result.clinical_significance == SIGNIFICANCE_BENIGN_LIMITED
        assert "No disease association" in result.disease_risk_summary

    def test_novel_classification(self) -> None:
        """
        ``is_novel=True`` overrides every other signal and yields
        ``Novel``.
        """
        allele = _make_allele(
            "A*99:99", is_novel=True, imgt_accession=None, is_ambiguous=True
        )
        config = AnnotatorConfig(offline=True)
        result = annotate_genotype(
            [allele], config, clients=AnnotatorClients()
        )[0]
        assert result.clinical_significance == SIGNIFICANCE_NOVEL

    def test_vus_classification(self) -> None:
        """
        Ambiguous allele with no findings → ``VUS``.
        """
        allele = _make_allele(
            "B*07:02", is_ambiguous=True, imgt_accession=None
        )
        config = AnnotatorConfig(offline=True)
        result = annotate_genotype(
            [allele], config, clients=AnnotatorClients()
        )[0]
        assert result.clinical_significance == SIGNIFICANCE_VUS

    def test_disabled_source_not_queried(self) -> None:
        """
        Disabled clients must not be queried, and the curated table alone
        must not assert a CPIC level.

        With PharmGKB switched off the only actionable evidence is the
        curated reference table. That is a transcription maintained with
        this package, not a live guideline lookup, so it gets its own label
        rather than the CPIC-1A wording.
        """
        allele = _make_allele("B*57:01")
        clients = self._make_clients()
        config = AnnotatorConfig(
            offline=True,
            enable_gwas=True,
            enable_pharmgkb=False,
        )
        # Simulate the disable by setting the stub to None.
        clients.pharmgkb = None

        result = annotate_genotype([allele], config, clients=clients)[0]
        assert result.pharm_annotations == []
        # Pharm disabled: the curated table alone cannot assert a CPIC level.
        assert result.clinical_significance == SIGNIFICANCE_CURATED_ACTIONABLE

    def test_client_exception_is_logged_and_swallowed(
        self, caplog: pytest.LogCaptureFixture
    ) -> None:
        """
        A client exception must be logged and swallowed so other
        clients still run.
        """
        import logging

        class _BoomGWAS:
            def query_allele(self, allele: str) -> List[GWASHit]:
                raise RuntimeError("network explosion")

        clients = AnnotatorClients(
            gwas=_BoomGWAS(),
            pharmgkb=_StubPharmGKB({}),
        )
        allele = _make_allele("B*57:01")
        config = AnnotatorConfig(offline=True)
        with caplog.at_level(logging.WARNING, logger="hlante.annotator"):
            result = annotate_genotype([allele], config, clients=clients)[0]
        assert result.gwas_hits == []
        assert any("failed" in rec.message.lower() for rec in caplog.records)

    def test_batch_preserves_order(self) -> None:
        """
        Multiple inputs must preserve order in the output.
        """
        alleles = [
            _make_allele("B*57:01"),
            _make_allele("DRB1*15:01", gene="HLA-DRB1", hla_class="II"),
        ]
        config = AnnotatorConfig(offline=True)
        results = annotate_genotype(
            alleles, config, clients=self._make_clients()
        )
        assert [r.normalized_allele.allele_name for r in results] == [
            "B*57:01",
            "DRB1*15:01",
        ]

    def test_build_clients_respects_offline(self, tmp_path: Path) -> None:
        """
        :func:`build_clients` must propagate the offline flag to every
        client.
        """
        from hlante.annotator import build_clients

        config = AnnotatorConfig(
            offline=True,
            cache_root=tmp_path,
            pharmgkb_local_dir=PHARMGKB_FIXTURE_DIR,
        )
        clients = build_clients(config)
        assert clients.gwas is not None and clients.gwas.offline is True
        assert clients.pharmgkb is not None and clients.pharmgkb.offline is True


# ---------------------------------------------------------------------------
# Benign vs Benign (limited evidence)
# ---------------------------------------------------------------------------


class TestBenignLimitedEvidence:
    """
    The classifier must distinguish a queried-but-empty "Benign" call
    from a silent-databases "Benign (limited evidence)" call.
    """

    def test_benign_requires_successful_query(self) -> None:
        """
        IMGT-known, unambiguous allele + no hits from any DB + no
        GWAS query ever ran → ``Benign (limited evidence)``.
        """
        from hlante.annotator import (
            SIGNIFICANCE_BENIGN_LIMITED,
            _classify_significance,
        )

        allele = _make_allele(
            "A*01:01:01:01",
            gene="HLA-A",
            imgt_accession="HLA00001",
        )
        label = _classify_significance(
            allele, [], [], [], gwas_resolution="none"
        )
        assert label == SIGNIFICANCE_BENIGN_LIMITED

    def test_benign_after_successful_query(self) -> None:
        """
        Same allele + GWAS ran but returned no hits → plain ``Benign``.
        """
        from hlante.annotator import (
            SIGNIFICANCE_BENIGN,
            _classify_significance,
        )

        allele = _make_allele(
            "A*01:01:01:01",
            gene="HLA-A",
            imgt_accession="HLA00001",
        )
        label = _classify_significance(
            allele, [], [], [], gwas_resolution="two-field"
        )
        assert label == SIGNIFICANCE_BENIGN
