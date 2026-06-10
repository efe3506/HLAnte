# Changelog

This file follows the [Keep a Changelog](https://keepachangelog.com/en/1.1.0/)
format and the project uses [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Click-based CLI (`annotate`, `validate`, `db-update`, `version`).
- Parsers for ARCAS-HLA, T1K, HLA-HD, and OptiType outputs.
- IPD-IMGT/HLA nomenclature normalizer with G / P group handling.
- Resolution-aware fallback lookup (6 → 4 → 2 field) for GWAS,
  PharmGKB, and AFND.
- GWAS Catalog bulk-TSV ingestion with HLA subset caching.
- PharmGKB clinical-annotation parser with
  ``clinical_ann_evidence.tsv`` PMID join (CPIC Level 1A/1B by default).
- AFND allele-frequency client with a universal geographic taxonomy
  (EUR / AFR / ASN (alias EAS) / SAS / MID / AMR / OCE / global). No
  country-level aliasing.
- Confidence score (0.0–1.0) combining novelty, frequency,
  resolution, and ambiguity signals.
- TSV / CSV / Markdown / JSON reporters with a research-use
  disclaimer header.
- Bioconda-ready `environment.yml`, PEP 621 `pyproject.toml`.

### Changed
- All user-facing text and developer-facing text are in English; no
  Turkish strings remain in the codebase.
- ``AFNDClient`` population matching is purely geographic; the
  previous country-specific alias was removed.

## [0.1.0] - 2026-04-21

### Added
- Initial public release.
