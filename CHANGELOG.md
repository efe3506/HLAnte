# Changelog

This file follows the [Keep a Changelog](https://keepachangelog.com/en/1.1.0/)
format and the project uses [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Changed
- **Caller-reported quality is now surfaced, and stopped being mislabelled.**
  `HLAGenotype.quality_score` held three different quantities depending on the
  tool — T1K abundance, the OptiType solution objective, or nothing — and was
  discarded before reaching any output. It is replaced by
  `caller_quality1` / `caller_quality2`, populated only from T1K's native
  per-allele quality columns, and surfaced as the `caller_allele_quality`
  column. OptiType inputs carrying more than one enumerated solution now emit
  a warning: only the top-ranked solution is annotated.
- **The composite confidence score is renamed and reframed.** `confidence_score`
  / `confidence_tier` / `confidence_rationale` become `input_quality_score` /
  `input_quality_tier` / `input_quality_rationale`, and the tiers
  `HIGH` / `MODERATE` / `LOW` become `detailed` / `partial` / `limited`. The
  score combines reference-data completeness with characteristics of the
  submitted call; it is not a measure of genotype accuracy and not a posterior
  probability, and the previous naming invited exactly that reading. Values are
  unchanged — only the names. Committed benchmark tables are relabelled
  accordingly.
- **Resolution is now reported as a field count, not a digit count.** Current
  HLA nomenclature admits at most four colon-separated fields, so the previous
  labels (`2-field`, `4-field`, `6-field`, `8-field`) conflated fields with
  digits. The `resolution`, `gwas_annotation_resolution` and
  `confidence_rationale` columns now read `one-field` … `four-field`, and
  `NormalizedAllele.resolution_level` is 1–4 rather than 2/4/6/8. Word forms
  were chosen so a filter written against v0.1.0 fails loudly instead of
  silently selecting the wrong resolution.
- `--resolution` now takes a field count (1–4). The digit-scale values `6` and
  `8` are still accepted and translated with a warning; passing `2` or `4`
  warns that the scale changed, because those values are valid on both scales
  but no longer mean the same thing.

### Fixed
- `db-update --db imgt --imgt-ref 3.64.0` failed with HTTP 404: the ANHIG
  mirror names its release branches without separators (`3640`). Release
  numbers are now translated automatically, so the documented command pins the
  release as intended.
- Allele normalisation aborted with an error naming an internal Python
  function; the message now gives the `hlante db-update --db imgt` command.
- The integration test fixture requested separate stderr capture in a way that
  only worked on Click 8.3+, so a fresh install on Python 3.9 (Click 8.1)
  reported three failures.

### Removed
- The bundled NMDP allele-frequency table. NMDP data are licensed by
  NMDP/Be The Match and may not be redistributed; the NMDP source is now
  inactive unless the user supplies their own extract.

## [0.1.0] - 2026-04-21

### Added
- Initial public release.
- Click-based CLI (`annotate`, `validate`, `db-update`, `version`).
- Parsers for ARCAS-HLA, T1K, HLA-HD, and OptiType outputs.
- IPD-IMGT/HLA nomenclature normalizer with G / P group handling.
- Resolution-aware fallback lookup for GWAS, PharmGKB, and AFND.
- GWAS Catalog bulk-TSV ingestion with HLA subset caching.
- PharmGKB clinical-annotation parser with
  ``clinical_ann_evidence.tsv`` PMID join (CPIC Level 1A/1B by default).
- AFND allele-frequency client with a universal geographic taxonomy
  (EUR / AFR / EAS (alias ASN) / SAS / MID / AMR / OCE / global). No
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
