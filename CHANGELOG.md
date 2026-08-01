# Changelog

This file follows the [Keep a Changelog](https://keepachangelog.com/en/1.1.0/)
format and the project uses [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- **The AFND snapshot can be pinned, and is now identified in every report.**
  The mirror was fetched from its moving `main` branch and left nothing behind:
  no ref, no checksum, no entry in the report's `db_versions`. A published
  figure that depends on population frequencies could not be traced to the
  table that produced it. `db-update --db afnd --afnd-ref <sha>` now pins a
  commit, a `version.json` beside the installed table records the source URL,
  the ref, the acquisition date and the file's SHA-256, and report headers
  carry an `afnd` entry naming the ref and digest.

### Fixed
- **`--imgt-ref` did not replace an installed release, and `version.json`
  claimed it had.** Asking for a release other than the one on disk left the
  cached files in place — the download was skipped because the files existed —
  while the metadata was rewritten with the requested ref. A snapshot recorded
  as `ref: 3650` could therefore hold 3.64.0 data together with 3.64.0
  checksums, so the record used to identify a historical annotation described a
  snapshot that was never installed. A ref that differs from the recorded one
  now refreshes the files without needing `--force`, and the acquisition
  timestamp is preserved when nothing is fetched. An installation carrying no
  `version.json` is re-fetched rather than described.

## [0.2.0] - 2026-07-29

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
- **Pipe-joined columns dropped missing entries instead of reserving the slot,
  so positionally aligned columns fell out of step.** Fifteen columns joined
  their values with a helper that skips empties while the columns they are
  meant to align with did not, which silently shifted every value after the
  gap. The columns now all reserve the slot with `NA`. Measured on the
  2,692-sample 1000 Genomes cohort:
  - `gwas_odds_ratios` — a GWAS Catalog record with no effect size shortened
    the list, so every later odds ratio attached to the wrong trait. 1,426
    rows affected; the JSON and Markdown outputs were already correct.
  - `gwas_annotation_scope` and `gwas_index_siblings` — 1,527 rows attributed
    the second allele's value to the first.
  - `imgt_accession`, `hla_serotype`, `protein_group` — 356 values attributed
    to the wrong allele. An exactly-matched allele has an accession but no
    G-group and a prefix-matched one has the reverse, so within a single row
    two of these columns could describe different alleles.
  - `gwas_traits`, `gwas_p_values`, `gwas_pmids`, `gwas_annotation_resolution`,
    `pharm_drugs`, `pharm_evidence`, `pharm_cpic_action`, `pharm_pmids`,
    `disease_risk_summary`, `drug_response_summary`, `clinical_significance`
    and `input_quality_tier` — no misalignment observed in this cohort, but
    they carried the same latent defect and are corrected with the rest.
- **Every GWAS p-value in the TSV read `0.0000`.** The column was written as
  fixed-point with four decimals, but only associations at or below the
  genome-wide threshold of 5×10⁻⁸ are retained, so no p-value could ever be
  represented — a p of 4×10⁻²⁴⁶ and one of 2×10⁻⁸ were indistinguishable.
  p-values are now written in scientific notation (`9.00e-13`) and remain
  parseable as plain floats. Only the TSV was affected; the JSON always
  carried the true value.
- **`db-update --db pharmgkb` could not download anything.** PharmGKB retired
  the `api.pharmgkb.org` host — it no longer resolves — so every attempt failed
  with a name-resolution error. The bulk files are now fetched from
  `s3.pgkb.org`, where PharmGKB serves them.
- **The CPIC Level 1A label was unreachable.** `clinical_annotations.tsv`
  carries no CPIC column, so the guideline link was always empty and the
  `Actionable pharmacogenomic risk (CPIC 1A — avoid)` label could never be
  emitted for any allele — HLA-B*57:01/abacavir included. The link is now
  joined from the `Guideline Annotation` rows of `clinical_ann_evidence.tsv`,
  and DPWG guidelines filed under the same evidence type are excluded, so only
  a genuine CPIC guideline backs a CPIC assertion.
- `hlante version` reported the GWAS query cache instead of the GWAS Catalog
  bulk dump, so a freshly installed dump still read `not installed`. The status
  line now reads the dump directory, and `version` takes `--gwas-cache-dir` in
  place of `--cache-dir`.

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
