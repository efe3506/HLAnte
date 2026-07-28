# Limitations

## Research use only

HLAnte is a **research annotation tool**, not a clinical decision support system.
Outputs must not be used to guide patient care, prescribing, or diagnostic decisions
without independent validation by qualified medical professionals.

## Scope limitations

- **Allele-level only**: annotation is per-allele, not haplotype- or diplotype-level.
  DQ2 (DQA1\*05:01 + DQB1\*02:01) and DQ8 heterodimer combinations are not computed.
- **Batch mode**: the CLI accepts a file or directory per invocation. Multi-sample
  cohort processing requires a wrapper script.
- **KIR genotyping**: out of scope.
- **Transplant matching**: permanently out of scope.

## Database limitations

- **GWAS Catalog**: data quality varies. Effect sizes from quantitative-trait GWAS
  (β coefficients) are stored alongside disease OR values; always check the
  `effect_size_warning` flag.
- **PharmGKB**: PharmGKB coverage is strongest for European-ancestry populations.
- **AFND**: no standard bulk download endpoint exists. HLAnte ships a curated
  frequency table derived from Gonzalez-Galarza 2020. Custom population coverage
  may require manual TSV preparation (`--afnd-dir`).
- **Built-in curated table**: 40 manually curated entries; only well-replicated
  associations with published PMID are included. Coverage is intentionally
  conservative — absence of an entry does not mean the allele has no association.
- **GWAS Catalog EMBL-EBI terms of use**: data obtained from GWAS Catalog must
  comply with EMBL-EBI terms of use (non-commercial use).

## Population frequency limitations

Coverage is uneven across populations. Confidence scores may be systematically
underestimated for alleles common in underrepresented populations (e.g., AFR, OCE, MID).

## Reproducibility

GWAS Catalog and PharmGKB databases update continuously. Report headers include
a `gwas_cache_date` field recording when the local bulk dump was downloaded.
For reproducible results, archive the local database copies alongside outputs.

## Caller confidence and alternative genotypes

HLAnte annotates one genotype per locus per sample: the genotype the typing
tool reported. What the tool knew about its own confidence is largely not
carried through, and this bounds how the annotations should be read.

| Tool | Provided by the tool | Retained by HLAnte |
|------|----------------------|--------------------|
| arcasHLA | no per-allele quality or likelihood in the genotype JSON | — |
| T1K | per-allele abundance and quality (native layout) | quality, in `caller_allele_quality` |
| HLA-HD | no numeric confidence in the result file | — |
| OptiType | a solution-level objective, and optionally several enumerated optimal solutions | neither; a warning is emitted when more than one solution is present |

Consequences:

- **No posterior probabilities are available.** No supported caller reports a
  posterior probability for a genotype, so neither HLAnte nor its
  `input_quality_score` can express one.
- **Alternative genotypes are not evaluated.** When a caller enumerates
  several optimal solutions, only the top-ranked one is annotated. Actionable
  annotations can differ between plausible alternative genotypes, so a report
  reflects one solution rather than the full set the caller considered.
- **`caller_allele_quality` is on the tool's own scale.** It is comparable
  only between T1K runs, and is not comparable with `input_quality_score`.

Where a genotype is uncertain, the primary typing output remains the
authoritative record.
