# Test fixtures

This directory contains sample outputs used by the HLAnte unit tests.

Expected files
--------------

| File                                 | Source tool   |
|--------------------------------------|---------------|
| `sample.genotype.json`               | ARCAS-HLA     |
| `sample_nested.genotype.json`        | ARCAS-HLA (nested schema) |
| `sample_t1k_genotype.tsv`            | T1K           |
| `sample_final.result.txt`            | HLA-HD        |
| `sample_result.tsv`                  | OptiType      |
| `malformed.json`                     | ARCAS-HLA (invalid JSON, for error-path tests) |
| `imgt_mini/`                         | IPD-IMGT/HLA mini copy |
| `pharmgkb/`                          | PharmGKB mini copy |
| `gwas/`                              | GWAS Catalog mini copy |
| `afnd/`                              | AFND mini copy |

Use synthetic / example alleles rather than real patient data.
