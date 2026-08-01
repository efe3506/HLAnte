# HLAnte 1000 Genomes Annotation Benchmark

Generated: 2026-08-01  
HLAnte version: 0.2.0  IMGT-HLA: IPD-IMGT/HLA 3.64.0  
N samples: 2693

## Scope

This benchmark evaluates HLAnte's annotation pipeline using Sanger-validated HLA types from the 1000 Genomes Project (Abi-Rached et al. 2018) as input. HLAnte is an annotation tool, not a typing tool; typing accuracy is out of scope.

## Parser Fidelity

| Tool | Alleles parsed | Expected | Success rate |
|------|----------------|----------|--------------|
| arcashla  |          26301 |    26332 |        99.9% |
| t1k       |          26294 |    26332 |        99.9% |
| hlahd     |          26294 |    26332 |        99.9% |
| optitype  |          16124 |    16152 |        99.8% |

## Normalisation Success (IMGT-HLA recognised; exact or prefix match)

Note: two-field input alleles have no exact IMGT accession (ambiguous by design) but are counted as recognised when a prefix match exists (`is_novel=False`).  Novel alleles are those with no IMGT prefix.

| Tool | Allele count | IMGT recognised | Rate |
|------|--------------|-----------------|------|
| arcashla  |        26301 |           26300 | 100.0% |
| t1k       |        26294 |           26293 | 100.0% |
| hlahd     |        26294 |           26293 | 100.0% |
| optitype  |        16124 |           16123 | 100.0% |

## CPIC Level 1A Pharmacogenomic Recall

Results shown for arcashla (representative; other tools similar).

| Drug-allele | Carriers | Hits | Recall |
|-------------|----------|------|--------|
| B*57:01 → abacavir       |      162 |  162 | 100.0% |
| B*58:01 → allopurinol    |      176 |  176 | 100.0% |
| B*15:02 → carbamazepine  |       98 |   98 | 100.0% |
| A*31:01 → carbamazepine  |      131 |  131 | 100.0% |

## GWAS / Curated Disease Recall

| Allele → Trait | Carriers | Hits | Recall |
|----------------|----------|------|--------|
| DRB1*03:01 → SLE / T1D              |      327 |  327 | 100.0% |
| DRB1*04:01 → rheumatoid arthritis   |      107 |  107 | 100.0% |
| DRB1*15:01 → multiple sclerosis     |      357 |  357 | 100.0% |
| B*27:05 → ankylosing spondylitis |       66 |   66 | 100.0% |
| B*51:01 → Behcet disease         |      259 |  259 | 100.0% |
| DQB1*02:01 → celiac disease         |      289 |  289 | 100.0% |
| DQB1*06:02 → narcolepsy             |      390 |  390 | 100.0% |

## Population-Stratified Performance (arcashla)

| Super-pop | N | AFND coverage | Mean input quality | limited tier % |
|-----------|---|---------------|-----------------|------------|
| EUR       | 529 |        100.0% |           0.894 |       0.1% |
| AFR       | 712 |        100.0% |           0.880 |       0.5% |
| EAS       | 536 |         99.9% |           0.892 |       0.2% |
| SAS       | 543 |         99.8% |           0.886 |       0.6% |
| AMR       | 372 |         99.8% |           0.883 |       1.0% |

## Input-Quality Tier Distribution (arcashla)

Note: this run used --input-source validated, because the 1000 Genomes types are Sanger-validated reference data. The ambiguity penalty (×0.75) is suppressed under that mode, leaving the two-field resolution penalty (×0.90), so most calls land in the detailed tier.

| Tier | Count | Percentage |
|------|-------|------------|
| detailed | 24546 |      93.3% |
| partial  |  1632 |       6.2% |
| limited  |   123 |       0.5% |
| NA       |     0 |       0.0% |

---

> **Scope note**: This benchmark validates annotation pipeline behaviour, not typing accuracy. HLAnte does not perform HLA typing.
