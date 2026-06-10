# HLAnte 1000 Genomes Annotation Benchmark

Generated: 2026-05-04  
HLAnte version: 0.1.0  IMGT-HLA: IPD-IMGT/HLA 3.64.0  
N samples: 2692

## Scope

This benchmark evaluates HLAnte's annotation pipeline using Sanger-validated HLA types from the 1000 Genomes Project (Abi-Rached et al. 2018) as input. HLAnte is an annotation tool, not a typing tool; typing accuracy is out of scope.

## Parser Fidelity

| Tool | Alleles parsed | Expected | Success rate |
|------|----------------|----------|--------------|
| arcashla  |          26308 |    26332 |        99.9% |
| t1k       |          26284 |    26332 |        99.8% |
| hlahd     |          26284 |    26332 |        99.8% |
| optitype  |          15868 |    16152 |        98.2% |

## Normalisation Success (IMGT-HLA recognised; exact or prefix match)

Note: 2-field input alleles have no exact IMGT accession (ambiguous by design) but are counted as recognised when a prefix match exists (`is_novel=False`).  Novel alleles are those with no IMGT prefix.

| Tool | Allele count | IMGT recognised | Rate |
|------|--------------|-----------------|------|
| arcashla  |        26308 |           26304 | 100.0% |
| t1k       |        26284 |           26280 | 100.0% |
| hlahd     |        26284 |           26280 | 100.0% |
| optitype  |        15868 |           15864 | 100.0% |

## CPIC Level 1A Pharmacogenomic Recall

Results shown for arcashla (representative; other tools similar).

| Drug-allele | Carriers | Hits | Recall |
|-------------|----------|------|--------|
| B*57:01 → abacavir       |      162 |  162 | 100.0% |
| B*58:01 → allopurinol    |      176 |  176 | 100.0% |
| B*15:02 → carbamazepine  |       98 |   98 | 100.0% |
| A*31:01 → carbamazepine  |      131 |  131 | 100.0% |

## GWAS / Curated Disease Recall

arcasHLA (representative; 100% recall across all sentinels):

| Allele → Trait | Carriers | Hits | Recall |
|----------------|----------|------|--------|
| DRB1*03:01 → SLE / T1D              |      327 |  327 | 100.0% |
| DRB1*04:01 → rheumatoid arthritis   |      107 |  107 | 100.0% |
| DRB1*15:01 → multiple sclerosis     |      357 |  357 | 100.0% |
| B*27:05 → ankylosing spondylitis |       66 |   66 | 100.0% |
| B*51:01 → Behcet disease         |      259 |  259 | 100.0% |
| DQB1*02:01 → celiac disease         |      289 |  289 | 100.0% |
| DQB1*06:02 → narcolepsy             |      390 |  390 | 100.0% |

T1K and HLA-HD show 99.7% recall (356/357) for DRB1\*15:01 → multiple sclerosis.
The single missed carrier is HG01284: both T1K and HLA-HD fixture files for this
sample contain only A, B, C, and DQB1 loci — DRB1 was not typed by those tools
for that sample. arcasHLA reports DRB1\*15:01 for the same sample, achieving 100%
recall. This is a typing-tool coverage gap, not an HLAnte annotation failure.

## Population-Stratified Performance (arcashla)

| Super-pop | N | AFND coverage | Mean confidence | LOW tier % |
|-----------|---|---------------|-----------------|------------|
| EUR       | 529 |         96.7% |           0.895 |       0.1% |
| AFR       | 712 |         85.0% |           0.880 |       0.0% |
| EAS       | 536 |         91.2% |           0.888 |       0.0% |
| SAS       | 543 |         87.4% |           0.883 |       0.0% |
| AMR       | 372 |         87.4% |           0.883 |       0.0% |

## Confidence Tier Distribution (arcashla)

Note: This benchmark uses --input-source validated because the 1000 Genomes types are Sanger-validated reference data. Under this mode, the ambiguity penalty (×0.75) is suppressed for 2-field inputs since the call is exactly correct at its reported resolution; only the resolution penalty (×0.90 for 4-field-equivalent matching) applies. This produces predominantly HIGH-tier scores. For typing-tool inputs (--input-source typing_tool, default), 2-field outputs would receive both penalties and score in the LOW range.

| Tier | Count | Percentage |
|------|-------|------------|
| HIGH     | 23480 |      89.3% |
| MODERATE |  2824 |      10.7% |
| LOW      |     4 |       0.0% |
| NA       |     0 |       0.0% |

---

> **Scope note**: This benchmark validates annotation pipeline behaviour, not typing accuracy. HLAnte does not perform HLA typing.
