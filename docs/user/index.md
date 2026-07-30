# HLAnte

HLAnte is a command-line tool that annotates HLA alleles with clinical evidence.

It ingests outputs from four HLA genotyping tools (ARCAS-HLA, T1K, HLA-HD, OptiType),
normalises alleles against a local IPD-IMGT/HLA database, and annotates each allele
with pharmacogenomic evidence (PharmGKB), disease associations (GWAS Catalog, built-in
curated table), and population allele frequencies (AFND).

## Documentation

- [Step-by-step tutorial](TUTORIAL.md) — a guided first session for users new to
  Python and the command line
- [Quickstart](quickstart.md) — install, first annotation, common flags
- [Worked examples](EXAMPLES.md) — runnable recipes: one per typing tool, batch
  and mixed-tool cohorts, pharmacogenomic and disease screening of the TSV,
  population-stratified runs, and reference-release pinning
- [Inputs](inputs.md) — supported typing-tool formats and how they are parsed
- [Outputs](outputs.md) — the TSV/JSON/Markdown schema and every column
- [TSV Interpretation Guide](TSV_INTERPRETATION_GUIDE.md) — a column-by-column
  walkthrough for biologists and clinicians reading the report
- [Input-quality tiers](tier_definitions.md) — how the input-quality score is derived
- [Updating databases](update_db.md) — refreshing IMGT, PharmGKB, GWAS, AFND
- [Limitations](limitations.md) and [Validation](validation.md)
- [API reference](api_reference.md)

!!! note
    Documentation is being actively developed. New users should start with the
    [step-by-step tutorial](TUTORIAL.md); see [Quickstart](quickstart.md) for the
    condensed version.
