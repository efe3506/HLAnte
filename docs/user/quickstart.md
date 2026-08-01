# Quickstart

## Installation

HLAnte is installed from this repository, not from PyPI — `pip install hlante`
will not work, because the name is not published. Clone first, then install:

```bash
git clone https://github.com/efe3506/HLAnte
cd HLAnte
pip install .
```

Add `-e ".[dev]"` in place of `.` only if you intend to run the test suite.

## Step 1 — Download the databases

HLAnte needs local copies of IPD-IMGT/HLA, PharmGKB, GWAS Catalog, and AFND.
Run this once after installation (requires internet access, ~200 MB total):

```bash
hlante db-update
```

To update a single database:

```bash
hlante db-update --db imgt
hlante db-update --db pharmgkb
hlante db-update --db gwas
```

Check what is installed:

```bash
hlante version
```

## Step 2 — Annotate a file

### arcasHLA

```bash
hlante annotate -i HG00096.genotype.json -t arcashla
```

### T1K

```bash
hlante annotate -i HG00096_genotype.tsv -t t1k
```

### HLA-HD

```bash
hlante annotate -i HG00096_final.result.txt -t hlahd
```

### OptiType

```bash
hlante annotate -i HG00096_result.tsv -t optitype
```

## Step 3 — Read the output

By default HLAnte writes three files to `./hlante_output/`:

```
hlante_output/
├── hlante_report.tsv       # main annotation table (Excel-friendly)
├── hlante_report.json      # full provenance, nested
└── hlante_report.md        # human-readable summary
```

To specify an output location and prefix:

```bash
hlante annotate -i sample.genotype.json -t arcashla \
    --output-dir results/ \
    --prefix sample_HG00096
```

## Common options

| Flag | Example | What it does |
|------|---------|--------------|
| `--format` | `--format tsv` | Write only TSV (choices: `tsv`, `json`, `markdown`, `all`) |
| `--population` | `--population EUR` | AFND population for allele frequency (`EUR`, `AFR`, `EAS`, `SAS`, `MID`, `AMR`, `OCE`, `global`; `ASN` = alias for `EAS`) |
| `--offline` | `--offline` | Use only local caches; no HTTP calls |
| `--threads` | `--threads 8` | Parallelise normalisation (default: 4) |
| `--overwrite` | `--overwrite` | Overwrite existing output files |

Full option list:

```bash
hlante annotate --help
```

## Annotate a whole directory

When `--input` is a directory, HLAnte globs for the appropriate file
extension automatically:

```bash
hlante annotate -i arcashla_outputs/ -t arcashla --output-dir results/
```

## Validate before annotating

Check that your files are parseable without running the full pipeline:

```bash
hlante validate -i sample.genotype.json -t arcashla
```

## Research use only

All output files include the disclaimer:

> RESEARCH USE ONLY. Nothing in this output constitutes a clinical
> diagnosis, medical advice, or pharmacogenomic recommendation.
