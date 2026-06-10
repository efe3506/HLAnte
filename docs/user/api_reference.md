# Python API Reference

HLAnte can be used directly from Python as well as from the command line.
The public API consists of four functions covering the same
parse → normalize → annotate → report pipeline exposed by the CLI.

!!! warning
    The Python API is not yet stable. Function signatures may change
    between minor versions until v1.0.

---

## Pipeline functions

### `parse_hla_output`

```python
from hlante.parser import parse_hla_output

genotypes = parse_hla_output(path, tool)
```

Parse one HLA typing tool output file.

| Parameter | Type | Description |
|-----------|------|-------------|
| `path` | `Path` or `str` | Path to the input file |
| `tool` | `str` | One of `"arcashla"`, `"t1k"`, `"hlahd"`, `"optitype"` |

Returns a list of `HLAGenotype` objects (one per locus call).

Raises `HLAnteParseError` on malformed input.

---

### `batch_normalize`

```python
from hlante.normalizer import batch_normalize, load_imgt_db

imgt_db = load_imgt_db()
normalized = batch_normalize(genotypes, imgt_db=imgt_db, max_workers=4)
```

Normalize a list of `HLAGenotype` objects against the IPD-IMGT/HLA database.

| Parameter | Type | Description |
|-----------|------|-------------|
| `genotypes` | list of `HLAGenotype` | Output from `parse_hla_output` |
| `imgt_db` | dict | Loaded IMGT database (`load_imgt_db()`) |
| `max_workers` | int | Thread count (default: 4) |

Returns a list of `NormalizedAllele` objects.

---

### `annotate_genotype`

```python
from hlante.annotator import annotate_genotype, AnnotatorConfig

config = AnnotatorConfig(
    offline=False,
    population_group="EUR",
)
annotated = annotate_genotype(normalized, config)
```

Annotate a list of `NormalizedAllele` objects with GWAS, PharmGKB, AFND,
and curated disease evidence.

| Parameter | Type | Description |
|-----------|------|-------------|
| `normalized` | list of `NormalizedAllele` | Output from `batch_normalize` |
| `config` | `AnnotatorConfig` | Annotation settings |

Returns a list of `AnnotatedHLA` objects.

**`AnnotatorConfig` options:**

| Attribute | Default | Description |
|-----------|---------|-------------|
| `offline` | `False` | Use only local caches |
| `population_group` | `"global"` | AFND population group |
| `enable_gwas` | `True` | Query GWAS Catalog |
| `enable_pharmgkb` | `True` | Query PharmGKB |
| `enable_afnd` | `True` | Query AFND |
| `input_source` | `InputSource.TYPING_TOOL` | Provenance tag |

---

### `generate_all`

```python
from hlante.reporter import generate_all, ReportContext
from pathlib import Path

context = ReportContext()
paths = generate_all(annotated, Path("output/"), prefix="report", context=context)
# paths == {"tsv": Path("output/report.tsv"),
#            "json": Path("output/report.json"),
#            "markdown": Path("output/report.md")}
```

Write TSV, JSON, and Markdown reports to a directory.

---

## End-to-end example

```python
from pathlib import Path
from hlante.parser import parse_hla_output
from hlante.normalizer import batch_normalize, load_imgt_db
from hlante.annotator import annotate_genotype, AnnotatorConfig
from hlante.reporter import generate_all, ReportContext

# 1. Parse
genotypes = parse_hla_output(Path("HG00096.genotype.json"), "arcashla")

# 2. Normalize
imgt_db = load_imgt_db()
normalized = batch_normalize(genotypes, imgt_db=imgt_db)

# 3. Annotate
config = AnnotatorConfig(offline=True, population_group="global")
annotated = annotate_genotype(normalized, config)

# 4. Report
context = ReportContext()
paths = generate_all(annotated, Path("output/"), context=context)
print("TSV written to:", paths["tsv"])
```

---

## Key data classes

### `AnnotatedHLA`

The central result object. Key attributes:

| Attribute | Type | Description |
|-----------|------|-------------|
| `normalized_allele` | `NormalizedAllele` | Normalized allele info |
| `gwas_hits` | list of `GWASHit` | GWAS Catalog records |
| `pharm_annotations` | list of `PharmAnnotation` | PharmGKB records |
| `disease_entries` | list of `DiseaseEntry` | Built-in curated records |
| `clinical_significance` | str | Evidence-strength label |
| `confidence_score` | float | Score in `[0, 1]` |
| `confidence_tier` | str | `HIGH` / `MODERATE` / `LOW` |
| `allele_frequency` | float or None | AFND frequency |

---

## Exceptions

| Exception | Module | Raised when |
|-----------|--------|-------------|
| `HLAnteParseError` | `hlante.parser` | Malformed input file |
| `UnsupportedToolError` | `hlante.parser` | Unknown tool name |
| `IMGTDatabaseMissingError` | `hlante.normalizer` | IMGT DB not installed |
| `HLAReportError` | `hlante.reporter` | Report generation failure |
| `OutputFileExistsError` | `hlante.reporter` | Output file exists and `--overwrite` not set |
