# HLAnte

**HLA genotype → disease / drug research annotation aid (research use only)**

HLAnte parses HLA typing tool outputs (ARCAS-HLA, T1K, HLA-HD, OptiType)
and produces a unified research-annotation report by querying
IPD-IMGT/HLA, GWAS Catalog, PharmGKB (CPIC 1A/1B), AFND, and a built-in
curated HLA–disease/drug table. Outputs are emitted in TSV, Markdown,
and JSON formats for downstream analysis and API consumption.

**Research use only.** HLAnte does not implement ACMG/AMP criteria.
The evidence-strength labels surfaced in every output are descriptive,
not diagnostic. See the disclaimer embedded in every report.

## Features

- Multi-tool HLA typing parser (ARCAS-HLA, T1K, HLA-HD, OptiType)
- IPD-IMGT/HLA-compliant allele normalization (one to four fields,
  G and P groups)
- Resolution-aware fallback for GWAS / PharmGKB / AFND lookups
- Built-in connectors for GWAS Catalog, PharmGKB (CPIC 1A/1B), AFND, and
  a curated HLA–disease/drug table
- Population-aware allele-frequency scoring via AFND with a universal
  geographic taxonomy (EUR / AFR / EAS / SAS / MID / AMR / OCE / global)
- Input-quality score combining novelty, rarity, resolution, and
  ambiguity signals
- TSV / CSV machine-readable output + Markdown research-annotation summary + JSON
- Batch mode, local database update, configurable log levels
- Research-oriented disclaimer embedded in every report

## Installation

### Conda (Bioconda environment)

```bash
git clone https://github.com/efe3506/HLAnte.git
cd HLAnte
conda env create -f environment.yml
conda activate hlante
```

### pip (from source)

HLAnte is installed from this repository, not from PyPI — `pip install hlante`
will not work. Clone first, then install:

```bash
git clone https://github.com/efe3506/HLAnte.git
cd HLAnte
pip install .
```

For development (editable install with linters and test runner):

```bash
pip install -e ".[dev]"
```

### Required one-time database download

HLAnte normalises every allele call against a local IPD-IMGT/HLA release, so
this download is **required before the first annotation run** (~10 MB):

```bash
hlante db-update --db imgt
```

Without it, `hlante annotate` stops with an error and writes no report. To
reproduce a specific release, pin it: `hlante db-update --db imgt --imgt-ref 3.64.0`.
The AFND mirror pins the same way: `hlante db-update --db afnd --afnd-ref <commit>`.
Both record the ref and a SHA-256 in `version.json`.
The GWAS Catalog, PharmGKB, and AFND downloads are optional — see
[INSTALL.md](INSTALL.md#5-database-setup).

## Usage

New to Python or the command line? Work through
[`docs/user/TUTORIAL.md`](docs/user/TUTORIAL.md) instead — it walks through
install, the one-time database download, a first annotation and how to read
the report, one command at a time.

### Single-sample annotation

arcasHLA emits JSON, so the input file is the genotype JSON produced by
`arcasHLA genotype` (a ready-made example ships in `tests/fixtures/`):

```bash
hlante annotate \
  -i tests/fixtures/sample.genotype.json \
  -t arcashla \
  -o results/ \
  -p EUR \
  --format tsv
```

This writes `results/hlante_report.tsv`. With `--format all` (the default) you
additionally get `hlante_report.md` and `hlante_report.json`. Use `--prefix` to
change the `hlante_report` stem.

### Directory batch

Point ``-i`` at a directory; HLAnte discovers files matching the chosen
tool's default glob patterns and produces a single combined report.

```bash
hlante annotate \
  -i cohort/ \
  -t t1k \
  -o results/ \
  -p global \
  --format all
```

For multi-sample cohort workflows including the recommended metadata
manifest schema, ancestry-aware per-group annotation, and worked
examples, see [`docs/COHORT_METADATA.md`](docs/COHORT_METADATA.md).

### Refreshing local databases

```bash
hlante db-update --db imgt
hlante db-update --db pharmgkb
hlante db-update --db gwas
# AFND has no stable bulk endpoint; place the TSV manually at
# ~/.hlante/afnd/afnd_frequencies.tsv (or use --afnd-url for a mirror).
```

## Commands

| Command      | Description                                               |
|--------------|-----------------------------------------------------------|
| `annotate`   | Run the full annotation pipeline on a file or directory   |
| `validate`   | Sanity-check input format and report file/sample counts   |
| `db-update`  | Refresh local copies of IMGT, PharmGKB, GWAS, AFND        |
| `version`    | Print HLAnte and installed database versions               |

## Population groups

The ``--population`` flag accepts the following universal geographic
codes (matching AFND's Population Group column via keyword substring):

| Code     | Keywords                                          |
|----------|---------------------------------------------------|
| `EUR`    | European, Caucasian, White, *regional European*   |
| `AFR`    | African, Sub-Saharan, West / East African, Black  |
| `EAS`    | Asian, East Asian, Southeast Asian (`ASN` is accepted as an alias) |
| `SAS`    | South Asian                                       |
| `MID`    | Middle Eastern, Arab, North African, West Asian   |
| `AMR`    | American, Hispanic, Latino, Mestizo               |
| `OCE`    | Oceanian, Pacific Islander                        |
| `global` | Aggregate across every population (sample-weighted) |

No country names appear in the taxonomy; a row labelled e.g.
``"Japan pop 1"`` with ``Population Group = "East Asian"`` matches
``EAS`` via the Population Group substring, not via the country name.

## Data sources

- [IPD-IMGT/HLA](https://www.ebi.ac.uk/ipd/imgt/hla/)
- [GWAS Catalog](https://www.ebi.ac.uk/gwas/)
- [PharmGKB](https://www.pharmgkb.org/)
- [AFND](http://www.allelefrequencies.net/)

## Development

```bash
pytest
pytest -m integration   # regression tests against real local dumps
(cd src && mypy hlante)  # strict type check (src layout)
ruff check src/hlante
```

## Disclaimer

HLAnte is intended for **research purposes only** and must not be used
for clinical decision-making. Treatment decisions must be made by a
qualified clinician.

## About the name

HLAnte is a portmanteau of **HLA** and ***andante*** — the musical
tempo marking for a moderate, walking pace.

## License

MIT License — see the `LICENSE` file for details.
