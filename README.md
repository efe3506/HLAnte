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
- IPD-IMGT/HLA-compliant allele normalization (2/4/6/8 digit,
  G and P groups)
- Resolution-aware fallback for GWAS / PharmGKB / AFND lookups
- Built-in connectors for GWAS Catalog, PharmGKB (CPIC 1A/1B), AFND, and
  a curated HLA–disease/drug table
- Population-aware allele-frequency scoring via AFND with a universal
  geographic taxonomy (EUR / AFR / EAS / SAS / MID / AMR / OCE / global)
- Confidence score combining novelty, rarity, resolution, and
  ambiguity signals
- TSV / CSV machine-readable output + Markdown research-annotation summary + JSON
- Batch mode, local database update, configurable log levels
- Research-oriented disclaimer embedded in every report

## Installation

### Conda (Bioconda environment)

```bash
conda env create -f environment.yml
conda activate hlante
```

### Developer / user install (pip)

```bash
git clone https://github.com/efe3506/HLAnte.git
cd HLAnte
pip install -e .
```

For development (includes linters and test runner):

```bash
pip install -e ".[dev]"
```

## Usage

### Single-sample annotation

```bash
hlante annotate \
  -i sample.genotype.tsv \
  -t arcas-hla \
  -o results/ \
  -p EUR \
  --format tsv
```

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
