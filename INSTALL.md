# HLAnte — Installation Guide

> **Version:** 0.2.0 | **Python:** ≥3.9 | **Platform:** Linux, macOS

---

## Table of Contents

1. [System Requirements](#1-system-requirements)
2. [Quick Install (Automated Script)](#2-quick-install-automated-script)
3. [Manual Install — Conda](#3-manual-install--conda)
4. [Manual Install — pip + venv](#4-manual-install--pip--venv)
5. [Database Setup](#5-database-setup)
6. [Verifying the Installation](#6-verifying-the-installation)
7. [Running the Benchmark](#7-running-the-benchmark)
8. [Troubleshooting](#8-troubleshooting)

---

## 1. System Requirements

| Requirement | Minimum | Recommended |
|-----------|---------|---------|
| Python | 3.9 | 3.11 |
| RAM | 2 GB | 8 GB (for large cohorts) |
| Disk | 500 MB | 2 GB (with all databases) |
| CPU | 1 core | 4 cores |
| Operating system | Ubuntu 20.04+ / Rocky 8+ / macOS 12+ | Ubuntu 22.04 |

**Prerequisites:**
```bash
# Ubuntu/Debian
sudo apt-get update && sudo apt-get install -y git python3 python3-pip python3-venv

# Rocky Linux / CentOS
sudo dnf install -y git python3 python3-pip
```

---

## 2. Quick Install (Automated Script)

```bash
# Clone the repository
git clone https://github.com/efe3506/HLAnte.git
cd HLAnte

# Install the package only (with conda)
bash install.sh --conda

# Install the package + databases
bash install.sh --conda --with-dbs

# Install with pip (if conda is not available)
bash install.sh --pip

# Install with the development dependencies
bash install.sh --conda --dev
```

---

## 3. Manual Install — Conda

### 3.1 Install Miniconda (if not present)

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -p $HOME/miniconda3
eval "$($HOME/miniconda3/bin/conda shell.bash hook)"
conda init bash && source ~/.bashrc
```

### 3.2 The HLAnte Environment

```bash
git clone https://github.com/efe3506/HLAnte.git
cd HLAnte

# Create the conda environment
conda env create -f environment.yml

# Activate the environment
conda activate hlante

# Install the package (development mode)
pip install -e .

# Verify
hlante version
```

### 3.3 Updating the Conda Environment

```bash
conda env update -n hlante -f environment.yml --prune
pip install -e .
```

---

## 4. Manual Install — pip + venv

```bash
git clone https://github.com/efe3506/HLAnte.git
cd HLAnte

# Check the Python version
python3 --version   # must be ≥3.9

# Create a virtual environment
python3 -m venv .venv
source .venv/bin/activate   # Linux/macOS
# .venv\Scripts\activate    # Windows (unsupported; for reference only)

# Install dependencies
pip install --upgrade pip
pip install -r requirements.txt

# Install HLAnte
pip install -e .

# Verify
hlante version
```

---

## 5. Database Setup

HLAnte uses 5 databases. Each is cached under `~/.hlante/`.

### 5.1 IPD-IMGT/HLA (required, ~10 MB)

```bash
hlante db-update --db imgt
```

**Source:** GitHub mirror (ANHIG/IMGTHLA)
**Time:** ~30 seconds
**Update frequency:** quarterly

### 5.2 PharmGKB (recommended, ~50 MB)

```bash
hlante db-update --db pharmgkb
```

**Source:** pharmgkb.org/downloads
**Time:** ~1 minute
**Update frequency:** monthly

### 5.3 GWAS Catalog (recommended, ~59 MB)

```bash
hlante db-update --db gwas
```

**Source:** EBI FTP
**Time:** 1–5 minutes (depending on connection speed)
**Update frequency:** continuous (monthly updates recommended)

### 5.4 AFND (built-in fallback) and NMDP (user-supplied)

The package ships a small built-in AFND fallback table
(`hlante/db/afnd_builtin.tsv`, 7 loci × 5 population groups) so that HLAnte
runs out of the box. For full coverage, use `hlante db-update --db afnd` or
supply your own extract.

**NMDP frequency data are not redistributed with HLAnte.** The resource is
licensed by NMDP/Be The Match and its terms do not permit redistribution, so no
NMDP table is bundled and the NMDP source stays inactive unless you supply your
own extract.

```bash
# Optional: your own AFND TSV
mkdir -p ~/.hlante/afnd/
cp /path/to/afnd_frequencies.tsv ~/.hlante/afnd/

# Optional: NMDP data you obtained yourself from https://frequency.nmdp.org/
# under that resource's terms of use
mkdir -p ~/.hlante/nmdp/
cp /path/to/nmdp_frequencies.tsv ~/.hlante/nmdp/
```

### 5.5 Install All at Once

```bash
hlante db-update --db all
```

---

## 6. Verifying the Installation

### 6.1 Version Check

```bash
hlante version
# Expected output:
# HLAnte v0.2.0
#   IPD-IMGT/HLA : 3.64.0  (~/.hlante/imgt_hla)
#   PharmGKB     : installed  (~/.hlante/pharmgkb)
#   GWAS Catalog : installed  (~/.hlante/gwas)
#   AFND         : installed  (~/.hlante/afnd)
```

### 6.2 Unit Tests

```bash
# Standard tests (excluding integration, ~1 second)
pytest tests/ -q -m "not integration and not qa"

# Expected: 381 passed, 8 deselected

# Integration tests (require local databases)
pytest tests/ -q -m integration

# Full QA panel (slow; writes real database output)
pytest tests/test_qa_full_panel.py -s -m qa
```

### 6.3 Example Annotation

> **Prerequisite:** step 5.1 (`hlante db-update --db imgt`) must have completed.
> Allele normalisation requires a local IPD-IMGT/HLA release; without it
> `hlante annotate` exits with an error and writes no report. Note that
> `--offline` means "make no network calls", not "run without databases".

```bash
# Quick check with a test fixture (parses only; no database needed)
hlante validate -i tests/fixtures/sample.genotype.json -t arcashla

# Real annotation, using the local IPD-IMGT/HLA release and the
# built-in AFND fallback table; no network access
hlante annotate \
    -i tests/fixtures/sample.genotype.json \
    -t arcashla \
    -o /tmp/hlante_test \
    --offline \
    --no-gwas \
    --no-pharmgkb

cat /tmp/hlante_test/hlante_report.tsv
```

Expected files in `/tmp/hlante_test/`: `hlante_report.tsv`, `hlante_report.md`,
and `hlante_report.json` (the `--format all` default; the run above writes all
three because `--format` was not restricted).

---

## 7. Running the Benchmark

HLAnte ships with a fully reproducible benchmark against the 1000 Genomes
Project Sanger-validated HLA genotypes, spanning 5 super-populations and 4
input-tool formats. The complete, step-by-step instructions (fixture
generation, the annotation benchmark, the inter-format concordance check, and
how to interpret every output file) are in
[`benchmarks/README.md`](benchmarks/README.md).

In short:

```bash
# 1. Make sure the IMGT database is present (pin the published release)
hlante db-update --db imgt --imgt-ref 3.64.0

# 2. Generate the per-tool fixtures from the ground-truth table
python scripts/benchmark/convert_1000g_to_tool_formats.py \
    --input  benchmarks/1000g_HLA_types.tsv \
    --output benchmarks/fixtures_1000g/

# 3. Run the annotation benchmark
python scripts/benchmark/run_annotation_benchmark.py \
    --fixtures-dir benchmarks/fixtures_1000g/ \
    --ground-truth benchmarks/1000g_HLA_types.tsv \
    --output-dir   benchmarks/1000g_total_results_v2/ \
    --input-source validated \
    --threads 8
```

The committed reference outputs in `benchmarks/1000g_total_results/` correspond
to the published run; see `benchmarks/README.md` for the full flag reference and
the reproducibility notes.

---

## 8. Troubleshooting

### `hlante: command not found`

```bash
# PATH issues with a pip installation
export PATH="$HOME/.local/bin:$PATH"
# or
python -m hlante.cli --help
```

### `IMGTDatabaseMissingError`

```bash
hlante db-update --db imgt
```

### GWAS download is very slow

```bash
# Show a progress bar
hlante --log-level DEBUG db-update --db gwas
```

### `ModuleNotFoundError: No module named 'numpy'`

```bash
pip install numpy
# or
conda install numpy
```

### Conda `PackageNotInstalledError`

```bash
conda env remove -n hlante
conda env create -f environment.yml
```

### `afnd_builtin.tsv` error in the TSV output

`package_data` must be declared in `pyproject.toml`. Check:

```bash
python -c "from hlante.db.afnd import BUILTIN_AFND_TSV; print(BUILTIN_AFND_TSV.is_file())"
# should print True
```

If it does not:
```bash
pip install -e . --force-reinstall
```

---

## Directory Structure (After Installation)

```
~/.hlante/
├── imgt_hla/
│   ├── Allelelist.txt
│   ├── hla_nom_g.txt
│   └── hla_nom_p.txt
├── pharmgkb/
│   └── clinical_annotations.tsv
├── gwas/
│   └── *.json           (indexed GWAS data)
├── afnd/
│   └── afnd_frequencies.tsv  (optional)
└── nmdp/
    └── nmdp_frequencies.tsv  (optional)
```
