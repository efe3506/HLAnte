#!/usr/bin/env bash
# =============================================================================
# install.sh — HLAnte installation script
#
# Usage:
#   bash install.sh [--conda | --pip] [--with-dbs] [--dev]
#
# Options:
#   --conda     Create a conda environment (default)
#   --pip       Install directly with pip (no conda)
#   --with-dbs  Also download the databases (IMGT, PharmGKB, GWAS)
#   --dev       Also install the development dependencies
# =============================================================================

set -euo pipefail

CONDA_MODE=true
WITH_DBS=false
DEV_MODE=false

while [[ $# -gt 0 ]]; do
  case "$1" in
    --conda)    CONDA_MODE=true;  shift ;;
    --pip)      CONDA_MODE=false; shift ;;
    --with-dbs) WITH_DBS=true;    shift ;;
    --dev)      DEV_MODE=true;    shift ;;
    -h|--help)
      echo "Usage: bash install.sh [--conda|--pip] [--with-dbs] [--dev]"
      exit 0
      ;;
    *) echo "Unknown option: $1"; exit 1 ;;
  esac
done

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
echo "============================================================"
echo "HLAnte installation script"
echo "  Directory: ${REPO_DIR}"
echo "  Mode: $( [[ "${CONDA_MODE}" == "true" ]] && echo "conda" || echo "pip" )"
echo "  Databases: ${WITH_DBS}"
echo "  Development: ${DEV_MODE}"
echo "============================================================"

# ---------------------------------------------------------------------------
# Conda installation
# ---------------------------------------------------------------------------
if [[ "${CONDA_MODE}" == "true" ]]; then
  if ! command -v conda &>/dev/null; then
    echo "ERROR: conda not found. Please install Miniconda/Anaconda:"
    echo "  https://docs.conda.io/en/latest/miniconda.html"
    exit 1
  fi

  echo ""
  echo "[1/3] Creating conda environment (hlante)..."

  # Create the environment if it does not exist
  if conda info --envs | grep -q "^hlante "; then
    echo "  'hlante' environment already exists, updating..."
    conda env update -n hlante -f "${REPO_DIR}/environment.yml" --prune
  else
    conda env create -f "${REPO_DIR}/environment.yml"
  fi

  if [[ "${DEV_MODE}" == "true" ]]; then
    echo "  Installing development dependencies..."
    conda run -n hlante pip install -r "${REPO_DIR}/requirements-dev.txt"
  fi

  echo ""
  echo "[2/3] Installing the package..."
  conda run -n hlante pip install -e "${REPO_DIR}"

  echo ""
  echo "[3/3] Verifying the installation..."
  conda run -n hlante hlante version
  if [[ "${DEV_MODE}" == "true" ]]; then
    conda run -n hlante python -m pytest "${REPO_DIR}/tests" -q --tb=short \
      -m "not integration and not qa"
  else
    echo "  (test suite skipped — rerun with --dev to install pytest and run it)"
  fi

  # The databases are downloaded further down, outside this branch, so record
  # how to reach the installed console script from the current shell.
  HLANTE_CMD=(conda run -n hlante hlante)

  echo ""
  echo "============================================================"
  echo "✓ Installation complete"
  echo ""
  echo "Activate the environment:"
  echo "  conda activate hlante"
  echo "============================================================"

# ---------------------------------------------------------------------------
# Pip installation
# ---------------------------------------------------------------------------
else
  if ! command -v python3 &>/dev/null; then
    echo "ERROR: python3 not found."
    exit 1
  fi

  PYTHON_VERSION=$(python3 -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')")
  echo "Python version: ${PYTHON_VERSION}"

  echo ""
  echo "[1/3] Creating virtual environment..."
  python3 -m venv "${REPO_DIR}/.venv"
  source "${REPO_DIR}/.venv/bin/activate"

  echo ""
  echo "[2/3] Installing dependencies..."
  pip install --upgrade pip

  if [[ "${DEV_MODE}" == "true" ]]; then
    pip install -r "${REPO_DIR}/requirements-dev.txt"
  else
    pip install -r "${REPO_DIR}/requirements.txt"
  fi

  pip install -e "${REPO_DIR}"

  echo ""
  echo "[3/3] Verifying the installation..."
  hlante version
  if [[ "${DEV_MODE}" == "true" ]]; then
    python -m pytest "${REPO_DIR}/tests" -q --tb=short -m "not integration and not qa"
  else
    echo "  (test suite skipped — rerun with --dev to install pytest and run it)"
  fi

  HLANTE_CMD=("${REPO_DIR}/.venv/bin/hlante")

  echo ""
  echo "============================================================"
  echo "✓ Installation complete"
  echo ""
  echo "Activate it in each new shell:"
  echo "  source ${REPO_DIR}/.venv/bin/activate"
  echo "============================================================"
fi

# ---------------------------------------------------------------------------
# Database download (optional)
# ---------------------------------------------------------------------------
if [[ "${WITH_DBS}" == "true" ]]; then
  echo ""
  echo "Downloading databases (this may take a few minutes)..."
  echo ""

  # Neither branch leaves the installed environment active in this shell, so
  # call the console script through the path recorded above rather than
  # relying on it being on PATH.
  echo "[DB 1/3] Downloading IPD-IMGT/HLA..."
  "${HLANTE_CMD[@]}" db-update --db imgt

  echo "[DB 2/3] Downloading PharmGKB..."
  "${HLANTE_CMD[@]}" db-update --db pharmgkb

  echo "[DB 3/3] Downloading GWAS Catalog (~59 MB)..."
  "${HLANTE_CMD[@]}" db-update --db gwas

  echo ""
  "${HLANTE_CMD[@]}" version
fi
