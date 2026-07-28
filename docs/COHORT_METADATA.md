# Cohort Metadata Guide

HLAnte processes one sample per input file and derives ``sample_id``
from the file name. For cohort analyses (N ≥ 10 samples with
phenotype / covariate joins) you will usually need to **pair HLAnte's
per-sample output with an external metadata manifest** and join them
post-hoc.

This document describes the recommended metadata layout, the HLAnte
workflow that consumes it, and two worked examples.

---

## 1. Recommended manifest schema

Maintain a single TSV file next to your cohort directory. Minimum
columns:

| Column              | Type   | Required | Description                                                                                   |
|---------------------|--------|----------|-----------------------------------------------------------------------------------------------|
| `sample_id`         | string | **yes**  | Must match the HLAnte-derived `sample_id` (file stem, see §3).                                |
| `file_name`         | string | **yes**  | Basename of the tool output file (e.g., `subject_042.genotype.json`).                         |
| `tool`              | string | **yes**  | One of `arcashla`, `t1k`, `hlahd`, `optitype`.                                                |
| `population_group`  | string | **yes**  | AFND group code: `EUR` / `AFR` / `EAS` / `SAS` / `MID` / `AMR` / `OCE` / `global` (`ASN` = alias for `EAS`). |
| `phenotype`         | string | optional | Primary study group (e.g., `case`, `control`, disease name).                                  |
| `sex`               | string | optional | `F`, `M`, or free text.                                                                       |
| `age`               | int    | optional | Age in years at sample collection.                                                            |
| `ancestry_self_reported` | string | optional | Free text — kept for transparency even if `population_group` is coarser.                 |
| `sequencing_batch`  | string | optional | Plate / flow-cell / batch identifier (for batch-effect analysis).                             |
| `sequencing_depth`  | float  | optional | Mean on-target coverage, useful for quality filtering.                                        |
| `sample_collection_date` | date | optional | ISO-8601 (`YYYY-MM-DD`).                                                                  |
| `notes`             | string | optional | Any free-text comment.                                                                        |

### Example manifest (10 samples)

`cohort_manifest.tsv`:

```tsv
sample_id	file_name	tool	population_group	phenotype	sex	age	sequencing_batch	sequencing_depth	notes
S001	S001.genotype.json	arcashla	EUR	case	F	42	B01	98.5	
S002	S002.genotype.json	arcashla	EUR	case	M	55	B01	102.3	
S003	S003.genotype.json	arcashla	EUR	control	F	39	B01	95.8	
S004	S004.genotype.json	arcashla	MID	case	M	47	B02	88.2	
S005	S005.genotype.json	arcashla	MID	control	F	52	B02	91.7	
S006	S006.genotype.json	arcashla	EAS	case	M	61	B02	89.4	
S007	S007.genotype.json	arcashla	EAS	control	F	34	B03	94.6	
S008	S008.genotype.json	arcashla	AFR	case	M	48	B03	87.1	withdrawn-consent-review-pending
S009	S009.genotype.json	arcashla	SAS	control	F	44	B03	96.2	
S010	S010.genotype.json	arcashla	global	case	M	50	B03	93.8	ancestry-unknown
```

---

## 2. `sample_id` derivation rules

HLAnte derives `sample_id` from the file name **stem**, with the
following tool-specific trimming:

| Tool       | Input filename                          | Derived `sample_id` |
|------------|-----------------------------------------|---------------------|
| ARCAS-HLA  | `subject_042.genotype.json`             | `subject_042`       |
| T1K        | `subject_042_genotype.tsv`              | `subject_042_genotype` (whole stem — rename to `subject_042.tsv` if cleaner IDs are needed) |
| HLA-HD     | `subject_042_final.result.txt`          | `subject_042`       |
| OptiType   | `subject_042_result.tsv`                | `subject_042`       |

For a consistent identifier across tools, **rename files to a single
convention before running HLAnte**:

```bash
# Example: normalize all ARCAS filenames to <sample_id>.genotype.json
for f in cohort/raw/*.json; do
  id=$(basename "$f" .genotype.json)
  cp "$f" "cohort/staged/${id}.genotype.json"
done
```

---

## 3. Workflow

### Step 1 — Run HLAnte on the cohort directory

HLAnte discovers every tool-matching file under the directory and
produces one unified report:

```bash
hlante annotate \
  -i cohort/staged/ \
  -t arcashla \
  --population global \
  -o cohort/reports/ \
  --format all
```

Outputs:

- `cohort/reports/hlante_report.tsv` — one row per `(sample_id, locus)`.
- `cohort/reports/hlante_report.md` — per-sample Markdown summaries.
- `cohort/reports/hlante_report.json` — nested JSON for API / pipeline use.

### Step 2 — Join with the manifest

The TSV is directly join-friendly. With **pandas** (Python):

```python
import pandas as pd

ann = pd.read_csv(
    "cohort/reports/hlante_report.tsv",
    sep="\t",
    comment="#",           # Skip the HLAnte metadata header lines
)
meta = pd.read_csv("cohort/cohort_manifest.tsv", sep="\t")

joined = ann.merge(meta, on="sample_id", how="left")

# Example downstream analysis: case/control counts per risk flag
counts = (
    joined.groupby(["phenotype", "clinical_significance"])
    .size()
    .unstack(fill_value=0)
)
print(counts)
```

With **R / tidyverse**:

```r
library(readr); library(dplyr)

ann <- read_tsv("cohort/reports/hlante_report.tsv", comment = "#")
meta <- read_tsv("cohort/cohort_manifest.tsv")

joined <- ann %>% left_join(meta, by = "sample_id")
joined %>% count(phenotype, clinical_significance)
```

### Step 3 — Per-population annotation (optional)

If your cohort spans several population groups and you want
population-specific allele frequencies in the report, **split the
cohort by group and run HLAnte once per group**:

```bash
# Simple split based on the manifest
for grp in EUR MID EAS AFR SAS AMR OCE; do
  mkdir -p cohort/staged/$grp
  awk -v g=$grp -F'\t' 'NR>1 && $4==g {print $2}' cohort/cohort_manifest.tsv | \
    while read fname; do
      cp "cohort/staged/$fname" "cohort/staged/$grp/"
    done

  if [ -n "$(ls -A cohort/staged/$grp 2>/dev/null)" ]; then
    hlante annotate \
      -i cohort/staged/$grp/ \
      -t arcashla \
      --population $grp \
      -o cohort/reports/$grp/ \
      --format tsv
  fi
done

# Concatenate group-level TSVs
head -1 cohort/reports/EUR/hlante_report.tsv > cohort/reports/all.tsv
for grp in EUR MID EAS AFR SAS AMR OCE; do
  tail -n +7 cohort/reports/$grp/hlante_report.tsv 2>/dev/null >> cohort/reports/all.tsv
done
```

(`head -1` / `tail -n +7` accounts for the six-line HLAnte metadata
header; adjust if header length changes.)

---

## 4. Worked example A — 50-sample case/control study

Goal: compare PharmGKB 1A drug-allele hit rates between cases and
controls.

```bash
# 1. Run
hlante annotate -i study/ -t t1k --population EUR \
  -o out/ --format tsv --overwrite

# 2. Join & tabulate in Python
python3 - <<'PY'
import pandas as pd
ann = pd.read_csv("out/hlante_report.tsv", sep="\t", comment="#")
meta = pd.read_csv("study/manifest.tsv", sep="\t")
joined = ann.merge(meta, on="sample_id")

# Proportion of samples carrying any 1A PharmGKB evidence
has_1a = joined["pharm_evidence"].fillna("NA").str.contains("1A")
by_phenotype = (
    joined.assign(has_1a=has_1a)
          .groupby("phenotype")["has_1a"]
          .mean()
)
print(by_phenotype)
PY
```

---

## 5. Worked example B — Multi-cohort ancestry-aware analysis

Goal: compute HLA-B*57:01 carrier frequency per population group,
with input-quality-weighted aggregation.

Given the manifest from §1 and its companion cohort directory,
the script below runs one `hlante annotate` call per group, joins
the outputs, and reports population-specific statistics:

```python
import pandas as pd
import subprocess
import pathlib

manifest = pd.read_csv("cohort/cohort_manifest.tsv", sep="\t")
results = []

for grp, grp_meta in manifest.groupby("population_group"):
    stage = pathlib.Path(f"cohort/staged/{grp}")
    stage.mkdir(parents=True, exist_ok=True)
    for fname in grp_meta["file_name"]:
        src = pathlib.Path("cohort/staged") / fname
        if src.is_file():
            (stage / fname).write_bytes(src.read_bytes())

    out = pathlib.Path(f"cohort/reports/{grp}")
    subprocess.run([
        "hlante", "annotate",
        "-i", str(stage), "-t", "arcashla",
        "--population", grp,
        "-o", str(out), "--format", "tsv",
        "--overwrite",
    ], check=True)

    ann = pd.read_csv(out / "hlante_report.tsv", sep="\t", comment="#")
    ann["population_group"] = grp
    results.append(ann)

combined = pd.concat(results, ignore_index=True)
combined = combined.merge(manifest, on=["sample_id", "population_group"])

# B*57:01 carrier rate per group
b5701 = combined["allele1"].eq("B*57:01") | combined["allele2"].eq("B*57:01")
rate_by_group = (
    combined.assign(has_b5701=b5701)
            .groupby("population_group")
            .agg(carrier_rate=("has_b5701", "mean"),
                 n_samples=("sample_id", "nunique"),
                 mean_confidence=("input_quality_score", "mean"))
)
print(rate_by_group)
```

---

## 6. Reproducibility checklist

Before submitting a cohort paper that uses HLAnte, verify you have
captured:

- [ ] `cohort_manifest.tsv` with one row per sample.
- [ ] HLAnte version: `hlante --version` output.
- [ ] Reference DB versions: `hlante version` output.
- [ ] The `ReportContext` header lines (lines 1–6 of every HLAnte TSV)
      embed the run timestamp, DB versions, and disclaimer.
- [ ] Raw tool outputs (typing results) archived with SHA-256 sums.
- [ ] The exact command lines used (`hlante annotate ...`) recorded in
      a `commands.log` file.
- [ ] Any manual AFND TSV you bundled, with its source URL and
      download date.

HLAnte embeds the DB versions into every report header so post-hoc
audit is straightforward — a reader can confirm which GWAS /
PharmGKB / IPD-IMGT/HLA release produced each finding.

---

## 7. Privacy and consent

HLA genotypes can **uniquely identify individuals** and expose
disease-susceptibility information. Before sharing any cohort output:

- Ensure your study protocol has IRB / ethics board approval for
  HLA-based analysis.
- Where possible, publish only **aggregate counts** (carrier rates per
  group) rather than sample-level annotations.
- If you must publish sample-level data, confirm your participant
  consent documents permit HLA data release.
- Use pseudonymous `sample_id` values that cannot be back-referenced
  to direct identifiers.

This guide does not substitute legal or ethical review.
