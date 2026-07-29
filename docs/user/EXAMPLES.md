# Worked Examples

A collection of copy-and-paste recipes for laboratory and research use.
Each recipe was executed against the shipped test fixtures and the output
shown below is the output that was actually produced.

This document assumes HLAnte is **already installed** and that the
IPD-IMGT/HLA database has been downloaded. If that is not yet the case,
start with [INSTALL.md](../../INSTALL.md) and the
[Tutorial](TUTORIAL.md); the [Quickstart](quickstart.md) covers the
common flags. Column meanings are documented in
[Outputs](outputs.md) and, at length, in the
[TSV Interpretation Guide](TSV_INTERPRETATION_GUIDE.md).

---

## 0. Conventions used here

- Commands are run from the **repository root** unless a full path is
  given, so that `tests/fixtures/...` resolves.
- Reports are written to scratch directories under `/tmp/ex_out/`.
  Substitute your own output directory.
- `hlante` refers to the console script created by
  `pip install .`.
- Log lines (`INFO`/`WARNING`) are omitted from the pasted output,
  except where the log line *is* the point of the recipe; those are
  shown with their `<timestamp> INFO <module>:` prefix stripped.
- Long report rows are shown with a subset of columns
  (`cut -f...`) so they fit on the page. The real TSV has 40 columns.
- Tabular output is pasted after `| column -t -s$'\t'`, which aligns the
  tab-separated fields for reading. The files themselves are strictly
  tab-separated.

### One-time prerequisite

```bash
hlante db-update --db imgt
```

Annotation stops with an explicit error and writes no report if this
has not been run. The GWAS Catalog and AFND downloads are optional but
improve coverage:

```bash
hlante db-update --db gwas
hlante db-update --db afnd
```

### A note on the PharmGKB recipes

PharmGKB is **not redistributed** with HLAnte. The cohort runs in §2.1
and §4 — whose reports §3 then filters — pass
`--pharmgkb-dir tests/fixtures/pharmgkb`, a miniature seven-row extract
shipped for testing, so that the examples are reproducible immediately
after cloning. In production, run `hlante db-update --db pharmgkb`
once and **drop the `--pharmgkb-dir` flag** — HLAnte then reads
`~/.hlante/pharmgkb` automatically. Do not treat the fixture extract as
a substitute for the real resource; it carries seven illustrative rows,
not the PharmGKB clinical annotation set.

---

## 1. One recipe per typing tool

The four supported typers emit four different file formats. HLAnte
takes each one natively; you only have to name the tool with `-t`.
Before annotating an unfamiliar file, `hlante validate` confirms that
the parser recognises it.

### 1.1 arcasHLA

**What arcasHLA gives you** — `arcasHLA genotype` writes a JSON file
per sample, conventionally `<sample>.genotype.json`. Alleles are
usually reported as G groups:

```json
{
  "HLA-A": ["A*02:01:01G", "A*24:02:01G"],
  "HLA-B": ["B*07:02:01G", "B*15:01:01G"],
  "HLA-C": ["C*07:01:01G", "C*07:02:01G"],
  "HLA-DRB1": ["DRB1*15:01:01G", "DRB1*04:01:01G"],
  "HLA-DQB1": ["DQB1*06:02:01G", "DQB1*03:02:01G"]
}
```

(`tests/fixtures/sample.genotype.json`. A nested variant with a
top-level `"alleles"` key — as in
`tests/fixtures/sample_nested.genotype.json` — is also accepted.)

**Check it parses, then annotate:**

```bash
hlante validate -i tests/fixtures/sample.genotype.json -t arcashla

hlante annotate \
  -i tests/fixtures/sample.genotype.json \
  -t arcashla \
  -o /tmp/ex_out/tools/arcashla \
  --prefix arcashla \
  -p EUR \
  --format tsv --overwrite
```

`validate` prints:

```
✓ Input is valid
  Tool: arcashla
  Files: 1
  Samples: 1
  Locus calls: 5
  Unique loci: HLA-A, HLA-B, HLA-C, HLA-DQB1, HLA-DRB1
```

`annotate` writes `/tmp/ex_out/tools/arcashla/arcashla.tsv`. Columns
1–4, 6 and 9:

```
sample_id  locus     allele1         allele2         tool      hla_serotype
sample     HLA-A     A*02:01:01G     A*24:02:01G     arcashla  A2|A24
sample     HLA-B     B*07:02:01G     B*15:01:01G     arcashla  B7|B15
sample     HLA-C     C*07:01:01G     C*07:02:01G     arcashla  Cw7|Cw7
sample     HLA-DRB1  DRB1*15:01:01G  DRB1*04:01:01G  arcashla  DR15|DR4
sample     HLA-DQB1  DQB1*06:02:01G  DQB1*03:02:01G  arcashla  DQ6|DQ8
```

G-group input is carried through unchanged; the WHO serotype in
`hla_serotype` is derived per allele, which is often the fastest way to
reconcile an NGS report against a historical serological record.

### 1.2 T1K

**What T1K gives you** — `T1K --genotype` writes a TSV, one row per
gene, with a per-allele quality score. Untyped second alleles appear as
a bare `*`:

```
gene	allele1	allele2	score1	score2
HLA-A	A*02:01	A*24:02	98.5	97.2
HLA-B	B*07:02	B*15:01	96.1	95.8
HLA-C	C*07:01	*	94.3	0.0
HLA-DRB1	DRB1*15:01	DRB1*04:01	92.7	91.4
```

(`tests/fixtures/sample_t1k_genotype.tsv`)

```bash
hlante annotate \
  -i tests/fixtures/sample_t1k_genotype.tsv \
  -t t1k \
  -o /tmp/ex_out/tools/t1k \
  --prefix t1k \
  -p EUR \
  --format tsv --overwrite
```

```
sample_id            locus     allele1     allele2     tool  hla_serotype
sample_t1k_genotype  HLA-A     A*02:01     A*24:02     t1k   A2|A24
sample_t1k_genotype  HLA-B     B*07:02     B*15:01     t1k   B7|B15
sample_t1k_genotype  HLA-C     C*07:01     NA          t1k   Cw7
sample_t1k_genotype  HLA-DRB1  DRB1*15:01  DRB1*04:01  t1k   DR15|DR4
```

Two things to note:

- The `*` null token became `NA`, and *HLA-C* is reported as a
  single-allele locus. A single-allele locus is **hemizygous or
  not-fully-reported, not homozygous** — HLAnte never infers
  homozygosity from a missing second call.
- `sample_id` is the whole file stem, `sample_t1k_genotype`. HLAnte
  trims `.genotype`, `_final.result` and `_result` suffixes but not
  `_genotype`; rename T1K outputs to `<sample_id>.tsv` if you want
  clean identifiers. See
  [COHORT_METADATA.md §2](../COHORT_METADATA.md) for the full
  `sample_id` derivation table.

### 1.3 HLA-HD

**What HLA-HD gives you** — a tab-separated
`<sample>_final.result.txt` with one row per locus. Loci that failed to
type carry `-`:

```
HLA-A	A*02:01:01	A*24:02:01
HLA-B	B*07:02:01	B*15:01:01
HLA-C	C*07:01:01	-
HLA-DRB1	DRB1*15:01:01	DRB1*04:01:01
HLA-DQB1	DQB1*06:02:01	DQB1*03:02:01
HLA-DPB1	-	-
```

(`tests/fixtures/sample_final.result.txt`)

```bash
hlante annotate \
  -i tests/fixtures/sample_final.result.txt \
  -t hlahd \
  -o /tmp/ex_out/tools/hlahd \
  --prefix hlahd \
  -p EUR \
  --format tsv --overwrite
```

```
sample_id  locus     allele1        allele2        tool   hla_serotype
sample     HLA-A     A*02:01:01     A*24:02:01     hlahd  A2|A24
sample     HLA-B     B*07:02:01     B*15:01:01     hlahd  B7|B15
sample     HLA-C     C*07:01:01     NA             hlahd  Cw7
sample     HLA-DRB1  DRB1*15:01:01  DRB1*04:01:01  hlahd  DR15|DR4
sample     HLA-DQB1  DQB1*06:02:01  DQB1*03:02:01  hlahd  DQ6|DQ8
```

The fully untyped *HLA-DPB1* row is dropped rather than emitted with
two `NA` alleles. **A locus that is absent from the report was not
assessed — that is not the same as "no risk allele found".** The
Markdown report states this explicitly in a "Loci not typed" banner.

### 1.4 OptiType

**What OptiType gives you** — a single-row TSV of class I calls plus
read support. OptiType reports *HLA-A*, *HLA-B* and *HLA-C* only:

```
	A1	A2	B1	B2	C1	C2	Reads	Objective
0	A*02:01	A*24:02	B*07:02	B*15:01	C*07:01	C*07:02	1523	1456.78
```

(`tests/fixtures/sample_result.tsv`)

```bash
hlante annotate \
  -i tests/fixtures/sample_result.tsv \
  -t optitype \
  -o /tmp/ex_out/tools/optitype \
  --prefix optitype \
  -p EUR \
  --format tsv --overwrite
```

```
sample_id  locus  allele1  allele2  tool      hla_serotype
sample     HLA-A  A*02:01  A*24:02  optitype  A2|A24
sample     HLA-B  B*07:02  B*15:01  optitype  B7|B15
sample     HLA-C  C*07:01  C*07:02  optitype  Cw7|Cw7
```

Because OptiType is class I only, class II associations (celiac,
type 1 diabetes, rheumatoid arthritis, multiple sclerosis) simply
cannot be screened from an OptiType run. Use arcasHLA, T1K or HLA-HD
when class II matters.

### 1.5 Laboratory-validated typing

If the allele list came from an accredited SSP, SSO or SBT assay rather
than from a read-based typer, declare it with `--input-source
validated`. This suppresses the ambiguity penalty that HLAnte otherwise
applies to two-field (four-digit) calls:

```bash
hlante annotate -i /tmp/ex_out/lab/LAB01.genotype.json -t arcashla \
  -o /tmp/ex_out/src_typing_tool --format tsv --overwrite \
  -p EUR --input-source typing_tool

hlante annotate -i /tmp/ex_out/lab/LAB01.genotype.json -t arcashla \
  -o /tmp/ex_out/src_validated --format tsv --overwrite \
  -p EUR --input-source validated
```

Same input, the two settings side by side (columns 31–32):

```
-- --input-source typing_tool (default) --
sample_id  locus     allele1     allele2     input_quality_score  input_quality_tier
LAB01      HLA-A     A*01:01     A*02:01     0.6750|0.6750     LOW|LOW
LAB01      HLA-B     B*57:01     B*08:01     0.6750|0.6750     LOW|LOW
LAB01      HLA-DRB1  DRB1*03:01  DRB1*15:01  0.6750|0.6750     LOW|LOW

-- --input-source validated --
sample_id  locus     allele1     allele2     input_quality_score  input_quality_tier
LAB01      HLA-A     A*01:01     A*02:01     0.9000|0.9000     HIGH|HIGH
LAB01      HLA-B     B*57:01     B*08:01     0.9000|0.9000     HIGH|HIGH
LAB01      HLA-DRB1  DRB1*03:01  DRB1*15:01  0.9000|0.9000     HIGH|HIGH
```

The input-quality score summarises how completely the **allele call
itself** is characterised — novelty, resolution, ambiguity, and
population-frequency support. It is a heuristic descriptor of the input
and of annotation completeness, **not** an estimate of genotyping
accuracy and **not** a probability that any associated risk is real. A
`limited` tier never down-weights an actionable pharmacogenomic finding.
The naming of this score and of the resolution labels is being revised;
see [TSV_INTERPRETATION_GUIDE.md](TSV_INTERPRETATION_GUIDE.md) for the
current definition.

---

## 2. Batch and cohort recipes

### 2.1 A directory of results into one report

Point `-i` at a directory. HLAnte discovers the files matching the
tool's glob, parses them all, and writes **one** combined report:

| Tool       | Globs searched                    |
|------------|-----------------------------------|
| `arcashla` | `*.json`                          |
| `t1k`      | `*_genotype.tsv`, `*.tsv`         |
| `hlahd`    | `*final.result.txt`, `*.txt`      |
| `optitype` | `*_result.tsv`, `*.tsv`           |

Using the six-sample demonstration cohort from
[Appendix A](#appendix-a--the-demonstration-cohort):

```bash
hlante annotate \
  -i /tmp/ex_out/cohort \
  -t arcashla \
  -o /tmp/ex_out/cohort_report \
  -p EUR \
  --format tsv --overwrite \
  --pharmgkb-dir tests/fixtures/pharmgkb
```

```
  Tool: arcashla | File count: 6
✓ Parse complete: 30 locus call(s) (6 sample(s))
✓ Normalization complete: 60 allele(s)
✓ Annotation complete: 60 record(s) processed.
✓ Report(s) written:
  tsv: /tmp/ex_out/cohort_report/hlante_report.tsv
```

The result is one TSV with one row per `(sample_id, locus)` — 30 data
rows here. For the recommended metadata manifest, `sample_id`
conventions, phenotype joins and per-group runs, see
[COHORT_METADATA.md](../COHORT_METADATA.md); that document is the
reference for cohort study design and is not repeated here.

### 2.2 A cohort typed with more than one tool

`-t` takes exactly one tool per invocation, so a mixed cohort is run
once per tool and the reports are concatenated. Stage the files by
tool:

```
/tmp/ex_out/mixed/arcashla/PT01.genotype.json
/tmp/ex_out/mixed/arcashla/PT02.genotype.json
/tmp/ex_out/mixed/hlahd/PT04_final.result.txt
/tmp/ex_out/mixed/optitype/PT05_result.tsv
/tmp/ex_out/mixed/t1k/PT03_genotype.tsv
```

Then loop, and combine:

```bash
for t in arcashla t1k hlahd optitype; do
  hlante annotate -i /tmp/ex_out/mixed/$t -t $t \
    -o /tmp/ex_out/mixed_reports/$t --prefix $t \
    -p EUR --format tsv --overwrite
done

# One header row, then every data row from every per-tool report.
grep -v '^#' /tmp/ex_out/mixed_reports/arcashla/arcashla.tsv | head -1 \
  > /tmp/ex_out/mixed_reports/combined.tsv
for t in arcashla t1k hlahd optitype; do
  grep -v '^#' /tmp/ex_out/mixed_reports/$t/$t.tsv | tail -n +2 \
    >> /tmp/ex_out/mixed_reports/combined.tsv
done
```

```
sample_id      locus     allele1         allele2         tool
PT01           HLA-A     A*02:01:01G     A*24:02:01G     arcashla
PT01           HLA-B     B*07:02:01G     B*15:01:01G     arcashla
PT01           HLA-C     C*07:01:01G     C*07:02:01G     arcashla
PT01           HLA-DRB1  DRB1*15:01:01G  DRB1*04:01:01G  arcashla
PT01           HLA-DQB1  DQB1*06:02:01G  DQB1*03:02:01G  arcashla
PT02           HLA-A     A*01:01:01G     A*02:01:01G     arcashla
PT02           HLA-B     B*08:01:01G     B*44:02:01G     arcashla
PT02           HLA-C     C*07:01:01G     C*05:01:01G     arcashla
PT03_genotype  HLA-A     A*02:01         A*24:02         t1k
PT03_genotype  HLA-B     B*07:02         B*15:01         t1k
PT03_genotype  HLA-C     C*07:01         NA              t1k
PT03_genotype  HLA-DRB1  DRB1*15:01      DRB1*04:01      t1k
PT04           HLA-A     A*02:01:01      A*24:02:01      hlahd
PT04           HLA-B     B*07:02:01      B*15:01:01      hlahd
PT04           HLA-C     C*07:01:01      NA              hlahd
PT04           HLA-DRB1  DRB1*15:01:01   DRB1*04:01:01   hlahd
PT04           HLA-DQB1  DQB1*06:02:01   DQB1*03:02:01   hlahd
PT05           HLA-A     A*02:01         A*24:02         optitype
PT05           HLA-B     B*07:02         B*15:01         optitype
PT05           HLA-C     C*07:01         C*07:02         optitype
```

Use `grep -v '^#'` rather than a fixed `tail -n +N`: the `#` metadata
block grows or shrinks depending on which reference databases are
installed (seven lines in the runs above, five when only IPD-IMGT/HLA
is present).

The `tool` column is retained precisely so that tool-of-origin can be
carried into the analysis and, if necessary, adjusted for — resolution
and locus coverage differ systematically between typers.

### 2.3 Surviving one bad file in a batch

By default an unparseable file is skipped, named in the summary, and
the rest of the cohort is still annotated:

```bash
hlante annotate -i /tmp/ex_out/dirty -t arcashla \
  -o /tmp/ex_out/dirty_out --format tsv --overwrite
```

```
  Tool: arcashla | File count: 3
⚠ Skipped 1 unparseable file(s):
⚠   - /tmp/ex_out/dirty/PT03.genotype.json: ARCAS-HLA JSON parse error (…): Expecting property name enclosed in double quotes: line 1 column 3 (char 2)
✓ Parse complete: 8 locus call(s) (2 sample(s))
✓ Normalization complete: 16 allele(s)
✓ Annotation complete: 16 record(s) processed.
✓ Report(s) written:
  tsv: /tmp/ex_out/dirty_out/hlante_report.tsv
```

Exit status is `0`. In a validated pipeline you usually want the
opposite — add `--strict` and the whole run aborts on the first
unparseable file:

```
  Tool: arcashla | File count: 3
ERROR: Parse error: ARCAS-HLA JSON parse error (…): Expecting property name enclosed in double quotes: line 1 column 3 (char 2)
```

Exit status is `1`, so `set -e` in a wrapper script will stop there.

---

## 3. Downstream interpretation

All recipes in this section read
`/tmp/ex_out/cohort_report/hlante_report.tsv` produced in §2.1.

Two rules govern every filter:

1. **Strip the metadata block first** (`grep -v '^#'`, or
   `comment="#"` in pandas).
2. **Per-allele columns are positional.** In `clinical_significance`,
   `disease_risk_summary`, `drug_response_summary`, `allele_frequency`,
   `input_quality_score` and `input_quality_tier`, the first `|`-separated
   slot belongs to `allele1` and the second to `allele2`. A locus
   reported with one allele carries one slot. Columns named
   `*_allele1` / `*_allele2` are already split for you; the unsuffixed
   `gwas_traits` and `pharm_drugs` columns are aggregates over both
   alleles and should not be used for attribution.

### 3.1 Actionable pharmacogenomic findings

**Shell triage** — every row carrying a CPIC 1A avoid-level allele,
with the drug and the CPIC action verb:

```bash
awk -F'\t' '$1 !~ /^#/ && $28 ~ /Actionable pharmacogenomic risk/ \
  {print $1"\t"$2"\t"$3"\t"$4"\t"$20"\t"$23"\t"$24}' \
  /tmp/ex_out/cohort_report/hlante_report.tsv
```

```
SUBJ001	HLA-B	B*57:01:01	B*08:01:01	abacavir	1A	Contraindicated (do not use)
SUBJ002	HLA-B	B*58:01:01	B*15:02:01	allopurinol|carbamazepine	1A|1A	Contraindicated — use alternative|Avoid — contraindicated if carbamazepine-naïve
```

Columns are 28 `clinical_significance`, 20 `pharm_drugs`,
23 `pharm_evidence`, 24 `pharm_cpic_action`.

**pandas, per allele** — the row-level view above cannot tell you
*which* of the two alleles triggered the flag. Reshaping to one row per
allele does:

```python
import pandas as pd

REPORT = "/tmp/ex_out/cohort_report/hlante_report.tsv"
ann = pd.read_csv(REPORT, sep="\t", comment="#", dtype=str)

def per_allele(df):
    """One row per (sample, locus, allele). Slot 0 = allele1, slot 1 = allele2."""
    frames = []
    for slot in (0, 1):
        frames.append(pd.DataFrame({
            "sample_id": df["sample_id"],
            "locus": df["locus"],
            "allele": df[f"allele{slot + 1}"],
            "pharm_drugs": df[f"pharm_drugs_allele{slot + 1}"],
            "gwas_traits": df[f"gwas_traits_allele{slot + 1}"],
            "significance": df["clinical_significance"].str.split("|").str[slot],
            "frequency": df["allele_frequency"].str.split("|").str[slot],
            "tier": df["input_quality_tier"].str.split("|").str[slot],
        }))
    long = pd.concat(frames, ignore_index=True)
    return long[long["allele"].ne("NA") & long["allele"].notna()]

long = per_allele(ann)
actionable = long[long["significance"].str.contains("Actionable pharmacogenomic", na=False)]
print(actionable[["sample_id", "locus", "allele", "pharm_drugs", "frequency", "tier"]]
      .to_string(index=False))
```

```
sample_id locus     allele   pharm_drugs frequency     tier
  SUBJ001 HLA-B B*57:01:01      abacavir  0.031041 MODERATE
  SUBJ002 HLA-B B*58:01:01   allopurinol  0.012617 MODERATE
  SUBJ002 HLA-B B*15:02:01 carbamazepine  0.000348      LOW
```

`per_allele()` is reused in §3.2 and §4.

Widen the net by also matching `Strong pharmacogenomic risk
association` (the CPIC 1A alleles that are not SJS/TEN-level, such as
A\*31:01 and B\*13:01):

```python
mask = long["significance"].str.contains(
    "Actionable pharmacogenomic|Strong pharmacogenomic", na=False)
```

The label glossary is in [tier_definitions.md](tier_definitions.md).

### 3.2 Disease-association screening

**Shell** — every curated disease association in the cohort, attributed
to the allele that carries it:

```bash
grep -v '^#' /tmp/ex_out/cohort_report/hlante_report.tsv \
| awk -F'\t' 'NR>1 {n=split($26,s,"|");
    for(i=1;i<=n;i++) if (s[i] ~ /Curated/) {
      a=(i==1?$3:$4); m=s[i]; sub(/.*Curated /,"",m); sub(/ \[.*/,"",m);
      print $1"\t"$2"\t"a"\t"m}}'
```

```
SUBJ001  HLA-B     B*57:01:01     pathogenic: Abacavir hypersensitivity
SUBJ001  HLA-B     B*08:01:01     risk factor: Myasthenia gravis
SUBJ001  HLA-C     C*06:02:01     risk factor: Psoriasis
SUBJ001  HLA-DRB1  DRB1*15:01:01  risk factor: Multiple sclerosis
SUBJ001  HLA-DRB1  DRB1*03:01:01  risk factor: Type 1 diabetes mellitus
SUBJ001  HLA-DQB1  DQB1*06:02:01  risk factor: Narcolepsy with cataplexy
SUBJ001  HLA-DQB1  DQB1*02:01:01  risk factor: Celiac disease
SUBJ002  HLA-B     B*58:01:01     pathogenic: Allopurinol-induced SJS/TEN
SUBJ002  HLA-B     B*15:02:01     pathogenic: Carbamazepine-induced SJS/TEN
SUBJ002  HLA-DRB1  DRB1*09:01:02  risk factor: Type 1 diabetes mellitus
SUBJ003  HLA-B     B*27:05:02     risk factor: Ankylosing spondylitis
SUBJ004  HLA-A     A*31:01:02     likely pathogenic: Carbamazepine-induced hypersensitivity
SUBJ004  HLA-B     B*51:01:01     risk factor: Behcet disease
SUBJ004  HLA-DRB1  DRB1*04:01:01  risk factor: Type 1 diabetes mellitus
SUBJ004  HLA-DRB1  DRB1*04:04:01  risk factor: Rheumatoid arthritis
SUBJ004  HLA-DQB1  DQB1*03:02:01  risk factor: Celiac disease
SUBJ006  HLA-B     B*13:01:01     likely pathogenic: Dapsone hypersensitivity syndrome
SUBJ006  HLA-DRB1  DRB1*04:05:01  risk factor: Type 1 diabetes mellitus
SUBJ006  HLA-DRB1  DRB1*15:01:01  risk factor: Multiple sclerosis
SUBJ006  HLA-DQB1  DQB1*06:02:01  risk factor: Narcolepsy with cataplexy
```

To screen for one disease, filter the same column:

```bash
grep -v '^#' /tmp/ex_out/cohort_report/hlante_report.tsv \
| awk -F'\t' 'tolower($26) ~ /ankylosing spondylitis/ {print $1"\t"$2"\t"$3"\t"$4}'
```

**pandas — cohort carrier counts.** Inside one allele's summary,
individual entries are separated by `"; "`; between alleles by `|`:

```python
import pandas as pd

REPORT = "/tmp/ex_out/cohort_report/hlante_report.tsv"
ann = pd.read_csv(REPORT, sep="\t", comment="#", dtype=str)

frames = []
for slot in (0, 1):
    frames.append(pd.DataFrame({
        "sample_id": ann["sample_id"],
        "locus": ann["locus"],
        "allele": ann[f"allele{slot + 1}"],
        "disease": ann["disease_risk_summary"].str.split("|").str[slot],
    }))
long = pd.concat(frames, ignore_index=True)
long = long[long["allele"].ne("NA") & long["allele"].notna()]
long = long.assign(disease=long["disease"].str.split("; ")).explode("disease")

curated = long[long["disease"].str.startswith("Curated", na=False)].copy()
curated["entry"] = curated["disease"].str.replace(r"^Curated [a-z ]+: ", "", regex=True)
curated["entry"] = curated["entry"].str.replace(r" \[.*", "", regex=True)

carriers = (curated.groupby("entry")["sample_id"].nunique()
            .sort_values(ascending=False)
            .rename("n_carriers")
            .to_frame()
            .assign(carrier_pct=lambda d:
                    (100 * d["n_carriers"] / ann["sample_id"].nunique()).round(1)))
print(carriers.to_string())
```

```
                                        n_carriers  carrier_pct
entry
Type 1 diabetes mellitus                         4         66.7
Celiac disease                                   2         33.3
Multiple sclerosis                               2         33.3
Narcolepsy with cataplexy                        2         33.3
Abacavir hypersensitivity                        1         16.7
Allopurinol-induced SJS/TEN                      1         16.7
Ankylosing spondylitis                           1         16.7
Behcet disease                                   1         16.7
Carbamazepine-induced SJS/TEN                    1         16.7
Carbamazepine-induced hypersensitivity           1         16.7
Dapsone hypersensitivity syndrome                1         16.7
Myasthenia gravis                                1         16.7
Psoriasis                                        1         16.7
Rheumatoid arthritis                             1         16.7
```

To join these counts onto phenotype or ancestry metadata — and for the
R/tidyverse equivalent of the join — see
[COHORT_METADATA.md §3, Step 2](../COHORT_METADATA.md).

Three interpretation caveats apply to every table above:

- **`disease_risk_summary` mixes two sources.** Entries prefixed
  `Curated ...` come from the built-in, manually curated table; entries
  prefixed `Strong/Moderate/Reported/Inverse association` come from the
  GWAS Catalog, whose HLA rows are frequently matched at subtype or
  locus scope rather than to your exact allele. Check
  `gwas_annotation_scope` (column 18) and `gwas_fallback_expansion`
  (column 19) before quoting a GWAS effect size for a specific allele.
- **Alleles are annotated independently.** HLAnte does not compute
  DQ2.5/DQ8 heterodimers or cis/trans configuration, so a "Celiac
  disease" flag on DQB1\*02:01 is an allele-level signal, not a DQ2.5
  call.
- **Absence is not negative.** An allele with no hit is unqueried or
  uncovered, not established as risk-free.

---

## 4. Population-stratified annotation

`--population` (`-p`) selects the AFND population group used for the
`allele_frequency` column and for the frequency term of the confidence
heuristic. Codes are `EUR`, `AFR`, `EAS` (alias `ASN`), `SAS`, `MID`,
`AMR`, `OCE` and `global`.

Run the same cohort once per group, keeping the reports apart with
`--prefix`:

```bash
for pop in EUR AFR EAS; do
  hlante annotate \
    -i /tmp/ex_out/cohort \
    -t arcashla \
    -o /tmp/ex_out/pop/$pop \
    --prefix cohort_$pop \
    -p $pop \
    --format tsv --overwrite \
    --pharmgkb-dir tests/fixtures/pharmgkb
done
```

Collect the per-group frequencies into one table:

```python
import pandas as pd

frames = []
for pop in ("EUR", "AFR", "EAS"):
    df = pd.read_csv(f"/tmp/ex_out/pop/{pop}/cohort_{pop}.tsv",
                     sep="\t", comment="#", dtype=str)
    for slot in (0, 1):
        frames.append(pd.DataFrame({
            "population": pop,
            "allele": df[f"allele{slot + 1}"],
            "frequency": pd.to_numeric(
                df["allele_frequency"].str.split("|").str[slot], errors="coerce"),
        }))

long = pd.concat(frames, ignore_index=True)
long = long[long["allele"].ne("NA") & long["allele"].notna()]

wide = (long.drop_duplicates(["population", "allele"])
            .pivot(index="allele", columns="population", values="frequency"))
print(wide.loc[["B*15:02:01", "B*58:01:01", "B*57:01:01", "B*13:01:01",
                "A*31:01:02", "DRB1*04:05:01"]].to_string())
```

```
population          AFR       EAS       EUR
allele
B*15:02:01     0.006039  0.075574  0.000348
B*58:01:01     0.047715  0.066517  0.012617
B*57:01:01     0.010182  0.011378  0.031041
B*13:01:01     0.004827  0.058298  0.001272
A*31:01:02     0.003800  0.017252  0.030817
DRB1*04:05:01  0.019923  0.072100  0.012579
```

This is the number that decides whether a pre-prescription screening
programme is worth running. B\*15:02 is ~220-fold more frequent in the
EAS group than in the EUR group in these AFND data, which is why
carbamazepine screening is guideline-recommended in some populations
and not in others; B\*57:01 runs the other way.

Practical points:

- `--population` changes only the frequency lookup and the frequency
  term of the confidence heuristic. GWAS, PharmGKB and curated-table
  annotations are **not** filtered by ancestry — an association
  established in one population is still reported for a sample from
  another. The curated table records the study population in its
  `[... ; East Asian; ...]` bracket; read it before generalising.
- `global` aggregates every AFND population weighted by sample size.
  Use it when self-reported ancestry is unavailable, and record that
  choice.
- AFND coverage is uneven; a missing frequency is reported as `NA`
  with the reason code `freq_unknown` in `input_quality_rationale`, not as
  zero.
- For splitting a cohort by an ancestry column in a metadata manifest,
  the loop is written out in
  [COHORT_METADATA.md §3, Step 3](../COHORT_METADATA.md).

---

## 5. Reproducibility: pinning the reference release

### 5.1 Pin IPD-IMGT/HLA

`hlante db-update --db imgt` tracks the moving `Latest` branch of the
ANHIG/IMGTHLA mirror. For an analysis you intend to publish, pin an
explicit release with `--imgt-ref` and keep it in its own directory:

```bash
hlante db-update --db imgt --imgt-ref 3640 --imgt-dir /tmp/ex_out/imgt_3640
```

```
Downloading IMGT file: https://raw.githubusercontent.com/ANHIG/IMGTHLA/3640/Allelelist.txt → /tmp/ex_out/imgt_3640/Allelelist.txt
Downloading IMGT file: https://raw.githubusercontent.com/ANHIG/IMGTHLA/3640/wmda/hla_nom_g.txt → /tmp/ex_out/imgt_3640/hla_nom_g.txt
Downloading IMGT file: https://raw.githubusercontent.com/ANHIG/IMGTHLA/3640/wmda/hla_nom_p.txt → /tmp/ex_out/imgt_3640/hla_nom_p.txt
IPD-IMGT/HLA download complete (version=IPD-IMGT/HLA 3.64.0) → /tmp/ex_out/imgt_3640
✓ IPD-IMGT/HLA updated → /tmp/ex_out/imgt_3640
```

> **`--imgt-ref` takes a git ref of the ANHIG/IMGTHLA repository, not a
> release number.** Release 3.64.0 lives on the branch `3640` (the tag
> `v3.64.0-alpha` also resolves). Passing `3.64.0` fails with
> `HTTP Error 404: Not Found`. Check the branch and tag list at
> <https://github.com/ANHIG/IMGTHLA> before pinning.

The download records what it fetched, with per-file SHA-256 checksums:

```bash
cat /tmp/ex_out/imgt_3640/version.json
```

```json
{
  "version": "IPD-IMGT/HLA 3.64.0",
  "downloaded_at": "2026-07-28T09:16:18.017507+00:00",
  "source_base": "https://raw.githubusercontent.com/ANHIG/IMGTHLA/3640",
  "ref": "3640",
  "sha256": {
    "Allelelist.txt": "d36d917890e482def53c640fed50f5e18b556b4269ae627d2fe5bd8253eee40f",
    "hla_nom_g.txt": "10e8838d154f093d5b2286f1f74cf974bb16c2b1934bc6d900044e3bccd62094",
    "hla_nom_p.txt": "c82ecb4fbce00be3c79dcf1f6987fcec7878c60ef1b391b3d5cc22a694cf26b6"
  }
}
```

### 5.2 Annotate against the pinned release

```bash
hlante annotate -i /tmp/ex_out/cohort/SUBJ001.genotype.json -t arcashla \
  --imgt-db-path /tmp/ex_out/imgt_3640 \
  -o /tmp/ex_out/pinned --format all --overwrite -p EUR
```

```
IPD-IMGT/HLA loaded: 46005 allele(s), 777 G-group(s), 1554 P-group(s) (version=IPD-IMGT/HLA 3.64.0)
✓ Report(s) written:
  tsv: /tmp/ex_out/pinned/hlante_report.tsv
  markdown: /tmp/ex_out/pinned/hlante_report.md
  json: /tmp/ex_out/pinned/hlante_report.json
```

(The default `~/.hlante/imgt_hla` in this environment holds release
3.65.0 with 46 652 alleles; the pinned directory gives 46 005. The
counts are a quick check that the pin took effect.)

### 5.3 Where the provenance is recorded

**In the TSV**, as the `#` block before the header row:

```bash
grep '^#' /tmp/ex_out/pinned/hlante_report.tsv
```

```
# HLAnte 0.2.0
# Generated at: 2026-07-28T09:16:37.140220+00:00
# gwas_cache_date version: 2026-07-28
# imgt version: IPD-IMGT/HLA 3.64.0
# pharmgkb version: local
# input_source: typing_tool
# Disclaimer: RESEARCH USE ONLY. …
```

**In the JSON**, under `meta` — same versions, plus the exact command
line that produced the file:

```json
{
  "hlante_version": "0.2.0",
  "generated_at": "2026-07-28T09:16:37.140220+00:00",
  "db_versions": {
    "imgt": "IPD-IMGT/HLA 3.64.0",
    "pharmgkb": "local",
    "gwas_cache_date": "2026-07-28"
  },
  "cli_invocation": "hlante annotate -i /tmp/ex_out/cohort/SUBJ001.genotype.json -t arcashla --imgt-db-path /tmp/ex_out/imgt_3640 -o /tmp/ex_out/pinned --format all --overwrite -p EUR",
  "disclaimer": "RESEARCH USE ONLY. …"
}
```

`cli_invocation` is the field to quote in a methods section: it is the
verbatim command line, so the typing tool, the population group, the
output format and the pinned database path are all recoverable from one
string.

Two caveats on the header block. `pharmgkb version: local` records only
that the local-directory code path was used — it is not a PharmGKB
release identifier, and it appears even when the PharmGKB directory is
empty. `gwas_cache_date` is the modification date of the local GWAS
dump, not a GWAS Catalog release number. Archive the dumps themselves
if the analysis must be re-run exactly.

**In the Markdown report**, as the header bullet list:

```
# HLAnte Research Annotation Report

- **HLAnte version**: 0.2.0
- **Generated at**: 2026-07-28T09:16:37.140220+00:00
- **gwas_cache_date version**: 2026-07-28
- **imgt version**: IPD-IMGT/HLA 3.64.0
- **pharmgkb version**: local
```

### 5.4 Reproducibility checklist for a manuscript

- Pin IPD-IMGT/HLA with `--imgt-ref` and archive `version.json`
  (it carries the SHA-256 of all three reference files).
- Archive the GWAS Catalog and PharmGKB dumps used, or at minimum
  record `gwas_cache_date` from the report header. Both resources
  change continuously and are **not** versioned by release number.
- Keep the JSON report: `cli_invocation` reconstructs the run, and the
  JSON retains full p-value precision (the TSV rounds to four decimal
  places, so p-values below 10⁻⁴ appear as `0.0000`).
- Record the AFND table in use — the bundled fallback and the
  downloaded mirror give different frequencies.
- The full checklist, including manifest and consent items, is in
  [COHORT_METADATA.md §6](../COHORT_METADATA.md).

---

## Appendix A — the demonstration cohort

Recipes §2.1, §3 and §4 use six synthetic arcasHLA files carrying
alleles chosen to exercise the pharmacogenomic and disease paths.
Reproduce them with:

```bash
mkdir -p /tmp/ex_out/cohort && cd /tmp/ex_out/cohort
python - <<'PY'
import json, pathlib
cohort = {
 "SUBJ001": {"HLA-A": ["A*01:01:01","A*02:01:01"], "HLA-B": ["B*57:01:01","B*08:01:01"],
             "HLA-C": ["C*06:02:01","C*07:01:01"], "HLA-DRB1": ["DRB1*15:01:01","DRB1*03:01:01"],
             "HLA-DQB1": ["DQB1*06:02:01","DQB1*02:01:01"]},
 "SUBJ002": {"HLA-A": ["A*11:01:01","A*24:02:01"], "HLA-B": ["B*58:01:01","B*15:02:01"],
             "HLA-C": ["C*03:02:01","C*08:01:01"], "HLA-DRB1": ["DRB1*09:01:02","DRB1*12:02:01"],
             "HLA-DQB1": ["DQB1*03:03:02","DQB1*03:01:01"]},
 "SUBJ003": {"HLA-A": ["A*02:01:01","A*29:02:01"], "HLA-B": ["B*27:05:02","B*44:03:01"],
             "HLA-C": ["C*01:02:01","C*16:01:01"], "HLA-DRB1": ["DRB1*01:01:01","DRB1*07:01:01"],
             "HLA-DQB1": ["DQB1*05:01:01","DQB1*02:02:01"]},
 "SUBJ004": {"HLA-A": ["A*31:01:02","A*03:01:01"], "HLA-B": ["B*35:01:01","B*51:01:01"],
             "HLA-C": ["C*04:01:01","C*14:02:01"], "HLA-DRB1": ["DRB1*04:01:01","DRB1*04:04:01"],
             "HLA-DQB1": ["DQB1*03:02:01","DQB1*03:01:01"]},
 "SUBJ005": {"HLA-A": ["A*02:01:01","A*68:01:02"], "HLA-B": ["B*07:02:01","B*40:01:02"],
             "HLA-C": ["C*07:02:01","C*03:04:01"], "HLA-DRB1": ["DRB1*11:01:01","DRB1*13:01:01"],
             "HLA-DQB1": ["DQB1*03:01:01","DQB1*06:03:01"]},
 "SUBJ006": {"HLA-A": ["A*33:03:01","A*24:02:01"], "HLA-B": ["B*13:01:01","B*46:01:01"],
             "HLA-C": ["C*03:04:01","C*01:02:01"], "HLA-DRB1": ["DRB1*04:05:01","DRB1*15:01:01"],
             "HLA-DQB1": ["DQB1*04:01:01","DQB1*06:02:01"]},
}
for sid, g in cohort.items():
    pathlib.Path(f"{sid}.genotype.json").write_text(
        json.dumps({"sample": sid, "alleles": g}, indent=2) + "\n")
PY
```

The mixed-tool cohort in §2.2 and the deliberately broken batch in §2.3
are staged from the shipped fixtures:

```bash
FIX=tests/fixtures
mkdir -p /tmp/ex_out/mixed/arcashla /tmp/ex_out/mixed/t1k \
         /tmp/ex_out/mixed/hlahd /tmp/ex_out/mixed/optitype
cp "$FIX/sample.genotype.json"        /tmp/ex_out/mixed/arcashla/PT01.genotype.json
cp "$FIX/sample_nested.genotype.json" /tmp/ex_out/mixed/arcashla/PT02.genotype.json
cp "$FIX/sample_t1k_genotype.tsv"     /tmp/ex_out/mixed/t1k/PT03_genotype.tsv
cp "$FIX/sample_final.result.txt"     /tmp/ex_out/mixed/hlahd/PT04_final.result.txt
cp "$FIX/sample_result.tsv"           /tmp/ex_out/mixed/optitype/PT05_result.tsv

mkdir -p /tmp/ex_out/dirty
cp "$FIX/sample.genotype.json"        /tmp/ex_out/dirty/PT01.genotype.json
cp "$FIX/sample_nested.genotype.json" /tmp/ex_out/dirty/PT02.genotype.json
cp "$FIX/malformed.json"              /tmp/ex_out/dirty/PT03.genotype.json
```

The laboratory-validated example in §1.5:

```bash
mkdir -p /tmp/ex_out/lab
printf '%s\n' '{"sample":"LAB01","alleles":{"HLA-A":["A*01:01","A*02:01"],"HLA-B":["B*57:01","B*08:01"],"HLA-DRB1":["DRB1*03:01","DRB1*15:01"]}}' \
  > /tmp/ex_out/lab/LAB01.genotype.json
```

These are synthetic allele lists for demonstration. Do not use real
patient identifiers in example directories — HLA genotypes are
identifying data.

---

## Research use only

Every HLAnte report carries the disclaimer reproduced in full in
[outputs.md](outputs.md). Nothing in this document, and nothing in an
HLAnte report, constitutes a clinical diagnosis, medical advice or a
prescribing recommendation. The evidence-strength labels are annotation
descriptors, not ACMG/AMP classifications. Any clinical decision based
on an HLA allele must rest on a certified laboratory result interpreted
by a qualified clinician.
