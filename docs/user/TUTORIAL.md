# Step-by-Step Tutorial

A guided first session with HLAnte, written for HLA laboratory scientists who
use a terminal only occasionally and who have never used Python.

Nothing here assumes prior programming experience. Every command below is
shown together with the output it actually produced, so you can compare your
screen with the page as you go.

By the end you will have installed HLAnte, downloaded the one database it
requires, annotated a genotype, and read the result.

Allow about 20 minutes.

---

## How to read this page

- Lines you type are shown in a code block. Type the command and press Enter.
- The lines underneath, in the same block after a `# ---- output ----` marker,
  are what HLAnte printed. Your version numbers, dates and file paths will
  differ.
- Where output was long it has been shortened; a shortened section is marked
  with a `# ...` comment.
- HLAnte prints progress messages (`INFO`, `WARNING`) while it works. They are
  normal. Section 10 explains the ones that matter.

---

## 1. What you need before you start

| Requirement | How to check |
|---|---|
| A terminal (macOS Terminal, or any Linux shell) | It is already open if you can see a `$` prompt |
| Python 3.9 or newer | `python3 --version` |
| `git` | `git --version` |
| An internet connection | Needed once, for the database download in section 3 |

```bash
python3 --version
git --version

# ---- output ----
Python 3.9.6
git version 2.50.1 (Apple Git-155)
```

If `python3 --version` reports anything below 3.9, install a newer Python
before continuing; see [INSTALL.md](../../INSTALL.md#1-system-requirements).

---

## 2. Install HLAnte

### 2.1 Download the source code

HLAnte is **not** distributed through PyPI, so `pip install hlante` will not
work (section 10 shows the error). Take a copy of the repository instead:

```bash
git clone https://github.com/efe3506/HLAnte.git
cd HLAnte
```

Everything in the rest of this tutorial is run from inside that `HLAnte`
folder unless stated otherwise.

### 2.2 Create a virtual environment

A virtual environment is a private folder that holds HLAnte and the Python
packages it needs, so that installing HLAnte cannot disturb any other software
on the machine.

```bash
python3 -m venv .venv
source .venv/bin/activate
```

After the second command your prompt gains a `(.venv)` prefix. That prefix is
the only sign that the environment is active. Open a new terminal window later
and you must run `source .venv/bin/activate` again before `hlante` will be
found.

### 2.3 Install

```bash
pip install .

# ---- output (shortened) ----
# ... one line per dependency as pip collects it ...
Building wheels for collected packages: hlante
  Building wheel for hlante (pyproject.toml): started
  Building wheel for hlante (pyproject.toml): finished with status 'done'
  Created wheel for hlante: filename=hlante-0.2.0-py3-none-any.whl size=95416 ...
Successfully built hlante
Installing collected packages: pytz, tzdata, tqdm, six, numpy, click, python-dateutil, pandas, hlante

Successfully installed click-8.1.8 hlante-0.2.0 numpy-2.0.2 pandas-2.3.3 python-dateutil-2.9.0.post0 pytz-2026.3.post1 six-1.17.0 tqdm-4.70.0 tzdata-2026.3
```

### 2.4 Check that it worked

```bash
hlante --log-level WARNING version

# ---- output ----
HLAnte v0.2.0
  IPD-IMGT/HLA : not installed
  PharmGKB     : not installed
  GWAS Catalog : not installed
  AFND         : not installed
```

`HLAnte v0.2.0` means the software is installed. The four `not installed`
lines are expected on a fresh machine — the next section fixes the first one.

`--log-level WARNING` simply suppresses routine progress messages; you can
leave it out.

Other ways to install (conda, the automated `install.sh` script) are described
in [INSTALL.md](../../INSTALL.md).

---

## 3. Download IPD-IMGT/HLA — required, once

HLAnte checks every allele call against a local copy of the IPD-IMGT/HLA
release. Without it, annotation stops immediately and no report is written.
This download is needed **once**, not once per run:

```bash
hlante db-update --db imgt

# ---- output (shortened) ----
# ... the trimmed lines report each file as it is fetched ...
2026-07-28 12:13:05,185 INFO    hlante.db.imgt: IPD-IMGT/HLA download complete (version=IPD-IMGT/HLA 3.65.0) → /home/user/.hlante/imgt_hla
✓ IPD-IMGT/HLA updated → /home/user/.hlante/imgt_hla
```

The files land in `~/.hlante/imgt_hla` and stay there. Confirm:

```bash
hlante --log-level WARNING version

# ---- output (shortened) ----
HLAnte v0.2.0
  IPD-IMGT/HLA : IPD-IMGT/HLA 3.65.0  (/home/user/.hlante/imgt_hla)
# ... the optional databases are listed on the following lines ...
```

The `IPD-IMGT/HLA` line must show a release number before you go on. If it
still says `not installed`, see section 10.

> **Reproducibility.** `db-update --db imgt` tracks the current release. To pin
> a specific one — which you should do for any analysis you intend to publish —
> add `--imgt-ref`, for example
> `hlante db-update --db imgt --imgt-ref 3.64.0`. The release actually used is
> recorded in the header of every report.

### Optional: extra evidence sources

HLAnte runs with IPD-IMGT/HLA alone; the built-in curated HLA–disease/drug
table and a small built-in allele-frequency table are compiled into the
package. Two further downloads widen the evidence base and are used in the
worked example below:

```bash
hlante db-update --db afnd    # full Allele Frequency Net Database extract, ~5 MB
hlante db-update --db gwas    # GWAS Catalog; large — about 0.7 GB on disk
```

PharmGKB, user-supplied NMDP frequency tables, and the disk requirements of
each source are covered in
[INSTALL.md § 5 Database Setup](../../INSTALL.md#5-database-setup).
If you skip these downloads, the corresponding report columns simply read `NA`
and everything else is unchanged.

---

## 4. Prepare an input file

HLAnte does not perform genotyping. It reads the output of a genotyping tool
you have already run, and you must tell it which tool produced the file:

| Genotyping tool | Typical file | `--tool` value |
|---|---|---|
| arcasHLA | `*.genotype.json` | `arcashla` |
| T1K | `*_genotype.tsv` | `t1k` |
| HLA-HD | `*_final.result.txt` | `hlahd` |
| OptiType | `*_result.tsv` | `optitype` |

The exact fields each parser reads are documented in [Inputs](inputs.md).

For this tutorial, create a small arcasHLA-style file so that you can follow
along without your own data. Copy the block below into a text editor and save
it as `example.genotype.json`, or paste the whole thing into the terminal:

```bash
cat > example.genotype.json <<'EOF'
{
  "HLA-A": ["A*01:01", "A*02:01"],
  "HLA-B": ["B*57:01", "B*08:01"],
  "HLA-C": ["C*06:02", "C*07:01"],
  "HLA-DRB1": ["DRB1*03:01", "DRB1*15:01"],
  "HLA-DQB1": ["DQB1*02:01", "DQB1*06:02"]
}
EOF
```

Ready-made example files for all four tools also ship with the repository, in
`tests/fixtures/`.

---

## 5. Check the file before you annotate

`hlante validate` reads the file, reports what it found, and stops. It is
quick, it touches no database, and it is the right first move with any new
file:

```bash
hlante validate -i example.genotype.json -t arcashla

# ---- output ----
✓ Input is valid
  Tool: arcashla
  Files: 1
  Samples: 1
  Locus calls: 5
  Unique loci: HLA-A, HLA-B, HLA-C, HLA-DQB1, HLA-DRB1
```

Read the counts, not just the tick. If a file you expect to contain 12 loci
reports 5, the parser did not see what you assumed it would.

---

## 6. Run your first annotation

```bash
hlante annotate -i example.genotype.json -t arcashla -o results

# ---- output (shortened) ----
  Tool: arcashla | File count: 1
✓ Parse complete: 5 locus call(s) (1 sample(s))
# ... INFO lines: IPD-IMGT/HLA loaded, curated table loaded, GWAS/AFND loaded ...
# ... one WARNING per allele for each database you have not downloaded ...
✓ Annotation complete: 10 record(s) processed.
✓ Report(s) written:
  tsv: results/hlante_report.tsv
  markdown: results/hlante_report.md
  json: results/hlante_report.json
```

The three flags are all you need:

| Flag | Meaning |
|---|---|
| `-i` | the input file (or a directory of input files) |
| `-t` | which genotyping tool produced it |
| `-o` | where to put the reports |

Two options you will want early:

```bash
# restrict the frequency lookup to a population, and write only the TSV
hlante annotate -i example.genotype.json -t arcashla -o results_eur -p EUR --format tsv
```

`-p` accepts `EUR`, `AFR`, `EAS`, `SAS`, `MID`, `AMR`, `OCE` or `global`
(the default). The full option list is in [Quickstart](quickstart.md), or run
`hlante annotate --help`.

---

## 7. Where the results are

```bash
ls -1 results/

# ---- output ----
hlante_report.json
hlante_report.md
hlante_report.tsv
```

All three files describe the same annotation. Pick by what you want to do:

| File | Open it with | Use it when |
|---|---|---|
| `hlante_report.tsv` | Excel, LibreOffice, R, Python | You want one row per sample × locus for filtering, sorting or statistics |
| `hlante_report.md` | Any text editor, or a Markdown viewer | You want to read the findings as prose, with the supporting PMIDs listed under each one |
| `hlante_report.json` | A script, or an API | You want the full nested record, including every source identifier, for downstream software |

Two things to know before you open the TSV in a spreadsheet:

- The file starts with comment lines beginning with `#`. They record the
  HLAnte version, the generation time and — importantly — the database
  releases used:

  ```
  # HLAnte 0.2.0
  # Generated at: 2026-07-28T09:18:09.308613+00:00
  # gwas_cache_date version: 2026-07-28
  # imgt version: IPD-IMGT/HLA 3.65.0
  ```

  Keep them with the file; they are what makes the run reproducible. Skip them
  when importing (in Excel's text-import dialogue, start at the row beginning
  `sample_id`).

- The last comment line is the research-use disclaimer. It is repeated in all
  three files.

The Markdown report also states plainly which loci were *not* typed. For the
example file above it reads:

> ⚠️ **Loci not typed (not assessed for HLA-linked risk):** HLA-DQA1, HLA-DPB1. A locus that was not typed is *indeterminate*, not negative — the absence of an alert does **not** mean the risk allele is absent.

---

## 8. Reading the key columns

The TSV has 40 columns. You do not need all of them on a first look. This
command shows six of them side by side:

```bash
grep -v '^#' results/hlante_report.tsv \
| awk -F'\t' -v OFS='\t' 'NR==1 {for (i=1; i<=NF; i++) c[$i]=i}
    {print $c["locus"], $c["allele1"], $c["allele2"],
           $c["clinical_significance"], $c["allele_frequency"], $c["input_quality_tier"]}' \
| column -t -s $'\t'

# ---- output ----
locus     allele1     allele2     clinical_significance                                                    allele_frequency   input_quality_tier
HLA-A     A*01:01     A*02:01     Suggestive risk factor|Suggestive risk factor                            0.064181|0.151974  limited|limited
HLA-B     B*57:01     B*08:01     Actionable pharmacogenomic risk (CPIC 1A — avoid)|Suggestive risk factor  0.020265|0.040392  limited|limited
HLA-C     C*06:02     C*07:01     Suggestive risk factor|Suggestive risk factor                            0.080365|0.090534  limited|limited
HLA-DRB1  DRB1*03:01  DRB1*15:01  Suggestive risk factor|Suggestive risk factor                            0.067246|0.080539  limited|limited
HLA-DQB1  DQB1*02:01  DQB1*06:02  Suggestive risk factor|Suggestive risk factor                            0.117904|0.078371  limited|limited
```

The `awk` block reads the header row into `c`, so the columns are selected by
name. That is deliberate: the schema has grown from 33 to 40 columns across
releases, and a fixed column number silently returns the wrong column instead
of failing.

(In a spreadsheet, hiding the columns you do not need achieves the same thing.)

**One row is one sample at one locus.** Both alleles of that locus are
reported on the same row.

**The `|` character separates allele 1 from allele 2.** Wherever you see a
pipe, the value on the left belongs to `allele1` and the value on the right to
`allele2`, in that order. The same rule aligns the GWAS columns with each
other: the third trait in `gwas_traits` goes with the third p-value, the third
odds ratio and the third PMID.

The columns worth knowing first:

| Column | What it tells you |
|---|---|
| `locus`, `allele1`, `allele2` | The genotype as HLAnte normalised it against IPD-IMGT/HLA |
| `disease_risk_summary` | One readable sentence per allele, combining GWAS Catalog hits and the built-in curated table |
| `drug_response_summary` | The same, for drug-response findings |
| `clinical_significance` | A single evidence-strength label per allele — **not** an ACMG/AMP class |
| `allele_frequency`, `allele_freq_population` | The AFND frequency, and which population it came from. Always read the two together |
| `input_quality_score`, `input_quality_tier`, `input_quality_rationale` | How completely the allele **call** is characterised |

Two cautions that the reports themselves repeat:

- `clinical_significance` labels are descriptions of evidence strength. They
  are not diagnostic categories, and HLAnte does not implement ACMG/AMP
  criteria.
- The input-quality score is an uncalibrated heuristic reflecting how well the
  allele call is characterised, **not** the certainty or correctness of any
  associated risk. A `limited` tier never means an actionable association may be
  ignored — in the table above every row is `limited` because the calls carry two
  fields and are therefore subtype-ambiguous, while the *HLA-B* row still
  carries the strongest pharmacogenomic flag HLAnte can raise.

Every column, every label and a set of worked TSV rows are documented in the
[TSV Interpretation Guide](TSV_INTERPRETATION_GUIDE.md) — § 3 for the column
reference, § 4 for the `clinical_significance` glossary and § 7 for the
frequently asked questions.

---

## 9. A worked interpretation: the *HLA-B* row

Take the row for *HLA-B* from the table above:

```
HLA-B     B*57:01     B*08:01     Actionable pharmacogenomic risk (CPIC 1A — avoid)|Suggestive risk factor  0.020265|0.040392  limited|limited
```

Now read the disease summary for the same row, one allele per line:

```bash
grep -v '^#' results/hlante_report.tsv \
| awk -F'\t' 'NR==1 {for (i=1; i<=NF; i++) c[$i]=i; next}
    $c["locus"]=="HLA-B" {print $c["disease_risk_summary"]}' | tr '|' '\n'

# ---- output ----
Strong association: drug-induced liver injury (OR=36.62); Curated pathogenic: Abacavir hypersensitivity [OR=960.0; Global; Strong evidence; Mallal 2002 Lancet; Hetherington 2002 Lancet; Mallal 2008 NEJM (PREDICT-1)]
Moderate association: sclerosing cholangitis (OR=2.82); Curated risk factor: Myasthenia gravis [OR=4.2; European; Moderate evidence; Gregersen 2012 Ann Neurol]
```

Reading it line by line:

1. **Allele 1 is HLA-B\*57:01**, and its `clinical_significance` is
   `Actionable pharmacogenomic risk (CPIC 1A — avoid)`. That is the strongest
   label HLAnte assigns, and it is reserved for allele–drug pairs recorded as
   *pathogenic* in the built-in curated table. The reason follows in
   `disease_risk_summary`: **abacavir hypersensitivity**, with the three
   primary references named (Mallal 2002, Hetherington 2002, and the PREDICT-1
   trial, Mallal 2008).

2. **The GWAS finding on the same allele is separate.** `Strong association:
   drug-induced liver injury (OR=36.62)` comes from the GWAS Catalog, not from
   the curated table. The Markdown report attaches the study accession and
   PMID to it, and flags the odds ratio for checking:

   ```
   - **HLA-B*57:01**: drug-induced liver injury (OR=36.62)
     - Source: GWAS Catalog (GCST007606), PMID: 30661239
     - ⚠️ Effect size note: OR=36.62 exceeds 10 — unusually large for a binary-disease odds ratio. Verify the source study (PMID: 30661239).
   ```

3. **Allele 2 is HLA-B\*08:01**, labelled `Suggestive risk factor` — a
   population-level signal (sclerosing cholangitis from GWAS; myasthenia
   gravis from the curated table), not an allele-specific action.

4. **The frequency column puts it in context.** `0.020265` is the global AFND
   frequency of HLA-B\*57:01, i.e. about 2%. Re-running with `-p EUR` returns
   `0.031041` for the same allele — always read `allele_frequency` together
   with `allele_freq_population`.

5. **Both alleles carry `input_quality_tier = limited`.** These calls are named with
   two fields, so more than one IPD-IMGT/HLA allele shares each name. That is a
   statement about the precision of the *call*, and it does not weaken the
   abacavir association.

**What this row is, and is not.** It is a pointer to published evidence for a
research-use annotation: this genotype carries an allele with a very
well-replicated hypersensitivity association. It is not a clinical result. A
laboratory report of HLA-B\*57:01 status for a patient must come from a
validated assay in an accredited laboratory and be interpreted by a qualified
clinician.

---

## 10. If something goes wrong

### `pip install hlante` fails

```bash
pip install hlante

# ---- output ----
ERROR: Could not find a version that satisfies the requirement hlante (from versions: none)
ERROR: No matching distribution found for hlante
```

**Why.** HLAnte is not published on PyPI, so there is nothing for pip to
download by name.

**Fix.** Install from the repository, as in section 2:

```bash
git clone https://github.com/efe3506/HLAnte.git
cd HLAnte
pip install .
```

### Annotation stops with an IPD-IMGT/HLA error and writes no report

```bash
hlante annotate -i example.genotype.json -t arcashla -o results

# ---- output ----
  Tool: arcashla | File count: 1
✓ Parse complete: 5 locus call(s) (1 sample(s))
ERROR: IPD-IMGT/HLA could not be loaded: IPD-IMGT/HLA Allelelist not found: /home/user/.hlante/imgt_hla/Allelelist.txt.
The IPD-IMGT/HLA release is a required, one-time download (~10 MB). Install it with:
    hlante db-update --db imgt
To reproduce a specific release, pin it with `--imgt-ref` (for example: `hlante db-update --db imgt --imgt-ref 3.64.0`).
```

**Why.** The one-time download in section 3 has not been run, so there is no
allele list to normalise against. The message names the directory HLAnte
looked in — `~/.hlante/imgt_hla`, printed as a full path.

**Fix.** Run the command the message names, then annotate again:

```bash
hlante db-update --db imgt
```

Adding `--offline` does **not** avoid this — it means "make no network calls",
not "run without databases", and produces exactly the same error.

### `hlante: command not found`

The virtual environment is not active in this terminal window. Go back to the
`HLAnte` folder and run `source .venv/bin/activate`. Further PATH remedies are
in [INSTALL.md § 8 Troubleshooting](../../INSTALL.md#8-troubleshooting).

### `ERROR: Output file already exists`

```
ERROR: Output file already exists: results/hlante_report.tsv. Pass overwrite=True (or --overwrite on the CLI) to replace it.
```

HLAnte will not silently replace a previous report. Either write to a new
directory with `-o`, or add `--overwrite`.

### `ERROR: Invalid input: ... row has missing columns`

```
ERROR: Invalid input: T1K row has missing columns (line 1): '{' (source: example.genotype.json)
```

The file is fine but `-t` names the wrong tool — here an arcasHLA JSON was
passed as `-t t1k`. Check the table in section 4 and re-run `hlante validate`.

### Warnings about a database you did not download

```
2026-07-28 12:18:09,252 WARNING hlante.annotator: PharmGKB B*57:01 query failed: PharmGKB TSV not found in /home/user/.hlante/pharmgkb. Call update() first or pass the correct local_dir.
```

Expected, and harmless. That source is skipped and the corresponding columns
read `NA`. Download it if you want the columns populated
([INSTALL.md § 5](../../INSTALL.md#5-database-setup)).

---

## 11. Where to go next

- [Quickstart](quickstart.md) — the same workflow in condensed form, plus batch
  annotation of a whole directory and the full flag table
- [Worked examples](EXAMPLES.md) — copy-and-paste recipes for each typing tool,
  cohorts, population stratification and pinned reference releases
- [TSV Interpretation Guide](TSV_INTERPRETATION_GUIDE.md) — every column,
  every label, worked rows and an FAQ
- [Inputs](inputs.md) and [Outputs](outputs.md) — the file formats in detail
- [INSTALL.md](../../INSTALL.md) — conda, the automated installer, database
  setup and troubleshooting
- [Limitations](limitations.md) — what HLAnte deliberately does not do

---

> **Research use only.** Nothing produced by HLAnte constitutes a clinical
> diagnosis, medical advice or a pharmacogenomic recommendation. Any clinical
> decision based on an HLA allele must rely on a certified laboratory result
> interpreted by a qualified clinician.
