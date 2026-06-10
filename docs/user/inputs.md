# Input Files

HLAnte accepts output files from four HLA genotyping tools.
Pass the tool name with `-t` / `--tool` and the file (or directory) with `-i` / `--input`.

---

## arcasHLA

**Flag:** `-t arcashla` (aliases: `-t arcas-hla`)

**Expected extension:** `*.genotype.json`

**Format:** JSON produced by `arcasHLA genotype`. The file must contain a
top-level `"genotype"` key mapping locus names (e.g. `"A"`, `"B"`, `"DRB1"`)
to a list of one or two allele strings.

```json
{
  "sample": "HG00096",
  "genotype": {
    "A": ["A*01:01:01:01", "A*02:01:01:01"],
    "B": ["B*08:01:01:01", "B*57:01:01:01"],
    "C": ["C*07:01:01:01", "C*06:02:01:01"],
    "DRB1": ["DRB1*03:01:01:01", "DRB1*07:01:01:01"],
    "DQB1": ["DQB1*02:01:01:01", "DQB1*03:03:02:01"]
  }
}
```

**Edge cases handled:**
- Trailing asterisk (`DRB1*04:92*`) — stripped automatically; a DEBUG log records each occurrence.
- Space-separated ambiguous pair (`B*49:01 50:01`) — primary call (first token) is used.

---

## T1K

**Flag:** `-t t1k`

**Expected extension:** `*_genotype.tsv`

**Format:** TSV produced by `T1K --genotype`. HLAnte reads the columns
`#gene_name`, `allele1`, `allele2`, and `abundance` (optional).

```
#gene_name	allele1	allele2	abundance
A	A*01:01:01:01	A*02:01:01:01	120.4
B	B*08:01:01:01	B*57:01:01:01	118.2
C	C*07:01:01:01	C*06:02:01:01	109.7
DRB1	DRB1*03:01:01:01	DRB1*07:01:01:01	98.3
```

---

## HLA-HD

**Flag:** `-t hlahd` (alias: `-t hla-hd`)

**Expected extension:** `*_final.result.txt`

**Format:** Tab-separated text produced by HLA-HD. HLAnte reads lines
that are not comment lines (`#`). Each line has locus name, allele1,
allele2 (and optionally a separator).

```
A	A*01:01:01:01	A*02:01:01:01
B	B*08:01:01:01	B*57:01:01:01
C	C*07:01:01:01	C*06:02:01:01
DRB1	DRB1*03:01:01:01	DRB1*07:01:01:01
DQB1	DQB1*02:01:01:01	DQB1*03:03:02:01
```

---

## OptiType

**Flag:** `-t optitype`

**Expected extension:** `*_result.tsv`

**Format:** TSV produced by OptiType. HLAnte reads the first result row.
OptiType reports **Class I loci only** (HLA-A, HLA-B, HLA-C).

```
	A1	A2	B1	B2	C1	C2	Reads	Objective
0	A*01:01	A*02:01	B*08:01	B*57:01	C*07:01	C*06:02	1420	1395.80
```

---

## Allele resolution

HLAnte accepts alleles at any resolution (2-field through 8-field).
The confidence score is penalised for lower-resolution calls:

| Resolution | Example | Penalty |
|------------|---------|---------|
| 8-field | `A*02:01:01:01` | none |
| 6-field | `A*02:01:01` | none |
| 4-field | `A*02:01` | none |
| 2-field | `A*02` | ×0.85 (typing_tool) or none (validated) |

Use `--input-source validated` when the input alleles are Sanger-validated
(e.g. 1000 Genomes types) to suppress the 2-field ambiguity penalty.

---

## Directory input

When `--input` is a directory, HLAnte globs for the tool-appropriate
file extension and processes all matching files in a single run:

```bash
hlante annotate -i /data/arcashla_results/ -t arcashla
```

---

## Supported loci

| Class I | Class II |
|---------|----------|
| HLA-A | HLA-DRB1 |
| HLA-B | HLA-DQB1 |
| HLA-C | HLA-DPB1 |
| | HLA-DQA1 |
| | HLA-DPA1 |

OptiType reports Class I only. For Class II, use arcasHLA, T1K, or HLA-HD.
