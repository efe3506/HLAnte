# HLAnte TSV Output — Biologist/Clinician Interpretation Guide

> **For research purposes.** HLAnte outputs are not designed or validated for use as a clinical
> decision support tool. Authorized expert opinion is required for clinical application.

---

## Table of Contents

1. [File Structure](#1-file-structure)
2. [Basic Reading Rules](#2-basic-reading-rules)
3. [Column Reference](#3-column-reference)
   - [Sample and Locus Information](#31-sample-and-locus-information-columns-1-6)
   - [IMGT / Molecular Classification](#32-imgt--molecular-classification-columns-7-10)
   - [GWAS Catalog Findings](#33-gwas-catalog-findings-columns-11-19)
   - [Pharmacogenomic Findings](#34-pharmacogenomic-findings-pharm_-columns-20-25)
   - [Summary and Clinical Significance Labels](#35-summary-and-clinical-significance-labels-columns-26-28)
   - [Allele Frequency](#36-allele-frequency-columns-29-30)
   - [Confidence Score](#37-confidence-score-columns-31-33)
4. [clinical_significance Label Glossary](#4-clinical_significance-label-glossary)
5. [Confidence Score Interpretation](#5-confidence-score-interpretation)
6. [Real Data Examples](#6-real-data-examples)
7. [Frequently Asked Questions](#7-frequently-asked-questions)

---

## 1. File Structure

Every HLAnte TSV file consists of three sections:

### 1.1 Metadata Block (lines beginning with `#`)

```
# HLAnte 0.1.0
# Generated at: 2026-05-04T10:35:18.827820+00:00
# gwas_cache_date version: 2026-04-28
# imgt version: IPD-IMGT/HLA 3.64.0
```

| Field | Meaning |
|---|---|
| `HLAnte 0.1.0` | The HLAnte version used |
| `Generated at` | Date/time the file was produced (UTC) |
| `gwas_cache_date` | Download date of the GWAS Catalog data — runs with different dates may yield different findings |
| `imgt version` | IPD-IMGT/HLA database version — allele normalization was performed according to this version |

These lines must be skipped before passing the file to a standard table parser.

### 1.2 Header Row

Lists the column names. Contains 33 columns.

### 1.3 Data Rows

Each row carries all findings for **one sample × one HLA locus**.
If a sample has 20+ loci, that sample is emitted as 20+ rows.

---

## 2. Basic Reading Rules

### 2.1 Delimiters

HLAnte TSV uses three different internal delimiters; each has exactly one meaning:

| Delimiter | Meaning | Where it appears |
|---|---|---|
| `\|` (pipe) | Lists values within a single cell; **order matters** (each value matches the value at the same position in the related other columns) | Almost every multi-value column |
| `,` (comma) | Joins multiple PMIDs of a single record | Only the `pharm_pmids` column |
| `;;` (double semicolon) | Separates the allele1 block from the allele2 block in the `confidence_rationale` column | Only `confidence_rationale` |

### 2.2 The `NA` Value

`NA` = "no data / not applicable". Empty cells are not used; a missing value is always written as `NA`.

### 2.3 Pipe-Matching Rule

The values in the `gwas_traits`, `gwas_p_values`, `gwas_odds_ratios`, and `gwas_pmids` columns are
**positionally aligned**:

```
gwas_traits:        sistemik lupus eritematoz | sklerozan kolanjit
gwas_p_values:      0.0000                   | 0.0000
gwas_odds_ratios:   1.9400                   | 2.8226
gwas_pmids:         28714469                 | 23603763
```

That is, the p-value of the 1st trait is the 1st p-value; the p-value of the 2nd trait is the 2nd p-value.

### 2.4 per-allele Columns

Columns whose names contain `_allele1` or `_allele2` belong to a specific allele.
The "legacy" columns whose names do not contain these, such as `gwas_traits` and `pharm_drugs`,
combine the findings of both alleles (may contain duplication).

---

## 3. Column Reference

### 3.1 Sample and Locus Information (Columns 1–6)

#### `sample_id` — Sample Identifier
The identifier of the patient/sample on which the analysis was performed. HLAnte takes this value directly from the input file.

#### `locus` — HLA Locus
Which HLA gene the result belongs to. Examples: `HLA-A`, `HLA-B`, `HLA-C`, `HLA-DRB1`, `HLA-DQB1`.

HLA genes are divided into two classes:
- **Class I** (`HLA-A`, `HLA-B`, `HLA-C`): Expressed on all nucleated cells; recognized by cytotoxic T cells. Critical for organ transplantation compatibility and drug hypersensitivity.
- **Class II** (`HLA-DRB1`, `HLA-DQA1`, `HLA-DQB1`, `HLA-DPA1`, `HLA-DPB1`, etc.): Expressed on antigen-presenting cells; especially important for autoimmune disease susceptibility.

#### `allele1` / `allele2` — The Two Alleles

Every individual has two alleles at each locus (one from the mother, one from the father). Standard HLA nomenclature:

```
B*57:01:01
│  │  │  └─ 4th field (synonymous codon difference)
│  │  └──── 3rd field (non-coding nucleotide difference)
│  └─────── 2nd field (amino acid change present)
└─────────── Gene (B locus)
```

As the number of fields increases, specificity increases. In practice, what is most important for the clinician is the **2nd field (serological equivalent)** and **3rd field (amino acid level)** information.

- If the same allele is written in both copies (**homozygous**): e.g. `B*57:01 / B*57:01`
- If there are different alleles (**heterozygous**): e.g. `B*57:01 / B*15:02`

#### `resolution` — Typing Resolution

| Value | How many fields | Example | Meaning |
|---|---|---|---|
| `2-field` | 2 | `A*02:01` | Serological group level; ambiguous at the amino acid level |
| `4-field` | 3–4 | `A*02:01:01` | Protein (amino acid) level; field missing |
| `6-field` | 5–6 | `B*57:01:01` | Full resolution at the exon level |
| `8-field` | 7–8 | `B*57:01:01:01` | Full genomic resolution |

Low resolution (`2-field`) means that the finding is ambiguous and lowers the `confidence_score`.

#### `tool` — HLA Typing Tool

| Value | Tool Name | Input Type |
|---|---|---|
| `arcashla` | arcasHLA | RNA-seq / WES |
| `t1k` | T1K | WGS / WES |
| `hlahd` | HLA-HD | WGS / WES |
| `optitype` | OptiType | WES / RNA-seq (Class I only) |

---

### 3.2 IMGT / Molecular Classification (Columns 7–10)

#### `imgt_accession` — IMGT Accession Code

The unique identifier in the IPD-IMGT/HLA database (e.g. `HLA22205`). If `NA`, this allele
is not defined in the IMGT version used (it may be a new or rare variant).

#### `hla_class` — HLA Class

`I` or `II`. Automatically derived from the locus.

#### `hla_serotype` — Serological Type

WHO serotype; corresponds to the older serological typing system (e.g. `B57`, `DR15`, `Cw7`).
Pipe-separated (`allele1|allele2`). If serological information is not available, `NA`.

#### `protein_group` — G/P Group

The G-group to which the allele belongs (alleles with the same exon 2–3 sequence) or the P-group.
E.g. `B*57:01:01G`. Can be used in transplantation compatibility assessment.

---

### 3.3 GWAS Catalog Findings (Columns 11–19)

> **Important:** GWAS findings are **population-level statistical associations**, not individual diagnosis
> or risk prediction. The appearance of an allele in GWAS does not prove that the person has or will
> develop that disease.

#### `gwas_traits_allele1` / `gwas_traits_allele2` — Allele-Specific GWAS Traits *(Preferred)*

For each allele, the list of traits for which a genome-wide significant (p < 5×10⁻⁸) association
has been detected in the GWAS Catalog. Pipe-separated. If `NA`, there are no findings for that allele.

Example:
```
gwas_traits_allele1: systemic lupus erythematosus|sclerosing cholangitis|tuberculosis susceptibility
gwas_traits_allele2: NA
```
→ The 1st allele (B*08:01) is associated with three different diseases in GWAS; no findings for the 2nd allele.

#### `gwas_traits` — Combined GWAS Traits *(Legacy; contains duplication)*

The union of the GWAS findings of both alleles. This column does not indicate which finding belongs
to which allele; it is used only for an aggregate overview.

#### `gwas_p_values` — p-Values

The statistical significance of the GWAS association. Positionally aligned with `gwas_traits`.

- Written with 4 decimal places in the TSV; values smaller than 1×10⁻⁴ appear as `0.0000`.
- Use the JSON output for the exact value.

#### `gwas_odds_ratios` — Odds Ratio / Beta Coefficients

| Value range | Interpretation |
|---|---|
| OR > 1 | Disease/trait more frequent in the presence of the allele; positive association |
| OR < 1 | Disease/trait less frequent in the presence of the allele; protective association |
| OR ≈ 1 | No practical effect |
| Very large OR (e.g. 70–112) | Normal for quantitative traits such as immunoglobulin measurements; requires careful assessment for disease |

For quantitative traits (e.g. "beta-2 microglobulin measurement"), a β coefficient may be reported and
is interpreted differently from an OR.

#### `gwas_pmids` — GWAS Publication PMIDs

The PubMed article identifier. You can look up the relevant study at
`https://pubmed.ncbi.nlm.nih.gov/<PMID>`.

#### `gwas_annotation_resolution` — GWAS Match Resolution

| Value | Meaning |
|---|---|
| `4-field` | The GWAS match was made at the full protein level |
| `2-field` | A match was found only at the serological group level |
| `none` | No GWAS record was found for this allele |

#### `gwas_annotation_scope` — GWAS Scope Breadth

| Value | Meaning |
|---|---|
| `allele` | Full allele match |
| `subtype` | Derived from the GWAS record of another subtype in the same protein group |
| `locus` | Match only at the gene level; a finding belonging to the entire locus, not the specific allele |

`locus` or `subtype` values indicate that the finding is attributed to a broader group rather than this
specific allele — interpretation should be more cautious.

#### `gwas_fallback_expansion` — Fallback Expansion

The number of alleles sharing the IMGT prefix used for the match.
- `1` → exact (specific) match
- `> 1` → a broad match covering more than one allele

---

### 3.4 Pharmacogenomic Findings (`pharm_*` Columns 20–25)

> **Clinical Importance:** The findings in this section are the most important outputs that can directly
> affect drug prescribing decisions. The `pharm_cpic_action` column in particular should be reviewed first.

#### `pharm_drugs_allele1` / `pharm_drugs_allele2` — Allele-Specific Drugs *(Preferred)*

The names of drugs that have a clinical annotation for that allele in PharmGKB. Pipe-separated.

#### `pharm_drugs` — Combined Drug List *(Legacy)*

The union of the drug findings of both alleles. It is not indicated which drug belongs to which allele.

#### `pharm_evidence` — PharmGKB Evidence Level

PharmGKB's own evidence grading system:

| Level | Description |
|---|---|
| `1A` | Single allele + CPIC Level A; strongest evidence |
| `1B` | Single allele + CPIC Level B |
| `2A` | Multiple alleles; strong evidence |
| `2B` | Multiple alleles; moderate evidence |
| `3 (low evidence)` | Limited case series or mechanistic study |

#### `pharm_cpic_action` — CPIC Recommended Clinical Action ⚠️

This column contains the standard CPIC actions directed at the clinician for the relevant allele-drug interaction.
Pipe-separated; positionally aligned with `pharm_drugs`.

| Value | Meaning |
|---|---|
| `Contraindicated` | This drug is contraindicated for this allele (e.g. risk of SJS/TEN) |
| `Test required before prescribing` | Testing is required before prescribing |
| `Use with caution — monitor` | Use with caution; close monitoring required |
| `Alternative therapy recommended` | Alternative therapy is recommended |
| `No specific action required` | No specific clinical action is required |
| `NA` | No CPIC action is defined for this drug |

**Example:**
```
allele1:         B*40:01:02
pharm_drugs:     carbamazepine
pharm_cpic_action: Contraindicated
```
→ In a patient carrying the B*40:01 allele, carbamazepine is contraindicated
(due to the risk of SJS/TEN/DRESS).

#### `pharm_pmids` — Pharmacogenomic PMIDs

Comma-separated (within a single record) and pipe-separated between records. Long lists
may be truncated: `"18256392,19001001 (+41 more)"` — use the PharmGKB site for the full list.

---

### 3.5 Summary and Clinical Significance Labels (Columns 26–28)

#### `disease_risk_summary` — Disease Association Summary

Human-readable disease information for each allele. Pipe-separated. Contents:

- Derived from a GWAS finding: `Strong association: sklerozan kolanjit (OR=2.82)`
- Derived from the built-in curated table: `Curated risk factor: SLE [OR=3.2; European; Fernando 2007]`
- If no findings: `No disease association reported`

#### `drug_response_summary` — Drug Response Summary

Readable drug response information for each allele. Pipe-separated. Example:
```
carbamazepine DRESS;SJS;TEN: 2A evidence
```

#### `clinical_significance` — Clinical Significance Label

A single summary label for each allele. **Also described separately below (Section 4).**

---

### 3.6 Allele Frequency (Columns 29–30)

#### `allele_frequency` — Population Allele Frequency

The frequency value taken from the AFND (Allele Frequency Net Database) database.
6 decimal places. E.g. `0.061500` → this allele is observed at a frequency of 6.15% in the relevant population.

If `NA`, there is no frequency record for this allele in AFND (it may be a rare or understudied allele).

#### `allele_freq_population` — Frequency Reference Population

| Code | Population |
|---|---|
| `EUR` | European |
| `AFR` | African |
| `EAS` | East/Southeast Asian (alias: `ASN`) |
| `SAS` | South Asian |
| `MID` | Middle East |
| `AMR` | Amerindian |
| `OCE` | Oceanian |
| `global` | Weighted average of all populations (default) |

Always read the frequency value together with this column to know which population it belongs to.

---

### 3.7 Confidence Score (Columns 31–33)

#### `confidence_score` — Confidence Score

A numerical score between 0.0 and 1.0. Each allele starts at 1.0; uncertainty factors
pull the score down. Pipe-separated (allele1|allele2).

#### `confidence_tier` — Confidence Tier

| Tier | Score range | Interpretation |
|---|---|---|
| `HIGH` | ≥ 0.85 | Well-characterized allele; ≥ 4-field; known frequency |
| `MODERATE` | 0.70 – 0.84 | Minor uncertainty (ambiguity or unknown frequency) |
| `LOW` | < 0.70 | Significant uncertainty (2-field, novel, or very rare allele) |

#### `confidence_rationale` — Reasons for Confidence Reduction

The allele1 and allele2 sections are separated by `;;`; within each section the reasons are listed with `|`.

| Code | Meaning |
|---|---|
| `standard` | No penalty applied |
| `novel_allele` | Not defined in IMGT |
| `rare_allele(freq=X)` | Frequency < 0.001 (very rare) |
| `uncommon_allele(freq=X)` | Frequency between 0.001–0.01 |
| `freq_unknown` | No frequency record in AFND |
| `low_resolution(2-field)` | 2-field resolution |
| `medium_resolution(4-field)` | 4-field resolution |
| `ambiguous` | Subtype ambiguity that the typing tool could not distinguish |

---

## 4. `clinical_significance` Label Glossary

> These labels are **evidence-strength descriptors**; they are not an ACMG/AMP classification.
> They cannot be used as a clinical pathology classification.

| Label | What it means |
|---|---|
| **`Actionable pharmacogenomic risk (CPIC 1A — avoid)`** | There is a "pathogenic" record for this allele in HLAnte's built-in curated table. Drug-allele pairs with CPIC Level A "AVOID" or SJS/TEN risk. *E.g.: B\*57:01 + abacavir, B\*58:01 + allopurinol, B\*15:02 + carbamazepine.* |
| **`Strong pharmacogenomic risk association`** | There is a "likely pathogenic" record in the curated table. Severe drug reaction; CPIC Level A. *E.g.: A\*31:01 + carbamazepine.* |
| **`Suggestive risk factor`** | There is strong PharmGKB evidence (1A/1B) or a genome-wide significant association in GWAS; but there is no pathogenic/likely-pathogenic record in the curated table. A population-level risk signal. |
| **`Inconclusive evidence`** | The allele is ambiguous or unknown in IMGT, and no validation source could be found. The result is insufficient to the point of being uninterpretable. |
| **`No reported risk`** | An allele known in IMGT with no ambiguity; queried in at least one database and no signal was returned. |
| **`Not assessed — insufficient coverage`** | An allele known in IMGT, but all queries returned zero — cannot distinguish between "clean" and "no data". |
| **`Not in IMGT`** | The allele name is valid in format but is not defined in this IMGT version. It may be a new or rare allele; cross-referencing cannot be performed. |

### Quick Decision Guide

```
What is the clinical_significance value?
│
├─ Actionable pharmacogenomic risk (CPIC 1A — avoid)  →  look at pharm_cpic_action; the drug-allele interaction is critical
├─ Strong pharmacogenomic risk association   →  look at pharm_cpic_action; a drug-specific action may be needed
├─ Suggestive risk factor       →  look at gwas_traits and drug_response_summary;
│                                   a population risk signal; use caution in individual interpretation
├─ Inconclusive evidence        →  confidence_tier is usually LOW; the finding is not reliable
├─ No reported risk             →  known allele, no signal in the available databases
├─ Not assessed — insufficient coverage → known allele, insufficient data — not negative evidence
└─ Not in IMGT                  →  novel allele; no cross-referencing can be performed
```

---

## 5. Confidence Score Interpretation

### Applied Penalty Multipliers

| Factor | Multiplier | Effect |
|---|---|---|
| Not defined in IMGT (novel) | × 0.30 | Heaviest penalty |
| Rare allele (frequency < 0.001) | × 0.50 | Heavy penalty |
| Uncommon allele (0.001–0.01) | × 0.80 | Moderate penalty |
| Frequency unknown (not in AFND) | × 0.85 | Light penalty |
| 2-field resolution | × 0.70 | Significant penalty |
| 4-field resolution | × 0.90 | Light penalty |
| Ambiguous allele (`typing_tool` source) | × 0.75 | Moderate penalty |

Penalties are applied **multiplicatively**. Example:

```
4-field resolution + frequency unknown + ambiguous:
1.0 × 0.90 × 0.85 × 0.75 = 0.5738 → LOW
```

### Effect of `input_source`

In analyses run with `--input-source validated` (e.g. 1000 Genomes validated data),
the "ambiguity" penalty is not applied. For this reason, the same allele:

- `typing_tool` (default): low score → `LOW`
- `validated`: higher score → `MODERATE`

**Check the `# input_source` line in the metadata block** — look at this information before
interpreting the scores.

---

## 6. Real Data Examples

### Example A — Pharmacogenomic Warning (HLA-B, T1K result)

```
sample_id:          HLA141_hla_genotype
locus:              HLA-B
allele1:            B*40:01:02
allele2:            B*15:09:01
resolution:         6-field
tool:               t1k
clinical_significance: Inconclusive evidence | Suggestive risk factor
pharm_drugs:        carbamazepine | oxcarbazepine
pharm_evidence:     2A | 3 (low evidence)
pharm_cpic_action:  Contraindicated | Contraindicated
drug_response_summary:
  carbamazepine DRESS;SJS;TEN: 2A evidence  |  No drug response reported
confidence_tier:    MODERATE | HIGH
```

**Interpretation:**
- **allele2 (B*15:09):** `HIGH` confidence, `Suggestive risk factor`. PharmGKB 2A evidence for
  carbamazepine-associated SJS/TEN risk → **Contraindicated**.
- **allele1 (B*40:01):** `MODERATE` confidence, `Inconclusive evidence`. Contraindicated for
  oxcarbazepine with low evidence — must be interpreted with caution because it came through 3 (low evidence).
- **Practical conclusion:** In this patient, the HLA-B status should be evaluated before carbamazepine is prescribed;
  the relevant clinical guideline and an expert should be consulted.

---

### Example B — GWAS Risk, No Drug Signal (HLA-B, arcasHLA result)

```
sample_id:          Tes503_Aligned
locus:              HLA-B
allele1:            B*08:01:01
allele2:            B*08:01:01
resolution:         6-field
clinical_significance: Suggestive risk factor | Suggestive risk factor
gwas_traits_allele1:   systemic lupus erythematosus | neuromyelitis optica | myositis | HIV-1 infection ...
gwas_annotation_scope: subtype | subtype
pharm_drugs:        infliximab | Antithyroid Preparations
pharm_evidence:     3 (low evidence) | 3 (low evidence)
pharm_cpic_action:  NA
confidence_tier:    MODERATE | MODERATE
confidence_rationale: ambiguous ;; ambiguous
```

**Interpretation:**
- Homozygous (both alleles the same: B*08:01).
- There is a GWAS association with multiple autoimmune diseases; however, the evidence is at the `subtype`
  level — not the specific subtype, but the GWAS data of the same protein group has been transferred.
- The drug associations are **3 (low evidence)** — low evidence level; no CPIC action defined.
- The `Suggestive risk factor` label can be evaluated in a monitoring or research context
  but does not directly require a clinical decision.
- Both alleles have the `ambiguous` penalty → the typing tool could not fully distinguish the subtype;
  confidence is `MODERATE`.

---

### Example C — Heterozygous, Allele-Specific Different Findings (HLA-C)

```
sample_id:          Tes503_Aligned
locus:              HLA-C
allele1:            C*07:744
allele2:            C*07:02:01
resolution:         4-field
clinical_significance: Not assessed — insufficient coverage | Suggestive risk factor
gwas_traits_allele1:   NA
gwas_traits_allele2:   Merkel cell polyomavirus seropositivity | HIV-1 infection | ...
allele_frequency:   NA | 0.081600
confidence_tier:    MODERATE | MODERATE
confidence_rationale: freq_unknown ;; ambiguous
```

**Interpretation:**
- **allele1 (C*07:744):** No frequency record in AFND (`freq_unknown`) → `Not assessed — insufficient coverage`.
  This does not mean "safe"; it merely means the data is insufficient.
- **allele2 (C*07:02:01):** Global frequency 8.2%; there are several GWAS associations → `Suggestive risk factor`.
- Using pipe-matching, reading `gwas_traits_allele1=NA` and `gwas_traits_allele2=...` clarifies which finding
  belongs to which allele.

---

### Example D — Insufficient Evidence / Novel Allele (HLA-A)

```
sample_id:          Tes503_Aligned
locus:              HLA-A
allele1:            A*03:02:01
allele2:            A*03:02:01
resolution:         6-field
imgt_accession:     NA
clinical_significance: Inconclusive evidence | Inconclusive evidence
allele_frequency:   NA
confidence_score:   0.6375 | 0.6375
confidence_tier:    LOW | LOW
confidence_rationale: freq_unknown|ambiguous ;; freq_unknown|ambiguous
```

**Interpretation:**
- `imgt_accession = NA` → This allele is not recorded in IMGT 3.64.0, or the IMGT accession code
  could not be matched. Cross-referencing cannot be performed.
- Frequency is also unknown; the ambiguity penalty was also added → score 0.6375 → `LOW`.
- The `Inconclusive evidence` label indicates that this row should not be used for clinical inference.

---

## 7. Frequently Asked Questions

### "All p-values appear as 0.0000 — is this correct?"

Yes. In the TSV format, p-values are written with 4 decimal places; since most GWAS findings
are p < 10⁻⁵ or smaller, they appear as `0.0000`. For the exact value, use the JSON output
(`.json` file) of the same run.

### "There are too many values in the gwas_traits column, which one is important?"

Review the `disease_risk_summary` column first — here HLAnte highlights findings with high OR and
from curated sources. Also, those with `gwas_annotation_scope = allele` are the most specific
matches; interpretation of those with `locus` scope should be more cautious.

### "At one locus allele1 = allele2 — is this homozygous?"

Yes. If HLA typing tools cannot distinguish two copies, they report the same allele twice.
This may be true homozygosity or a failure to distinguish caused by the tool.

### "What does it mean if pharm_cpic_action = NA but pharm_drugs is populated?"

The relevant allele-drug pair is recorded in PharmGKB but no standard action is defined in CPIC.
You need to evaluate the literature yourself by looking at the PharmGKB evidence level (`pharm_evidence`)
and the publication PMIDs (`pharm_pmids`).

### "Why are loci such as DMA, DMB, DOA, DOB always Inconclusive evidence?"

The HLA-DM and HLA-DO genes play a role in the antigen processing process; direct GWAS association
is limited and there is no record for these loci in our curated table. For these loci,
an `Inconclusive evidence` result is an expected output.

### "Can I ignore the rows where confidence_tier = LOW?"

Rows that are `LOW` should not be excluded; however, extra caution should be exercised when drawing
conclusions based on the findings of these rows. A `LOW` score usually corresponds to:
(a) low typing resolution, (b) the allele is rare or its frequency is unknown, or
(c) the allele is not defined in IMGT. For these rows, higher-resolution
typing or validation is generally recommended.

### "arcasHLA and T1K outputs gave different results for the same sample — which one should I trust?"

Tool selection varies according to clinical context and sequencing technology. To measure the
concordance between typing tools, you can use the `scripts/benchmark/inter_format_concordance.py`
script; in discordant findings, the results of both tools should be evaluated by an expert.

---

*This guide corresponds to the HLAnte 0.1.0 TSV schema.*
*The schema is defined by `hlante.reporter.TSV_COLUMNS` (33 columns).*
