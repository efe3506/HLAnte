# Output Files

HLAnte writes three output files per run. All three are written by default
(`--format all`); use `--format tsv`, `--format json`, or `--format markdown`
to write only one.

```
hlante_output/
├── hlante_report.tsv       # main annotation table
├── hlante_report.json      # full nested structure
└── hlante_report.md        # human-readable per-sample summary
```

---

## TSV output

One row per `(sample_id, locus)` pair. Homozygous or single-allele loci
have `NA` in the second-allele cells.

The file begins with a `#`-prefixed metadata block recording the HLAnte
version, generation timestamp, database versions, and research-use disclaimer.
Strip it before parsing: `grep -v '^#' report.tsv`.

### Column reference

| # | Column | Description |
|---|--------|-------------|
| 1 | `sample_id` | Input sample identifier |
| 2 | `locus` | Gene symbol (`HLA-A`, `HLA-B`, `HLA-DRB1`, …) |
| 3 | `allele1` | First allele (e.g. `B*57:01`) |
| 4 | `allele2` | Second allele, or `NA` |
| 5 | `resolution` | Reported resolution (`one-field` / `two-field` / `three-field` / `four-field`) |
| 6 | `gl_string` | Locus genotype as a GL String (Milius 2013; Mack 2023), e.g. `HLA-A*01:01+HLA-A*29:02`. Only `+` (gene copies) and `^` (loci, JSON sample level) are emitted; HLAnte models neither phase (`~`) nor genotype ambiguity (`\|`) |
| 7 | `tool` | Upstream HLA typer (`arcashla` / `t1k` / `hlahd` / `optitype`) |
| 8 | `imgt_accession` | `HLA#####` IMGT accession per allele, pipe-joined |
| 9 | `hla_class` | `I` or `II` |
| 10 | `hla_serotype` | WHO serotype (`DR3`, `DQ8`, `B57`, …); `NA` when not covered |
| 11 | `protein_group` | G or P group when known |
| 12 | `imgt_match_category` | How the call resolved against IPD-IMGT/HLA: `exact` (one listed allele), `prefix_unique` (one longer name), `prefix_multiple` (several), `g_group` / `p_group`, or `unmatched` |
| 13 | `imgt_match_candidates` | Number of IPD-IMGT/HLA alleles the call maps to — `1` for an exact hit, higher when the name is a prefix of several listed alleles |
| 14 | `gwas_traits` | GWAS trait names, pipe-joined (both alleles aggregated) |
| 15 | `gwas_traits_allele1` | GWAS traits for allele1 only |
| 16 | `gwas_traits_allele2` | GWAS traits for allele2 only |
| 17 | `gwas_p_values` | p-values aligned with `gwas_traits`, 4-decimal |
| 18 | `gwas_odds_ratios` | OR or β aligned with `gwas_traits`, 4-decimal |
| 19 | `gwas_pmids` | PubMed IDs aligned with `gwas_traits` |
| 20 | `gwas_annotation_resolution` | Resolution at which GWAS hit was found |
| 21 | `gwas_annotation_scope` | `allele` / `subtype` / `locus` (breadth of match) |
| 22 | `gwas_matched_allele` | Catalogue key that actually matched, after any truncation |
| 23 | `gwas_match_broadened` | `yes` when a record was found only after truncating the submitted allele, i.e. the association is reported for a less specific name; `no` when it matched as submitted |
| 24 | `gwas_index_siblings` | Keys in the GWAS index sharing the matched prefix (not IPD-IMGT/HLA alleles) |
| 25 | `pharm_drugs` | Drug names from PharmGKB, pipe-joined (both alleles) |
| 26 | `pharm_drugs_allele1` | PharmGKB drugs for allele1 only |
| 27 | `pharm_drugs_allele2` | PharmGKB drugs for allele2 only |
| 28 | `pharm_evidence` | PharmGKB evidence level (`1A` / `1B` / `2A`) per drug |
| 29 | `pharm_cpic_action` | CPIC action verb per drug (e.g. `Contraindicated`) |
| 30 | `pharm_pmids` | PMIDs per PharmGKB annotation |
| 31 | `disease_risk_summary` | One-line human-readable disease association summary per allele |
| 32 | `drug_response_summary` | One-line human-readable drug response summary per allele |
| 33 | `clinical_significance` | Evidence-strength label (see below) |
| 34 | `significance_basis` | Evidence layer behind `clinical_significance`: `cpic_guideline`, `database_association`, or `curated_table` |
| 35 | `allele_frequency` | AFND population-group frequency, 6-decimal |
| 36 | `allele_freq_population` | AFND group code (`EUR` / `AFR` / `EAS` / `SAS` / `MID` / `AMR` / `OCE` / `global`; `ASN` = alias for `EAS`) |
| 37 | `input_quality_score` | Heuristic score in `[0, 1]`, 4-decimal |
| 38 | `input_quality_tier` | `detailed` (≥0.85) / `partial` (0.70–0.84) / `limited` (<0.70) |
| 39 | `input_quality_rationale` | Pipe-delimited reason codes explaining the score |
| 40 | `caller_allele_quality` | Per-allele quality reported by the typing tool itself, per allele slot. Supplied only by T1K's native layout; arcasHLA and HLA-HD report none, and OptiType reports a solution-level objective rather than an allele quality. Distinct from `input_quality_score`, which HLAnte computes |

### `clinical_significance` labels

These are **evidence-strength descriptors**, not ACMG/AMP classifications.

| Label | Meaning |
|-------|---------|
| `Actionable pharmacogenomic risk (CPIC 1A — avoid)` | CPIC 1A AVOID or SJS/TEN-level drug reaction in built-in curated table |
| `Strong pharmacogenomic risk association` | Severe drug reaction (CPIC 1A), not SJS/TEN level |
| `Suggestive risk factor` | Strong GWAS or PharmGKB (1A/1B) signal without an actionable curated entry |
| `Inconclusive evidence` | Ambiguous or IMGT-unknown allele with no confirming source |
| `No reported risk` | IMGT-known, unambiguous, at least one DB lookup returned no signal |
| `Not assessed — insufficient coverage` | IMGT-known but every query returned zero hits |
| `Not in IMGT` | Well-formed allele name absent from loaded IPD-IMGT/HLA release |

### Separator conventions

| Separator | Meaning |
|-----------|---------|
| `\|` (pipe) | Joins multiple values within one cell (hits, alleles) |
| `,` (comma) | Joins PMIDs within a single PharmGKB record (`pharm_pmids` only) |
| `;;` | Separates allele1 from allele2 in `input_quality_rationale` |
| `NA` | Explicit missing-value token; empty cells never appear |

---

## JSON output

Fully nested structure with complete provenance. Use it when you need exact
p-values (the TSV truncates to 4 decimal places) or programmatic access to
individual annotation records.

Top-level structure:

```json
{
  "meta": {
    "hlante_version": "0.2.0",
    "generated": "2026-05-04T12:00:00Z",
    "db_versions": {"imgt": "3.64.0", "gwas_cache_date": "2026-04-20"},
    "input_source": "typing_tool",
    "disclaimer": "RESEARCH USE ONLY …"
  },
  "records": [
    {
      "sample_id": "HG00096",
      "locus": "HLA-B",
      "allele1": { … },
      "allele2": { … }
    }
  ]
}
```

Each allele object contains `gwas_hits`, `pharm_annotations`, `disease_entries`,
`input_quality_score`, `input_quality_tier`, `clinical_significance`, and
`allele_frequency`.

---

## Markdown output

Per-sample report formatted for human review. Includes a summary table,
GWAS association section, PharmGKB drug-response section, and curated
disease-entry section for each locus.

The report always opens with the research-use disclaimer.
