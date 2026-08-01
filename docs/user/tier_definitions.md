# Clinical Significance Labels

HLAnte assigns each annotated allele an overall `clinical_significance` label.
These are **evidence-strength descriptors**, not ACMG/AMP clinical variant
classifications.

| Label | Meaning |
|-------|---------|
| `Actionable pharmacogenomic risk (CPIC Level 1A)` | A PharmGKB Level 1A record with a CPIC guideline link. Curated-table entries without such a record get `Actionable pharmacogenomic risk (curated reference)` instead. |
| `Strong pharmacogenomic risk association` | Allele is in the built-in curated table with `likely pathogenic` significance (e.g., CPIC Level 1A A\*31:01, B\*13:01) |
| `Suggestive risk factor` | GWAS Catalog hit at p ≤ 5×10⁻⁸ or PharmGKB evidence level 1A / 1B |
| `Inconclusive evidence` | GWAS or PharmGKB hits present but below the above thresholds |
| `No reported risk` | Allele is IMGT-known and unambiguous; at least one database was queried and returned no hit |
| `Not assessed — insufficient coverage` | Allele is IMGT-known but no database query returned any result (offline mode, novel locus, or DB not installed) |
| `Not in IMGT` | Allele is absent from the local IPD-IMGT/HLA database |

## Important caveats

- These labels are **research annotations**, not clinical diagnoses or treatment
  recommendations.
- `No reported risk` does not mean an allele is benign — only that no hit was
  found in the configured databases. Absence of evidence is not evidence of
  absence.
- `Suggestive risk factor` encompasses population-level GWAS findings; individual
  risk depends on additional genetic and environmental factors.
