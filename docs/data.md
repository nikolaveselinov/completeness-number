# Data organization

## Canonical census

[`results/census/completeness_numbers.csv`](../results/census/completeness_numbers.csv)
is the authoritative per-characteristic table. It contains one row for each
of 1,352 pairs \((q,\mathfrak p)\) in 31 complete
\((q,\deg\mathfrak p)\)-families.

| Column | Meaning |
|---|---|
| `q` | Cardinality of the constant field |
| `p` | Monic irreducible characteristic |
| `deg_p` | Degree of \(\mathfrak p\) |
| `class_count` | Number of supersingular vertices |
| `U` | Proven cutoff above which every prime level is complete |
| `E` | Exact completeness number |
| `largest_incomplete_degree` | \(E-1\), left empty when no incomplete level occurs |
| `evidence_grade` | Type of validation used for the row |
| `source_record` | Published table from which the row is reconciled |

Two derived tables are stored beside it:

- [`family_summary.csv`](../results/census/family_summary.csv), with one row
  for each complete family; and
- [`e_distribution.csv`](../results/census/e_distribution.csv), with the
  overall distribution of \(E(\mathfrak p)\).

## Expanded level data

The directory [`results/detailed/`](../results/detailed) contains detailed
tables for 1,266 characteristics:

- `characteristics.csv` gives one row per characteristic;
- `degree_data.csv` aggregates complete and incomplete levels by degree; and
- `incomplete_levels.csv` records every incomplete prime level and its zero
  positions.

The binary degree-eight computation uses a compact certificate format and is
stored separately in
[`results/q2_degree8_compact/`](../results/q2_degree8_compact). Its tables
cover all 30 irreducible degree-eight characteristics.

The binary degree-nine computation records the exact value of
\(E(\mathfrak p)\) without materializing every lower-level Brandt matrix. Its
56-row table is
[`results/q2_degree9_e_only/characteristics_summary.csv`](../results/q2_degree9_e_only/characteristics_summary.csv).

Together these three sources contain exactly the 1,352 rows of the canonical
census.

## Independent certificate

[`certificates/q2_degree9/`](../certificates/q2_degree9) contains a
Sage-independent verification of the two exceptional binary degree-nine
characteristics. The retained JSON records prove positivity for all
degree-fourteen levels and exhibit the exact degree-thirteen obstructions.
`SHA256SUMS` binds the records and verification scripts.

## Consistency audit

Run

```bash
python3 scripts/audit_census.py
```

The audit checks uniqueness, family sizes, \(E\)-distributions, agreement
among all published tables, the certificate manifest, and the stated
degree-thirteen and degree-fourteen conclusions.
