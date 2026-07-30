# Canonical census

`completeness_numbers.csv` is the authoritative table of exact completeness
numbers in this repository. It contains 1,352 characteristics in 31 complete
\((q,\deg\mathfrak p)\)-families.

The companion files are:

- `family_summary.csv`, with one row per complete family; and
- `e_distribution.csv`, with the overall distribution of
  \(E(\mathfrak p)\).

The source tables and certificate structure are described in
[`../../docs/data.md`](../../docs/data.md). Run

```bash
python3 scripts/audit_census.py
```

from the repository root to reconcile the census with all published source
tables and the independent binary degree-nine certificate.
