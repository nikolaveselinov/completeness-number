# Binary degree-eight family

All 30 monic irreducible degree-eight characteristics over
\(\mathbf F_2\) were computed. The resulting distribution is

| \(E(\mathfrak p)\) | characteristics |
|---:|---:|
| 11 | 20 |
| 12 | 10 |

For each characteristic, the compact computation reconstructs the exact
2,537-level universe below \(U(\mathfrak p)=15\). Its unordered-pair
certificates decide every entry either by exhaustive Hom-space traversal or
by a successful all-target witness traversal. The theorem in the
forthcoming paper proves completeness at every prime level of degree at
least 15.

The files in this directory are:

- `characteristics.csv`: exact value for each characteristic;
- `degree_data.csv`: complete and incomplete levels grouped by degree; and
- `incomplete_levels.csv`: every incomplete level and its zero positions.

The archive format and validator are documented in
[`../../fast_q2/COMPACT_DEGREE8_SCHEMA.md`](../../fast_q2/COMPACT_DEGREE8_SCHEMA.md).
