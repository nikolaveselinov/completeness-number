# Compact exact-zero archive for q=2, degree(p)=8

Degree eight has:

- theorem cutoff `U = 15`;
- 85 supersingular vertices, all of weight one;
- 3,655 unordered vertex pairs;
- 2,537 prime levels below the cutoff;
- 9,272,735 unordered-pair/level decisions per characteristic.

The archive kind is `q2_degree8_compact_zero_v1`. It uses the same strict run
and per-level records documented in `COMPACT_ARCHIVE_SCHEMA.md`, with these
degree-eight specializations:

- every pair has an exhaustive run of dimension 17 and maximum degree 11,
  selecting all 411 levels of degrees 1 through 11;
- every pair has a degree-14 bounded run of dimension 23 selecting all 2,126
  levels of degrees 12 through 14;
- the high run either records `all_targets_witnessed`, or reaches all
  \(2^{23}-1\) nonzero vectors and records exact absent norms;
- `spectral_degrees` is empty because the exact pairwise inequality does not
  succeed below the global cutoff;
- weighted-symmetry transport is ordinary symmetry because all weights equal
  one.

`validate_compact_degree8.py` rebuilds the characteristic, context, level
universe, pair universe, run coverage, exact zero sets, degree summaries, and
`E(p)`. Its optional `--recompute` mode reruns the archived counter calls.
