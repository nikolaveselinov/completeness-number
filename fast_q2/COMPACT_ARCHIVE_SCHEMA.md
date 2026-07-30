# Compact exact-zero archive for \(q=2,\deg p=7\)

The degree-seven computation has \(43\) supersingular vertices, \(946\)
unordered vertex pairs, and \(4719\) prime levels below the exact cutoff
\(U(p)=16\).  Storing every Brandt entry would therefore obscure the only
information needed for \(E(p)\): which entries are zero.

`validate_compact_archive.py` validates the following JSON schema.  It
rebuilds the supersingular context and the level universe rather than trusting
their archived values.

## Top-level object

```json
{
  "archive_kind": "q2_compact_zero_v1",
  "schema_version": 1,
  "q": 2,
  "p": "T^7 + T + 1",
  "degree_p": 7,
  "theorem_bound_U": 16,
  "E": 12,
  "class_count": 43,
  "weights": [3, 1],
  "j_invariants": ["0", "..."],
  "pair_certificates": [],
  "total_vectors_visited": 123456,
  "levels": [],
  "degree_summary": {}
}
```

The illustrative `weights` array above is abbreviated; the real array has
one entry per vertex.  Vertex indices are meaningful only relative to the
archived `j_invariants`.  Both arrays must exactly equal the canonical
context rebuilt by the validator.

There must be exactly one `pair_certificates` record for each
\(0\le i\le j<43\):

```json
{
  "pair": [0, 17],
  "spectral_degrees": [15],
  "runs": []
}
```

Only the canonical unordered orientation `i <= j` is stored.  For an exact
raw Hom count \(r\),

\[
  b_{ij}=r/w_j,\qquad b_{ji}=r/w_i,
\]

and hence \(w_jb_{ij}=w_ib_{ji}\).  A zero or positivity certificate therefore
transports to the reversed ordered pair.  An early-stop count is only a lower
bound, so it transports positivity but must never be archived as an exact
Brandt entry.

Every degree in `spectral_degrees` is independently checked using the exact
pair-specific inequality

\[
  (2^e+1)^2 >
  4(w_iM-1)(w_jM-1)2^e,\qquad M=(2^7-1)/(2^2-1).
\]

All levels in such a degree receive certificate `spectral_positive`.  Run
records cover only the remaining degrees.  The validator permits a producer
to use fewer spectral degrees than are available, but never permits a claimed
spectral degree that fails the strict inequality.

## Run records

A run selects either whole level degrees:

```json
{
  "target_degrees": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
  "max_degree": 10,
  "dimension": 16,
  "iterations": 65535,
  "target_count": 225,
  "seen": 221,
  "invalid_norms": 0,
  "exhaustive": true,
  "outcome": "exhaustive",
  "zero_codes": [7, 11, 13, 19],
  "direct_sage_norm_cross_checks": 8,
  "seconds": 0.125
}
```

or an explicit, canonically sorted list:

```json
{
  "target_codes": [2053, 2065],
  "...": "the same remaining fields"
}
```

Exactly one of `target_degrees` and `target_codes` is allowed.  Polynomial
code bit \(k\) is the coefficient of \(T^k\).  Targets must be precisely
members of the expected level universe, so the characteristic \(p\) is never
a target.

For `outcome: "exhaustive"` the validator requires

```text
exhaustive = true
iterations = 2^dimension - 1
seen = target_count - len(zero_codes)
invalid_norms = 0
```

Every code in `zero_codes` then has certificate `exhaustive_zero`; every
other selected code has certificate `exhaustive_positive`.

For a successful witness pass:

```json
{
  "target_degrees": [11, 12, 13, 14],
  "max_degree": 14,
  "dimension": 24,
  "iterations": 3141592,
  "target_count": 2312,
  "seen": 2312,
  "invalid_norms": 0,
  "exhaustive": false,
  "outcome": "all_targets_witnessed",
  "zero_codes": [],
  "direct_sage_norm_cross_checks": 8,
  "seconds": 0.5
}
```

Every selected code then has certificate `all_targets_witnessed`.  A
non-exhaustive pass that has not witnessed every target proves nothing about
the absent targets and is rejected.  If a pass reaches the end of the Hom
space, it must instead be recorded as `outcome: "exhaustive"`.

For each unordered pair, run selectors must be mutually disjoint, must not
overlap a spectral degree, and together with the spectral degrees must cover
the complete expected level universe through degree \(15\).  Each run also
records a positive number of direct Sage norm cross-checks and a finite,
nonnegative elapsed time.  Thus there are exactly

```text
946 * 4719 = 4,464,174
```

pair-level decisions, even though they are represented by only a few thousand
run records.

## Per-level exact zero records

There must be exactly one level record for every monic irreducible
\(\ell\ne p\) with \(1\le\deg\ell<16\):

```json
{
  "ell": "T^3 + T + 1",
  "code": 11,
  "degree": 3,
  "complete": false,
  "zero_entries": [[0, 4], [4, 0]]
}
```

`zero_entries` is the sorted set of *ordered* pairs.  The validator derives
it independently from the unordered pair certificates, including both
orientations of every off-diagonal zero, and requires exact equality.  It
then recomputes `complete`, every entry of `degree_summary`, and

\[
E(p)=1+\max\{\deg\ell:\ell\text{ has a zero entry}\}.
\]

The exact irreducible counts in degrees \(1,\ldots,15\) are

```text
2, 1, 2, 3, 6, 9, 17, 30, 56, 99, 186, 335, 630, 1161, 2182.
```

The degree-seven value is \(17\), not \(18\), because the characteristic
itself is excluded.

## Structural validation versus exact replay

The default validator proves that the archive is internally complete and
that every claimed counter result has the necessary arithmetic shape.  Like
the existing matrix archive auditor, it necessarily treats archived counter
summaries as computational evidence.

Pass `--recompute` for an independent exact replay.  The validator rebuilds
each bounded Hom space and quadratic norm form, reruns the compiled counter
with the archived target selector, and checks the complete result.  This is
the strongest audit mode, but it has essentially the cost of the original
computation.

Neither a log hash nor `seen == target_count` by itself is a self-contained
mathematical witness.  If replay-free proof objects are wanted later, the
producer must store a coefficient vector for every witnessed positive target;
zero targets still require an exhaustive replay (or a separate theoretical
certificate).

## Producer pitfalls

- A bounded Hom basis at a high degree cannot be truncated by selecting some
  of its basis vectors.  The lower bounded-Hom space is the kernel of the
  high Ore-coefficient maps and is generally spanned by linear combinations
  of the high basis.
- If a high-degree quadratic form is reused, first compute a kernel basis
  \(r_a\), then restrict by
  \(D'_a=Q(r_a)\) and
  \(C'_{ab}=Q(r_a+r_b)+Q(r_a)+Q(r_b)\).  Every lane above the lower bound must
  vanish before packed values are trimmed.
- Partial early-stop counts are lower bounds.  They need not be divisible by
  automorphism orders and cannot be used in row-sum checks.
- The Brandt matrix is not ordinarily symmetric at the exceptional
  \(j=0\) vertex.  Only weighted symmetry is valid; zero patterns are
  symmetric because all weights are positive.
- A structural audit cannot turn an unverified counter summary into a proof.
  Use `--recompute` (or retain independently checkable witness vectors) for
  an audit that does not trust the producer's run summaries.
