# Exact counters in characteristic two

This directory contains two C++20 backends:

| Executable | Purpose |
|---|---|
| `q2_norm_counter` | Exact Gray-code traversal for the compact binary computations |
| `q2_walsh_counter` | Independent Boolean–Walsh implementation used for comparison tests |

Build both executables with

```bash
make -C fast_q2
```

`q2_norm_counter` receives the quadratic reduced-norm form constructed and
checked by SageMath. It traverses the bounded Hom space in binary-reflected
Gray order and groups all nonzero morphisms by their exact monic norm.

The counter supports two logically exact modes:

- exhaustive mode returns exact multiplicities, including certified zeros;
- witness mode may stop only after every requested irreducible norm has
  occurred, which certifies positivity for all requested targets.

If any requested target is absent, witness mode necessarily reaches the end
of the finite space.

The focused single-pair harness is `gray_counter.py`. For example:

```bash
./run_sage.sh fast_q2/gray_counter.py \
  --p 'T^6 + T + 1' \
  --source 0 \
  --target 0 \
  --max-degree 10 \
  --compare-python
```

The compact family drivers are
`drinfeld_complete/compact_q2.py`,
`drinfeld_complete/compact_q2_degree8.py`, and
`drinfeld_complete/e_only_q2.py`. Their archive formats are described in
[`COMPACT_ARCHIVE_SCHEMA.md`](COMPACT_ARCHIVE_SCHEMA.md) and
[`COMPACT_DEGREE8_SCHEMA.md`](COMPACT_DEGREE8_SCHEMA.md).

Validators in this directory reconstruct the expected prime-level universe,
check traversal sizes and spectral certificates, and recompute the reported
value of \(E(\mathfrak p)\).
