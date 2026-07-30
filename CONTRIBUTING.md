# Contributing

Corrections and improvements to the implementation, tests, documentation, or
data validation are welcome.

Before submitting a change, run

```bash
make check
```

Changes that use SageMath or alter a mathematical backend should also pass

```bash
make test
```

If the complete suite is too costly, report the precise tests that were run.

Computational changes should preserve exact arithmetic, deterministic
iteration, and the certificate checks described in
[`docs/algorithm.md`](docs/algorithm.md). A corrected numerical result should
include a regression test and enough provenance to reproduce the
recomputation. Generated binaries, scheduler logs, checkpoints, and large
primary JSON archives should not be committed.
