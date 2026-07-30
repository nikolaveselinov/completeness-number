# Reproducibility

## Environment

The released computations use SageMath 10.9, Python 3.13, and C++20. Create
the recorded environment with

```bash
conda env create -f environment.yml
conda activate drinfeld-completeness
```

`run_sage.sh` invokes Sage's Python interpreter and sets
`OPENBLAS_NUM_THREADS=1` and `OMP_NUM_THREADS=1`. Set
`DRINFELD_SAGE_PYTHON` only if the Sage-enabled interpreter is installed
under a nonstandard name.

## Verification levels

The following audit uses only the Python standard library:

```bash
make audit
```

It reconciles the 1,352-row census with its three published source tables and
checks the independent binary degree-nine certificate.

The portable release checks require Python, GNU Make, and a C++20 compiler:

```bash
make check
```

They run the census audit, compile all Python sources, validate shell syntax,
build the native counters, and execute the Sage-free tests.

The complete regression suite requires SageMath:

```bash
make test
```

It checks supersingular constructions, Brandt-matrix identities, exact Gray
traversals, compact archive validators, and the known low-degree examples.
Some tests exhaust large finite spaces and are substantially slower than the
portable checks.

## One characteristic

The standard backend supports \(q=2,3,4,8\):

```bash
./run_sage.sh compute.py \
  --q 2 \
  --p 'T^3 + T + 1' \
  --full-matrices \
  --output results/local/p_T3_T_1.json
```

Without `--full-matrices`, sufficiently high-degree positive entries may be
certified by the exact pair-specific spectral inequality. Zero entries are
always obtained by exhaustive Hom-space enumeration.

Each output archive records the input, cutoff \(U(\mathfrak p)\), exact value
\(E(\mathfrak p)\), supersingular vertices, all prime levels below the
cutoff, zero positions, and construction checks.

## Complete families

A family for the standard backend can be computed with

```bash
./run_sage.sh compute_family.py \
  --q 2 \
  --degree 5 \
  --output-dir results/local/q2_degree5

./run_sage.sh summarize_results.py results/local/q2_degree5
```

`compute_family.py` accepts `--index` for deterministic array-style
parallelism.

For an odd prime \(q\), build the native counter and run

```bash
make -C fast_prime

./run_sage.sh compute_prime_compact_family.py \
  --q 5 \
  --degree 3 \
  --index 0 \
  --output-dir results/local/q5_degree3
```

For the compact binary backends:

```bash
make -C fast_q2

./run_sage.sh compute_q2_degree7_compact.py \
  --index 0 \
  --output-dir results/local/q2_degree7

./run_sage.sh compute_q2_degree8_compact.py \
  --index 0 \
  --output-dir results/local/q2_degree8
```

The corresponding validators in `fast_q2/` rebuild the level universe,
recheck each spectral inequality, and verify the recorded traversal sizes.

The descending computation used for binary degree nine is

```bash
./run_sage.sh compute_q2_e_only_family.py \
  --degree 9 \
  --index 0 \
  --output-dir results/local/q2_degree9
```

## Determinism and retained evidence

All mathematical decisions use exact finite-field and integer arithmetic.
Finite sets are traversed in deterministic order, and primary archives are
written atomically. No completeness decision depends on random sampling or a
floating-point threshold.

The full primary JSON archives are not stored in this repository because
they are large and reproducible from the commands above. The released tree
contains the canonical result tables, the detailed incomplete-level tables,
and the compact independent certificate described in
[`data.md`](data.md).
