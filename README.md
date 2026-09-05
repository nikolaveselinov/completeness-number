# Completeness numbers for supersingular Drinfeld isogeny graphs

This repository is the computational supplement to
[*On Supersingular Isogeny Graphs of Drinfeld Modules*](https://arxiv.org/abs/2608.29812).
It contains the SageMath implementation and data used to compute completeness numbers for supersingular
rank-two Drinfeld modules over $A=\mathbb F_q[T]$.

For a finite characteristic $\mathfrak p$, the completeness number $E(\mathfrak p)$ is the least integer such that the supersingular isogeny graph in characteristic $\mathfrak p$ is complete for every prime $\mathfrak q\ne\mathfrak p$ with $\deg\mathfrak q\ge E(\mathfrak p)$.

All decisions are made with exact finite-field and integer arithmetic. The
program checks every prime level below the proven upper bound
$U(\mathfrak p)$; the theorem in the paper covers all larger
degrees.

## Data

The canonical table is
[`results/census/completeness_numbers.csv`](results/census/completeness_numbers.csv).
It contains 1,352 characteristics in 31 complete
$(q,\deg\mathfrak p)$-families:

| $q$ | degrees of $\mathfrak p$ | characteristics |
|---:|---:|---:|
| 2 | 1–9 | 127 |
| 3 | 1–5 | 80 |
| 4 | 1–4 | 90 |
| 5 | 1–4 | 205 |
| 7 | 1–3 | 140 |
| 8 | 1–3 | 204 |
| 11 | 1–3 | 506 |

The principal numerical observations are summarized in
[`docs/results.md`](docs/results.md). The
[`data guide`](docs/data.md) describes the canonical census, expanded
per-level tables, and the independent certificate for the binary
degree-nine family.

## Requirements

- SageMath 10.9;
- Python 3.13;
- a C++20 compiler; and
- GNU Make.

The recorded Conda environment can be created with

```bash
conda env create -f environment.yml
conda activate drinfeld-completeness
```

## Computing one characteristic

For example:

```bash
./run_sage.sh compute.py \
  --q 2 \
  --p 'T^3 + T + 1' \
  --full-matrices \
  --output results/local/p_T3_T_1.json
```

The output records the supersingular vertices, the cutoff
$U(\mathfrak p)$, the exact value $E(\mathfrak p)$, every prime level
below the cutoff, and the evidence for each Brandt entry. Omitting
`--full-matrices` permits exact pair-specific spectral positivity
certificates at the higher levels.

Further commands for complete families and compact backends are given in
[`docs/reproducibility.md`](docs/reproducibility.md).

## Verification

The census has a Sage-free consistency audit:

```bash
make audit
```

The portable release checks additionally compile the Python sources and
native counters and run the Sage-free tests:

```bash
make check
```

The full mathematical regression suite requires the Conda environment:

```bash
make test
```

The algorithms and the checks attached to each type of certificate are
described in [`docs/algorithm.md`](docs/algorithm.md).

## Repository structure

| Path | Contents |
|---|---|
| [`drinfeld_complete/`](drinfeld_complete) | SageMath implementation and compact computation backends |
| [`fast_q2/`](fast_q2) | Exact C++ counters for characteristic two |
| [`fast_prime/`](fast_prime) | Exact C++ counter for odd prime fields |
| [`results/`](results) | Canonical and expanded result tables |
| [`certificates/`](certificates) | Independently checkable retained certificates |
| [`tests/`](tests) | Regression and validation tests |
| [`docs/`](docs) | Mathematical and reproducibility documentation |

Large primary JSON archives and scheduler logs are not included.

## Citation and license

Citation metadata are provided in [`CITATION.cff`](CITATION.cff). The source
code, documentation, and included data are released under the
[MIT License](LICENSE).
