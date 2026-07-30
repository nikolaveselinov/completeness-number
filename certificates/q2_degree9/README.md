# Independent certificate for the binary degree-nine family

This directory verifies the two exceptional characteristics in the
\(q=2\), \(\deg\mathfrak p=9\) family:

\[
\begin{aligned}
\mathfrak p_0&=T^9+T^7+T^2+T+1,\\
\mathfrak p_1&=T^9+T^8+T^7+T^6+T^5+T^4+T^3+T+1.
\end{aligned}
\]

They satisfy \(\mathfrak p_1(T)=\mathfrak p_0(T+1)\). The retained records
establish that

- every degree-fourteen prime level is positive for all 14,706 unordered
  vertex pairs;
- the \(j=0\) self-entry vanishes for exactly 92 degree-thirteen levels for
  each characteristic; and
- the two obstruction sets correspond under \(T\mapsto T+1\).

Consequently \(E(\mathfrak p_0)=E(\mathfrak p_1)=14\).

## Files

- `raw/d9_degree14_diagonal_full.json` and
  `raw/d9_degree14_offdiag_full.json` contain the complete degree-fourteen
  pair records;
- `raw/d9_degree13_full.json` and
  `raw/d9_partner_degree13_witness.json` contain the degree-thirteen
  obstructions;
- the two `raw/d3_*.json` files reproduce the known
  degree-four-complete, degree-five-incomplete example;
- `d9_full_scan.py` performs the exact scan; and
- `d9_j0_verify.py` independently verifies the \(j=0\) calculation.

`SHA256SUMS` binds every retained script and record.

## Verification

The repository-level audit checks the manifest and the mathematical
conclusions above:

```bash
python3 scripts/audit_census.py
```

To recompute the records, first build the exact characteristic-two counter:

```bash
make -C fast_q2 q2_norm_counter
```

Then, for example, run

```bash
python3 certificates/q2_degree9/d9_full_scan.py \
  --degree-p 9 \
  --p-code 647 \
  --level-degree 14 \
  --counter fast_q2/q2_norm_counter \
  --workers 32 \
  --diagonal-only \
  --output /tmp/q2-d9-diagonal.json
```

The verifier is independent of Sage and of the main computation package. It
uses exact binary polynomial arithmetic and invokes the exhaustive
Gray-code counter.
