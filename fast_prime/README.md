# Exact compiled odd-prime counter

`prime_norm_counter.cpp` exhausts a bounded Hom-space quadratic norm form
over an odd prime field in reflected \(p\)-ary Gray order. It visits exactly
\(p^r-1\) nonzero vectors for a dimension-\(r\) space and returns exact
multiplicities for the requested monic prime norms.

The Python compact backend interpolates the complete quadratic form with
Sage, requires every coefficient to lie in the declared prime field, checks
deterministic multi-coordinate vectors against Sage's direct reduced norm,
and only then invokes this counter. A zero count is archived only after full
exhaustion.

Build with:

```bash
make -C fast_prime
```
