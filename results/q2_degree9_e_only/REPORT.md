# Binary degree-nine family

All 56 monic irreducible degree-nine characteristics over
\(\mathbf F_2\) were computed:

| \(E(\mathfrak p)\) | characteristics |
|---:|---:|
| 14 | 2 |
| 16 | 54 |

The two characteristics with \(E(\mathfrak p)=14\) are

\[
T^9+T^7+T^2+T+1
\quad\text{and}\quad
T^9+T^8+T^7+T^6+T^5+T^4+T^3+T+1.
\]

They are exchanged by \(T\mapsto T+1\). For each one, the \(j=0\)
self-entry vanishes at exactly 92 degree-thirteen prime levels, while every
degree-fourteen level is complete. The two zero sets correspond under the
same translation.

`characteristics_summary.csv` gives the complete polynomial-by-polynomial
list. The Sage-independent verification of the exceptional pair is retained
in [`../../certificates/q2_degree9/`](../../certificates/q2_degree9).

Theorem 1.2 in the forthcoming paper proves completeness for every prime
level of degree at least \(U(\mathfrak p)=20\).
