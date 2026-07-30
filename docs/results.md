# Computed completeness numbers

The canonical census contains every monic irreducible characteristic in the
families below. In the fourth column, the number in parentheses is the
multiplicity of the corresponding value of \(E(\mathfrak p)\).

| \(q\) | \(\deg\mathfrak p\) | characteristics | \(E(\mathfrak p)\)-distribution | \(U(\mathfrak p)\) |
|---:|---:|---:|---|---:|
| 2 | 1 | 2 | 1 (2) | 1 |
| 2 | 2 | 1 | 1 (1) | 1 |
| 2 | 3 | 2 | 6 (2) | 8 |
| 2 | 4 | 3 | 3 (1), 5 (2) | 6 |
| 2 | 5 | 6 | 10 (6) | 12 |
| 2 | 6 | 9 | 8 (4), 9 (5) | 11 |
| 2 | 7 | 18 | 12 (16), 14 (2) | 16 |
| 2 | 8 | 30 | 11 (20), 12 (10) | 15 |
| 2 | 9 | 56 | 14 (2), 16 (54) | 20 |
| 3 | 1 | 3 | 1 (3) | 1 |
| 3 | 2 | 3 | 1 (3) | 1 |
| 3 | 3 | 8 | 4 (8) | 6 |
| 3 | 4 | 18 | 4 (9), 5 (9) | 6 |
| 3 | 5 | 48 | 8 (48) | 10 |
| 4 | 1 | 4 | 1 (4) | 1 |
| 4 | 2 | 6 | 1 (6) | 1 |
| 4 | 3 | 20 | 4 (8), 6 (12) | 6 |
| 4 | 4 | 60 | 4 (12), 5 (48) | 5 |
| 5 | 1 | 5 | 1 (5) | 1 |
| 5 | 2 | 10 | 1 (10) | 1 |
| 5 | 3 | 40 | 4 (40) | 6 |
| 5 | 4 | 150 | 4 (130), 5 (20) | 5 |
| 7 | 1 | 7 | 1 (7) | 1 |
| 7 | 2 | 21 | 1 (21) | 1 |
| 7 | 3 | 112 | 4 (112) | 5 |
| 8 | 1 | 8 | 1 (8) | 1 |
| 8 | 2 | 28 | 1 (28) | 1 |
| 8 | 3 | 168 | 4 (168) | 5 |
| 11 | 1 | 11 | 1 (11) | 1 |
| 11 | 2 | 55 | 1 (55) | 1 |
| 11 | 3 | 440 | 4 (440) | 5 |

Thus the released data comprise 1,352 characteristics in 31 complete
families.

## Dependence on the characteristic

The value \(E(\mathfrak p)\) is not determined by
\(\deg\mathfrak p\). For example, among the three irreducible quartics over
\(\mathbf F_2\), one has \(E=3\) and two have \(E=5\). Further splits occur
for binary degrees six through nine, for ternary degree four, and for several
families over \(\mathbf F_4\) and \(\mathbf F_5\).

The complete polynomial-by-polynomial list is
[`results/census/completeness_numbers.csv`](../results/census/completeness_numbers.csv).

## Failure of the degree threshold \(\deg\mathfrak l>\deg\mathfrak p\)

Let

\[
\mathfrak p=T^3+T+1,\qquad
\mathfrak l=T^5+T^3+T^2+T+1
\]

over \(\mathbf F_2\). Both polynomials are irreducible and
\(\deg\mathfrak l>\deg\mathfrak p\), but in the canonical vertex ordering
beginning with \(j=0\),

\[
B_{\mathfrak p}(\mathfrak l)=
\begin{pmatrix}
0&18&15\\
6&11&16\\
5&16&12
\end{pmatrix}.
\]

The zero diagonal entry shows that the graph is not complete. Each row sums
to \(2^5+1=33\), and the matrix satisfies weighted symmetry with weights
\((3,1,1)\).

This example is also non-monotone in the level degree: every degree-four
prime level is complete, the displayed degree-five level is incomplete, and
every degree-six and degree-seven prime level is complete.

## Exceptional binary degree-nine pair

Exactly two irreducible binary degree-nine characteristics have \(E=14\);
the other 54 have \(E=16\). The exceptional pair is

\[
\begin{aligned}
\mathfrak p_0&=T^9+T^7+T^2+T+1,\\
\mathfrak p_1&=T^9+T^8+T^7+T^6+T^5+T^4+T^3+T+1,
\end{aligned}
\]

and \(\mathfrak p_1(T)=\mathfrak p_0(T+1)\). For each characteristic, the
\(j=0\) self-entry vanishes at exactly 92 degree-thirteen prime levels, while
every degree-fourteen level is complete. The two sets of obstructions
correspond under \(T\mapsto T+1\).

The independent records for this computation are in
[`certificates/q2_degree9/`](../certificates/q2_degree9).
