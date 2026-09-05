# Exact implementation of Appendix A

This note documents how the code implements Algorithm 1 in Appendix A of the
forthcoming paper *On Supersingular Isogeny Graphs of Drinfeld Modules*.
To avoid the paper's overloaded notation, \(q\) below is the size of the
constant field and \(\ell\) is a finite prime of \(A=\mathbf F_q[T]\).

## What is computed

For every monic irreducible \(\ell\ne p\), let
\(B_p(\ell)=(b_{ij}(\ell))\) be the Brandt matrix. By Definition 4.3,

\[
E(p)=1+\max\{\deg\ell:B_p(\ell)\text{ has a zero entry}\},
\]

with the maximum of the empty set defined as zero.  Theorem 5.5 identifies
the same zero test with the modular-polynomial remainder test:

\[
\overline{\Phi}_\ell(X,Y)
 \bmod (S_p(X),S_p(Y))=0
\quad\Longleftrightarrow\quad
b_{ij}(\ell)>0\ \text{for all }i,j.
\]

The implementation therefore follows Appendix A's finite scan exactly, but
computes the equivalent Brandt zero test directly instead of first
constructing the enormous universal polynomial \(\Phi_\ell\).

## 1. Proven finite cutoff

For \(d=\deg p\ge3\), the code computes

\[
C_p=
\begin{cases}
(q^d-q)/(q-1),&d\text{ odd},\\
(q^d-q^2)/(q^2-1),&d\text{ even},
\end{cases}
\]

and

\[
U(p)=\left\lfloor
2\log_q\!\left(C_p+\sqrt{C_p^2-1}\right)
\right\rfloor+1.
\]

No floating-point arithmetic is used.  The code finds the first integer
\(e\) satisfying

\[
(q^e+1)^2>4C_p^2q^e.
\]

Theorem 1.2 proves that every prime level of degree at least \(U(p)\) is
complete.  Thus the program enumerates the exact finite set

\[
\{\ell:\ell\text{ monic irreducible},\ \ell\ne p,\ 1\le\deg\ell<U(p)\}.
\]

For \(d=1,2\), the theorem gives \(E(p)=U(p)=1\).

## 2. Supersingular polynomial and vertices

Let \(\theta\) be the image of \(T\) in characteristic \(p\).  For the
generic normalized module

\[
\phi_T=\theta+\tau+x\tau^2,
\]

the code extracts the coefficient of \(\tau^d\) in \(\phi_p\), obtaining
Gekeler's Hasse/Deuring polynomial \(H_p(x)\).  If
\(m=\deg H_p\), it constructs

\[
S_p(J)=
\begin{cases}
J^mH_p(1/J)/H_p(0),&d\text{ even},\\
J^{m+1}H_p(1/J)/H_p(0),&d\text{ odd}.
\end{cases}
\]

An independent recurrence is evaluated as a construction certificate:

\[
P_0=P_1=1,\qquad
P_n=
J^{(q^{n-1}-(-1)^{n-1})/(q+1)}P_{n-1}
-(\theta^{q^{n-1}}-\theta)P_{n-2}.
\]

The expected result is \(P_d\) in even degree and \(JP_d\) in odd degree.
The code requires this recurrence to equal the Hasse construction.

The polynomial is also required to be squarefree, to split over
\(\mathbf F_{q^{2d}}\), and to have the theoretical class number

\[
h(p)=
\begin{cases}
(q^d-1)/(q^2-1),&d\text{ even},\\
(q^d-q)/(q^2-1)+1,&d\text{ odd}.
\end{cases}
\]

Every root \(j_i\) is used to construct a normalized module

\[
(\phi_i)_T=\theta+g_i\tau+\tau^2,\qquad g_i^{q+1}=j_i.
\]

For every vertex the program independently checks Sage's supersingularity
test and the identity

\[
(\phi_i)_p=\tau^{2d}.
\]

The latter identity also proves that all coefficients of every morphism
between these modules lie in \(\mathbf F_{q^{2d}}\).

## 3. Exact Brandt entries

For each ordered pair \((i,j)\) and required degree bound \(m\), Sage
constructs a complete \(\mathbf F_q\)-basis of

\[
\{u\in\mathrm{Hom}(\phi_i,\phi_j):\deg_\tau u\le m\}.
\]

Every vector in this finite space is exhausted.  The reduced norm is the
determinant of the corresponding motive matrix.  After monic normalization,
the vectors are grouped by their exact norm in \(\mathbf F_q[T]\).  Definition
4.1 then gives

\[
b_{ij}(\ell)=
\frac{\#\{u\in\mathrm{Hom}(\phi_i,\phi_j):N(u)=(\ell)\}}
     {|\mathrm{Aut}(\phi_j)|}.
\]

The implementation checks that every numerator is divisible by the target
automorphism order.  It also checks the expected bounded-Hom dimension,
reconstructs sample morphisms, and compares the packed quadratic norm with
Sage's direct norm.

This is a proof of absence as well as presence: a zero is reported only
after the entire relevant Hom space has been exhausted.

## 4. Exact finite-field Gray traversals

The reduced norm is a quadratic form on each bounded Hom space.  Computing it
from scratch for every vector would be unnecessarily expensive.

- For \(q=2\), the code interpolates diagonal and cross terms, packs finite
  field coefficients into bit blocks, and walks \(\mathbf F_2^m\) in binary
  reflected Gray order.
- For \(q=3\), it packs ternary coefficients into disjoint 1- and 2-bitplanes
  and walks \(\mathbf F_3^m\) in ternary reflected Gray order, updating the
  norm and all gradients exactly.
- For \(q=4\) and \(q=8\), it represents the splitting field absolutely over
  \(\mathbf F_2\), builds exact multiplication tables for every embedded
  constant-field scalar, and walks the Hom space in reflected radix-\(q\)
  order. The update
  \(Q(x+d e_k)-Q(x)=d^2Q(e_k)+dB(x,e_k)\) is evaluated with packed binary
  polynomial arithmetic.
- For odd prime fields, the compact backend sends the Sage-interpolated
  \(\mathbf F_q[T]\)-valued quadratic form to a compiled reflected
  \(q\)-ary counter. It maintains the norm and every gradient coefficient
  modulo \(q\), normalizes each nonzero norm to its monic polynomial code,
  and records exact target multiplicities.

If the bounded Hom space has dimension \(r\) over \(\mathbf F_q\), a full
traversal visits precisely \(q^r-1\) nonzero morphisms.  There is no random
sampling and no floating-point decision.

## 5. Pair-specific spectral certificates

High-degree positive entries may be skipped without weakening the result.
Set

\[
M=\frac{q^d-1}{q^2-1},\qquad
w_i=\frac{|\mathrm{Aut}(\phi_i)|}{q-1},\qquad
C_{ij}^2=(w_iM-1)(w_jM-1).
\]

Proposition 4.11 certifies \(b_{ij}(\ell)>0\) whenever

\[
(q^{\deg\ell}+1)^2>4C_{ij}^2q^{\deg\ell}.
\]

The code evaluates this strict inequality with exact rational/integer
arithmetic.  JSON entries marked `spectral_positive` are therefore proven
positive; entries marked `exhaustive_hom` were obtained by full Hom-space
enumeration.

## 6. End-to-end consistency checks

Whenever a complete Brandt matrix is materialized, the code requires

\[
\sum_j b_{ij}(\ell)=q^{\deg\ell}+1
\]

for every row and

\[
w_jb_{ij}(\ell)=w_ib_{ji}(\ell)
\]

for every ordered pair.  The archive auditor separately rebuilds the
supersingular context and exact irreducible-level set, recomputes all
spectral inequalities, and directly rechecks every stored full-matrix row
sum and weighted-symmetry identity.

The standard full-matrix backends support \(q=2,3,4,8\). The compact
odd-prime backend supports any odd prime field whose interpolated norm form
passes the required base-field and direct-Sage checks; the computed census
uses \(q=3,5,7,11\). Degrees one and two use Theorem 1.2 directly over any
field for which the supersingular context can be constructed.

## 7. Compact exact-zero mode in characteristic two

For \(\deg p=7\) and \(8\), compact backends compute the zero/positive
predicate required by Appendix A without materializing every positive Brandt
entry at every level. They enumerate only \(i\leq j\). The exact identity

\[
w_jb_{ij}(\ell)=w_ib_{ji}(\ell)
\]

and positivity of the weights prove
\(b_{ij}(\ell)=0\Longleftrightarrow b_{ji}(\ell)=0\).

Every non-spectral decision is either `exhaustive`, meaning all \(2^r-1\)
nonzero vectors in the complete bounded Hom space were visited, or
`all_targets_witnessed`, meaning every selected prime norm occurred before
the traversal ended. Partial occurrence counts are presence witnesses only,
never exact Brandt entries. Remaining high degrees are marked
`spectral_positive` only after the exact pairwise inequality succeeds.

Each packed norm form is checked on deterministic multi-coordinate vectors
against Sage's direct norm, exercising diagonal and cross terms. The compact
auditor independently rebuilds the context and level universe, recomputes
spectral inequalities, verifies traversal sizes and target coverage,
reconstructs every ordered zero set, and recomputes \(E(p)\). Optional
`--recompute` mode reruns every archived counter invocation.

At degree eight there are 85 vertices, 3,655 unordered pairs, 2,537 levels,
and no sub-cutoff pair-specific spectral certificate. Every pair therefore
has a complete traversal through degree 11 and a degree-14 all-target run for
degrees 12--14. The latter either witnesses every requested level or reaches
all \(2^{23}-1\) vectors and proves every remaining absence.

## 8. Compact exact-zero mode over odd prime fields

For each unordered pair, the odd-prime backend finds every degree below
\(U(p)\) not already covered by the exact pair-specific spectral inequality.
It builds the complete bounded Hom space at their largest degree and exhausts
all \(q^r-1\) vectors with the compiled counter. Exact multiplicities are
retained for all selected prime norms; zero is archived only when its count is
zero after that full traversal.

Before every compiled run, deterministic vectors containing both single and
multiple nonzero coordinates are evaluated twice: once in the interpolated
quadratic form and once through Sage's direct reduced norm. The compiled
counter itself is regression-tested against the independent ternary backend
on every target norm in a full Hom space. The strict archive validator then
rebuilds the context and level universe, recomputes spectral degrees,
verifies \(q^r-1\) traversal sizes and automorphism divisibility, reconstructs
ordered zero positions by weighted symmetry, and derives \(E(p)\) anew.
