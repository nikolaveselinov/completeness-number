"""Core exact algorithms for the completeness number E(p).

This module implements the finite scan in Appendix A of the accompanying
paper.  The expensive universal modular-polynomial primitive is replaced by
an equivalent characteristic-p Brandt computation:

* construct the supersingular polynomial exactly;
* enumerate every bounded morphism in each relevant Hom space;
* group morphisms by their exact ideal norm;
* read off all Brandt entries and hence the same zero/nonzero decision as the
  modular-polynomial remainder criterion.

The q=2, q=3, q=4, and q=8 paths evaluate the raw reduced norm as a quadratic
form and traverse the Hom space in binary, ternary, or extension-binary
q-ary Gray-code order.  This is exact, not probabilistic.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from itertools import product
from time import perf_counter
from typing import Any

from sage.all import (
    DrinfeldModule,
    GF,
    Hom,
    PolynomialRing,
)


@dataclass
class SupersingularContext:
    q: int
    A: Any
    T: Any
    p: Any
    degree: int
    residue_field: Any
    splitting_field: Any
    theta: Any
    hasse_polynomial: Any
    supersingular_polynomial: Any
    j_invariants: list[Any]
    modules: list[Any]
    automorphism_orders: list[int]
    checks: dict[str, Any]


def monic_irreducibles(A, degree: int) -> list[Any]:
    """Return all monic irreducibles of one exact positive degree."""
    if degree < 1:
        return []
    return [
        f
        for f in A.polynomials(of_degree=degree)
        if f.is_monic() and f.is_irreducible()
    ]


def theorem_C(q: int, degree: int) -> int:
    """Return C_p from Theorem 1.2."""
    if degree < 3:
        return 0
    if degree % 2:
        return (q**degree - q) // (q - 1)
    return (q**degree - q**2) // (q**2 - 1)


def theorem_bound(q: int, degree: int) -> int:
    """Compute U(p) exactly, without floating-point logarithms.

    This is floor(2 log_q(C + sqrt(C^2 - 1))) + 1.  Equivalently it is
    the first e for which (q^e + 1)^2 > 4 C^2 q^e.
    """
    if degree <= 2:
        return 1
    C = theorem_C(q, degree)
    e = 1
    while True:
        Q = q**e
        if (Q + 1) ** 2 > 4 * C**2 * Q:
            return e
        e += 1


def _expected_class_number(q: int, degree: int) -> int:
    if degree % 2 == 0:
        return (q**degree - 1) // (q**2 - 1)
    return (q**degree - q) // (q**2 - 1) + 1


def _supersingular_recurrence(coefficient_field, theta, degree: int, q: int):
    """Return the independent Gekeler recurrence for S_p(J)."""
    RJ = PolynomialRing(coefficient_field, "J")
    J = RJ.gen()
    P0 = RJ.one()
    P1 = RJ.one()
    for n in range(2, degree + 1):
        exponent = (q ** (n - 1) - (-1) ** (n - 1)) // (q + 1)
        bracket = theta ** (q ** (n - 1)) - theta
        P0, P1 = P1, J**exponent * P1 - bracket * P0
    return J * P1 if degree % 2 else P1


def build_supersingular_context(q: int, p) -> SupersingularContext:
    """Construct S_p, every supersingular j, and canonical modules.

    ``p`` must be a monic irreducible in GF(q)[T].  All modules are
    normalized as phi_T = theta + g*tau + tau^2.  Supersingularity then
    gives phi_p = tau^(2 deg p), so every Hom coefficient lies in
    F_{q^(2 deg p)}.
    """
    q = int(q)
    A = p.parent()
    Fq = A.base_ring()
    T = A.gen()
    p = A(p)
    if int(Fq.cardinality()) != q:
        raise ValueError("q does not match the coefficient field of p")
    if not p.is_monic() or not p.is_irreducible():
        raise ValueError("p must be monic and irreducible")
    degree = int(p.degree())

    if int(Fq.degree()) == 1:
        # Retain Sage's compact two-stage representation over prime fields.
        residue_field = Fq.extension(p, names="theta")
        theta0 = residue_field.gen()
        splitting_field = residue_field.extension(2, names="z")
        theta = splitting_field(theta0)
    else:
        # Polynomial-quotient extensions do not implement all finite-field
        # operations needed below.  Construct F_{q^(2d)} directly over F_q;
        # this also gives Sage the canonical embedding of the function-ring
        # constants, then choose an image of T by finding a root of p.
        splitting_field = Fq.extension(2 * degree, names="z")
        p_over_splitting_field = p.change_ring(splitting_field)
        characteristic_roots = p_over_splitting_field.roots(
            multiplicities=False
        )
        characteristic_roots.sort(key=repr)
        if len(characteristic_roots) != degree:
            raise ArithmeticError(
                f"expected {degree} roots of p in F_(q^(2d)), "
                f"found {len(characteristic_roots)}"
            )
        theta = characteristic_roots[0]
        theta0 = theta
        # H_p and S_p may harmlessly be constructed over the splitting field.
        residue_field = splitting_field

    Rx = PolynomialRing(residue_field, "x")
    x = Rx.gen()
    generic = DrinfeldModule(A, [theta0, 1, x])
    hasse = Rx(generic(p)[degree])
    if hasse[0] == 0:
        raise ArithmeticError("H_p(0) unexpectedly vanished")

    # J^m H_p(1/J) / H_p(0), with the j=0 factor in odd degree.
    RJ = PolynomialRing(residue_field, "J")
    J = RJ.gen()
    reverse = RJ(list(reversed(hasse.list()))) / hasse[0]
    supersingular = J * reverse if degree % 2 else reverse
    supersingular = RJ(supersingular)
    recurrence = _supersingular_recurrence(
        residue_field, theta0, degree, q
    )

    SK = supersingular.change_ring(splitting_field)
    roots = SK.roots(multiplicities=False)
    roots.sort(key=lambda a: (0 if not a else 1, repr(a)))
    expected_n = _expected_class_number(q, degree)
    if len(roots) != expected_n:
        raise ArithmeticError(
            f"expected {expected_n} supersingular roots, found {len(roots)}"
        )

    modules = []
    for j in roots:
        g = splitting_field(0) if not j else j.nth_root(q + 1)
        module = DrinfeldModule(A, [theta, g, 1])
        modules.append(module)

    tau = modules[0].ore_variable()
    pure_frobenius = tau ** (2 * degree)
    phi_p_checks = [module(p) == pure_frobenius for module in modules]
    supersingular_checks = [bool(module.is_supersingular()) for module in modules]

    # Avoid constructing J^(q^(2d)) explicitly.
    split_check = (
        (pow(J, q ** (2 * degree), supersingular) - J) % supersingular
    ).is_zero()
    expected_hasse_degree = (
        (q**degree - 1) // (q**2 - 1)
        if degree % 2 == 0
        else (q**degree - q) // (q**2 - 1)
    )
    checks = {
        "hasse_degree": int(hasse.degree()),
        "expected_hasse_degree": int(expected_hasse_degree),
        "supersingular_degree": int(supersingular.degree()),
        "expected_class_number": int(expected_n),
        "squarefree": bool(supersingular.gcd(supersingular.derivative()).degree() == 0),
        "splits_over_q_2d": bool(split_check),
        "recurrence_matches_hasse": bool(recurrence == supersingular),
        "all_phi_p_pure_frobenius": bool(all(phi_p_checks)),
        "all_sage_supersingular": bool(all(supersingular_checks)),
    }
    if not all(
        [
            checks["hasse_degree"] == checks["expected_hasse_degree"],
            checks["supersingular_degree"] == checks["expected_class_number"],
            checks["squarefree"],
            checks["splits_over_q_2d"],
            checks["recurrence_matches_hasse"],
            checks["all_phi_p_pure_frobenius"],
            checks["all_sage_supersingular"],
        ]
    ):
        raise ArithmeticError(f"supersingular construction failed checks: {checks}")

    automorphism_orders = [q**2 - 1 if not j else q - 1 for j in roots]
    return SupersingularContext(
        q=q,
        A=A,
        T=T,
        p=p,
        degree=degree,
        residue_field=residue_field,
        splitting_field=splitting_field,
        theta=theta,
        hasse_polynomial=hasse,
        supersingular_polynomial=supersingular,
        j_invariants=roots,
        modules=modules,
        automorphism_orders=automorphism_orders,
        checks=checks,
    )


def _pair_C_squared(ctx: SupersingularContext, i: int, j: int) -> Fraction:
    q = ctx.q
    d = ctx.degree
    M = Fraction(q**d - 1, q**2 - 1)
    wi = Fraction(ctx.automorphism_orders[i], q - 1)
    wj = Fraction(ctx.automorphism_orders[j], q - 1)
    return (wi * M - 1) * (wj * M - 1)


def pair_is_spectrally_positive(
    ctx: SupersingularContext, i: int, j: int, level_degree: int
) -> bool:
    """Certify b_ij(ell)>0 for every ell of this degree."""
    Q = ctx.q**level_degree
    C2 = _pair_C_squared(ctx, i, j)
    return Fraction((Q + 1) ** 2, 1) > 4 * C2 * Q


def _raw_norm(morphism):
    """Return the coherent (not monic-normalized) quadratic norm."""
    return morphism._motive_matrix().det()


def _raw_coeff_int(coefficient) -> int:
    if hasattr(coefficient, "backend"):
        coefficient = coefficient.backend()
    return int(coefficient)


def _pack_raw_q2(raw, max_degree: int, field_degree: int) -> int:
    """Pack a K[T] polynomial as concatenated polynomial-basis bit strings."""
    packed = 0
    for k, coefficient in enumerate(raw.list()):
        if k > max_degree:
            raise ArithmeticError("raw norm exceeded the requested degree")
        packed |= _raw_coeff_int(coefficient) << (k * field_degree)
    return packed


def _quadratic_form_q2(homset, basis, max_degree: int, field_degree: int):
    """Interpolate the raw reduced norm quadratic form in characteristic 2."""
    m = len(basis)
    diagonal = [
        _pack_raw_q2(_raw_norm(b), max_degree, field_degree) for b in basis
    ]
    cross = [[0] * m for _ in range(m)]
    for i in range(m):
        for j in range(i + 1, m):
            value = _raw_norm(basis[i] + basis[j])
            cross[i][j] = (
                _pack_raw_q2(value, max_degree, field_degree)
                ^ diagonal[i]
                ^ diagonal[j]
            )
            cross[j][i] = cross[i][j]
    return diagonal, cross


def _normalize_packed_q2(
    packed: int, max_degree: int, field_degree: int
) -> tuple[int, int] | None:
    """Return (degree, F2-polynomial bit code), or None for zero."""
    field_mask = (1 << field_degree) - 1
    degree = max_degree
    while degree >= 0:
        coefficient = (packed >> (degree * field_degree)) & field_mask
        if coefficient:
            break
        degree -= 1
    if degree < 0:
        return None
    leading = coefficient
    code = 0
    for k in range(degree + 1):
        coefficient = (packed >> (k * field_degree)) & field_mask
        if coefficient == leading:
            code |= 1 << k
        elif coefficient != 0:
            raise ArithmeticError("monic raw norm has a coefficient outside F_2")
    if not (code & (1 << degree)):
        raise ArithmeticError("normalization lost the leading coefficient")
    return degree, code


def enumerate_pair_norms_q2(
    ctx: SupersingularContext,
    source: int,
    target: int,
    max_degree: int,
    *,
    cross_check_samples: int = 8,
) -> tuple[dict[int, int], dict[str, Any]]:
    """Exhaustively count all monic norms for one ordered pair when q=2.

    Dictionary keys are bit encodings of F2[T] polynomials.
    """
    if ctx.q != 2:
        raise ValueError("the Gray-code backend is specific to q=2")
    homset = Hom(ctx.modules[source], ctx.modules[target])
    basis = homset.basis(degree=max_degree)
    dimension = len(basis)
    expected_dimension = (
        2 * (max_degree + 1) - (ctx.degree - 1)
        if max_degree >= ctx.degree - 2
        else None
    )
    if expected_dimension is not None and dimension != expected_dimension:
        raise ArithmeticError(
            f"bounded Hom dimension {dimension}, expected {expected_dimension}"
        )
    field_degree = 2 * ctx.degree
    diagonal, cross = _quadratic_form_q2(
        homset, basis, max_degree, field_degree
    )

    counts: dict[int, int] = {}
    witnesses: dict[int, int] = {}
    active = 0
    packed = 0
    exact_degree_counts = [0] * (max_degree + 1)
    valid = 0
    started = perf_counter()

    # Binary-reflected Gray code: transition n-1 -> n flips ctz(n).
    for n in range(1, 1 << dimension):
        flip = (n & -n).bit_length() - 1
        delta = diagonal[flip]
        others = active & ~(1 << flip)
        while others:
            low = others & -others
            j = low.bit_length() - 1
            delta ^= cross[flip][j]
            others ^= low
        packed ^= delta
        active ^= 1 << flip
        normalized = _normalize_packed_q2(packed, max_degree, field_degree)
        if normalized is None:
            raise ArithmeticError("nonzero morphism had zero raw norm")
        degree, code = normalized
        exact_degree_counts[degree] += 1
        counts[code] = counts.get(code, 0) + 1
        witnesses.setdefault(code, active)
        valid += 1

    if valid != (1 << dimension) - 1:
        raise ArithmeticError("Gray traversal did not cover the Hom space")

    # Reconstruct a few witnesses and compare against Sage's direct norm.
    checked = 0
    for code, vector_code in witnesses.items():
        if checked >= cross_check_samples:
            break
        ore = ctx.modules[source].ore_polring().zero()
        for k, b in enumerate(basis):
            if vector_code & (1 << k):
                ore += b.ore_polynomial()
        direct = homset(ore).norm().gen()
        direct_code = sum(int(direct[k]) << k for k in range(direct.degree() + 1))
        if direct_code != code:
            raise ArithmeticError(
                f"quadratic norm {code:b} disagrees with Sage {direct_code:b}"
            )
        checked += 1

    metadata = {
        "source": int(source),
        "target": int(target),
        "max_degree": int(max_degree),
        "dimension": int(dimension),
        "vectors_enumerated": int((1 << dimension) - 1),
        "exact_degree_vector_counts": [int(x) for x in exact_degree_counts],
        "quadratic_cross_checks": int(checked),
        "seconds": perf_counter() - started,
    }
    return counts, metadata


def _coefficient_backend(coefficient):
    if hasattr(coefficient, "backend"):
        return coefficient.backend()
    return coefficient


def _binary_polynomial_code(element) -> int:
    """Return polynomial-basis coordinates of a characteristic-two element."""
    element = _coefficient_backend(element)
    if not element:
        return 0
    polynomial = element.polynomial()
    return sum(int(value) << index for index, value in enumerate(polynomial.list()))


def _gf2_multiply_codes(
    left: int,
    right: int,
    modulus_code: int,
    field_degree: int,
) -> int:
    """Multiply two polynomial-basis codes in an absolute binary field."""
    product_code = 0
    multiplier = int(right)
    multiplicand = int(left)
    while multiplier:
        if multiplier & 1:
            product_code ^= multiplicand
        multiplier >>= 1
        multiplicand <<= 1
    for exponent in range(product_code.bit_length() - 1, field_degree - 1, -1):
        if product_code & (1 << exponent):
            product_code ^= modulus_code << (exponent - field_degree)
    return product_code


def _gf2_scalar_table(
    scalar_code: int,
    modulus_code: int,
    field_degree: int,
) -> list[int]:
    """Tabulate multiplication by one scalar in an absolute binary field."""
    basis_products = [
        _gf2_multiply_codes(
            1 << exponent, scalar_code, modulus_code, field_degree
        )
        for exponent in range(field_degree)
    ]
    table = [0] * (1 << field_degree)
    for value in range(1, len(table)):
        low = value & -value
        table[value] = (
            table[value ^ low] ^ basis_products[low.bit_length() - 1]
        )
    return table


def _multiply_packed_binary_scalar(
    packed: int,
    multiplication_table: list[int],
    max_degree: int,
    field_degree: int,
) -> int:
    """Multiply every coefficient block of a packed polynomial by a scalar."""
    field_mask = (1 << field_degree) - 1
    result = 0
    for degree in range(max_degree + 1):
        shift = degree * field_degree
        coefficient = (packed >> shift) & field_mask
        result |= multiplication_table[coefficient] << shift
    return result


def _pack_raw_binary_extension(raw, max_degree: int, field_degree: int) -> int:
    """Pack a polynomial over an absolute binary extension field."""
    packed = 0
    for degree, coefficient in enumerate(raw.list()):
        if degree > max_degree:
            raise ArithmeticError("raw norm exceeded the requested degree")
        code = _binary_polynomial_code(coefficient)
        if code >= 1 << field_degree:
            raise ArithmeticError("finite-field coefficient exceeded its basis")
        packed |= code << (degree * field_degree)
    return packed


def _raw_is_prime_field_q3(raw) -> bool:
    """Return whether every coefficient of ``raw`` is fixed by x |-> x^3."""
    return all(
        _coefficient_backend(coefficient) ** 3
        == _coefficient_backend(coefficient)
        for coefficient in raw.list()
    )


def _pack_raw_q3(
    raw, max_degree: int, field_degree: int
) -> tuple[int, int]:
    """Pack a K[T] polynomial into the 1- and 2-trit bitplanes.

    If the interpolated quadratic form is already defined over F_3, the
    effective ``field_degree`` is one and bit k represents the coefficient
    of T^k.  The general representation concatenates polynomial-basis blocks
    of length ``field_degree`` and is retained as an exact fallback.
    """
    ones = 0
    twos = 0
    for k, coefficient in enumerate(raw.list()):
        if k > max_degree:
            raise ArithmeticError("raw norm exceeded the requested degree")
        coefficient = _coefficient_backend(coefficient)
        if field_degree == 1:
            if coefficient**3 != coefficient:
                raise ArithmeticError(
                    "raw norm coefficient unexpectedly left F_3"
                )
            coordinates = [int(coefficient)]
        elif hasattr(coefficient, "polynomial"):
            coordinates = [
                int(value) for value in coefficient.polynomial().list()
            ]
        else:
            coordinates = [int(coefficient)]
        if len(coordinates) > field_degree:
            raise ArithmeticError("finite-field coefficient exceeded its basis")
        offset = k * field_degree
        for index, value in enumerate(coordinates):
            if value == 1:
                ones |= 1 << (offset + index)
            elif value == 2:
                twos |= 1 << (offset + index)
            elif value != 0:
                raise ArithmeticError("coefficient is not represented over F_3")
    return ones, twos


def _quadratic_form_q3(basis, max_degree: int, absolute_field_degree: int):
    """Interpolate and pack the reduced norm quadratic form for q=3."""
    dimension = len(basis)
    diagonal_raw = [_raw_norm(basis_element) for basis_element in basis]
    cross_raw = [[None] * dimension for _ in range(dimension)]
    raw_coefficients = list(diagonal_raw)
    for i in range(dimension):
        for j in range(i + 1, dimension):
            value = (
                _raw_norm(basis[i] + basis[j])
                - diagonal_raw[i]
                - diagonal_raw[j]
            )
            cross_raw[i][j] = value
            cross_raw[j][i] = value
            raw_coefficients.append(value)

    # Sage's Hom bases in the tested q=3 families yield a quadratic form over
    # F_3 itself.  Detect this rather than assume it, retaining a full
    # F_{3^(2d)} polynomial-basis representation as a correctness fallback.
    field_degree = (
        1
        if all(_raw_is_prime_field_q3(value) for value in raw_coefficients)
        else absolute_field_degree
    )
    diagonal = [
        _pack_raw_q3(value, max_degree, field_degree)
        for value in diagonal_raw
    ]
    zero = (0, 0)
    cross = [[zero] * dimension for _ in range(dimension)]
    for i in range(dimension):
        for j in range(i + 1, dimension):
            packed = _pack_raw_q3(
                cross_raw[i][j], max_degree, field_degree
            )
            cross[i][j] = packed
            cross[j][i] = packed
    return diagonal, cross, field_degree


def _add_packed_q3(
    left: tuple[int, int],
    right: tuple[int, int],
    lane_mask: int,
) -> tuple[int, int]:
    """Add two packed ternary vectors without carries between trits."""
    left_one, left_two = left
    right_one, right_two = right
    left_zero = lane_mask ^ (left_one | left_two)
    right_zero = lane_mask ^ (right_one | right_two)
    ones = (
        (left_zero & right_one)
        | (left_one & right_zero)
        | (left_two & right_two)
    )
    twos = (
        (left_zero & right_two)
        | (left_two & right_zero)
        | (left_one & right_one)
    )
    return ones, twos


def _negate_packed_q3(value: tuple[int, int]) -> tuple[int, int]:
    """Negate a packed ternary vector by interchanging 1 and 2."""
    return value[1], value[0]


def _normalize_packed_q3(
    packed: tuple[int, int], max_degree: int, field_degree: int
) -> tuple[int, tuple[int, int]] | None:
    """Return (degree, monic F3[T] bitplanes), or None for zero."""
    ones, twos = packed
    occupied = ones | twos
    if occupied == 0:
        return None
    if ones & twos:
        raise ArithmeticError("invalid packed ternary coefficient")

    if field_degree == 1:
        degree = occupied.bit_length() - 1
        if degree > max_degree:
            raise ArithmeticError("raw norm exceeded the requested degree")
        if (twos >> degree) & 1:
            ones, twos = twos, ones
        return degree, (ones, twos)

    field_mask = (1 << field_degree) - 1
    degree = (occupied.bit_length() - 1) // field_degree
    if degree > max_degree:
        raise ArithmeticError("raw norm exceeded the requested degree")
    shift = degree * field_degree
    leading = (
        (ones >> shift) & field_mask,
        (twos >> shift) & field_mask,
    )
    negative_leading = _negate_packed_q3(leading)
    normalized_ones = 0
    normalized_twos = 0
    for k in range(degree + 1):
        shift = k * field_degree
        coefficient = (
            (ones >> shift) & field_mask,
            (twos >> shift) & field_mask,
        )
        if coefficient == leading:
            normalized_ones |= 1 << k
        elif coefficient == negative_leading:
            normalized_twos |= 1 << k
        elif coefficient != (0, 0):
            raise ArithmeticError(
                "monic raw norm has a coefficient outside F_3"
            )
    if not ((normalized_ones >> degree) & 1):
        raise ArithmeticError("normalization lost the leading coefficient")
    return degree, (normalized_ones, normalized_twos)


def enumerate_pair_norms_q3(
    ctx: SupersingularContext,
    source: int,
    target: int,
    max_degree: int,
    *,
    cross_check_samples: int = 8,
) -> tuple[dict[tuple[int, int], int], dict[str, Any]]:
    """Exhaustively count all monic norms for one ordered pair when q=3."""
    if ctx.q != 3:
        raise ValueError("the ternary Gray-code backend is specific to q=3")
    homset = Hom(ctx.modules[source], ctx.modules[target])
    basis = homset.basis(degree=max_degree)
    dimension = len(basis)
    expected_dimension = (
        2 * (max_degree + 1) - (ctx.degree - 1)
        if max_degree >= ctx.degree - 2
        else None
    )
    if expected_dimension is not None and dimension != expected_dimension:
        raise ArithmeticError(
            f"bounded Hom dimension {dimension}, expected {expected_dimension}"
        )

    diagonal, cross, field_degree = _quadratic_form_q3(
        basis, max_degree, 2 * ctx.degree
    )
    lane_mask = (
        1 << ((max_degree + 1) * field_degree)
    ) - 1
    counts: dict[tuple[int, int], int] = {}
    witnesses: dict[tuple[int, int], tuple[int, ...]] = {}
    coefficients = [0] * dimension
    directions = [1] * dimension
    gradients = [(0, 0)] * dimension
    packed = (0, 0)
    exact_degree_counts = [0] * (max_degree + 1)
    vector_count = 3**dimension - 1
    started = perf_counter()

    # Reflected ternary Gray code.  The norm update is
    # Q(x+d e_k)-Q(x) = d grad_k(Q)(x) + d^2 Q(e_k), where d=+/-1.
    for _ in range(vector_count):
        flip = 0
        while coefficients[flip] + directions[flip] not in (0, 1, 2):
            directions[flip] = -directions[flip]
            flip += 1
        direction = directions[flip]
        signed_gradient = (
            gradients[flip]
            if direction == 1
            else _negate_packed_q3(gradients[flip])
        )
        packed = _add_packed_q3(
            packed,
            _add_packed_q3(
                signed_gradient, diagonal[flip], lane_mask
            ),
            lane_mask,
        )
        for index in range(dimension):
            if index == flip:
                continue
            signed_cross = (
                cross[index][flip]
                if direction == 1
                else _negate_packed_q3(cross[index][flip])
            )
            gradients[index] = _add_packed_q3(
                gradients[index], signed_cross, lane_mask
            )
        # 2*d is -1 for d=1 and +1 for d=-1 in F_3.
        gradient_delta = (
            _negate_packed_q3(diagonal[flip])
            if direction == 1
            else diagonal[flip]
        )
        gradients[flip] = _add_packed_q3(
            gradients[flip], gradient_delta, lane_mask
        )
        coefficients[flip] += direction

        normalized = _normalize_packed_q3(
            packed, max_degree, field_degree
        )
        if normalized is None:
            raise ArithmeticError("nonzero morphism had zero raw norm")
        degree, key = normalized
        exact_degree_counts[degree] += 1
        counts[key] = counts.get(key, 0) + 1
        if key not in witnesses:
            witnesses[key] = tuple(coefficients)

    if sum(counts.values()) != vector_count:
        raise ArithmeticError("ternary Gray traversal missed Hom vectors")

    checked = 0
    for key, vector in witnesses.items():
        if checked >= cross_check_samples:
            break
        ore = ctx.modules[source].ore_polring().zero()
        for coefficient, basis_element in zip(vector, basis):
            if coefficient:
                ore += coefficient * basis_element.ore_polynomial()
        direct = homset(ore).norm().gen()
        direct_key = _poly_key_q3(direct)
        if direct_key != key:
            raise ArithmeticError(
                f"quadratic norm {key} disagrees with Sage {direct_key}"
            )
        checked += 1

    metadata = {
        "source": int(source),
        "target": int(target),
        "max_degree": int(max_degree),
        "dimension": int(dimension),
        "vectors_enumerated": int(vector_count),
        "exact_degree_vector_counts": [int(x) for x in exact_degree_counts],
        "quadratic_cross_checks": int(checked),
        "packing_field_degree": int(field_degree),
        "seconds": perf_counter() - started,
    }
    return counts, metadata


def _quadratic_form_binary_constant(
    basis,
    max_degree: int,
    field_degree: int,
) -> tuple[list[int], list[list[int]]]:
    """Interpolate a packed norm quadratic form in characteristic two."""
    dimension = len(basis)
    diagonal_raw = [_raw_norm(basis_element) for basis_element in basis]
    diagonal = [
        _pack_raw_binary_extension(value, max_degree, field_degree)
        for value in diagonal_raw
    ]
    cross = [[0] * dimension for _ in range(dimension)]
    for i in range(dimension):
        for j in range(i + 1, dimension):
            value = (
                _raw_norm(basis[i] + basis[j])
                - diagonal_raw[i]
                - diagonal_raw[j]
            )
            packed = _pack_raw_binary_extension(
                value, max_degree, field_degree
            )
            cross[i][j] = packed
            cross[j][i] = packed
    return diagonal, cross


def _normalize_packed_binary_constant(
    packed: int,
    max_degree: int,
    field_degree: int,
    constant_degree: int,
    constant_raw_codes: list[int],
    multiplication_tables: dict[int, list[int]],
) -> tuple[int, int] | None:
    """Return (degree, packed F_q coefficient code), or None for zero."""
    field_mask = (1 << field_degree) - 1
    degree = max_degree
    while degree >= 0:
        leading = (packed >> (degree * field_degree)) & field_mask
        if leading:
            break
        degree -= 1
    if degree < 0:
        return None

    normalized_lookup = {0: 0}
    for constant_code, raw_code in enumerate(constant_raw_codes[1:], start=1):
        product = multiplication_tables[raw_code][leading]
        normalized_lookup[product] = constant_code
    q = 1 << constant_degree
    if len(normalized_lookup) != q:
        raise ArithmeticError("embedded constant-field values are not distinct")

    key = 0
    for exponent in range(degree + 1):
        coefficient = (
            packed >> (exponent * field_degree)
        ) & field_mask
        try:
            normalized = normalized_lookup[coefficient]
        except KeyError as error:
            raise ArithmeticError(
                "monic raw norm has a coefficient outside the constant field"
            ) from error
        key |= normalized << (constant_degree * exponent)
    if ((key >> (constant_degree * degree)) & (q - 1)) != 1:
        raise ArithmeticError("normalization lost the leading coefficient")
    return degree, key


def _poly_key_binary_constant(poly, q: int) -> int:
    constant_degree = q.bit_length() - 1
    if q != 1 << constant_degree:
        raise ValueError("constant-field size must be a power of two")
    key = 0
    for exponent in range(int(poly.degree()) + 1):
        coefficient_code = _binary_polynomial_code(poly[exponent])
        if coefficient_code >= q:
            raise ArithmeticError(
                "polynomial coefficient is not in the constant field"
            )
        key |= coefficient_code << (constant_degree * exponent)
    return key


def _poly_key_q4(poly) -> int:
    return _poly_key_binary_constant(poly, 4)


def _poly_key_q8(poly) -> int:
    return _poly_key_binary_constant(poly, 8)


def _binary_constant_backend(ctx: SupersingularContext) -> dict[str, Any]:
    """Build and cache exact absolute-binary scalar arithmetic for one context."""
    cached = getattr(ctx, "_binary_constant_backend_cache", None)
    if cached is not None:
        return cached

    constant_field = ctx.A.base_ring()
    constant_degree = int(constant_field.degree())
    q = 1 << constant_degree
    if int(constant_field.characteristic()) != 2 or ctx.q != q:
        raise ArithmeticError("constant field is not an absolute binary field")

    splitting_field = ctx.splitting_field
    if int(splitting_field.base_ring().cardinality()) != 2:
        raise ArithmeticError("backend requires an absolute binary splitting field")
    field_degree = int(splitting_field.degree())
    if int(splitting_field.cardinality()) != 1 << field_degree:
        raise ArithmeticError("unexpected splitting-field representation")
    modulus = splitting_field.modulus()
    modulus_code = sum(
        int(value) << index for index, value in enumerate(modulus.list())
    )
    if modulus_code.bit_length() - 1 != field_degree:
        raise ArithmeticError("splitting-field modulus has the wrong degree")

    generator = constant_field.gen()
    constant_values = []
    for code in range(q):
        value = constant_field.zero()
        for exponent in range(constant_degree):
            if code & (1 << exponent):
                value += generator**exponent
        constant_values.append(value)
    if len(set(constant_values)) != q:
        raise ArithmeticError("failed to enumerate the constant field")
    if [
        _binary_polynomial_code(value) for value in constant_values
    ] != list(range(q)):
        raise ArithmeticError(
            "unexpected constant-field polynomial-basis coordinates"
        )
    constant_raw_codes = [
        _binary_polynomial_code(splitting_field(value))
        for value in constant_values
    ]
    if len(set(constant_raw_codes)) != q:
        raise ArithmeticError("constant-field embedding is not injective")
    multiplication_tables = {
        scalar_code: _gf2_scalar_table(
            scalar_code, modulus_code, field_degree
        )
        for scalar_code in constant_raw_codes
    }
    cached = {
        "q": q,
        "constant_degree": constant_degree,
        "field_degree": field_degree,
        "modulus_code": modulus_code,
        "constant_values": constant_values,
        "constant_raw_codes": constant_raw_codes,
        "multiplication_tables": multiplication_tables,
    }
    setattr(ctx, "_binary_constant_backend_cache", cached)
    return cached


def _enumerate_pair_norms_binary_constant(
    ctx: SupersingularContext,
    source: int,
    target: int,
    max_degree: int,
    *,
    cross_check_samples: int = 8,
) -> tuple[dict[int, int], dict[str, Any]]:
    """Exhaustively count monic norms over a binary constant field."""
    backend = _binary_constant_backend(ctx)
    q = int(backend["q"])
    constant_degree = int(backend["constant_degree"])
    field_degree = int(backend["field_degree"])
    modulus_code = int(backend["modulus_code"])
    constant_values = backend["constant_values"]
    constant_raw_codes = backend["constant_raw_codes"]
    multiplication_tables = backend["multiplication_tables"]

    homset = Hom(ctx.modules[source], ctx.modules[target])
    basis = homset.basis(degree=max_degree)
    dimension = len(basis)
    expected_dimension = (
        2 * (max_degree + 1) - (ctx.degree - 1)
        if max_degree >= ctx.degree - 2
        else None
    )
    if expected_dimension is not None and dimension != expected_dimension:
        raise ArithmeticError(
            f"bounded Hom dimension {dimension}, expected {expected_dimension}"
        )

    diagonal, cross = _quadratic_form_binary_constant(
        basis, max_degree, field_degree
    )
    possible_deltas = {
        constant_raw_codes[left] ^ constant_raw_codes[right]
        for left in range(q)
        for right in range(q)
        if left != right
    }
    if not possible_deltas <= set(multiplication_tables):
        raise ArithmeticError("constant-field differences left the embedding")
    scaled_diagonal: dict[int, list[int]] = {}
    scaled_cross: dict[int, list[list[int]]] = {}
    for delta in possible_deltas:
        delta_squared = multiplication_tables[delta][delta]
        delta_squared_table = multiplication_tables.get(delta_squared)
        if delta_squared_table is None:
            delta_squared_table = _gf2_scalar_table(
                delta_squared, modulus_code, field_degree
            )
            multiplication_tables[delta_squared] = delta_squared_table
        scaled_diagonal[delta] = [
            _multiply_packed_binary_scalar(
                value,
                delta_squared_table,
                max_degree,
                field_degree,
            )
            for value in diagonal
        ]
        scaled_cross[delta] = [
            [
                _multiply_packed_binary_scalar(
                    value,
                    multiplication_tables[delta],
                    max_degree,
                    field_degree,
                )
                for value in row
            ]
            for row in cross
        ]

    counts: dict[int, int] = {}
    witnesses: dict[int, tuple[int, ...]] = {}
    positions = [0] * dimension
    directions = [1] * dimension
    gradients = [0] * dimension
    packed = 0
    exact_degree_counts = [0] * (max_degree + 1)
    vector_count = q**dimension - 1
    started = perf_counter()

    # Reflected base-q Gray code.  The update in characteristic two is
    # Q(x+d e_k)-Q(x) = d^2 Q(e_k) + d*B(x,e_k).
    for _ in range(vector_count):
        flip = 0
        while positions[flip] + directions[flip] not in range(q):
            directions[flip] = -directions[flip]
            flip += 1
        old_position = positions[flip]
        positions[flip] += directions[flip]
        new_position = positions[flip]
        delta = (
            constant_raw_codes[old_position]
            ^ constant_raw_codes[new_position]
        )
        packed ^= scaled_diagonal[delta][flip]
        packed ^= _multiply_packed_binary_scalar(
            gradients[flip],
            multiplication_tables[delta],
            max_degree,
            field_degree,
        )
        for index in range(dimension):
            if index != flip:
                gradients[index] ^= scaled_cross[delta][index][flip]

        normalized = _normalize_packed_binary_constant(
            packed,
            max_degree,
            field_degree,
            constant_degree,
            constant_raw_codes,
            multiplication_tables,
        )
        if normalized is None:
            raise ArithmeticError("nonzero morphism had zero raw norm")
        degree, key = normalized
        exact_degree_counts[degree] += 1
        counts[key] = counts.get(key, 0) + 1
        witnesses.setdefault(key, tuple(positions))

    if sum(counts.values()) != vector_count:
        raise ArithmeticError("base-q Gray traversal missed Hom vectors")

    checked = 0
    for key, vector in witnesses.items():
        if checked >= cross_check_samples:
            break
        ore = ctx.modules[source].ore_polring().zero()
        for position, basis_element in zip(vector, basis):
            if position:
                ore += (
                    constant_values[position]
                    * basis_element.ore_polynomial()
                )
        direct = homset(ore).norm().gen()
        direct_key = _poly_key_binary_constant(direct, q)
        if direct_key != key:
            raise ArithmeticError(
                f"quadratic norm {key} disagrees with Sage {direct_key}"
            )
        checked += 1

    metadata = {
        "source": int(source),
        "target": int(target),
        "max_degree": int(max_degree),
        "dimension": int(dimension),
        "vectors_enumerated": int(vector_count),
        "exact_degree_vector_counts": [int(x) for x in exact_degree_counts],
        "quadratic_cross_checks": int(checked),
        "gray_radix": int(q),
        "packing_absolute_binary_field_degree": int(field_degree),
        "seconds": perf_counter() - started,
    }
    return counts, metadata


def enumerate_pair_norms_q4(
    ctx: SupersingularContext,
    source: int,
    target: int,
    max_degree: int,
    *,
    cross_check_samples: int = 8,
) -> tuple[dict[int, int], dict[str, Any]]:
    """Exhaustively count all monic norms for one ordered pair when q=4."""
    if ctx.q != 4:
        raise ValueError("this exhaustive backend is specific to q=4")
    return _enumerate_pair_norms_binary_constant(
        ctx,
        source,
        target,
        max_degree,
        cross_check_samples=cross_check_samples,
    )


def enumerate_pair_norms_q8(
    ctx: SupersingularContext,
    source: int,
    target: int,
    max_degree: int,
    *,
    cross_check_samples: int = 8,
) -> tuple[dict[int, int], dict[str, Any]]:
    """Exhaustively count all monic norms for one ordered pair when q=8."""
    if ctx.q != 8:
        raise ValueError("this exhaustive backend is specific to q=8")
    return _enumerate_pair_norms_binary_constant(
        ctx,
        source,
        target,
        max_degree,
        cross_check_samples=cross_check_samples,
    )


def _poly_code_q2(poly) -> int:
    return sum(int(poly[k]) << k for k in range(poly.degree() + 1))


def _poly_key_q3(poly) -> tuple[int, int]:
    ones = 0
    twos = 0
    for k in range(poly.degree() + 1):
        coefficient = int(poly[k])
        if coefficient == 1:
            ones |= 1 << k
        elif coefficient == 2:
            twos |= 1 << k
        elif coefficient != 0:
            raise ArithmeticError("polynomial coefficient is not in F_3")
    return ones, twos


def _poly_from_code_q2(A, code: int):
    return A([(code >> k) & 1 for k in range(code.bit_length())])


def _serialize_matrix(matrix: list[list[int | None]]) -> list[list[int | None]]:
    return [
        [None if value is None else int(value) for value in row] for row in matrix
    ]


def _verify_full_brandt_matrix(
    ctx: SupersingularContext, ell, matrix: list[list[int | None]]
) -> dict[str, bool]:
    if any(value is None for row in matrix for value in row):
        return {"full": False, "row_sums": False, "weighted_symmetry": False}
    expected = ctx.q ** int(ell.degree()) + 1
    row_sums = all(sum(int(v) for v in row) == expected for row in matrix)
    symmetry = True
    n = len(matrix)
    weights = [
        order // (ctx.q - 1) for order in ctx.automorphism_orders
    ]
    for i in range(n):
        for j in range(n):
            if weights[j] * int(matrix[i][j]) != weights[i] * int(matrix[j][i]):
                symmetry = False
    return {
        "full": True,
        "row_sums": bool(row_sums),
        "weighted_symmetry": bool(symmetry),
    }


def compute_completeness(
    ctx: SupersingularContext,
    *,
    full_matrices: bool = False,
) -> dict[str, Any]:
    """Compute E(p) exactly.

    Degrees one and two use the theorem directly over any constant field for
    which the supersingular context can be constructed.  Higher degrees use
    the exhaustive q=2, q=3, q=4, or q=8 norm backend.

    If ``full_matrices`` is false, pair-specific spectral positivity
    certificates skip high-degree Hom enumerations.  Such entries are stored
    as ``None`` in a Brandt matrix, but their positivity is proved and zero
    positions remain conclusive.
    """
    q = ctx.q
    d = ctx.degree
    U = theorem_bound(q, d)
    if d <= 2:
        # U=1, so no level or Hom-space data are needed.  Delaying backend
        # selection makes the theorem-only E=1 archive available generally.
        enumerate_pair_norms = None
        polynomial_key = None
    elif ctx.q == 2:
        enumerate_pair_norms = enumerate_pair_norms_q2
        polynomial_key = _poly_code_q2
    elif ctx.q == 3:
        enumerate_pair_norms = enumerate_pair_norms_q3
        polynomial_key = _poly_key_q3
    elif ctx.q == 4:
        enumerate_pair_norms = enumerate_pair_norms_q4
        polynomial_key = _poly_key_q4
    elif ctx.q == 8:
        enumerate_pair_norms = enumerate_pair_norms_q8
        polynomial_key = _poly_key_q8
    else:
        raise NotImplementedError(
            "the current optimized exhaustive backend supports q=2, q=3, "
            "q=4, and q=8; "
            "the supersingular construction itself is general"
        )
    n = len(ctx.modules)
    weights = [order // (q - 1) for order in ctx.automorphism_orders]

    levels_by_degree: dict[int, list[Any]] = {}
    all_levels: list[Any] = []
    for degree in range(1, U):
        levels = [
            ell
            for ell in monic_irreducibles(ctx.A, degree)
            if ell != ctx.p
        ]
        levels_by_degree[degree] = levels
        all_levels.extend(levels)

    matrices: dict[str, list[list[int | None]]] = {
        str(ell): [[None] * n for _ in range(n)] for ell in all_levels
    }
    evidence: dict[str, list[list[str]]] = {
        str(ell): [[""] * n for _ in range(n)] for ell in all_levels
    }
    enumeration_metadata = []

    if U > 1:
        for i in range(n):
            for j in range(n):
                if full_matrices:
                    max_degree = U - 1
                else:
                    uncertified = [
                        degree
                        for degree in range(1, U)
                        if not pair_is_spectrally_positive(ctx, i, j, degree)
                    ]
                    max_degree = max(uncertified, default=0)

                if max_degree:
                    counts, metadata = enumerate_pair_norms(
                        ctx, i, j, max_degree
                    )
                    enumeration_metadata.append(metadata)
                else:
                    counts = {}

                aut = ctx.automorphism_orders[j]
                for ell in all_levels:
                    key = str(ell)
                    degree = int(ell.degree())
                    if degree <= max_degree:
                        raw = counts.get(polynomial_key(ell), 0)
                        if raw % aut:
                            raise ArithmeticError(
                                f"raw count {raw} not divisible by Aut target {aut}"
                            )
                        matrices[key][i][j] = raw // aut
                        evidence[key][i][j] = "exhaustive_hom"
                    else:
                        if not pair_is_spectrally_positive(ctx, i, j, degree):
                            raise ArithmeticError("an uncomputed entry lacks a certificate")
                        matrices[key][i][j] = None
                        evidence[key][i][j] = "spectral_positive"

    level_records = []
    largest_bad_degree = 0
    degree_summary = {}
    for degree, levels in levels_by_degree.items():
        complete_count = 0
        incomplete = []
        for ell in levels:
            key = str(ell)
            matrix = matrices[key]
            zero_entries = [
                [i, j]
                for i in range(n)
                for j in range(n)
                if matrix[i][j] == 0
            ]
            complete = len(zero_entries) == 0
            if complete:
                complete_count += 1
            else:
                largest_bad_degree = max(largest_bad_degree, degree)
                incomplete.append(key)
            checks = _verify_full_brandt_matrix(ctx, ell, matrix)
            if checks["full"] and not (
                checks["row_sums"] and checks["weighted_symmetry"]
            ):
                raise ArithmeticError(f"Brandt checks failed for {ell}: {checks}")
            level_records.append(
                {
                    "ell": key,
                    "degree": int(degree),
                    "complete": bool(complete),
                    "zero_entries": zero_entries,
                    "brandt_matrix": _serialize_matrix(matrix),
                    "entry_evidence": evidence[key],
                    "checks": checks,
                }
            )
        degree_summary[str(degree)] = {
            "level_count": len(levels),
            "complete_count": int(complete_count),
            "incomplete_count": int(len(levels) - complete_count),
            "incomplete_levels": incomplete,
        }

    E = largest_bad_degree + 1
    if d <= 2:
        E = 1
    result = {
        "schema_version": 1,
        "q": int(q),
        "p": str(ctx.p),
        "degree_p": int(d),
        "theorem_bound_U": int(U),
        "E": int(E),
        "cutoff_certificate": (
            f"Theorem 1.2 certifies every level degree >= {U}"
        ),
        "class_count": int(n),
        "weights": [int(w) for w in weights],
        "hasse_polynomial": str(ctx.hasse_polynomial),
        "supersingular_polynomial": str(ctx.supersingular_polynomial),
        "j_invariants": [str(j) for j in ctx.j_invariants],
        "construction_checks": ctx.checks,
        "full_matrices_requested": bool(full_matrices),
        "degree_summary": degree_summary,
        "levels": level_records,
        "enumerations": enumeration_metadata,
        "total_vectors_enumerated": int(
            sum(item["vectors_enumerated"] for item in enumeration_metadata)
        ),
    }
    return result
