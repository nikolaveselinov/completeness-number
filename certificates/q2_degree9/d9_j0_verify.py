#!/usr/bin/env python3
"""Independent Sage-free verification of the j=0 Brandt self-entry for q=2.

For odd ``d`` the supersingular class with ``j=0`` is represented by

    phi_T = theta + tau^2,

where theta is the image of T in F_2[T]/(p).  This script constructs the
bounded Hom space directly as the F_2-nullspace of

    u phi_T = phi_T u,

computes the determinant of right multiplication by u on the rank-two
Drinfeld motive, interpolates its exact quadratic norm form, and passes that
form to the repository's deterministic exhaustive Gray-code counter.

It deliberately does not import Sage or the main Python package.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import subprocess
from pathlib import Path
from time import perf_counter


def polynomial_degree(value: int) -> int:
    return value.bit_length() - 1


def polynomial_mod(value: int, modulus: int) -> int:
    modulus_degree = polynomial_degree(modulus)
    while value and polynomial_degree(value) >= modulus_degree:
        value ^= modulus << (polynomial_degree(value) - modulus_degree)
    return value


def polynomial_gcd(left: int, right: int) -> int:
    while right:
        left, right = right, polynomial_mod(left, right)
    return left


def polynomial_square_mod(value: int, modulus: int) -> int:
    squared = 0
    while value:
        low = value & -value
        squared ^= 1 << (2 * (low.bit_length() - 1))
        value ^= low
    return polynomial_mod(squared, modulus)


def prime_divisors(value: int) -> list[int]:
    result: list[int] = []
    candidate = 2
    while candidate * candidate <= value:
        if value % candidate == 0:
            result.append(candidate)
            while value % candidate == 0:
                value //= candidate
        candidate += 1
    if value > 1:
        result.append(value)
    return result


def is_irreducible_binary(modulus: int, degree: int) -> bool:
    if polynomial_degree(modulus) != degree or not (modulus & 1):
        return False
    x = 2
    frobenius = x
    test_exponents = {degree // prime for prime in prime_divisors(degree)}
    for exponent in range(1, degree + 1):
        frobenius = polynomial_square_mod(frobenius, modulus)
        if (
            exponent in test_exponents
            and polynomial_gcd(frobenius ^ x, modulus) != 1
        ):
            return False
    return frobenius == x


def monic_irreducibles_binary(degree: int) -> list[int]:
    return [
        code
        for code in range((1 << degree) | 1, 1 << (degree + 1), 2)
        if is_irreducible_binary(code, degree)
    ]


def polynomial_text(code: int, variable: str = "T") -> str:
    terms: list[str] = []
    for exponent in range(polynomial_degree(code), -1, -1):
        if not (code >> exponent) & 1:
            continue
        if exponent == 0:
            terms.append("1")
        elif exponent == 1:
            terms.append(variable)
        else:
            terms.append(f"{variable}^{exponent}")
    return "+".join(terms) if terms else "0"


class QuadraticTower:
    """F_{2^(2d)} = F_{2^d}[z]/(z^2+z+1), for odd d."""

    def __init__(self, degree: int, modulus: int):
        if degree % 2 == 0:
            raise ValueError("z^2+z+1 is used only for odd base degree")
        if not is_irreducible_binary(modulus, degree):
            raise ValueError("the characteristic polynomial is not irreducible")
        self.degree = degree
        self.field_degree = 2 * degree
        self.modulus = modulus
        self.base_mask = (1 << degree) - 1
        size = 1 << degree
        self.base_product = [
            [self._base_multiply(left, right) for right in range(size)]
            for left in range(size)
        ]

    def _base_multiply(self, left: int, right: int) -> int:
        result = 0
        while right:
            if right & 1:
                result ^= left
            right >>= 1
            left <<= 1
            if left & (1 << self.degree):
                left ^= self.modulus
        return result & self.base_mask

    def multiply(self, left: int, right: int) -> int:
        a = left & self.base_mask
        b = left >> self.degree
        c = right & self.base_mask
        d = right >> self.degree
        product = self.base_product
        constant = product[a][c] ^ product[b][d]
        z_coefficient = product[a][d] ^ product[b][c] ^ product[b][d]
        return constant | (z_coefficient << self.degree)

    def square(self, value: int) -> int:
        return self.multiply(value, value)

    def frobenius_power(self, value: int, exponent: int) -> int:
        for _ in range(exponent):
            value = self.square(value)
        return value


def rref_nullspace(rows: list[int], variable_count: int) -> list[int]:
    original = rows[:]
    rows = [row for row in rows if row]
    pivot_columns: list[int] = []
    rank = 0
    for column in range(variable_count):
        pivot_index = next(
            (
                index
                for index in range(rank, len(rows))
                if (rows[index] >> column) & 1
            ),
            None,
        )
        if pivot_index is None:
            continue
        rows[rank], rows[pivot_index] = rows[pivot_index], rows[rank]
        pivot_row = rows[rank]
        for index in range(len(rows)):
            if index != rank and ((rows[index] >> column) & 1):
                rows[index] ^= pivot_row
        pivot_columns.append(column)
        rank += 1
        if rank == len(rows):
            break
    rows = rows[:rank]
    pivot_set = set(pivot_columns)
    free_columns = [
        column for column in range(variable_count) if column not in pivot_set
    ]
    basis: list[int] = []
    for free in free_columns:
        vector = 1 << free
        for pivot, row in zip(pivot_columns, rows):
            if (row >> free) & 1:
                vector |= 1 << pivot
        basis.append(vector)
    if any(
        (row & vector).bit_count() % 2
        for vector in basis
        for row in original
    ):
        raise ArithmeticError("nullspace construction failed")
    return basis


def hom_basis(
    field: QuadraticTower,
    max_degree: int,
    source_g: int,
    target_g: int,
) -> list[list[int]]:
    """Return an F_2-basis for one bounded Hom space."""
    coefficient_bits = field.field_degree
    equation_count = (max_degree + 3) * coefficient_bits
    rows = [0] * equation_count
    theta = 2
    theta_power = theta
    source_g_power = source_g
    for tau_degree in range(max_degree + 1):
        if tau_degree:
            theta_power = field.square(theta_power)
            source_g_power = field.square(source_g_power)
        diagonal_factor = theta_power ^ theta
        for coordinate in range(coefficient_bits):
            element = 1 << coordinate
            column = tau_degree * coefficient_bits + coordinate
            element_squared = field.square(element)
            fourth_power = field.square(field.square(element))
            contributions = (
                (
                    tau_degree,
                    field.multiply(element, diagonal_factor),
                ),
                (
                    tau_degree + 1,
                    field.multiply(element, source_g_power)
                    ^ field.multiply(target_g, element_squared),
                ),
                (tau_degree + 2, element ^ fourth_power),
            )
            for equation_degree, value in contributions:
                while value:
                    low = value & -value
                    output_coordinate = low.bit_length() - 1
                    rows[
                        equation_degree * coefficient_bits + output_coordinate
                    ] ^= 1 << column
                    value ^= low
    variable_count = (max_degree + 1) * coefficient_bits
    bit_basis = rref_nullspace(rows, variable_count)
    expected_dimension = (
        2 * (max_degree + 1) - (field.degree - 1)
    )
    if len(bit_basis) != expected_dimension:
        raise ArithmeticError(
            f"Hom dimension {len(bit_basis)}, expected {expected_dimension}"
        )
    mask = (1 << coefficient_bits) - 1
    return [
        [
            (vector >> (tau_degree * coefficient_bits)) & mask
            for tau_degree in range(max_degree + 1)
        ]
        for vector in bit_basis
    ]


def polynomial_add(left: list[int], right: list[int]) -> list[int]:
    result = [0] * max(len(left), len(right))
    for index in range(len(result)):
        result[index] = (
            (left[index] if index < len(left) else 0)
            ^ (right[index] if index < len(right) else 0)
        )
    while len(result) > 1 and result[-1] == 0:
        result.pop()
    return result


def polynomial_scale(
    polynomial: list[int], scalar: int, field: QuadraticTower
) -> list[int]:
    return [field.multiply(coefficient, scalar) for coefficient in polynomial]


def polynomial_multiply(
    left: list[int], right: list[int], field: QuadraticTower
) -> list[int]:
    result = [0] * (len(left) + len(right) - 1)
    for i, left_coefficient in enumerate(left):
        if not left_coefficient:
            continue
        for j, right_coefficient in enumerate(right):
            if right_coefficient:
                result[i + j] ^= field.multiply(
                    left_coefficient, right_coefficient
                )
    while len(result) > 1 and result[-1] == 0:
        result.pop()
    return result


def motive_coordinates(
    field: QuadraticTower, max_tau_degree: int, source_g: int
) -> tuple[list[list[int]], list[list[int]]]:
    """Coordinates of tau^k in the motive basis (1,tau)."""
    theta = 2
    constant_coordinates = [[1], [0]]
    tau_coordinates = [[0], [1]]
    for tau_degree in range(1, max_tau_degree):
        a = constant_coordinates[tau_degree]
        b = tau_coordinates[tau_degree]
        a_squared = [field.square(coefficient) for coefficient in a]
        b_squared = [field.square(coefficient) for coefficient in b]
        next_a = [0] * (len(b_squared) + 1)
        for index, coefficient in enumerate(b_squared):
            next_a[index] ^= field.multiply(coefficient, theta)
            next_a[index + 1] ^= coefficient
        constant_coordinates.append(next_a)
        tau_coordinates.append(
            polynomial_add(
                a_squared,
                polynomial_scale(b_squared, source_g, field),
            )
        )
    return constant_coordinates, tau_coordinates


def raw_motive_norm(
    coefficients: list[int],
    field: QuadraticTower,
    constant_coordinates: list[list[int]],
    tau_coordinates: list[list[int]],
) -> list[int]:
    """Determinant of right multiplication by an endomorphism."""
    u_constant = [0]
    u_tau = [0]
    tau_u_constant = [0]
    tau_u_tau = [0]
    for tau_degree, coefficient in enumerate(coefficients):
        if not coefficient:
            continue
        u_constant = polynomial_add(
            u_constant,
            polynomial_scale(
                constant_coordinates[tau_degree], coefficient, field
            ),
        )
        u_tau = polynomial_add(
            u_tau,
            polynomial_scale(tau_coordinates[tau_degree], coefficient, field),
        )
        squared = field.square(coefficient)
        tau_u_constant = polynomial_add(
            tau_u_constant,
            polynomial_scale(
                constant_coordinates[tau_degree + 1], squared, field
            ),
        )
        tau_u_tau = polynomial_add(
            tau_u_tau,
            polynomial_scale(
                tau_coordinates[tau_degree + 1], squared, field
            ),
        )
    return polynomial_add(
        polynomial_multiply(u_constant, tau_u_tau, field),
        polynomial_multiply(tau_u_constant, u_tau, field),
    )


def pack_polynomial(
    polynomial: list[int], max_degree: int, field_degree: int
) -> int:
    if len(polynomial) - 1 > max_degree:
        raise ArithmeticError("raw norm exceeded the requested degree")
    return sum(
        coefficient << (index * field_degree)
        for index, coefficient in enumerate(polynomial)
    )


def normalized_code(
    packed: int, max_degree: int, field_degree: int
) -> int:
    mask = (1 << field_degree) - 1
    degree = max_degree
    while degree >= 0:
        leading = (packed >> (degree * field_degree)) & mask
        if leading:
            break
        degree -= 1
    if degree < 0:
        raise ArithmeticError("nonzero Hom vector has zero norm")
    code = 0
    for index in range(degree + 1):
        coefficient = (packed >> (index * field_degree)) & mask
        if coefficient == leading:
            code |= 1 << index
        elif coefficient:
            raise ArithmeticError("raw norm is not coherent over F_2")
    return code


def interpolate_norm_form(
    field: QuadraticTower,
    basis: list[list[int]],
    max_degree: int,
    source_g: int,
) -> tuple[list[int], list[list[int]]]:
    coordinates = motive_coordinates(field, max_degree + 2, source_g)
    diagonal = [
        pack_polynomial(
            raw_motive_norm(vector, field, *coordinates),
            max_degree,
            field.field_degree,
        )
        for vector in basis
    ]
    dimension = len(basis)
    cross = [[0] * dimension for _ in range(dimension)]
    for left in range(dimension):
        for right in range(left + 1, dimension):
            vector = [
                a ^ b for a, b in zip(basis[left], basis[right])
            ]
            packed = pack_polynomial(
                raw_motive_norm(vector, field, *coordinates),
                max_degree,
                field.field_degree,
            )
            cross[left][right] = (
                packed ^ diagonal[left] ^ diagonal[right]
            )
            cross[right][left] = cross[left][right]

    # Exercise every diagonal and cross term before invoking compiled code.
    for index, packed in enumerate(diagonal):
        normalized_code(packed, max_degree, field.field_degree)
        if index + 1 < dimension:
            normalized_code(
                packed
                ^ diagonal[index + 1]
                ^ cross[index][index + 1],
                max_degree,
                field.field_degree,
            )
    return diagonal, cross


def supersingular_roots(field: QuadraticTower) -> list[int]:
    """Evaluate Gekeler's recurrence on every element of F_{2^(2d)}."""
    import numpy as np

    degree = field.degree
    size = 1 << field.field_degree
    base_mask = field.base_mask
    product = np.asarray(field.base_product, dtype=np.uint16)

    def multiply(left, right):
        left = np.asarray(left, dtype=np.uint32)
        right = np.asarray(right, dtype=np.uint32)
        a = left & base_mask
        b = left >> degree
        c = right & base_mask
        d = right >> degree
        constant = product[a, c] ^ product[b, d]
        z_coefficient = product[a, d] ^ product[b, c] ^ product[b, d]
        return constant.astype(np.uint32) | (
            z_coefficient.astype(np.uint32) << degree
        )

    def power(values, exponent: int):
        result = np.ones_like(values, dtype=np.uint32)
        base = values
        while exponent:
            if exponent & 1:
                result = multiply(result, base)
            exponent >>= 1
            if exponent:
                base = multiply(base, base)
        return result

    values = np.arange(size, dtype=np.uint32)
    previous = np.ones(size, dtype=np.uint32)
    current = np.ones(size, dtype=np.uint32)
    theta = 2
    theta_power = theta
    for recurrence_index in range(2, degree + 1):
        theta_power = field.square(theta_power)
        exponent = (
            2 ** (recurrence_index - 1)
            - (-1) ** (recurrence_index - 1)
        ) // 3
        bracket = theta_power ^ theta
        following = multiply(power(values, exponent), current) ^ multiply(
            bracket, previous
        )
        previous, current = current, following
    supersingular_values = multiply(values, current)
    roots = np.flatnonzero(supersingular_values == 0).astype(int).tolist()
    expected = (2**degree - 2) // 3 + 1
    if len(roots) != expected:
        raise ArithmeticError(
            f"found {len(roots)} supersingular roots, expected {expected}"
        )
    return roots


def cube_root_representatives(
    field: QuadraticTower, roots: list[int]
) -> dict[int, int]:
    root_set = set(roots)
    representatives: dict[int, int] = {}
    for value in range(1 << field.field_degree):
        cube = field.multiply(field.square(value), value)
        if cube in root_set and cube not in representatives:
            representatives[cube] = value
            if len(representatives) == len(roots):
                break
    missing = root_set - representatives.keys()
    if missing:
        raise ArithmeticError(
            f"{len(missing)} supersingular j-invariants have no cube root"
        )
    return representatives


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def packed_words(value: int, count: int) -> list[str]:
    mask = (1 << 64) - 1
    return [hex((value >> (64 * index)) & mask) for index in range(count)]


def run_counter(
    binary: Path,
    diagonal: list[int],
    cross: list[list[int]],
    *,
    max_degree: int,
    field_degree: int,
    target_codes: list[int],
    stop_when_seen: bool = False,
) -> dict[str, object]:
    dimension = len(diagonal)
    word_count = math.ceil((max_degree + 1) * field_degree / 64)
    tokens = [
        str(dimension),
        str(max_degree),
        str(field_degree),
        str(word_count),
        str(len(target_codes)),
        str(int(stop_when_seen)),
    ]
    for value in diagonal:
        tokens.extend(packed_words(value, word_count))
    for left in range(dimension):
        for right in range(left + 1, dimension):
            tokens.extend(packed_words(cross[left][right], word_count))
    tokens.extend(str(code) for code in target_codes)
    completed = subprocess.run(
        [str(binary)],
        input=" ".join(tokens),
        text=True,
        capture_output=True,
        check=True,
    )
    lines = completed.stdout.splitlines()
    header = lines[0].split()
    if header[0] != "SUMMARY" or len(header) != 7:
        raise RuntimeError(f"unexpected counter output: {lines[0]!r}")
    counts: dict[int, int] = {}
    for line in lines[1:]:
        label, code, count = line.split()
        if label != "COUNT":
            raise RuntimeError(f"unexpected counter output: {line!r}")
        counts[int(code)] = int(count)
    return {
        "iterations": int(header[1]),
        "exhaustive": bool(int(header[2])),
        "seen": int(header[3]),
        "target_count": int(header[4]),
        "invalid_norms": int(header[5]),
        "counter_seconds": float(header[6]),
        "counts": counts,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--degree-p", type=int, required=True)
    parser.add_argument(
        "--p-code",
        type=lambda value: int(value, 0),
        required=True,
        help="binary bit code, including the monic leading bit",
    )
    parser.add_argument("--level-degree", type=int, required=True)
    parser.add_argument("--counter", type=Path, required=True)
    parser.add_argument(
        "--all-diagonals",
        action="store_true",
        help="scan every supersingular diagonal and stop at the first zero",
    )
    args = parser.parse_args()

    started = perf_counter()
    field = QuadraticTower(args.degree_p, args.p_code)
    targets = monic_irreducibles_binary(args.level_degree)
    if args.all_diagonals:
        roots = supersingular_roots(field)
        cube_roots = cube_root_representatives(field, roots)
        diagonal_classes = [
            (index, j, cube_roots[j]) for index, j in enumerate(roots)
        ]
    else:
        diagonal_classes = [(0, 0, 0)]

    selected: dict[str, object] | None = None
    scanned: list[dict[str, object]] = []
    for position, (class_index, j_invariant, g) in enumerate(
        diagonal_classes, start=1
    ):
        pair_started = perf_counter()
        basis = hom_basis(
            field,
            args.level_degree,
            source_g=g,
            target_g=g,
        )
        diagonal, cross = interpolate_norm_form(
            field,
            basis,
            args.level_degree,
            source_g=g,
        )
        result = run_counter(
            args.counter.resolve(),
            diagonal,
            cross,
            max_degree=args.level_degree,
            field_degree=field.field_degree,
            target_codes=targets,
        )
        counts = result.pop("counts")
        zero_codes = sorted(
            code for code in targets if counts.get(code, 0) == 0
        )
        entry = {
            "class_index": class_index,
            "j_code": j_invariant,
            "g_code": g,
            "hom_dimension": len(basis),
            **result,
            "zero_codes": zero_codes,
            "zero_levels": [polynomial_text(code) for code in zero_codes],
            "zero_multiplicities": [counts[code] for code in zero_codes],
            "seconds": perf_counter() - pair_started,
        }
        scanned.append(entry)
        print(
            f"diagonal {position}/{len(diagonal_classes)}: "
            f"j={j_invariant}, zeros={len(zero_codes)}, "
            f"seconds={entry['seconds']:.3f}",
            flush=True,
        )
        if zero_codes:
            selected = entry
            break

    report = {
        "method": "independent_sage_free_motive_determinant",
        "q": 2,
        "p": polynomial_text(args.p_code),
        "p_code": args.p_code,
        "degree_p": args.degree_p,
        "level_degree": args.level_degree,
        "target_count": len(targets),
        "counter_binary": str(args.counter.resolve()),
        "counter_binary_sha256": sha256(args.counter.resolve()),
        "all_diagonals_requested": args.all_diagonals,
        "diagonals_available": len(diagonal_classes),
        "diagonals_scanned": len(scanned),
        "first_obstruction": selected,
        "scan": scanned,
        "total_seconds": perf_counter() - started,
    }
    if any(not entry["exhaustive"] for entry in scanned):
        raise ArithmeticError("a counter did not traverse the full Hom space")
    if any(entry["invalid_norms"] for entry in scanned):
        raise ArithmeticError("a counter found an invalid raw norm")
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
