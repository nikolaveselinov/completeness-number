"""Compact exact-zero computation over odd prime constant fields.

This backend retains the Appendix A predicate rather than materializing every
positive Brandt entry. For each unordered vertex pair it exhausts the complete
bounded Hom space through the largest degree not covered by the exact
pair-specific spectral inequality. Weighted Brandt symmetry transports only
zero versus positivity to the reverse orientation.
"""

from __future__ import annotations

import hashlib
import subprocess
from dataclasses import dataclass
from pathlib import Path
from time import perf_counter
from typing import Any

from sage.all import Hom, is_prime

from .core import (
    SupersingularContext,
    _coefficient_backend,
    _raw_norm,
    monic_irreducibles,
    pair_is_spectrally_positive,
    theorem_bound,
)


REPOSITORY = Path(__file__).resolve().parents[1]
COUNTER_SOURCE = REPOSITORY / "fast_prime" / "prime_norm_counter.cpp"
COUNTER_BINARY = REPOSITORY / "fast_prime" / "prime_norm_counter"
COMPACT_KIND = "odd_prime_compact_zero_v1"


@dataclass(frozen=True)
class PrimeCounterResult:
    iterations: int
    target_count: int
    invalid_norms: int
    seconds: float
    counts: dict[int, int]


def compile_prime_counter() -> None:
    """Build the deterministic helper when absent or older than its source."""
    if (
        COUNTER_BINARY.exists()
        and COUNTER_BINARY.stat().st_mtime >= COUNTER_SOURCE.stat().st_mtime
    ):
        return
    subprocess.run(
        [
            "g++",
            "-O3",
            "-std=c++20",
            "-DNDEBUG",
            str(COUNTER_SOURCE),
            "-o",
            str(COUNTER_BINARY),
        ],
        check=True,
    )


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def polynomial_code(poly, prime: int) -> int:
    """Encode an F_p[T] polynomial in base p, low coefficient first."""
    place = 1
    code = 0
    for coefficient in poly.list():
        code += int(coefficient) * place
        place *= prime
    return int(code)


def _prime_coefficients(raw, prime: int, max_degree: int) -> list[int]:
    """Extract a padded F_p[T] coefficient vector, rejecting field escape."""
    result = [0] * (max_degree + 1)
    for degree, coefficient in enumerate(raw.list()):
        if degree > max_degree:
            raise ArithmeticError("raw norm exceeded the requested degree")
        coefficient = _coefficient_backend(coefficient)
        if coefficient**prime != coefficient:
            raise ArithmeticError(
                "interpolated norm form is not defined over the prime field"
            )
        value = int(coefficient)
        if not 0 <= value < prime:
            raise ArithmeticError("invalid prime-field coefficient encoding")
        result[degree] = value
    return result


def _normalize_coefficients(
    coefficients: list[int], prime: int
) -> tuple[int, int] | None:
    degree = len(coefficients) - 1
    while degree >= 0 and coefficients[degree] == 0:
        degree -= 1
    if degree < 0:
        return None
    inverse = pow(coefficients[degree], -1, prime)
    code = 0
    place = 1
    for value in coefficients[: degree + 1]:
        code += ((value * inverse) % prime) * place
        place *= prime
    return degree, code


def _evaluate_form(
    diagonal: list[list[int]],
    cross: list[list[list[int]]],
    vector: tuple[int, ...],
    prime: int,
) -> list[int]:
    value = [0] * len(diagonal[0])
    for i, coefficient in enumerate(vector):
        if coefficient:
            square = coefficient * coefficient
            for degree, entry in enumerate(diagonal[i]):
                value[degree] = (value[degree] + square * entry) % prime
            for j in range(i):
                product = coefficient * vector[j]
                if product:
                    for degree, entry in enumerate(cross[j][i]):
                        value[degree] = (
                            value[degree] + product * entry
                        ) % prime
    return value


def _sample_vectors(
    dimension: int, prime: int, count: int = 8
) -> list[tuple[int, ...]]:
    """Return deterministic samples that exercise diagonal and cross terms."""
    candidates: list[tuple[int, ...]] = []
    if dimension > 1:
        candidates.extend(
            [
                tuple(1 for _ in range(dimension)),
                tuple((index % (prime - 1)) + 1 for index in range(dimension)),
                tuple(1 if index % 2 == 0 else 0 for index in range(dimension)),
            ]
        )
    for index in range(dimension):
        vector = [0] * dimension
        vector[index] = (index % (prime - 1)) + 1
        if dimension > 1:
            vector[(index + 1) % dimension] = 1
        candidates.append(tuple(vector))
        unit = [0] * dimension
        unit[index] = 1
        candidates.append(tuple(unit))
    result: list[tuple[int, ...]] = []
    for vector in candidates:
        if any(vector) and vector not in result:
            result.append(vector)
        if len(result) == count:
            break
    if dimension > 1 and not any(
        sum(value != 0 for value in vector) > 1 for vector in result
    ):
        raise ArithmeticError("direct sample failed to exercise cross terms")
    return result


def interpolate_prime_form(
    ctx: SupersingularContext,
    source: int,
    target: int,
    max_degree: int,
) -> dict[str, Any]:
    """Interpolate a complete prime-field norm form and cross-check it."""
    started = perf_counter()
    prime = ctx.q
    homset = Hom(ctx.modules[source], ctx.modules[target])
    basis = homset.basis(degree=max_degree)
    dimension = len(basis)
    expected_dimension = 2 * (max_degree + 1) - (ctx.degree - 1)
    if dimension != expected_dimension:
        raise ArithmeticError(
            f"bounded Hom dimension {dimension}, expected "
            f"{expected_dimension} for pair {(source, target)}"
        )

    diagonal_raw = [_raw_norm(element) for element in basis]
    diagonal = [
        _prime_coefficients(value, prime, max_degree)
        for value in diagonal_raw
    ]
    zero = [0] * (max_degree + 1)
    cross = [
        [list(zero) for _ in range(dimension)]
        for _ in range(dimension)
    ]
    for i in range(dimension):
        for j in range(i + 1, dimension):
            raw = (
                _raw_norm(basis[i] + basis[j])
                - diagonal_raw[i]
                - diagonal_raw[j]
            )
            packed = _prime_coefficients(raw, prime, max_degree)
            cross[i][j] = packed
            cross[j][i] = packed

    checked = 0
    for vector in _sample_vectors(dimension, prime):
        normalized = _normalize_coefficients(
            _evaluate_form(diagonal, cross, vector, prime), prime
        )
        if normalized is None:
            raise ArithmeticError("nonzero Hom vector has zero packed norm")
        _, code = normalized
        ore = ctx.modules[source].ore_polring().zero()
        for coefficient, element in zip(vector, basis):
            if coefficient:
                ore += coefficient * element.ore_polynomial()
        direct = homset(ore).norm().gen()
        direct_code = polynomial_code(direct, prime)
        if code != direct_code:
            raise ArithmeticError(
                f"packed norm {code} disagrees with direct Sage norm "
                f"{direct_code} for pair {(source, target)}"
            )
        checked += 1

    return {
        "diagonal": diagonal,
        "cross": cross,
        "dimension": dimension,
        "direct_sage_norm_cross_checks": checked,
        "interpolation_seconds": perf_counter() - started,
    }


def run_prime_counter(
    prime: int,
    form: dict[str, Any],
    max_degree: int,
    target_codes: list[int],
) -> PrimeCounterResult:
    compile_prime_counter()
    dimension = int(form["dimension"])
    tokens = [
        str(prime),
        str(dimension),
        str(max_degree),
        str(len(target_codes)),
    ]
    for value in form["diagonal"]:
        tokens.extend(str(int(coefficient)) for coefficient in value)
    for i in range(dimension):
        for j in range(i + 1, dimension):
            tokens.extend(
                str(int(coefficient)) for coefficient in form["cross"][i][j]
            )
    tokens.extend(str(int(code)) for code in target_codes)
    completed = subprocess.run(
        [str(COUNTER_BINARY)],
        input=" ".join(tokens),
        text=True,
        capture_output=True,
        check=True,
    )
    lines = completed.stdout.splitlines()
    header = lines[0].split()
    if header[0] != "SUMMARY" or len(header) != 5:
        raise RuntimeError(f"unexpected counter output: {lines[0]!r}")
    counts: dict[int, int] = {}
    for line in lines[1:]:
        label, code, count = line.split()
        if label != "COUNT":
            raise RuntimeError(f"unexpected counter output: {line!r}")
        counts[int(code)] = int(count)
    return PrimeCounterResult(
        iterations=int(header[1]),
        target_count=int(header[2]),
        invalid_norms=int(header[3]),
        seconds=float(header[4]),
        counts=counts,
    )


def compute_compact_odd_prime(
    ctx: SupersingularContext,
    *,
    progress=None,
) -> dict[str, Any]:
    """Compute an exact compact completeness archive for one odd-prime p."""
    prime = ctx.q
    if not is_prime(prime) or prime == 2:
        raise ValueError("this compact backend requires an odd prime field")
    if ctx.degree < 3:
        raise ValueError("use the theorem directly for degree one or two")
    if not COUNTER_SOURCE.exists():
        raise FileNotFoundError(f"missing exact counter source: {COUNTER_SOURCE}")
    compile_prime_counter()

    cutoff = theorem_bound(prime, ctx.degree)
    levels_by_degree = {
        degree: [
            ell
            for ell in monic_irreducibles(ctx.A, degree)
            if ell != ctx.p
        ]
        for degree in range(1, cutoff)
    }
    level_by_code = {
        polynomial_code(ell, prime): ell
        for levels in levels_by_degree.values()
        for ell in levels
    }
    codes_by_degree = {
        degree: sorted(polynomial_code(ell, prime) for ell in levels)
        for degree, levels in levels_by_degree.items()
    }
    if len(level_by_code) != sum(map(len, levels_by_degree.values())):
        raise ArithmeticError("prime-polynomial encoding collision")

    vertex_count = len(ctx.modules)
    zero_pairs_by_code: dict[int, list[tuple[int, int]]] = {
        code: [] for code in level_by_code
    }
    pair_certificates: list[dict[str, Any]] = []
    total_iterations = 0
    total_interpolation_seconds = 0.0
    total_counter_seconds = 0.0
    pair_total = vertex_count * (vertex_count + 1) // 2
    pair_number = 0

    for source in range(vertex_count):
        for target in range(source, vertex_count):
            pair_number += 1
            non_spectral = [
                degree
                for degree in range(1, cutoff)
                if not pair_is_spectrally_positive(
                    ctx, source, target, degree
                )
            ]
            spectral = [
                degree
                for degree in range(1, cutoff)
                if degree not in non_spectral
            ]
            run = None
            if non_spectral:
                max_degree = max(non_spectral)
                target_codes = [
                    code
                    for degree in non_spectral
                    for code in codes_by_degree[degree]
                ]
                form = interpolate_prime_form(
                    ctx, source, target, max_degree
                )
                counter = run_prime_counter(
                    prime, form, max_degree, target_codes
                )
                expected_iterations = prime ** int(form["dimension"]) - 1
                if counter.iterations != expected_iterations:
                    raise ArithmeticError(
                        f"counter stopped after {counter.iterations}, "
                        f"expected {expected_iterations}"
                    )
                if counter.invalid_norms:
                    raise ArithmeticError("compiled counter found zero norms")
                if counter.target_count != len(target_codes):
                    raise ArithmeticError("counter target count changed")
                zero_codes = sorted(
                    code
                    for code in target_codes
                    if counter.counts.get(code, 0) == 0
                )
                for code, count in counter.counts.items():
                    if count % ctx.automorphism_orders[target]:
                        raise ArithmeticError(
                            f"count {count} is not divisible by target Aut"
                        )
                    if count % ctx.automorphism_orders[source]:
                        raise ArithmeticError(
                            f"count {count} is not divisible by source Aut"
                        )
                for code in zero_codes:
                    zero_pairs_by_code[code].append((source, target))
                run = {
                    "target_degrees": [int(value) for value in non_spectral],
                    "max_degree": int(max_degree),
                    "dimension": int(form["dimension"]),
                    "iterations": int(counter.iterations),
                    "target_count": int(counter.target_count),
                    "invalid_norms": int(counter.invalid_norms),
                    "exhaustive": True,
                    "zero_codes": [int(code) for code in zero_codes],
                    "counts": [
                        [int(code), int(counter.counts.get(code, 0))]
                        for code in target_codes
                    ],
                    "direct_sage_norm_cross_checks": int(
                        form["direct_sage_norm_cross_checks"]
                    ),
                    "interpolation_seconds": float(
                        form["interpolation_seconds"]
                    ),
                    "counter_seconds": float(counter.seconds),
                }
                total_iterations += counter.iterations
                total_interpolation_seconds += form["interpolation_seconds"]
                total_counter_seconds += counter.seconds

            pair_certificates.append(
                {
                    "pair": [int(source), int(target)],
                    "spectral_degrees": [int(value) for value in spectral],
                    "run": run,
                }
            )
            if progress is not None and (
                pair_number == 1
                or pair_number == pair_total
                or pair_number % 25 == 0
            ):
                progress(
                    f"pair {pair_number}/{pair_total} ({source},{target}); "
                    f"vectors={total_iterations}"
                )

    levels: list[dict[str, Any]] = []
    degree_summary: dict[str, dict[str, Any]] = {}
    largest_bad_degree = 0
    for degree in range(1, cutoff):
        incomplete: list[str] = []
        for code in codes_by_degree[degree]:
            ell = level_by_code[code]
            ordered: list[list[int]] = []
            for source, target in zero_pairs_by_code[code]:
                ordered.append([source, target])
                if source != target:
                    ordered.append([target, source])
            ordered.sort()
            complete = not ordered
            if not complete:
                incomplete.append(str(ell))
                largest_bad_degree = max(largest_bad_degree, degree)
            levels.append(
                {
                    "ell": str(ell),
                    "code": int(code),
                    "degree": int(degree),
                    "complete": bool(complete),
                    "zero_entries": ordered,
                }
            )
        degree_summary[str(degree)] = {
            "level_count": len(codes_by_degree[degree]),
            "complete_count": (
                len(codes_by_degree[degree]) - len(incomplete)
            ),
            "incomplete_count": len(incomplete),
            "incomplete_levels": incomplete,
        }

    return {
        "archive_kind": COMPACT_KIND,
        "schema_version": 1,
        "q": int(prime),
        "p": str(ctx.p),
        "degree_p": int(ctx.degree),
        "theorem_bound_U": int(cutoff),
        "E": int(largest_bad_degree + 1),
        "class_count": int(vertex_count),
        "weights": [
            int(order // (prime - 1))
            for order in ctx.automorphism_orders
        ],
        "automorphism_orders": [
            int(order) for order in ctx.automorphism_orders
        ],
        "j_invariants": [str(value) for value in ctx.j_invariants],
        "construction_checks": ctx.checks,
        "engine": {
            "name": "compiled_odd_prime_gray_zero_archive",
            "counter_source": str(COUNTER_SOURCE.relative_to(REPOSITORY)),
            "counter_source_sha256": _sha256(COUNTER_SOURCE),
            "counter_binary_sha256": _sha256(COUNTER_BINARY),
            "compiler_contract": "g++ -O3 -std=c++20 -DNDEBUG",
            "unordered_pair_identity": (
                "w_j*b_ij(ell)=w_i*b_ji(ell); positive weights make "
                "zero/positivity orientation-invariant"
            ),
            "zero_semantics": (
                "every non-spectral zero follows a complete p^r-1 "
                "bounded-Hom traversal"
            ),
            "total_interpolation_seconds": total_interpolation_seconds,
            "total_counter_seconds": total_counter_seconds,
        },
        "pair_certificates": pair_certificates,
        "levels": levels,
        "degree_summary": degree_summary,
        "total_vectors_visited": int(total_iterations),
    }
