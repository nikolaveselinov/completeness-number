"""Compact exact-zero computation for q=2, degree-seven characteristics.

The standard backend materializes ordered Brandt matrices.  In degree seven
that repeats both orientations of 946 unordered vertex pairs at 4,719 prime
levels.  This module records exactly the information needed for completeness:
the zero positions, together with exhaustive, witnessed-positive, spectral,
and weighted-symmetry certificates.

No early-stop multiplicity is exposed as an exact Brandt entry.
"""

from __future__ import annotations

import hashlib
from collections import defaultdict
from pathlib import Path
from time import perf_counter
from typing import Any

from sage.all import Hom

from fast_q2.gray_counter import CounterResult, run_counter

from .core import (
    SupersingularContext,
    _normalize_packed_q2,
    _poly_code_q2,
    _quadratic_form_q2,
    monic_irreducibles,
    pair_is_spectrally_positive,
    theorem_bound,
)


REPOSITORY = Path(__file__).resolve().parents[1]
COUNTER_SOURCE = REPOSITORY / "fast_q2" / "q2_norm_counter.cpp"
COUNTER_BINARY = REPOSITORY / "fast_q2" / "q2_norm_counter"
COMPACT_KIND = "q2_compact_zero_v1"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _evaluate_form(
    diagonal: list[int],
    cross: list[list[int]],
    vector: int,
) -> int:
    """Evaluate a characteristic-two packed quadratic form."""
    active: list[int] = []
    value = 0
    remaining = vector
    while remaining:
        low = remaining & -remaining
        index = low.bit_length() - 1
        value ^= diagonal[index]
        for earlier in active:
            value ^= cross[earlier][index]
        active.append(index)
        remaining ^= low
    return value


def _deterministic_sample_vectors(dimension: int, count: int) -> list[int]:
    """Return distinct nonzero coordinate masks without random state."""
    limit = (1 << dimension) - 1
    # Put multi-coordinate vectors first: cross terms are as important as
    # diagonal terms, and a short sample must exercise both.
    candidates = [
        limit,
        sum(1 << k for k in range(0, dimension, 2)),
        sum(1 << k for k in range(1, dimension, 2)),
    ]
    for k in range(dimension):
        if dimension > 1:
            candidates.append((1 << k) | (1 << ((k + 1) % dimension)))
        candidates.append(1 << k)
    state = 0x9E3779B97F4A7C15
    while len(candidates) < count * 3 + 8:
        state = (
            state * 6364136223846793005 + 1442695040888963407
        ) & ((1 << 64) - 1)
        candidates.append((state & limit) or 1)
    result: list[int] = []
    for vector in candidates:
        if vector and vector not in result:
            result.append(vector)
        if len(result) == min(count, limit):
            break
    if dimension > 1 and not any(vector.bit_count() > 1 for vector in result):
        raise ArithmeticError("direct norm sample did not exercise a cross term")
    return result


def _direct_sage_cross_checks(
    ctx: SupersingularContext,
    source: int,
    target: int,
    homset,
    basis,
    diagonal: list[int],
    cross: list[list[int]],
    max_degree: int,
    *,
    samples: int = 8,
) -> int:
    """Compare deterministic packed evaluations with Sage's direct norm."""
    field_degree = 2 * ctx.degree
    checked = 0
    for vector in _deterministic_sample_vectors(len(basis), samples):
        packed = _evaluate_form(diagonal, cross, vector)
        normalized = _normalize_packed_q2(
            packed, max_degree, field_degree
        )
        if normalized is None:
            raise ArithmeticError("nonzero Hom vector has zero packed norm")
        _, packed_code = normalized

        ore = ctx.modules[source].ore_polring().zero()
        for index, morphism in enumerate(basis):
            if vector & (1 << index):
                ore += morphism.ore_polynomial()
        direct = homset(ore).norm().gen()
        direct_code = sum(
            int(direct[k]) << k for k in range(int(direct.degree()) + 1)
        )
        if direct_code != packed_code:
            raise ArithmeticError(
                f"packed norm {packed_code:b} disagrees with direct Sage "
                f"norm {direct_code:b} for pair {(source, target)}"
            )
        checked += 1
    return checked


def _interpolate_form(
    ctx: SupersingularContext,
    source: int,
    target: int,
    max_degree: int,
) -> dict[str, Any]:
    """Construct and directly cross-check one complete bounded-Hom form."""
    started = perf_counter()
    homset = Hom(ctx.modules[source], ctx.modules[target])
    basis = homset.basis(degree=max_degree)
    dimension = len(basis)
    expected_dimension = 2 * (max_degree + 1) - (ctx.degree - 1)
    if dimension != expected_dimension:
        raise ArithmeticError(
            f"bounded Hom dimension {dimension}, expected "
            f"{expected_dimension} for pair {(source, target)}"
        )
    diagonal, cross = _quadratic_form_q2(
        homset, basis, max_degree, 2 * ctx.degree
    )
    checked = _direct_sage_cross_checks(
        ctx,
        source,
        target,
        homset,
        basis,
        diagonal,
        cross,
        max_degree,
    )
    return {
        "diagonal": diagonal,
        "cross": cross,
        "dimension": dimension,
        "direct_sage_norm_cross_checks": checked,
        "interpolation_seconds": perf_counter() - started,
    }


def _counter_run(
    form: dict[str, Any],
    *,
    ctx: SupersingularContext,
    max_degree: int,
    target_codes: list[int],
    target_degrees: list[int],
    early: bool,
) -> tuple[dict[str, Any], CounterResult]:
    result = run_counter(
        form["diagonal"],
        form["cross"],
        max_degree=max_degree,
        field_degree=2 * ctx.degree,
        target_codes=target_codes,
        stop_when_seen=early,
    )
    if result.invalid_norms:
        raise ArithmeticError(
            f"compiled counter found {result.invalid_norms} invalid norms"
        )
    if result.exhaustive:
        outcome = "exhaustive"
        zeros = sorted(
            code for code in target_codes if result.counts.get(code, 0) == 0
        )
        if result.seen != len(target_codes) - len(zeros):
            raise ArithmeticError("exhaustive seen count is inconsistent")
    else:
        if result.seen != len(target_codes):
            raise ArithmeticError(
                "an early counter run stopped without witnessing every target"
            )
        if any(result.counts.get(code, 0) == 0 for code in target_codes):
            raise ArithmeticError("all-target witness run contains an unseen target")
        outcome = "all_targets_witnessed"
        zeros = []

    record = {
        "target_degrees": [int(degree) for degree in target_degrees],
        "max_degree": int(max_degree),
        "dimension": int(form["dimension"]),
        "iterations": int(result.iterations),
        "target_count": int(result.target_count),
        "seen": int(result.seen),
        "invalid_norms": int(result.invalid_norms),
        "exhaustive": bool(result.exhaustive),
        "outcome": outcome,
        "zero_codes": [int(code) for code in zeros],
        "direct_sage_norm_cross_checks": int(
            form["direct_sage_norm_cross_checks"]
        ),
        "seconds": float(form["interpolation_seconds"] + result.seconds),
    }
    return record, result


def compute_compact_q2_degree7(
    ctx: SupersingularContext,
    *,
    low_exhaustive_degree: int = 10,
    progress=None,
) -> dict[str, Any]:
    """Compute a conclusive compact zero archive for one degree-seven p.

    For each unordered pair, degrees through ``low_exhaustive_degree`` are
    handled in one exhaustive bounded-Hom traversal.  Each remaining
    non-spectral degree gets its own all-target witness/exhaustion decision.
    The reverse ordered pair is transported only at the level of zero versus
    positivity via weighted Brandt symmetry.
    """
    if ctx.q != 2 or ctx.degree != 7:
        raise ValueError("the compact backend requires q=2 and degree(p)=7")
    U = theorem_bound(2, 7)
    if U != 16:
        raise ArithmeticError(f"degree-seven cutoff changed unexpectedly: {U}")
    if not 1 <= low_exhaustive_degree < U:
        raise ValueError("low exhaustive degree must lie in 1..U-1")
    if not COUNTER_SOURCE.exists() or not COUNTER_BINARY.exists():
        raise FileNotFoundError(
            "build the deterministic helper first with `make -C fast_q2`"
        )

    levels_by_degree = {
        degree: [
            ell
            for ell in monic_irreducibles(ctx.A, degree)
            if ell != ctx.p
        ]
        for degree in range(1, U)
    }
    level_by_code = {
        int(_poly_code_q2(ell)): ell
        for levels in levels_by_degree.values()
        for ell in levels
    }
    if len(level_by_code) != 4719:
        raise ArithmeticError(
            f"degree-seven scan has {len(level_by_code)} levels, expected 4719"
        )
    codes_by_degree = {
        degree: sorted(int(_poly_code_q2(ell)) for ell in levels)
        for degree, levels in levels_by_degree.items()
    }

    vertex_count = len(ctx.modules)
    if vertex_count != 43:
        raise ArithmeticError(f"degree-seven class count is {vertex_count}, expected 43")
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
                for degree in range(1, U)
                if not pair_is_spectrally_positive(
                    ctx, source, target, degree
                )
            ]
            spectral = [
                degree for degree in range(1, U) if degree not in non_spectral
            ]
            runs: list[dict[str, Any]] = []
            form_cache: dict[int, dict[str, Any]] = {}

            low_degrees = [
                degree
                for degree in non_spectral
                if degree <= low_exhaustive_degree
            ]
            if low_degrees:
                low_bound = max(low_degrees)
                form = _interpolate_form(ctx, source, target, low_bound)
                form_cache[low_bound] = form
                codes = [
                    code
                    for degree in low_degrees
                    for code in codes_by_degree[degree]
                ]
                run, counter = _counter_run(
                    form,
                    ctx=ctx,
                    max_degree=low_bound,
                    target_codes=codes,
                    target_degrees=low_degrees,
                    early=False,
                )
                runs.append(run)
                total_iterations += counter.iterations
                total_interpolation_seconds += form["interpolation_seconds"]
                total_counter_seconds += counter.seconds
                for code in run["zero_codes"]:
                    zero_pairs_by_code[code].append((source, target))

            high_degrees = [
                degree
                for degree in non_spectral
                if degree > low_exhaustive_degree
            ]
            for degree in high_degrees:
                # The unique exceptional self-pair is the only case in
                # which a degree-11 target is absent.  Giving self-pairs the
                # exact degree bound prevents one low missing target from
                # forcing a 2^26 traversal.  Other pairs share their largest
                # required form, avoiding thousands of Sage interpolations.
                if source == target and ctx.automorphism_orders[source] > 1:
                    form_bound = degree
                else:
                    form_bound = max(high_degrees)
                if form_bound not in form_cache:
                    form_cache[form_bound] = _interpolate_form(
                        ctx, source, target, form_bound
                    )
                    total_interpolation_seconds += form_cache[form_bound][
                        "interpolation_seconds"
                    ]
                form = form_cache[form_bound]
                run, counter = _counter_run(
                    form,
                    ctx=ctx,
                    max_degree=form_bound,
                    target_codes=codes_by_degree[degree],
                    target_degrees=[degree],
                    early=True,
                )
                runs.append(run)
                total_iterations += counter.iterations
                total_counter_seconds += counter.seconds
                for code in run["zero_codes"]:
                    zero_pairs_by_code[code].append((source, target))

            pair_certificates.append(
                {
                    "pair": [int(source), int(target)],
                    "spectral_degrees": [int(degree) for degree in spectral],
                    "runs": runs,
                }
            )
            if progress is not None and (
                pair_number == 1
                or pair_number == pair_total
                or pair_number % 25 == 0
            ):
                progress(
                    f"pair {pair_number}/{pair_total} ({source},{target}); "
                    f"iterations={total_iterations}"
                )

    levels = []
    degree_summary = {}
    largest_bad_degree = 0
    for degree in range(1, U):
        incomplete: list[str] = []
        for code in codes_by_degree[degree]:
            ell = level_by_code[code]
            unordered = zero_pairs_by_code[code]
            ordered = []
            for source, target in unordered:
                ordered.append([int(source), int(target)])
                if source != target:
                    ordered.append([int(target), int(source)])
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
            "complete_count": len(codes_by_degree[degree]) - len(incomplete),
            "incomplete_count": len(incomplete),
            "incomplete_levels": incomplete,
        }

    return {
        "archive_kind": COMPACT_KIND,
        "schema_version": 1,
        "q": 2,
        "p": str(ctx.p),
        "degree_p": 7,
        "theorem_bound_U": U,
        "E": largest_bad_degree + 1,
        "class_count": vertex_count,
        "weights": [int(order) for order in ctx.automorphism_orders],
        "j_invariants": [str(value) for value in ctx.j_invariants],
        "construction_checks": ctx.checks,
        "engine": {
            "name": "compiled_q2_gray_zero_archive",
            "counter_source": str(COUNTER_SOURCE.relative_to(REPOSITORY)),
            "counter_source_sha256": _sha256(COUNTER_SOURCE),
            "counter_binary_sha256": _sha256(COUNTER_BINARY),
            "compiler_contract": "g++ -O3 -std=c++20 -DNDEBUG",
            "unordered_pair_identity": (
                "w_j*b_ij(ell)=w_i*b_ji(ell); positive weights make "
                "zero/positivity orientation-invariant"
            ),
            "early_stop_semantics": (
                "partial counts are presence witnesses only and are never "
                "archived as exact Brandt entries"
            ),
            "low_exhaustive_degree": low_exhaustive_degree,
            "total_interpolation_seconds": total_interpolation_seconds,
            "total_counter_seconds": total_counter_seconds,
        },
        "pair_certificates": pair_certificates,
        "levels": levels,
        "degree_summary": degree_summary,
        "total_vectors_visited": int(total_iterations),
    }
