"""Compact exact-zero computation for q=2, degree-eight characteristics."""

from __future__ import annotations

from typing import Any

from .compact_q2 import (
    COUNTER_BINARY,
    COUNTER_SOURCE,
    REPOSITORY,
    _counter_run,
    _interpolate_form,
    _sha256,
)
from .core import (
    SupersingularContext,
    _poly_code_q2,
    monic_irreducibles,
    pair_is_spectrally_positive,
    theorem_bound,
)


COMPACT_KIND_DEGREE8 = "q2_degree8_compact_zero_v1"


def compute_compact_q2_degree8(
    ctx: SupersingularContext,
    *,
    low_exhaustive_degree: int = 11,
    progress=None,
) -> dict[str, Any]:
    """Compute a conclusive compact archive for one degree-eight p.

    Every unordered pair has two exact runs:

    * a full traversal through degree 11, deciding all levels of degrees
      1--11, and
    * one degree-14 bounded traversal targeting all levels of degrees 12--14.
      It either witnesses every target or exhausts the whole Hom space and
      thereby proves every missing target is zero.

    The degree-eight family has no pair-specific spectral positivity below
    the global theorem cutoff, and all vertex weights are one.
    """
    if ctx.q != 2 or ctx.degree != 8:
        raise ValueError("the degree-eight compact backend requires q=2, deg(p)=8")
    U = theorem_bound(2, 8)
    if U != 15:
        raise ArithmeticError(f"degree-eight cutoff changed unexpectedly: {U}")
    if low_exhaustive_degree != 11:
        raise ValueError("the validated degree-eight exhaustive cutoff is 11")
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
    codes_by_degree = {
        degree: sorted(int(_poly_code_q2(ell)) for ell in levels)
        for degree, levels in levels_by_degree.items()
    }
    level_by_code = {
        int(_poly_code_q2(ell)): ell
        for levels in levels_by_degree.values()
        for ell in levels
    }
    if len(level_by_code) != 2537:
        raise ArithmeticError(
            f"degree-eight scan has {len(level_by_code)} levels, expected 2537"
        )
    if len(ctx.modules) != 85:
        raise ArithmeticError(
            f"degree-eight class count is {len(ctx.modules)}, expected 85"
        )
    if any(int(order) != 1 for order in ctx.automorphism_orders):
        raise ArithmeticError("degree-eight compact backend expected unit weights")
    if any(
        pair_is_spectrally_positive(ctx, 0, 0, degree)
        for degree in range(1, U)
    ):
        raise ArithmeticError("unexpected sub-cutoff spectral degree")

    vertex_count = len(ctx.modules)
    pair_total = vertex_count * (vertex_count + 1) // 2
    low_degrees = list(range(1, low_exhaustive_degree + 1))
    high_degrees = list(range(low_exhaustive_degree + 1, U))
    low_codes = [
        code for degree in low_degrees for code in codes_by_degree[degree]
    ]
    high_codes = [
        code for degree in high_degrees for code in codes_by_degree[degree]
    ]
    if len(low_codes) != 411 or len(high_codes) != 2126:
        raise ArithmeticError("degree-eight target partition regression")

    zero_pairs_by_code: dict[int, list[tuple[int, int]]] = {
        code: [] for code in level_by_code
    }
    pair_certificates = []
    total_iterations = 0
    total_interpolation_seconds = 0.0
    total_counter_seconds = 0.0
    pair_number = 0

    for source in range(vertex_count):
        for target in range(source, vertex_count):
            pair_number += 1
            low_form = _interpolate_form(
                ctx, source, target, low_exhaustive_degree
            )
            low_run, low_counter = _counter_run(
                low_form,
                ctx=ctx,
                max_degree=low_exhaustive_degree,
                target_codes=low_codes,
                target_degrees=low_degrees,
                early=False,
            )
            high_form = _interpolate_form(ctx, source, target, U - 1)
            high_run, high_counter = _counter_run(
                high_form,
                ctx=ctx,
                max_degree=U - 1,
                target_codes=high_codes,
                target_degrees=high_degrees,
                early=True,
            )
            runs = [low_run, high_run]
            for run in runs:
                for code in run["zero_codes"]:
                    zero_pairs_by_code[code].append((source, target))
            pair_certificates.append(
                {
                    "pair": [int(source), int(target)],
                    "spectral_degrees": [],
                    "runs": runs,
                }
            )
            total_iterations += low_counter.iterations + high_counter.iterations
            total_interpolation_seconds += (
                low_form["interpolation_seconds"]
                + high_form["interpolation_seconds"]
            )
            total_counter_seconds += (
                low_counter.seconds + high_counter.seconds
            )
            if progress is not None and (
                pair_number == 1
                or pair_number == pair_total
                or pair_number % 50 == 0
            ):
                progress(
                    f"pair {pair_number}/{pair_total} ({source},{target}); "
                    f"iterations={total_iterations}"
                )

    levels = []
    degree_summary = {}
    largest_bad_degree = 0
    for degree in range(1, U):
        incomplete = []
        for code in codes_by_degree[degree]:
            ell = level_by_code[code]
            ordered = []
            for source, target in zero_pairs_by_code[code]:
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
        "archive_kind": COMPACT_KIND_DEGREE8,
        "schema_version": 1,
        "q": 2,
        "p": str(ctx.p),
        "degree_p": 8,
        "theorem_bound_U": U,
        "E": largest_bad_degree + 1,
        "class_count": vertex_count,
        "weights": [int(order) for order in ctx.automorphism_orders],
        "j_invariants": [str(value) for value in ctx.j_invariants],
        "construction_checks": ctx.checks,
        "engine": {
            "name": "compiled_q2_degree8_gray_zero_archive",
            "counter_source": str(COUNTER_SOURCE.relative_to(REPOSITORY)),
            "counter_source_sha256": _sha256(COUNTER_SOURCE),
            "counter_binary_sha256": _sha256(COUNTER_BINARY),
            "compiler_contract": "g++ -O3 -std=c++20 -DNDEBUG",
            "unordered_pair_identity": (
                "all weights are one, so b_ij(ell)=b_ji(ell); "
                "zero/positivity is orientation-invariant"
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
