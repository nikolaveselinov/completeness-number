#!/usr/bin/env python3
"""Cross-check the isolated C++ Walsh counter against both exact backends."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from time import perf_counter

from sage.all import GF, PolynomialRing

REPOSITORY = Path(__file__).resolve().parents[1]
if str(REPOSITORY) not in sys.path:
    sys.path.insert(0, str(REPOSITORY))

from drinfeld_complete import build_supersingular_context, monic_irreducibles
from fast_q2.compare_walsh_gray import (
    _read_only_binary_state,
    run_read_only_gray_counter,
)
from fast_q2.gray_counter import interpolate_pair
from fast_q2.q2_walsh_cpp import run_walsh_counter
from fast_q2.walsh_norm_counter import walsh_multiplicities


def compare(args: argparse.Namespace) -> dict:
    A = PolynomialRing(GF(2), "T")
    characteristics = monic_irreducibles(A, args.degree_p)
    if not 0 <= args.characteristic_index < len(characteristics):
        raise ValueError(
            "characteristic index must satisfy "
            f"0 <= index < {len(characteristics)}"
        )
    characteristic = characteristics[args.characteristic_index]

    context_started = perf_counter()
    ctx = build_supersingular_context(2, characteristic)
    context_seconds = perf_counter() - context_started
    if not 0 <= args.source < len(ctx.modules):
        raise ValueError("source index is outside the supersingular class set")
    if not 0 <= args.target < len(ctx.modules):
        raise ValueError("target index is outside the supersingular class set")

    interpolation_started = perf_counter()
    diagonal, cross = interpolate_pair(
        ctx, args.source, args.target, args.max_degree
    )
    interpolation_seconds = perf_counter() - interpolation_started
    field_degree = 2 * args.degree_p

    python_result = walsh_multiplicities(
        diagonal,
        cross,
        max_degree=args.max_degree,
        field_degree=field_degree,
    )
    cpp_wall_started = perf_counter()
    cpp_result = run_walsh_counter(
        diagonal,
        cross,
        max_degree=args.max_degree,
        field_degree=field_degree,
    )
    cpp_wall_seconds = perf_counter() - cpp_wall_started

    binary_before = _read_only_binary_state()
    target_codes = list(range(1, 1 << (args.max_degree + 1)))
    gray_wall_started = perf_counter()
    gray_result = run_read_only_gray_counter(
        diagonal,
        cross,
        max_degree=args.max_degree,
        field_degree=field_degree,
        target_codes=target_codes,
    )
    gray_wall_seconds = perf_counter() - gray_wall_started
    binary_after = _read_only_binary_state()
    if binary_after != binary_before:
        raise RuntimeError("the primary Gray-counter binary changed")
    if not gray_result.exhaustive or gray_result.invalid_norms:
        raise ArithmeticError("the primary Gray traversal was not valid")

    python_cpp_mismatches = [
        code
        for code, (python_count, cpp_count) in enumerate(
            zip(python_result.counts, cpp_result.counts)
        )
        if python_count != cpp_count
    ]
    cpp_gray_mismatches = [
        code
        for code in target_codes
        if cpp_result.count(code) != gray_result.counts[code]
    ]
    if cpp_result.count(0):
        cpp_gray_mismatches.insert(0, 0)

    return {
        "characteristic": str(characteristic),
        "characteristic_degree": args.degree_p,
        "characteristic_index": args.characteristic_index,
        "class_count": len(ctx.modules),
        "source": args.source,
        "target": args.target,
        "max_degree": args.max_degree,
        "hom_dimension": len(diagonal),
        "output_code_count": len(cpp_result.counts),
        "context_seconds": context_seconds,
        "interpolation_seconds": interpolation_seconds,
        "python_walsh": {
            "seconds": python_result.seconds,
            "nonzero_character_sums":
                python_result.nonzero_character_sums,
        },
        "cpp_walsh": {
            "reported_seconds": cpp_result.seconds,
            "wall_seconds": cpp_wall_seconds,
            "nonzero_character_sums": cpp_result.nonzero_character_sums,
        },
        "gray": {
            "reported_seconds": gray_result.seconds,
            "wall_seconds": gray_wall_seconds,
            "iterations": gray_result.iterations,
        },
        "python_cpp_match": not python_cpp_mismatches,
        "cpp_gray_match": not cpp_gray_mismatches,
        "python_cpp_mismatch_count": len(python_cpp_mismatches),
        "cpp_gray_mismatch_count": len(cpp_gray_mismatches),
        "first_python_cpp_mismatches": python_cpp_mismatches[:20],
        "first_cpp_gray_mismatches": cpp_gray_mismatches[:20],
        "vectors_represented": sum(cpp_result.counts),
        "gray_binary_sha256": binary_before[2],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--degree-p", type=int, required=True)
    parser.add_argument("--characteristic-index", type=int, default=0)
    parser.add_argument("--source", type=int, default=0)
    parser.add_argument("--target", type=int, default=1)
    parser.add_argument("--max-degree", type=int, required=True)
    args = parser.parse_args()
    print(json.dumps(compare(args), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
