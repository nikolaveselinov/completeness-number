#!/usr/bin/env python3
"""Standalone driver for the exact compiled binary reduced-norm counter."""

from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from time import perf_counter

from sage.all import GF, Hom, PolynomialRing

REPOSITORY = Path(__file__).resolve().parents[1]
if str(REPOSITORY) not in sys.path:
    sys.path.insert(0, str(REPOSITORY))

from drinfeld_complete import build_supersingular_context, monic_irreducibles
from drinfeld_complete.core import (
    _poly_code_q2,
    _quadratic_form_q2,
    enumerate_pair_norms_q2,
)


HERE = Path(__file__).resolve().parent
SOURCE = HERE / "q2_norm_counter.cpp"
BINARY = HERE / "q2_norm_counter"


@dataclass(frozen=True)
class CounterResult:
    iterations: int
    exhaustive: bool
    seen: int
    target_count: int
    invalid_norms: int
    seconds: float
    counts: dict[int, int]


def compile_counter() -> None:
    """Build the helper if absent or older than its source."""
    if BINARY.exists() and BINARY.stat().st_mtime >= SOURCE.stat().st_mtime:
        return
    subprocess.run(
        [
            "g++",
            "-O3",
            "-std=c++20",
            "-DNDEBUG",
            str(SOURCE),
            "-o",
            str(BINARY),
        ],
        check=True,
    )


def _words(value: int, count: int) -> list[str]:
    mask = (1 << 64) - 1
    return [hex((int(value) >> (64 * k)) & mask) for k in range(count)]


def run_counter(
    diagonal: list[int],
    cross: list[list[int]],
    *,
    max_degree: int,
    field_degree: int,
    target_codes: list[int],
    stop_when_seen: bool,
) -> CounterResult:
    """Run the compiled counter on one already-interpolated norm form."""
    compile_counter()
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
        tokens.extend(_words(value, word_count))
    for i in range(dimension):
        for j in range(i + 1, dimension):
            tokens.extend(_words(cross[i][j], word_count))
    tokens.extend(str(code) for code in target_codes)
    completed = subprocess.run(
        [str(BINARY)],
        input=" ".join(tokens),
        text=True,
        capture_output=True,
        check=True,
    )
    lines = completed.stdout.splitlines()
    header = lines[0].split()
    if header[0] != "SUMMARY" or len(header) != 7:
        raise RuntimeError(f"unexpected helper output: {lines[0]!r}")
    counts = {}
    for line in lines[1:]:
        label, code, count = line.split()
        if label != "COUNT":
            raise RuntimeError(f"unexpected helper output: {line!r}")
        counts[int(code)] = int(count)
    return CounterResult(
        iterations=int(header[1]),
        exhaustive=bool(int(header[2])),
        seen=int(header[3]),
        target_count=int(header[4]),
        invalid_norms=int(header[5]),
        seconds=float(header[6]),
        counts=counts,
    )


def interpolate_pair(ctx, source: int, target: int, max_degree: int):
    homset = Hom(ctx.modules[source], ctx.modules[target])
    basis = homset.basis(degree=max_degree)
    expected = 2 * (max_degree + 1) - (ctx.degree - 1)
    if len(basis) != expected:
        raise ArithmeticError(
            f"bounded Hom dimension {len(basis)}, expected {expected}"
        )
    diagonal, cross = _quadratic_form_q2(
        homset, basis, max_degree, 2 * ctx.degree
    )
    return diagonal, cross


def parse_context(p_text: str):
    Fq = GF(2)
    A = PolynomialRing(Fq, "T")
    T = A.gen()
    p = A(eval(p_text.replace("^", "**"), {"T": T}))
    if not p.is_irreducible():
        raise ValueError(f"{p} is not irreducible")
    return build_supersingular_context(2, p)


def target_codes(ctx, minimum_degree: int, maximum_degree: int) -> list[int]:
    return [
        _poly_code_q2(ell)
        for degree in range(minimum_degree, maximum_degree + 1)
        for ell in monic_irreducibles(ctx.A, degree)
        if ell != ctx.p
    ]


def benchmark(args) -> dict:
    ctx = parse_context(args.p)
    codes = target_codes(ctx, args.min_target_degree, args.max_degree)
    interpolation_started = perf_counter()
    diagonal, cross = interpolate_pair(
        ctx, args.source, args.target, args.max_degree
    )
    interpolation_seconds = perf_counter() - interpolation_started
    fast = run_counter(
        diagonal,
        cross,
        max_degree=args.max_degree,
        field_degree=2 * ctx.degree,
        target_codes=codes,
        stop_when_seen=args.early,
    )
    report = {
        "p": str(ctx.p),
        "degree_p": ctx.degree,
        "class_count": len(ctx.modules),
        "source": args.source,
        "target": args.target,
        "max_degree": args.max_degree,
        "min_target_degree": args.min_target_degree,
        "dimension": len(diagonal),
        "target_count": len(codes),
        "interpolation_seconds": interpolation_seconds,
        "counter": {
            "iterations": fast.iterations,
            "exhaustive": fast.exhaustive,
            "seen": fast.seen,
            "seconds": fast.seconds,
            "million_vectors_per_second": (
                fast.iterations / fast.seconds / 1_000_000
            ),
        },
        "zero_codes": sorted(
            code for code, count in fast.counts.items() if count == 0
        ),
    }

    if args.compare_python:
        if args.early:
            raise ValueError("--compare-python requires a full traversal")
        python_started = perf_counter()
        expected, metadata = enumerate_pair_norms_q2(
            ctx, args.source, args.target, args.max_degree
        )
        python_seconds = perf_counter() - python_started
        discrepancies = {
            code: [expected.get(code, 0), fast.counts.get(code, 0)]
            for code in codes
            if expected.get(code, 0) != fast.counts.get(code, 0)
        }
        report["python_comparison"] = {
            "seconds": python_seconds,
            "metadata_seconds": metadata["seconds"],
            "discrepancies": discrepancies,
            "passed": not discrepancies,
        }
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p", required=True, help="irreducible characteristic")
    parser.add_argument("--source", type=int, default=0)
    parser.add_argument("--target", type=int, default=0)
    parser.add_argument("--max-degree", type=int, required=True)
    parser.add_argument("--min-target-degree", type=int, default=1)
    parser.add_argument(
        "--early",
        action="store_true",
        help="stop once every requested irreducible norm has a witness",
    )
    parser.add_argument(
        "--compare-python",
        action="store_true",
        help="compare target multiplicities with direct Sage enumeration",
    )
    args = parser.parse_args()
    print(json.dumps(benchmark(args), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
