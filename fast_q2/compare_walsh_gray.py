#!/usr/bin/env python3
"""Benchmark the experimental Walsh counter against the existing q=2 helper.

This command refuses to build the compiled helper.  The checked-in binary must
already exist and be at least as new as its source, which keeps comparisons
read-only with respect to the primary counter.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import subprocess
import sys
from pathlib import Path
from time import perf_counter

from sage.all import GF, PolynomialRing

REPOSITORY = Path(__file__).resolve().parents[1]
if str(REPOSITORY) not in sys.path:
    sys.path.insert(0, str(REPOSITORY))

from drinfeld_complete import build_supersingular_context, monic_irreducibles
from fast_q2.gray_counter import (
    BINARY,
    SOURCE,
    CounterResult,
    _words,
    interpolate_pair,
)
from fast_q2.walsh_norm_counter import walsh_multiplicities


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _read_only_binary_state() -> tuple[int, int, str]:
    if not BINARY.exists():
        raise FileNotFoundError(
            f"{BINARY} is absent; build it separately before this comparison"
        )
    if BINARY.stat().st_mtime < SOURCE.stat().st_mtime:
        raise RuntimeError(
            f"{BINARY} is older than {SOURCE}; build it separately first"
        )
    return BINARY.stat().st_mtime_ns, BINARY.stat().st_size, _sha256(BINARY)


def run_read_only_gray_counter(
    diagonal: list[int],
    cross: list[list[int]],
    *,
    max_degree: int,
    field_degree: int,
    target_codes: list[int],
    stop_when_seen: bool = False,
) -> CounterResult:
    """Invoke the existing helper without ever calling its build function."""
    binary_before = _read_only_binary_state()
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
    for left in range(dimension):
        for right in range(left + 1, dimension):
            tokens.extend(_words(cross[left][right], word_count))
    tokens.extend(str(code) for code in target_codes)

    completed = subprocess.run(
        [str(BINARY)],
        input=" ".join(tokens),
        text=True,
        capture_output=True,
        check=True,
    )
    binary_after = _read_only_binary_state()
    if binary_after != binary_before:
        raise RuntimeError("the primary Gray-counter binary changed")

    lines = completed.stdout.splitlines()
    if not lines:
        raise RuntimeError("the primary Gray counter returned no output")
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

    walsh = walsh_multiplicities(
        diagonal,
        cross,
        max_degree=args.max_degree,
        field_degree=2 * args.degree_p,
    )

    binary_before = _read_only_binary_state()
    target_codes = list(range(1, 1 << (args.max_degree + 1)))
    gray_wall_started = perf_counter()
    gray = run_read_only_gray_counter(
        diagonal,
        cross,
        max_degree=args.max_degree,
        field_degree=2 * args.degree_p,
        target_codes=target_codes,
    )
    gray_wall_seconds = perf_counter() - gray_wall_started
    binary_after = _read_only_binary_state()
    if binary_after != binary_before:
        raise RuntimeError("the primary Gray-counter binary changed")
    if not gray.exhaustive or gray.invalid_norms:
        raise ArithmeticError("the primary Gray traversal was not valid")

    mismatches = [
        {
            "code": code,
            "walsh": walsh.count(code),
            "gray": gray.counts[code],
        }
        for code in target_codes
        if walsh.count(code) != gray.counts[code]
    ]
    if walsh.count(0):
        mismatches.insert(
            0,
            {
                "code": 0,
                "walsh": walsh.count(0),
                "gray": 0,
            },
        )

    return {
        "characteristic": str(characteristic),
        "characteristic_degree": args.degree_p,
        "characteristic_index": args.characteristic_index,
        "class_count": len(ctx.modules),
        "source": args.source,
        "target": args.target,
        "max_degree": args.max_degree,
        "hom_dimension": len(diagonal),
        "output_code_count": 1 << (args.max_degree + 1),
        "context_seconds": context_seconds,
        "interpolation_seconds": interpolation_seconds,
        "walsh": {
            "seconds": walsh.seconds,
            "character_sums_evaluated": walsh.character_sums_evaluated,
            "nonzero_character_sums": walsh.nonzero_character_sums,
            "vectors_represented": sum(walsh.counts),
        },
        "gray": {
            "reported_seconds": gray.seconds,
            "wall_seconds": gray_wall_seconds,
            "iterations": gray.iterations,
            "vectors_represented": sum(gray.counts.values()),
        },
        "all_multiplicities_match": not mismatches,
        "mismatch_count": len(mismatches),
        "first_mismatches": mismatches[:20],
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
