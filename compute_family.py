#!/usr/bin/env sage-python
"""Compute every characteristic of one fixed degree.

The optional ``--index`` argument is intended for Slurm array jobs.  Indices
are zero based and follow Sage's deterministic ordering of monic irreducibles.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
from typing import Sequence

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")

from sage.all import GF, PolynomialRing

from drinfeld_complete import (
    __version__,
    build_supersingular_context,
    compute_completeness,
    monic_irreducibles,
)


def polynomial_slug(poly) -> str:
    """Return the filename convention used by the result archive."""
    pieces = []
    for exponent in range(int(poly.degree()), -1, -1):
        coefficient = poly[exponent]
        if not coefficient:
            continue
        if coefficient != 1:
            try:
                coefficient_text = str(int(coefficient))
            except TypeError:
                coefficient_text = (
                    str(coefficient)
                    .replace(" ", "")
                    .replace("+", "p")
                    .replace("^", "")
                    .replace("*", "")
                )
            pieces.append(coefficient_text)
        if exponent == 1:
            pieces.append("T")
        elif exponent > 1:
            pieces.append(f"T{exponent}")
        elif coefficient == 1:
            pieces.append("1")
    return "_".join(pieces)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compute E(p) for every monic irreducible p of one degree."
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )
    parser.add_argument("--q", type=int, required=True)
    parser.add_argument("--degree", type=int, required=True)
    parser.add_argument(
        "--index",
        type=int,
        help="compute only this zero-based irreducible (or use SLURM_ARRAY_TASK_ID)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        help="default: results/local/q<Q>_degree<D>",
    )
    parser.add_argument("--full-matrices", action="store_true")
    parser.add_argument(
        "--force",
        action="store_true",
        help="replace a result that already exists",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    Fq = GF(args.q, name="a")
    A = PolynomialRing(Fq, "T")
    characteristics = monic_irreducibles(A, args.degree)
    index = args.index
    if index is None and "SLURM_ARRAY_TASK_ID" in os.environ:
        index = int(os.environ["SLURM_ARRAY_TASK_ID"])
    if index is not None:
        if not 0 <= index < len(characteristics):
            raise IndexError(
                f"index {index} is outside 0..{len(characteristics) - 1}"
            )
        characteristics = [characteristics[index]]

    output_dir = args.output_dir or Path(
        f"results/local/q{args.q}_degree{args.degree}"
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    for p in characteristics:
        path = output_dir / f"p_{polynomial_slug(p)}.json"
        if path.exists() and not args.force:
            print(f"exists: {path}")
            continue
        ctx = build_supersingular_context(args.q, p)
        result = compute_completeness(
            ctx,
            full_matrices=args.full_matrices,
        )
        # Publish only complete JSON; unique temporary paths also prevent
        # collisions between concurrent writers.
        temporary = path.with_suffix(path.suffix + f".{os.getpid()}.tmp")
        try:
            temporary.write_text(
                json.dumps(result, indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
            )
            temporary.replace(path)
        finally:
            temporary.unlink(missing_ok=True)
        print(
            f"q={args.q}, p={p}, U={result['theorem_bound_U']}, "
            f"E={result['E']} -> {path}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
