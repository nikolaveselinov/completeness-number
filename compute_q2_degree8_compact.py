#!/usr/bin/env sage-python
"""Compute one or all compact exact q=2, degree-eight archives."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")

from sage.all import GF, PolynomialRing

from compute_family import polynomial_slug
from drinfeld_complete import (
    build_supersingular_context,
    compute_compact_q2_degree8,
    monic_irreducibles,
)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--index",
        type=int,
        help="zero-based irreducible index (defaults to SLURM_ARRAY_TASK_ID)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/local/q2_degree8"),
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="replace an existing compact archive",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    A = PolynomialRing(GF(2), "T")
    characteristics = monic_irreducibles(A, 8)
    if len(characteristics) != 30:
        raise ArithmeticError(
            f"expected 30 degree-eight characteristics, got {len(characteristics)}"
        )
    index = args.index
    if index is None and "SLURM_ARRAY_TASK_ID" in os.environ:
        index = int(os.environ["SLURM_ARRAY_TASK_ID"])
    if index is not None:
        if not 0 <= index < len(characteristics):
            raise IndexError(f"index {index} is outside 0..29")
        characteristics = [characteristics[index]]

    args.output_dir.mkdir(parents=True, exist_ok=True)
    for p in characteristics:
        path = args.output_dir / f"p_{polynomial_slug(p)}.compact.json"
        if path.exists() and not args.force:
            print(f"exists: {path}", flush=True)
            continue
        print(f"constructing q=2, p={p}", flush=True)
        ctx = build_supersingular_context(2, p)
        result = compute_compact_q2_degree8(
            ctx,
            progress=lambda message: print(f"{p}: {message}", flush=True),
        )
        temporary = path.with_suffix(path.suffix + ".tmp")
        temporary.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
        temporary.replace(path)
        print(
            f"completed q=2, p={p}, U={result['theorem_bound_U']}, "
            f"E={result['E']} -> {path}",
            flush=True,
        )


if __name__ == "__main__":
    main()
