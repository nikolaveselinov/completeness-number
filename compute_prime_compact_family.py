#!/usr/bin/env sage-python
"""Compute exact compact archives over an odd prime constant field."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")

from sage.all import GF, PolynomialRing, is_prime

from compute_family import polynomial_slug
from drinfeld_complete.compact_prime import compute_compact_odd_prime
from drinfeld_complete.core import (
    build_supersingular_context,
    monic_irreducibles,
)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compute an odd-prime exact compact characteristic family."
    )
    parser.add_argument("--q", type=int, required=True)
    parser.add_argument("--degree", type=int, required=True)
    parser.add_argument(
        "--index",
        type=int,
        help="zero-based characteristic index (or SLURM_ARRAY_TASK_ID)",
    )
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def main():
    args = parse_args()
    if not is_prime(args.q) or args.q == 2:
        raise ValueError("--q must be an odd prime")
    if args.degree < 3:
        raise ValueError("--degree must be at least three")
    field = GF(args.q)
    ring = PolynomialRing(field, "T")
    characteristics = monic_irreducibles(ring, args.degree)
    index = args.index
    if index is None and "SLURM_ARRAY_TASK_ID" in os.environ:
        index = int(os.environ["SLURM_ARRAY_TASK_ID"])
    if index is not None:
        if not 0 <= index < len(characteristics):
            raise IndexError(
                f"index {index} outside 0..{len(characteristics) - 1}"
            )
        characteristics = [characteristics[index]]

    output_dir = args.output_dir or Path(
        f"results/local/q{args.q}_degree{args.degree}_compact"
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    for characteristic in characteristics:
        path = (
            output_dir
            / f"p_{polynomial_slug(characteristic)}.compact.json"
        )
        if path.exists() and not args.force:
            print(f"exists: {path}")
            continue
        context = build_supersingular_context(args.q, characteristic)
        result = compute_compact_odd_prime(
            context,
            progress=lambda message: print(
                f"p={characteristic}: {message}", flush=True
            ),
        )
        temporary = path.with_suffix(path.suffix + ".tmp")
        temporary.write_text(
            json.dumps(result, indent=2, sort_keys=True) + "\n"
        )
        temporary.replace(path)
        print(
            f"q={args.q}, p={characteristic}, "
            f"U={result['theorem_bound_U']}, E={result['E']} -> {path}",
            flush=True,
        )


if __name__ == "__main__":
    main()
