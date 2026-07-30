#!/usr/bin/env sage-python
"""Compute one or all exact descending q=2 E(p) certificates."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")

from sage.all import GF, PolynomialRing

from compute_family import polynomial_slug
from drinfeld_complete import build_supersingular_context, monic_irreducibles
from drinfeld_complete.e_only_q2 import compute_e_only_q2


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--degree", type=int, required=True)
    parser.add_argument(
        "--index",
        type=int,
        help="zero-based irreducible index (defaults to SLURM_ARRAY_TASK_ID)",
    )
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--checkpoint-every", type=int, default=250)
    return parser.parse_args()


def atomic_write(path: Path, record: dict) -> None:
    temporary = path.with_name(path.name + f".{os.getpid()}.tmp")
    temporary.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def main() -> None:
    args = parse_args()
    if args.degree < 3:
        raise ValueError("--degree must be at least 3")
    A = PolynomialRing(GF(2), "T")
    characteristics = monic_irreducibles(A, args.degree)
    index = args.index
    if index is None and "SLURM_ARRAY_TASK_ID" in os.environ:
        index = int(os.environ["SLURM_ARRAY_TASK_ID"])
    if index is not None:
        if not 0 <= index < len(characteristics):
            raise IndexError(
                f"index {index} is outside 0..{len(characteristics) - 1}"
            )
        selected = [(index, characteristics[index])]
    else:
        selected = list(enumerate(characteristics))

    output_dir = args.output_dir or Path(
        f"results/local/q2_degree{args.degree}_e_only"
    )
    checkpoint_dir = output_dir / ".checkpoints"
    output_dir.mkdir(parents=True, exist_ok=True)
    checkpoint_dir.mkdir(parents=True, exist_ok=True)
    for original_index, p in selected:
        slug = polynomial_slug(p)
        output = output_dir / f"p_{slug}.e-only.json"
        if output.exists():
            print(f"exists: {output}", flush=True)
            continue
        checkpoint = checkpoint_dir / f"p_{slug}.checkpoint.json"
        print(
            f"constructing q=2, degree={args.degree}, index={original_index}, "
            f"p={p}",
            flush=True,
        )
        ctx = build_supersingular_context(2, p)
        record = compute_e_only_q2(
            ctx,
            checkpoint_path=checkpoint,
            checkpoint_every=args.checkpoint_every,
            progress=lambda message: print(f"{p}: {message}", flush=True),
        )
        atomic_write(output, record)
        print(
            f"completed q=2, p={p}, U={record['theorem_bound_U']}, "
            f"E={record['E']} -> {output}",
            flush=True,
        )


if __name__ == "__main__":
    main()
