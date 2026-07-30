#!/usr/bin/env sage-python
"""Command-line driver for exact Drinfeld completeness computations."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
from typing import Sequence

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")

from sage.all import GF, PolynomialRing, sage_eval

from drinfeld_complete import (
    __version__,
    build_supersingular_context,
    compute_completeness,
)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compute the exact completeness number E(p).",
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )
    parser.add_argument("--q", type=int, required=True, help="base-field size")
    parser.add_argument(
        "--p",
        required=True,
        help="monic irreducible characteristic, e.g. 'T^3 + T + 1'",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="JSON output path (stdout if omitted)",
    )
    parser.add_argument(
        "--full-matrices",
        action="store_true",
        help="enumerate every entry through U-1 instead of using pair certificates",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    Fq = GF(args.q, name="a")
    A = PolynomialRing(Fq, "T")
    T = A.gen()
    p = A(
        sage_eval(
            args.p.replace("^", "**"),
            locals={"T": T, "a": Fq.gen()},
        )
    )
    ctx = build_supersingular_context(args.q, p)
    result = compute_completeness(ctx, full_matrices=args.full_matrices)
    payload = json.dumps(result, indent=2, sort_keys=True)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        temporary = args.output.with_suffix(
            args.output.suffix + f".{os.getpid()}.tmp"
        )
        try:
            temporary.write_text(payload + "\n", encoding="utf-8")
            temporary.replace(args.output)
        finally:
            temporary.unlink(missing_ok=True)
        print(
            f"q={args.q}, p={p}, U={result['theorem_bound_U']}, "
            f"E={result['E']} -> {args.output}"
        )
    else:
        print(payload)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
