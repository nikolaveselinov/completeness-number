#!/usr/bin/env python3
"""Parallel exact unordered-pair scan using the Sage-free q=2 backend."""

from __future__ import annotations

import argparse
import json
import multiprocessing
import os
from pathlib import Path
from time import perf_counter

from d9_j0_verify import (
    QuadraticTower,
    cube_root_representatives,
    hom_basis,
    interpolate_norm_form,
    monic_irreducibles_binary,
    polynomial_text,
    run_counter,
    sha256,
    supersingular_roots,
)


WORKER_CONTEXT: dict[str, object] = {}


def scan_pair(pair: tuple[int, int]) -> dict[str, object]:
    field = WORKER_CONTEXT["field"]
    roots = WORKER_CONTEXT["roots"]
    module_parameters = WORKER_CONTEXT["module_parameters"]
    target_codes = WORKER_CONTEXT["target_codes"]
    level_degree = WORKER_CONTEXT["level_degree"]
    counter = WORKER_CONTEXT["counter"]
    if not isinstance(field, QuadraticTower):
        raise TypeError("worker field was not initialized")
    source, target = pair
    source_g = module_parameters[source]
    target_g = module_parameters[target]
    started = perf_counter()
    basis = hom_basis(
        field,
        level_degree,
        source_g=source_g,
        target_g=target_g,
    )
    diagonal, cross = interpolate_norm_form(
        field,
        basis,
        level_degree,
        source_g=source_g,
    )
    result = run_counter(
        counter,
        diagonal,
        cross,
        max_degree=level_degree,
        field_degree=field.field_degree,
        target_codes=target_codes,
        stop_when_seen=True,
    )
    counts = result.pop("counts")
    zero_codes = sorted(
        code for code in target_codes if counts.get(code, 0) == 0
    )
    return {
        "pair": [source, target],
        "source_j_code": roots[source],
        "target_j_code": roots[target],
        "source_g_code": source_g,
        "target_g_code": target_g,
        "hom_dimension": len(basis),
        **result,
        "zero_codes": zero_codes,
        "zero_levels": [polynomial_text(code) for code in zero_codes],
        "seconds": perf_counter() - started,
    }


def atomic_json(path: Path, value: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".{os.getpid()}.tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--degree-p", type=int, required=True)
    parser.add_argument(
        "--p-code", type=lambda value: int(value, 0), required=True
    )
    parser.add_argument("--level-degree", type=int, required=True)
    parser.add_argument("--counter", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--output", type=Path, required=True)
    scope = parser.add_mutually_exclusive_group()
    scope.add_argument(
        "--off-diagonal-only",
        action="store_true",
        help="skip diagonal pairs already certified separately",
    )
    scope.add_argument(
        "--diagonal-only",
        action="store_true",
        help="scan only the diagonal pairs",
    )
    parser.add_argument(
        "--max-pairs",
        type=int,
        help="development-only prefix limit; omitted for a conclusive scan",
    )
    args = parser.parse_args()
    if args.workers < 1:
        raise ValueError("--workers must be positive")

    started = perf_counter()
    field = QuadraticTower(args.degree_p, args.p_code)
    roots = supersingular_roots(field)
    cube_roots = cube_root_representatives(field, roots)
    module_parameters = [cube_roots[j] for j in roots]
    target_codes = monic_irreducibles_binary(args.level_degree)
    if args.diagonal_only:
        pairs = [(index, index) for index in range(len(roots))]
    else:
        pairs = [
            (source, target)
            for source in range(len(roots))
            for target in range(
                source + int(args.off_diagonal_only), len(roots)
            )
        ]
    full_pair_count = len(pairs)
    if args.max_pairs is not None:
        pairs = pairs[: args.max_pairs]

    WORKER_CONTEXT.update(
        {
            "field": field,
            "roots": roots,
            "module_parameters": module_parameters,
            "target_codes": target_codes,
            "level_degree": args.level_degree,
            "counter": args.counter.resolve(),
        }
    )

    receipts: list[dict[str, object]] = []
    obstruction: dict[str, object] | None = None
    worker_count = min(args.workers, len(pairs))
    context = multiprocessing.get_context("fork")
    with context.Pool(worker_count) as pool:
        iterator = pool.imap_unordered(scan_pair, pairs, chunksize=1)
        for completed, receipt in enumerate(iterator, start=1):
            if (
                not receipt["exhaustive"]
                and receipt["seen"] != receipt["target_count"]
            ):
                raise ArithmeticError(
                    f"pair {receipt['pair']} stopped without all targets"
                )
            if receipt["invalid_norms"]:
                raise ArithmeticError(
                    f"pair {receipt['pair']} had invalid norms"
                )
            receipts.append(receipt)
            if receipt["zero_codes"]:
                obstruction = receipt
                pool.terminate()
                break
            if completed == 1 or completed % 100 == 0:
                elapsed = perf_counter() - started
                rate = completed / elapsed
                print(
                    f"pairs {completed}/{len(pairs)}; "
                    f"{rate:.2f} pairs/s; elapsed={elapsed:.1f}s",
                    flush=True,
                )

    receipts.sort(key=lambda receipt: receipt["pair"])
    prefix_limited = args.max_pairs is not None
    conclusive_complete = (
        obstruction is None
        and not prefix_limited
        and len(receipts) == full_pair_count
    )
    report: dict[str, object] = {
        "schema": "sage_free_q2_pair_scan_v1",
        "method": "direct_hom_nullspace_and_motive_determinant",
        "q": 2,
        "p": polynomial_text(args.p_code),
        "p_code": args.p_code,
        "degree_p": args.degree_p,
        "level_degree": args.level_degree,
        "class_count": len(roots),
        "target_count": len(target_codes),
        "off_diagonal_only": args.off_diagonal_only,
        "diagonal_only": args.diagonal_only,
        "expected_pair_count": full_pair_count,
        "requested_pair_count": len(pairs),
        "completed_pair_count": len(receipts),
        "prefix_limited": prefix_limited,
        "complete_for_scanned_scope": conclusive_complete,
        "first_obstruction": obstruction,
        "counter_binary": str(args.counter.resolve()),
        "counter_binary_sha256": sha256(args.counter.resolve()),
        "workers": worker_count,
        "total_iterations": sum(
            int(receipt["iterations"]) for receipt in receipts
        ),
        "total_counter_seconds": sum(
            float(receipt["counter_seconds"]) for receipt in receipts
        ),
        "total_pair_seconds": sum(
            float(receipt["seconds"]) for receipt in receipts
        ),
        "wall_seconds": perf_counter() - started,
        "pair_receipts": receipts,
    }
    atomic_json(args.output, report)
    print(
        json.dumps(
            {
                key: report[key]
                for key in (
                    "p",
                    "level_degree",
                    "class_count",
                    "target_count",
                    "completed_pair_count",
                    "complete_for_scanned_scope",
                    "first_obstruction",
                    "wall_seconds",
                )
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
