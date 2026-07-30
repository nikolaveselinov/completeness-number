"""Exact descending computation of only the completeness number for q=2.

The full compact archives decide every pair at every level below the theorem
cutoff.  To determine E(p), considerably less is required: prove every degree
above one obstruction positive, then exhibit one exhaustively certified zero
at the obstruction degree.  This module records precisely that certificate.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
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


E_ONLY_KIND = "q2_e_only_v1"
CHECKPOINT_KIND = "q2_e_only_checkpoint_v1"


def e_only_pair_order(
    ctx: SupersingularContext, level_degree: int
) -> list[tuple[int, int]]:
    """Return every non-spectral pair in deterministic search order.

    Diagonal pairs are tried first because they are frequent extremal
    witnesses in the known data.  The list remains exhaustive: if no
    diagonal obstruction exists, every off-diagonal pair is still checked.
    """
    pairs = [
        (source, target)
        for source in range(len(ctx.modules))
        for target in range(source, len(ctx.modules))
        if not pair_is_spectrally_positive(
            ctx, source, target, level_degree
        )
    ]
    pairs.sort(key=lambda pair: (pair[0] != pair[1], pair[0], pair[1]))
    return pairs


def _atomic_json(path: Path, value: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".{os.getpid()}.tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def _level_data(
    ctx: SupersingularContext, degree: int
) -> tuple[list[Any], list[int], dict[int, Any]]:
    levels = [
        ell
        for ell in monic_irreducibles(ctx.A, degree)
        if ell != ctx.p
    ]
    levels.sort(key=lambda ell: int(_poly_code_q2(ell)))
    codes = [int(_poly_code_q2(ell)) for ell in levels]
    return levels, codes, dict(zip(codes, levels))


def _new_checkpoint(ctx: SupersingularContext, cutoff: int) -> dict[str, Any]:
    return {
        "checkpoint_kind": CHECKPOINT_KIND,
        "schema_version": 1,
        "q": 2,
        "p": str(ctx.p),
        "degree_p": int(ctx.degree),
        "theorem_bound_U": int(cutoff),
        "counter_source_sha256": _sha256(COUNTER_SOURCE),
        "counter_binary_sha256": _sha256(COUNTER_BINARY),
        "positive_degrees": [],
        "current_degree": None,
        "current_pair_runs": [],
        "witness": None,
        "total_vectors_visited": 0,
        "total_interpolation_seconds": 0.0,
        "total_counter_seconds": 0.0,
    }


def _expected_next_degree(state: dict[str, Any]) -> int:
    return int(state["theorem_bound_U"]) - 1 - len(
        state["positive_degrees"]
    )


def _check_positive_run(run: dict[str, Any], degree: int) -> None:
    if run.get("target_degrees") != [degree]:
        raise ValueError("checkpoint run targets the wrong degree")
    if run.get("zero_codes"):
        raise ValueError("checkpoint positivity run contains a zero")
    if run.get("outcome") not in {"all_targets_witnessed", "exhaustive"}:
        raise ValueError("checkpoint contains an unsupported run outcome")


def _load_checkpoint(
    path: Path, ctx: SupersingularContext, cutoff: int
) -> dict[str, Any]:
    if not path.exists():
        return _new_checkpoint(ctx, cutoff)
    state = json.loads(path.read_text())
    expected_header = {
        "checkpoint_kind": CHECKPOINT_KIND,
        "schema_version": 1,
        "q": 2,
        "p": str(ctx.p),
        "degree_p": int(ctx.degree),
        "theorem_bound_U": int(cutoff),
        "counter_source_sha256": _sha256(COUNTER_SOURCE),
        "counter_binary_sha256": _sha256(COUNTER_BINARY),
    }
    for key, expected in expected_header.items():
        if state.get(key) != expected:
            raise ValueError(
                f"checkpoint {key}={state.get(key)!r}, expected {expected!r}"
            )
    positive = state.get("positive_degrees")
    if not isinstance(positive, list):
        raise ValueError("checkpoint positive_degrees must be a list")
    expected_degrees = list(
        range(cutoff - 1, cutoff - 1 - len(positive), -1)
    )
    if [item.get("degree") for item in positive] != expected_degrees:
        raise ValueError("checkpoint positive degrees are not a descending prefix")
    for item in positive:
        degree = int(item["degree"])
        expected_pairs = e_only_pair_order(ctx, degree)
        runs = item.get("pair_runs")
        if not isinstance(runs, list):
            raise ValueError("checkpoint pair_runs must be a list")
        archived_pairs = [tuple(entry["pair"]) for entry in runs]
        if archived_pairs != expected_pairs:
            raise ValueError(
                f"checkpoint degree {degree} does not cover every non-spectral pair"
            )
        for entry in runs:
            _check_positive_run(entry["run"], degree)
    current_degree = state.get("current_degree")
    current_runs = state.get("current_pair_runs")
    if current_degree is None:
        if current_runs != []:
            raise ValueError("checkpoint has runs but no current degree")
    else:
        if int(current_degree) != _expected_next_degree(state):
            raise ValueError("checkpoint current degree is not the next degree")
        if not isinstance(current_runs, list):
            raise ValueError("checkpoint current_pair_runs must be a list")
        expected_pairs = e_only_pair_order(ctx, int(current_degree))
        archived_pairs = [tuple(entry["pair"]) for entry in current_runs]
        if archived_pairs != expected_pairs[: len(archived_pairs)]:
            raise ValueError("checkpoint current pairs are not a canonical prefix")
        for entry in current_runs:
            _check_positive_run(entry["run"], int(current_degree))
    return state


def _final_record(
    ctx: SupersingularContext,
    state: dict[str, Any],
) -> dict[str, Any]:
    witness = state["witness"]
    if not isinstance(witness, dict):
        raise ValueError("cannot finalize an E-only archive without a witness")
    return {
        "archive_kind": E_ONLY_KIND,
        "schema_version": 1,
        "q": 2,
        "p": str(ctx.p),
        "degree_p": int(ctx.degree),
        "theorem_bound_U": int(state["theorem_bound_U"]),
        "E": int(witness["degree"]) + 1,
        "class_count": len(ctx.modules),
        "weights": [int(order) for order in ctx.automorphism_orders],
        "j_invariants": [str(value) for value in ctx.j_invariants],
        "construction_checks": ctx.checks,
        "positive_degrees": state["positive_degrees"],
        "witness": witness,
        "total_vectors_visited": int(state["total_vectors_visited"]),
        "engine": {
            "name": "compiled_q2_descending_e_only",
            "counter_source": str(COUNTER_SOURCE.relative_to(REPOSITORY)),
            "counter_source_sha256": state["counter_source_sha256"],
            "counter_binary_sha256": state["counter_binary_sha256"],
            "compiler_contract": "g++ -O3 -std=c++20 -DNDEBUG",
            "pair_order": (
                "all non-spectral diagonals first, then all remaining "
                "non-spectral unordered pairs in lexicographic order"
            ),
            "certificate_semantics": (
                "every degree above E-1 is positive for every pair by an "
                "exact spectral inequality or an all-target/exhaustive run; "
                "degree E-1 has an exhaustive zero or a row-sum obstruction"
            ),
            "total_interpolation_seconds": float(
                state["total_interpolation_seconds"]
            ),
            "total_counter_seconds": float(state["total_counter_seconds"]),
        },
    }


def compute_e_only_q2(
    ctx: SupersingularContext,
    *,
    checkpoint_path: Path,
    checkpoint_every: int = 250,
    progress=None,
) -> dict[str, Any]:
    """Compute a conclusive E(p) certificate in descending degree order."""
    if ctx.q != 2 or ctx.degree < 3:
        raise ValueError("the descending E-only backend requires q=2, deg(p)>=3")
    if checkpoint_every < 1:
        raise ValueError("checkpoint_every must be positive")
    if not COUNTER_SOURCE.exists() or not COUNTER_BINARY.exists():
        raise FileNotFoundError(
            "build the deterministic helper first with `make -C fast_q2`"
        )
    cutoff = theorem_bound(2, ctx.degree)
    state = _load_checkpoint(checkpoint_path, ctx, cutoff)
    if state.get("witness") is not None:
        return _final_record(ctx, state)

    pair_total = len(ctx.modules) * (len(ctx.modules) + 1) // 2
    while _expected_next_degree(state) >= 1:
        degree = _expected_next_degree(state)
        # A positive integer row with h entries has sum at least h.  If the
        # exact Brandt row sum is smaller, every level of this degree is
        # necessarily incomplete.
        row_sum = 2**degree + 1
        if row_sum < len(ctx.modules):
            state["witness"] = {
                "kind": "row_sum_obstruction",
                "degree": int(degree),
                "row_sum": int(row_sum),
                "class_count": len(ctx.modules),
            }
            state["current_degree"] = None
            state["current_pair_runs"] = []
            _atomic_json(checkpoint_path, state)
            return _final_record(ctx, state)

        levels, target_codes, level_by_code = _level_data(ctx, degree)
        if not levels:
            raise ArithmeticError(f"degree {degree} contains no target level")
        pairs = e_only_pair_order(ctx, degree)
        spectral_pair_count = pair_total - len(pairs)

        if state["current_degree"] is None:
            state["current_degree"] = int(degree)
            state["current_pair_runs"] = []
            _atomic_json(checkpoint_path, state)
        elif int(state["current_degree"]) != degree:
            raise ValueError("checkpoint current degree disagrees with descent")
        pair_runs = state["current_pair_runs"]
        start = len(pair_runs)
        if progress is not None:
            progress(
                f"degree {degree}: {len(pairs)} non-spectral pairs, "
                f"{spectral_pair_count} spectral pairs; resuming at {start}"
            )

        for pair_index in range(start, len(pairs)):
            source, target = pairs[pair_index]
            form = _interpolate_form(ctx, source, target, degree)
            run, counter = _counter_run(
                form,
                ctx=ctx,
                max_degree=degree,
                target_codes=target_codes,
                target_degrees=[degree],
                early=True,
            )
            state["total_vectors_visited"] += int(counter.iterations)
            state["total_interpolation_seconds"] += float(
                form["interpolation_seconds"]
            )
            state["total_counter_seconds"] += float(counter.seconds)
            if run["zero_codes"]:
                code = int(run["zero_codes"][0])
                state["witness"] = {
                    "kind": "exhaustive_zero",
                    "degree": int(degree),
                    "pair": [int(source), int(target)],
                    "code": code,
                    "ell": str(level_by_code[code]),
                    "run": run,
                }
                _atomic_json(checkpoint_path, state)
                if progress is not None:
                    progress(
                        f"degree {degree}: exact zero at pair "
                        f"({source},{target}), ell={level_by_code[code]}"
                    )
                return _final_record(ctx, state)
            pair_runs.append(
                {
                    "pair": [int(source), int(target)],
                    "run": run,
                }
            )
            if (
                (pair_index + 1) % checkpoint_every == 0
                or pair_index + 1 == len(pairs)
            ):
                _atomic_json(checkpoint_path, state)
            if progress is not None and (
                pair_index == start
                or pair_index + 1 == len(pairs)
                or (pair_index + 1) % 25 == 0
            ):
                progress(
                    f"degree {degree}: pair {pair_index + 1}/{len(pairs)} "
                    f"({source},{target}); "
                    f"vectors={state['total_vectors_visited']}"
                )

        state["positive_degrees"].append(
            {
                "degree": int(degree),
                "target_count": len(target_codes),
                "spectral_pair_count": int(spectral_pair_count),
                "pair_runs": pair_runs,
            }
        )
        state["current_degree"] = None
        state["current_pair_runs"] = []
        _atomic_json(checkpoint_path, state)
        if progress is not None:
            progress(f"degree {degree}: every level and pair is positive")

    raise ArithmeticError("descending scan reached degree zero without a witness")
