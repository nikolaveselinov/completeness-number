#!/usr/bin/env sage-python
"""Strictly validate one descending q=2 E(p)-only archive."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPOSITORY = Path(__file__).resolve().parents[1]
if str(REPOSITORY) not in sys.path:
    sys.path.insert(0, str(REPOSITORY))

from drinfeld_complete import build_supersingular_context, theorem_bound
from drinfeld_complete.compact_q2 import (
    COUNTER_BINARY,
    COUNTER_SOURCE,
    REPOSITORY,
    _sha256,
)
from drinfeld_complete.e_only_q2 import E_ONLY_KIND, e_only_pair_order
from fast_q2.validate_compact_archive import (
    ArchiveValidationError,
    _exact_int,
    _expected_levels,
    _nonnegative_number,
    _parse_characteristic,
    _replay_run,
    _require,
    _strict_keys,
    _validate_run,
)


def _levels_of_degree(expected_levels, degree):
    levels = [item for item in expected_levels if item["degree"] == degree]
    return levels, {item["code"]: item for item in levels}


def _validate_engine(engine_value):
    engine = _strict_keys(
        engine_value,
        required={
            "name",
            "counter_source",
            "counter_source_sha256",
            "counter_binary_sha256",
            "compiler_contract",
            "pair_order",
            "certificate_semantics",
            "total_interpolation_seconds",
            "total_counter_seconds",
        },
        optional=None,
        location="engine",
    )
    expected_text = {
        "name": "compiled_q2_descending_e_only",
        "counter_source": str(COUNTER_SOURCE.relative_to(REPOSITORY)),
        "counter_source_sha256": _sha256(COUNTER_SOURCE),
        "counter_binary_sha256": _sha256(COUNTER_BINARY),
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
    }
    for key, expected in expected_text.items():
        _require(
            engine[key] == expected,
            f"engine.{key}",
            f"expected {expected!r}",
        )
    _nonnegative_number(
        engine["total_interpolation_seconds"],
        "engine.total_interpolation_seconds",
    )
    _nonnegative_number(
        engine["total_counter_seconds"],
        "engine.total_counter_seconds",
    )


def validate_e_only_archive(record_value, *, recompute: bool = False):
    record = _strict_keys(
        record_value,
        required={
            "archive_kind",
            "schema_version",
            "q",
            "p",
            "degree_p",
            "theorem_bound_U",
            "E",
            "class_count",
            "weights",
            "j_invariants",
            "construction_checks",
            "positive_degrees",
            "witness",
            "total_vectors_visited",
            "engine",
        },
        optional=None,
        location="archive",
    )
    _require(
        record["archive_kind"] == E_ONLY_KIND,
        "archive_kind",
        f"expected {E_ONLY_KIND}",
    )
    _require(
        _exact_int(record["schema_version"], "schema_version") == 1,
        "schema_version",
        "expected 1",
    )
    _require(_exact_int(record["q"], "q") == 2, "q", "expected q=2")
    A, p = _parse_characteristic(record["p"])
    degree_p = int(p.degree())
    _require(
        _exact_int(record["degree_p"], "degree_p") == degree_p,
        "degree_p",
        "does not match p",
    )
    cutoff = theorem_bound(2, degree_p)
    _require(
        _exact_int(record["theorem_bound_U"], "theorem_bound_U") == cutoff,
        "theorem_bound_U",
        f"expected {cutoff}",
    )
    expected_levels, _ = _expected_levels(A, p, cutoff)
    ctx = build_supersingular_context(2, p)
    vertex_count = len(ctx.modules)
    _require(
        _exact_int(record["class_count"], "class_count") == vertex_count,
        "class_count",
        f"expected {vertex_count}",
    )
    _require(
        record["weights"]
        == [int(order) for order in ctx.automorphism_orders],
        "weights",
        "does not match rebuilt context",
    )
    _require(
        record["j_invariants"] == [str(value) for value in ctx.j_invariants],
        "j_invariants",
        "canonical vertex ordering differs",
    )
    _require(
        record["construction_checks"] == ctx.checks,
        "construction_checks",
        "does not match rebuilt context",
    )
    _validate_engine(record["engine"])

    witness = record["witness"]
    if not isinstance(witness, dict):
        raise ArchiveValidationError("witness: expected a JSON object")
    witness_degree = _exact_int(witness.get("degree"), "witness.degree")
    _require(
        1 <= witness_degree < cutoff,
        "witness.degree",
        "must lie in 1..U-1",
    )
    _require(
        _exact_int(record["E"], "E") == witness_degree + 1,
        "E",
        "must equal witness degree plus one",
    )

    positive = record["positive_degrees"]
    if not isinstance(positive, list):
        raise ArchiveValidationError(
            "positive_degrees: expected a JSON array"
        )
    expected_positive_degrees = list(
        range(cutoff - 1, witness_degree, -1)
    )
    _require(
        [item.get("degree") for item in positive]
        == expected_positive_degrees,
        "positive_degrees",
        "must cover every degree strictly above the witness",
    )
    # `total_vectors_visited` is a work counter, not merely a sum of the
    # retained certificate runs.  At the eventual witness degree the search
    # may first test positive pairs; those disposable attempts need not be
    # archived because they prove nothing needed for E(p), but their vector
    # visits remain part of the honest work total.
    certificate_iterations = 0
    replayed_runs = 0
    pair_total = vertex_count * (vertex_count + 1) // 2
    for degree_record in positive:
        degree = int(degree_record["degree"])
        required_pairs = e_only_pair_order(ctx, degree)
        runs = degree_record.get("pair_runs")
        if not isinstance(runs, list):
            raise ArchiveValidationError(
                f"positive degree {degree}: pair_runs must be a list"
            )
        archived_pairs = [tuple(entry["pair"]) for entry in runs]
        _require(
            archived_pairs == required_pairs,
            f"positive degree {degree}.pair_runs",
            "must cover every non-spectral pair in canonical order",
        )
        spectral_count = pair_total - len(required_pairs)
        _require(
            _exact_int(
                degree_record.get("spectral_pair_count"),
                f"positive degree {degree}.spectral_pair_count",
            )
            == spectral_count,
            f"positive degree {degree}.spectral_pair_count",
            f"expected {spectral_count}",
        )
        levels, by_code = _levels_of_degree(expected_levels, degree)
        _require(
            _exact_int(
                degree_record.get("target_count"),
                f"positive degree {degree}.target_count",
            )
            == len(levels),
            f"positive degree {degree}.target_count",
            f"expected {len(levels)}",
        )
        for index, entry in enumerate(runs):
            pair = tuple(entry["pair"])
            validated = _validate_run(
                entry["run"],
                pair=pair,
                index=index,
                degree_p=degree_p,
                expected_levels=levels,
                expected_by_code=by_code,
            )
            _require(
                validated["max_degree"] == degree,
                validated["location"],
                "E-only run must use the exact target degree bound",
            )
            _require(
                not validated["zero_codes"],
                validated["location"],
                "a certified positive degree cannot contain a zero",
            )
            certificate_iterations += validated["iterations"]
            if recompute:
                _replay_run(ctx, pair, validated)
                replayed_runs += 1

    witness_kind = witness.get("kind")
    if witness_kind == "row_sum_obstruction":
        expected_keys = {
            "kind",
            "degree",
            "row_sum",
            "class_count",
        }
        _require(
            set(witness) == expected_keys,
            "witness",
            f"expected keys {sorted(expected_keys)}",
        )
        row_sum = 2**witness_degree + 1
        _require(
            _exact_int(witness["row_sum"], "witness.row_sum") == row_sum,
            "witness.row_sum",
            f"expected {row_sum}",
        )
        _require(
            _exact_int(witness["class_count"], "witness.class_count")
            == vertex_count,
            "witness.class_count",
            f"expected {vertex_count}",
        )
        _require(
            row_sum < vertex_count,
            "witness",
            "row sum is not small enough to force a zero",
        )
    elif witness_kind == "exhaustive_zero":
        expected_keys = {
            "kind",
            "degree",
            "pair",
            "code",
            "ell",
            "run",
        }
        _require(
            set(witness) == expected_keys,
            "witness",
            f"expected keys {sorted(expected_keys)}",
        )
        pair = tuple(witness["pair"])
        _require(
            pair in e_only_pair_order(ctx, witness_degree),
            "witness.pair",
            "pair is spectral, malformed, or outside the context",
        )
        levels, by_code = _levels_of_degree(
            expected_levels, witness_degree
        )
        validated = _validate_run(
            witness["run"],
            pair=pair,
            index=0,
            degree_p=degree_p,
            expected_levels=levels,
            expected_by_code=by_code,
        )
        _require(
            validated["outcome"] == "exhaustive",
            "witness.run",
            "a zero requires complete Hom-space exhaustion",
        )
        _require(
            validated["max_degree"] == witness_degree,
            "witness.run",
            "must use the exact witness degree bound",
        )
        code = _exact_int(witness["code"], "witness.code")
        _require(
            code in validated["zero_codes"],
            "witness.code",
            "is not one of the exhaustively absent targets",
        )
        _require(
            witness["ell"] == by_code[code]["ell"],
            "witness.ell",
            f"expected {by_code[code]['ell']}",
        )
        certificate_iterations += validated["iterations"]
        witness_pairs = e_only_pair_order(ctx, witness_degree)
        discarded_witness_runs = witness_pairs.index(pair)
        witness_run_capacity = (1 << validated["dimension"]) - 1
        if recompute:
            _replay_run(ctx, pair, validated)
            replayed_runs += 1
    else:
        raise ArchiveValidationError(
            f"witness.kind: unsupported value {witness_kind!r}"
        )

    work_iterations = _exact_int(
        record["total_vectors_visited"], "total_vectors_visited"
    )
    if witness_kind == "row_sum_obstruction":
        discarded_witness_runs = 0
        minimum_work = maximum_work = certificate_iterations
    else:
        # Every earlier pair in the canonical witness-degree search produced
        # a positive run.  The archive intentionally discards these runs.
        # Each used at least one and at most the full Hom-space traversal.
        minimum_work = certificate_iterations + discarded_witness_runs
        maximum_work = (
            certificate_iterations
            + discarded_witness_runs * witness_run_capacity
        )
    _require(
        minimum_work <= work_iterations <= maximum_work,
        "total_vectors_visited",
        f"must lie in {minimum_work}..{maximum_work}",
    )
    return {
        "q": 2,
        "p": str(p),
        "degree_p": degree_p,
        "theorem_bound_U": cutoff,
        "E": witness_degree + 1,
        "class_count": vertex_count,
        "positive_degrees": len(positive),
        "witness_kind": witness_kind,
        "witness_degree": witness_degree,
        "total_iterations": work_iterations,
        "certificate_iterations": certificate_iterations,
        "discarded_witness_runs": discarded_witness_runs,
        "recompute": recompute,
        "replayed_runs": replayed_runs,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    parser.add_argument("--recompute", action="store_true")
    args = parser.parse_args()
    try:
        record = json.loads(args.archive.read_text())
        summary = validate_e_only_archive(record, recompute=args.recompute)
    except (OSError, json.JSONDecodeError, ArchiveValidationError) as error:
        raise SystemExit(f"validation failed: {error}") from error
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
