#!/usr/bin/env sage-python
"""Validate compact exact-zero archives for supported q=2 families.

The structural mode rebuilds all context and level data and checks that every
unordered vertex-pair/level decision has a valid archived counter
certificate.  ``--recompute`` additionally rebuilds and reruns every bounded
Hom-space counter invocation.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

from sage.all import GF, PolynomialRing, sage_eval

REPOSITORY = Path(__file__).resolve().parents[1]
if str(REPOSITORY) not in sys.path:
    sys.path.insert(0, str(REPOSITORY))

from drinfeld_complete import (
    build_supersingular_context,
    monic_irreducibles,
    theorem_bound,
)
from drinfeld_complete.core import _poly_code_q2, pair_is_spectrally_positive


ARCHIVE_KIND = "q2_compact_zero_v1"
ARCHIVE_KIND_DEGREE8 = "q2_degree8_compact_zero_v1"
SCHEMA_VERSION = 1
ARCHIVE_SPECS = {
    ARCHIVE_KIND: {
        "degree_p": 7,
        "cutoff": 16,
        "level_count": 4719,
        "vertex_count": 43,
        "unordered_pair_count": 946,
    },
    ARCHIVE_KIND_DEGREE8: {
        "degree_p": 8,
        "cutoff": 15,
        "level_count": 2537,
        "vertex_count": 85,
        "unordered_pair_count": 3655,
    },
}


class ArchiveValidationError(ValueError):
    """A compact archive failed an exact validation condition."""


def _fail(location: str, message: str) -> None:
    raise ArchiveValidationError(f"{location}: {message}")


def _require(condition: bool, location: str, message: str) -> None:
    if not condition:
        _fail(location, message)


def _exact_int(value: Any, location: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        _fail(location, "expected a JSON integer")
    return value


def _exact_bool(value: Any, location: str) -> bool:
    if not isinstance(value, bool):
        _fail(location, "expected a JSON boolean")
    return value


def _nonnegative_number(value: Any, location: str) -> float:
    if (
        isinstance(value, bool)
        or not isinstance(value, (int, float))
        or not math.isfinite(value)
        or value < 0
    ):
        _fail(location, "expected a finite nonnegative JSON number")
    return float(value)


def _strict_keys(
    obj: Any,
    *,
    required: set[str],
    optional: set[str] | None,
    location: str,
) -> dict[str, Any]:
    if not isinstance(obj, dict):
        _fail(location, "expected a JSON object")
    missing = required - set(obj)
    unknown = set(obj) - required - (optional or set())
    if missing:
        _fail(location, f"missing keys {sorted(missing)}")
    if unknown:
        _fail(location, f"unknown keys {sorted(unknown)}")
    return obj


def _parse_characteristic(text: Any):
    if not isinstance(text, str):
        _fail("p", "expected a polynomial string")
    F2 = GF(2)
    A = PolynomialRing(F2, "T")
    T = A.gen()
    try:
        p = A(sage_eval(text.replace("^", "**"), locals={"T": T}))
    except Exception as error:
        _fail("p", f"could not parse polynomial: {error}")
    _require(str(p) == text, "p", f"noncanonical spelling; expected {p}")
    _require(p.is_monic(), "p", "characteristic is not monic")
    _require(p.is_irreducible(), "p", "characteristic is not irreducible")
    return A, p


def _expected_levels(A, p, cutoff: int) -> tuple[list[dict[str, Any]], dict[int, dict]]:
    levels: list[dict[str, Any]] = []
    by_code: dict[int, dict[str, Any]] = {}
    for degree in range(1, cutoff):
        for ell in monic_irreducibles(A, degree):
            if ell == p:
                continue
            code = int(_poly_code_q2(ell))
            item = {"ell": str(ell), "code": code, "degree": degree}
            _require(code not in by_code, "level universe", f"duplicate code {code}")
            levels.append(item)
            by_code[code] = item
    levels.sort(key=lambda item: (item["degree"], item["code"]))
    return levels, by_code


def _sorted_unique_ints(values: Any, location: str) -> list[int]:
    if not isinstance(values, list):
        _fail(location, "expected a JSON array")
    result = [_exact_int(value, f"{location}[{index}]") for index, value in enumerate(values)]
    _require(result == sorted(set(result)), location, "must be sorted and duplicate-free")
    return result


def _pair(value: Any, location: str, vertex_count: int) -> tuple[int, int]:
    if not isinstance(value, list) or len(value) != 2:
        _fail(location, "expected [source, target]")
    i = _exact_int(value[0], f"{location}[0]")
    j = _exact_int(value[1], f"{location}[1]")
    _require(0 <= i <= j < vertex_count, location, "pair must satisfy 0 <= i <= j < h(p)")
    return i, j


def _ordered_pairs(value: Any, location: str, vertex_count: int) -> list[tuple[int, int]]:
    if not isinstance(value, list):
        _fail(location, "expected a JSON array")
    pairs: list[tuple[int, int]] = []
    for index, item in enumerate(value):
        item_location = f"{location}[{index}]"
        if not isinstance(item, list) or len(item) != 2:
            _fail(item_location, "expected [source, target]")
        i = _exact_int(item[0], f"{item_location}[0]")
        j = _exact_int(item[1], f"{item_location}[1]")
        _require(
            0 <= i < vertex_count and 0 <= j < vertex_count,
            item_location,
            "vertex index outside archive context",
        )
        pairs.append((i, j))
    _require(pairs == sorted(set(pairs)), location, "must be sorted and duplicate-free")
    return pairs


def _run_targets(
    run: dict[str, Any],
    *,
    expected_levels: list[dict[str, Any]],
    expected_by_code: dict[int, dict[str, Any]],
    location: str,
) -> list[int]:
    has_degrees = "target_degrees" in run
    has_codes = "target_codes" in run
    _require(
        has_degrees != has_codes,
        location,
        "exactly one of target_degrees and target_codes is required",
    )
    if has_degrees:
        degrees = _sorted_unique_ints(run["target_degrees"], f"{location}.target_degrees")
        _require(degrees, f"{location}.target_degrees", "selector cannot be empty")
        codes = [
            level["code"]
            for level in expected_levels
            if level["degree"] in set(degrees)
        ]
        actual_degrees = sorted({expected_by_code[code]["degree"] for code in codes})
        _require(
            actual_degrees == degrees,
            f"{location}.target_degrees",
            "contains a degree with no expected level",
        )
        return codes

    codes = _sorted_unique_ints(run["target_codes"], f"{location}.target_codes")
    _require(codes, f"{location}.target_codes", "selector cannot be empty")
    unknown = [code for code in codes if code not in expected_by_code]
    _require(not unknown, f"{location}.target_codes", f"unknown level codes {unknown}")
    return codes


def _validate_run(
    run_value: Any,
    *,
    pair: tuple[int, int],
    index: int,
    degree_p: int,
    expected_levels: list[dict[str, Any]],
    expected_by_code: dict[int, dict[str, Any]],
) -> dict[str, Any]:
    location = f"pair_certificates[{pair}].runs[{index}]"
    run = _strict_keys(
        run_value,
        required={
            "max_degree",
            "dimension",
            "iterations",
            "target_count",
            "seen",
            "invalid_norms",
            "exhaustive",
            "outcome",
            "zero_codes",
            "direct_sage_norm_cross_checks",
            "seconds",
        },
        optional={"target_degrees", "target_codes"},
        location=location,
    )
    codes = _run_targets(
        run,
        expected_levels=expected_levels,
        expected_by_code=expected_by_code,
        location=location,
    )
    max_degree = _exact_int(run["max_degree"], f"{location}.max_degree")
    largest_required_degree = max(
        level["degree"] for level in expected_levels
    )
    _require(
        max_degree >= max(expected_by_code[code]["degree"] for code in codes),
        f"{location}.max_degree",
        "is smaller than a selected target degree",
    )
    _require(
        max_degree <= largest_required_degree,
        f"{location}.max_degree",
        f"must not exceed the finite scan bound {largest_required_degree}",
    )
    _require(
        max_degree >= degree_p - 2,
        f"{location}.max_degree",
        "compact runs must lie in the stable bounded-Hom dimension range",
    )
    dimension = _exact_int(run["dimension"], f"{location}.dimension")
    expected_dimension = 2 * (max_degree + 1) - (degree_p - 1)
    _require(
        dimension == expected_dimension and dimension > 0,
        f"{location}.dimension",
        f"expected bounded-Hom dimension {expected_dimension}",
    )
    _require(
        dimension < 63,
        f"{location}.dimension",
        "counter traversal size does not fit its uint64 contract",
    )
    total = (1 << dimension) - 1
    iterations = _exact_int(run["iterations"], f"{location}.iterations")
    _require(0 < iterations <= total, f"{location}.iterations", f"must lie in 1..{total}")
    target_count = _exact_int(run["target_count"], f"{location}.target_count")
    _require(
        target_count == len(codes),
        f"{location}.target_count",
        f"selector contains {len(codes)} levels",
    )
    seen = _exact_int(run["seen"], f"{location}.seen")
    invalid = _exact_int(run["invalid_norms"], f"{location}.invalid_norms")
    _require(invalid == 0, f"{location}.invalid_norms", "must be zero")
    direct_checks = _exact_int(
        run["direct_sage_norm_cross_checks"],
        f"{location}.direct_sage_norm_cross_checks",
    )
    _require(
        0 < direct_checks <= total,
        f"{location}.direct_sage_norm_cross_checks",
        f"must lie in 1..{total}",
    )
    _nonnegative_number(run["seconds"], f"{location}.seconds")
    exhaustive = _exact_bool(run["exhaustive"], f"{location}.exhaustive")
    outcome = run["outcome"]
    _require(
        outcome in {"exhaustive", "all_targets_witnessed"},
        f"{location}.outcome",
        "unsupported certificate outcome",
    )
    zero_codes = _sorted_unique_ints(run["zero_codes"], f"{location}.zero_codes")
    _require(
        set(zero_codes) <= set(codes),
        f"{location}.zero_codes",
        "contains a level outside this run's selector",
    )

    status_by_code: dict[int, str]
    if outcome == "exhaustive":
        _require(exhaustive, f"{location}.exhaustive", "must be true for exhaustive outcome")
        _require(
            iterations == total,
            f"{location}.iterations",
            "an exhaustive traversal must visit every nonzero vector",
        )
        _require(
            seen == target_count - len(zero_codes),
            f"{location}.seen",
            "must equal target_count minus the exact zero count",
        )
        status_by_code = {
            code: ("exhaustive_zero" if code in set(zero_codes) else "exhaustive_positive")
            for code in codes
        }
    else:
        _require(
            not exhaustive,
            f"{location}.exhaustive",
            "a completed traversal must use outcome 'exhaustive'",
        )
        _require(
            iterations < total,
            f"{location}.iterations",
            "all_targets_witnessed must be an actual early stop",
        )
        _require(seen == target_count, f"{location}.seen", "not every target was witnessed")
        _require(not zero_codes, f"{location}.zero_codes", "witnessed targets cannot be zero")
        status_by_code = {code: "all_targets_witnessed" for code in codes}

    return {
        "codes": codes,
        "status_by_code": status_by_code,
        "max_degree": max_degree,
        "dimension": dimension,
        "iterations": iterations,
        "exhaustive": exhaustive,
        "seen": seen,
        "target_count": target_count,
        "zero_codes": zero_codes,
        "outcome": outcome,
        "direct_sage_norm_cross_checks": direct_checks,
        "raw": run,
        "location": location,
    }


def _replay_run(ctx, pair: tuple[int, int], validated_run: dict[str, Any]) -> None:
    # Imported lazily so structural validation never launches the helper.
    # The read-only wrapper refuses a missing/stale binary and fingerprints it
    # before and after; replay must never build or replace a live executable.
    from fast_q2.compare_walsh_gray import run_read_only_gray_counter
    from fast_q2.gray_counter import interpolate_pair

    i, j = pair
    diagonal, cross = interpolate_pair(ctx, i, j, validated_run["max_degree"])
    replay = run_read_only_gray_counter(
        diagonal,
        cross,
        max_degree=validated_run["max_degree"],
        field_degree=2 * ctx.degree,
        target_codes=validated_run["codes"],
        stop_when_seen=validated_run["outcome"] == "all_targets_witnessed",
    )
    location = validated_run["location"]
    _require(replay.iterations == validated_run["iterations"], location, "replay iteration mismatch")
    _require(replay.exhaustive == validated_run["exhaustive"], location, "replay exhaustion mismatch")
    _require(replay.seen == validated_run["seen"], location, "replay seen-count mismatch")
    _require(replay.target_count == validated_run["target_count"], location, "replay target-count mismatch")
    _require(replay.invalid_norms == 0, location, "replay encountered invalid norms")
    replay_zeros = sorted(code for code in validated_run["codes"] if replay.counts.get(code, 0) == 0)
    _require(
        replay_zeros == validated_run["zero_codes"],
        location,
        f"replay zero set differs: {replay_zeros}",
    )
    if validated_run["outcome"] == "all_targets_witnessed":
        _require(
            all(replay.counts.get(code, 0) > 0 for code in validated_run["codes"]),
            location,
            "replay did not witness every target",
        )


def _validate_degree_summary(
    archived: Any,
    *,
    expected_levels: list[dict[str, Any]],
    zero_entries_by_code: dict[int, list[tuple[int, int]]],
    cutoff: int,
) -> None:
    if not isinstance(archived, dict):
        _fail("degree_summary", "expected a JSON object")
    _require(
        set(archived) == {str(degree) for degree in range(1, cutoff)},
        "degree_summary",
        "degree keys must be exactly 1 through U-1",
    )
    levels_by_degree: dict[int, list[dict[str, Any]]] = defaultdict(list)
    for level in expected_levels:
        levels_by_degree[level["degree"]].append(level)
    for degree in range(1, cutoff):
        location = f"degree_summary.{degree}"
        summary = _strict_keys(
            archived[str(degree)],
            required={
                "level_count",
                "complete_count",
                "incomplete_count",
                "incomplete_levels",
            },
            optional=None,
            location=location,
        )
        levels = levels_by_degree[degree]
        incomplete = [
            level["ell"] for level in levels if zero_entries_by_code[level["code"]]
        ]
        _require(
            _exact_int(summary["level_count"], f"{location}.level_count") == len(levels),
            f"{location}.level_count",
            f"expected {len(levels)}",
        )
        _require(
            _exact_int(summary["incomplete_count"], f"{location}.incomplete_count")
            == len(incomplete),
            f"{location}.incomplete_count",
            f"expected {len(incomplete)}",
        )
        _require(
            _exact_int(summary["complete_count"], f"{location}.complete_count")
            == len(levels) - len(incomplete),
            f"{location}.complete_count",
            f"expected {len(levels) - len(incomplete)}",
        )
        _require(
            summary["incomplete_levels"] == incomplete,
            f"{location}.incomplete_levels",
            f"expected {incomplete}",
        )


def validate_archive(
    record_value: Any,
    *,
    recompute: bool = False,
) -> dict[str, Any]:
    """Validate a decoded compact archive and return independently derived data."""
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
            "pair_certificates",
            "levels",
            "degree_summary",
        },
        optional={"construction_checks", "engine", "total_vectors_visited"},
        location="archive",
    )
    _require(
        record["archive_kind"] in ARCHIVE_SPECS,
        "archive_kind",
        f"expected one of {sorted(ARCHIVE_SPECS)}",
    )
    spec = ARCHIVE_SPECS[record["archive_kind"]]
    _require(
        _exact_int(record["schema_version"], "schema_version") == SCHEMA_VERSION,
        "schema_version",
        f"expected {SCHEMA_VERSION}",
    )
    _require(_exact_int(record["q"], "q") == 2, "q", "compact backend is q=2 only")
    A, p = _parse_characteristic(record["p"])
    degree_p = int(p.degree())
    _require(
        degree_p == spec["degree_p"],
        "p",
        f"archive kind requires degree {spec['degree_p']}",
    )
    _require(_exact_int(record["degree_p"], "degree_p") == degree_p, "degree_p", "does not match p")

    cutoff = theorem_bound(2, degree_p)
    _require(
        cutoff == spec["cutoff"],
        "theorem_bound_U",
        f"internal cutoff regression: expected {spec['cutoff']}",
    )
    _require(
        _exact_int(record["theorem_bound_U"], "theorem_bound_U") == cutoff,
        "theorem_bound_U",
        f"exact strict-inequality cutoff is {cutoff}",
    )
    expected_levels, expected_by_code = _expected_levels(A, p, cutoff)
    _require(
        len(expected_levels) == spec["level_count"],
        "level universe",
        f"expected {spec['level_count']} levels",
    )

    ctx = build_supersingular_context(2, p)
    vertex_count = len(ctx.modules)
    _require(
        vertex_count == spec["vertex_count"],
        "supersingular context",
        f"expected {spec['vertex_count']} vertices",
    )
    _require(
        _exact_int(record["class_count"], "class_count") == vertex_count,
        "class_count",
        f"expected {vertex_count}",
    )
    expected_weights = [int(order) for order in ctx.automorphism_orders]
    _require(record["weights"] == expected_weights, "weights", f"expected {expected_weights}")
    expected_j = [str(value) for value in ctx.j_invariants]
    _require(record["j_invariants"] == expected_j, "j_invariants", "canonical vertex ordering differs")
    if "construction_checks" in record:
        _require(
            record["construction_checks"] == ctx.checks,
            "construction_checks",
            "does not equal the independently rebuilt context checks",
        )

    pair_values = record["pair_certificates"]
    if not isinstance(pair_values, list):
        _fail("pair_certificates", "expected a JSON array")
    expected_pairs = [
        (i, j) for i in range(vertex_count) for j in range(i, vertex_count)
    ]
    _require(
        len(expected_pairs) == spec["unordered_pair_count"],
        "pair universe",
        f"expected {spec['unordered_pair_count']} unordered pairs",
    )
    archived_pairs: list[tuple[int, int]] = []
    zero_unordered_by_code: dict[int, list[tuple[int, int]]] = {
        level["code"]: [] for level in expected_levels
    }
    certificate_counts: Counter[str] = Counter()
    replayed_runs = 0
    total_iterations = 0

    for pair_index, pair_value in enumerate(pair_values):
        pair_location = f"pair_certificates[{pair_index}]"
        pair_record = _strict_keys(
            pair_value,
            required={"pair", "spectral_degrees", "runs"},
            optional=None,
            location=pair_location,
        )
        pair = _pair(pair_record["pair"], f"{pair_location}.pair", vertex_count)
        archived_pairs.append(pair)
        spectral_degrees = _sorted_unique_ints(
            pair_record["spectral_degrees"],
            f"{pair_location}.spectral_degrees",
        )
        _require(
            all(1 <= degree < cutoff for degree in spectral_degrees),
            f"{pair_location}.spectral_degrees",
            "degrees must lie in the finite scan range 1..U-1",
        )
        for degree in spectral_degrees:
            _require(
                pair_is_spectrally_positive(ctx, pair[0], pair[1], degree),
                f"{pair_location}.spectral_degrees",
                f"degree {degree} fails the exact pair-specific inequality",
            )
        runs = pair_record["runs"]
        if not isinstance(runs, list):
            _fail(f"{pair_location}.runs", "expected a JSON array")
        status_by_code: dict[int, str] = {
            level["code"]: "spectral_positive"
            for level in expected_levels
            if level["degree"] in set(spectral_degrees)
        }
        for run_index, run in enumerate(runs):
            validated_run = _validate_run(
                run,
                pair=pair,
                index=run_index,
                degree_p=degree_p,
                expected_levels=expected_levels,
                expected_by_code=expected_by_code,
            )
            overlap = set(status_by_code) & set(validated_run["status_by_code"])
            _require(
                not overlap,
                validated_run["location"],
                f"selector overlaps earlier runs at codes {sorted(overlap)}",
            )
            status_by_code.update(validated_run["status_by_code"])
            total_iterations += validated_run["iterations"]
            if recompute:
                _replay_run(ctx, pair, validated_run)
                replayed_runs += 1

        missing = set(expected_by_code) - set(status_by_code)
        extra = set(status_by_code) - set(expected_by_code)
        _require(not missing, f"{pair_location}.runs", f"missing {len(missing)} target decisions")
        _require(not extra, f"{pair_location}.runs", f"contains {len(extra)} unknown decisions")
        for code, status in status_by_code.items():
            certificate_counts[status] += 1
            if status == "exhaustive_zero":
                zero_unordered_by_code[code].append(pair)

    _require(
        archived_pairs == expected_pairs,
        "pair_certificates",
        "records must be the exact lexicographically sorted unordered-pair universe",
    )
    _require(
        sum(certificate_counts.values()) == len(expected_pairs) * len(expected_levels),
        "pair certificates",
        "coverage count regression",
    )

    zero_entries_by_code: dict[int, list[tuple[int, int]]] = {}
    for code, pairs in zero_unordered_by_code.items():
        ordered: list[tuple[int, int]] = []
        for i, j in pairs:
            # This is the only transport used by the compact archive.  Since
            # all weights are positive, weighted Brandt symmetry transports
            # zero/nonzero exactly even when entries themselves differ.
            if i == j:
                ordered.append((i, i))
            else:
                ordered.extend([(i, j), (j, i)])
        zero_entries_by_code[code] = sorted(ordered)

    level_values = record["levels"]
    if not isinstance(level_values, list):
        _fail("levels", "expected a JSON array")
    _require(
        len(level_values) == len(expected_levels),
        "levels",
        f"expected exactly {len(expected_levels)} records",
    )
    archived_level_keys: list[tuple[int, int]] = []
    for index, (level_value, expected) in enumerate(zip(level_values, expected_levels)):
        location = f"levels[{index}]"
        level = _strict_keys(
            level_value,
            required={"ell", "code", "degree", "complete", "zero_entries"},
            optional=None,
            location=location,
        )
        code = _exact_int(level["code"], f"{location}.code")
        degree = _exact_int(level["degree"], f"{location}.degree")
        archived_level_keys.append((degree, code))
        _require(level["ell"] == expected["ell"], f"{location}.ell", f"expected {expected['ell']}")
        _require(code == expected["code"], f"{location}.code", f"expected {expected['code']}")
        _require(degree == expected["degree"], f"{location}.degree", f"expected {expected['degree']}")
        archived_zeros = _ordered_pairs(
            level["zero_entries"], f"{location}.zero_entries", vertex_count
        )
        derived_zeros = zero_entries_by_code[code]
        _require(
            archived_zeros == derived_zeros,
            f"{location}.zero_entries",
            f"certificate-derived set has {len(derived_zeros)} entries",
        )
        zero_set = set(derived_zeros)
        _require(
            all(
                any((source, target) not in zero_set for target in range(vertex_count))
                for source in range(vertex_count)
            ),
            f"{location}.zero_entries",
            "contains an impossible all-zero Brandt row",
        )
        complete = _exact_bool(level["complete"], f"{location}.complete")
        _require(complete == (not derived_zeros), f"{location}.complete", "does not match exact zero set")

    expected_level_keys = [
        (level["degree"], level["code"]) for level in expected_levels
    ]
    _require(
        archived_level_keys == expected_level_keys,
        "levels",
        "must be the exact canonically sorted irreducible-level universe",
    )
    _validate_degree_summary(
        record["degree_summary"],
        expected_levels=expected_levels,
        zero_entries_by_code=zero_entries_by_code,
        cutoff=cutoff,
    )

    bad_degrees = [
        level["degree"]
        for level in expected_levels
        if zero_entries_by_code[level["code"]]
    ]
    expected_E = 1 + max(bad_degrees, default=0)
    _require(
        _exact_int(record["E"], "E") == expected_E,
        "E",
        f"certificate-derived value is {expected_E}",
    )
    _require(expected_E <= cutoff, "E", "exceeds the proven theorem cutoff")

    if "total_vectors_visited" in record:
        _require(
            _exact_int(record["total_vectors_visited"], "total_vectors_visited")
            == total_iterations,
            "total_vectors_visited",
            f"expected {total_iterations}",
        )

    return {
        "q": 2,
        "p": str(p),
        "degree_p": degree_p,
        "theorem_bound_U": cutoff,
        "E": expected_E,
        "class_count": vertex_count,
        "level_count": len(expected_levels),
        "unordered_pair_count": len(expected_pairs),
        "pair_level_decisions": len(expected_pairs) * len(expected_levels),
        "certificate_counts": dict(sorted(certificate_counts.items())),
        "incomplete_level_count": len(bad_degrees),
        "largest_incomplete_degree": max(bad_degrees, default=None),
        "total_iterations": total_iterations,
        "recompute": recompute,
        "replayed_runs": replayed_runs,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    parser.add_argument(
        "--recompute",
        action="store_true",
        help="rebuild and replay every archived counter invocation",
    )
    args = parser.parse_args()
    try:
        record = json.loads(args.archive.read_text())
        summary = validate_archive(record, recompute=args.recompute)
    except (OSError, json.JSONDecodeError, ArchiveValidationError) as error:
        raise SystemExit(f"validation failed: {error}") from error
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
