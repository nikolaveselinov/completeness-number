#!/usr/bin/env python3
"""Audit the published completeness-number census and retained certificate."""

from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path


REPOSITORY = Path(__file__).resolve().parents[1]
CENSUS = REPOSITORY / "results/census/completeness_numbers.csv"
FAMILY_SUMMARY = REPOSITORY / "results/census/family_summary.csv"
DETAILED = REPOSITORY / "results/detailed/characteristics.csv"
DEGREE_EIGHT = (
    REPOSITORY / "results/q2_degree8_compact/characteristics.csv"
)
DEGREE_NINE = (
    REPOSITORY
    / "results/q2_degree9_e_only/characteristics_summary.csv"
)
INDEPENDENT = REPOSITORY / "certificates/q2_degree9/raw"
INDEPENDENT_ROOT = INDEPENDENT.parent


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def key(row: dict[str, str]) -> tuple[int, str]:
    return int(row["q"]), row["p"]


def integer(row: dict[str, str], field: str) -> int | None:
    value = row.get(field, "")
    return None if value in ("", None) else int(value)


def compare_rows(
    canonical: dict[tuple[int, str], dict[str, str]],
    rows: list[dict[str, str]],
    source: str,
) -> set[tuple[int, str]]:
    fields = (
        "deg_p",
        "class_count",
        "U",
        "E",
        "largest_incomplete_degree",
    )
    seen: set[tuple[int, str]] = set()
    for row in rows:
        row_key = key(row)
        if row_key in seen:
            raise AssertionError(f"duplicate {row_key} in {source}")
        seen.add(row_key)
        expected = canonical.get(row_key)
        if expected is None:
            raise AssertionError(f"{source} contains uncatalogued row {row_key}")
        for field in fields:
            if integer(row, field) != integer(expected, field):
                raise AssertionError(
                    f"{source} disagrees at {row_key}, field {field}: "
                    f"{row.get(field)!r} != {expected.get(field)!r}"
                )
    return seen


def multiply_binary(left: int, right: int) -> int:
    result = 0
    while right:
        if right & 1:
            result ^= left
        left <<= 1
        right >>= 1
    return result


def translate_binary(code: int) -> int:
    """Return the code of f(T+1) for a binary polynomial f(T)."""
    translated = 0
    power = 1
    exponent = 0
    while code:
        if code & 1:
            translated ^= power
        code >>= 1
        exponent += 1
        if code:
            power = multiply_binary(power, 0b11)
    return translated


def load_json(name: str) -> dict:
    return json.loads((INDEPENDENT / name).read_text(encoding="utf-8"))


def audit_manifest() -> int:
    manifest = INDEPENDENT_ROOT / "SHA256SUMS"
    entries = 0
    for line in manifest.read_text(encoding="utf-8").splitlines():
        expected, relative = line.split("  ", 1)
        path = INDEPENDENT_ROOT / relative
        actual = hashlib.sha256(path.read_bytes()).hexdigest()
        if actual != expected:
            raise AssertionError(f"checksum mismatch for {relative}")
        entries += 1
    if entries != 9:
        raise AssertionError(f"expected 9 degree-nine manifest entries, found {entries}")
    return entries


def audit_independent_completion() -> dict[str, int]:
    degree14_diagonal = load_json("d9_degree14_diagonal_full.json")
    degree14_off_diagonal = load_json("d9_degree14_offdiag_full.json")
    degree13 = load_json("d9_degree13_full.json")
    degree13_partner = load_json("d9_partner_degree13_witness.json")
    regression4 = load_json("d3_degree4.json")
    regression5 = load_json("d3_degree5.json")

    if (
        degree14_diagonal["level_degree"] != 14
        or degree14_diagonal["completed_pair_count"] != 171
        or degree14_diagonal["expected_pair_count"] != 171
        or not degree14_diagonal["complete_for_scanned_scope"]
    ):
        raise AssertionError("degree-14 diagonal certificate is incomplete")
    if (
        degree14_off_diagonal["level_degree"] != 14
        or degree14_off_diagonal["completed_pair_count"] != 14535
        or degree14_off_diagonal["expected_pair_count"] != 14535
        or not degree14_off_diagonal["complete_for_scanned_scope"]
    ):
        raise AssertionError("degree-14 off-diagonal certificate is incomplete")

    obstruction = degree13["first_obstruction"]
    partner_obstruction = degree13_partner["first_obstruction"]
    if (
        degree13["level_degree"] != 13
        or obstruction is None
        or obstruction["pair"] != [0, 0]
        or not obstruction["exhaustive"]
        or len(obstruction["zero_codes"]) != 92
    ):
        raise AssertionError("first degree-13 obstruction is malformed")
    if (
        degree13_partner["level_degree"] != 13
        or partner_obstruction is None
        or partner_obstruction["pair"] != [0, 0]
        or not partner_obstruction["exhaustive"]
        or len(partner_obstruction["zero_codes"]) != 92
    ):
        raise AssertionError("translated degree-13 obstruction is malformed")
    if translate_binary(degree13["p_code"]) != degree13_partner["p_code"]:
        raise AssertionError("the two exceptional characteristics are not translations")
    translated_zeros = {
        translate_binary(code) for code in obstruction["zero_codes"]
    }
    if translated_zeros != set(partner_obstruction["zero_codes"]):
        raise AssertionError("the degree-13 zero sets do not translate exactly")

    if (
        not regression4["complete_for_scanned_scope"]
        or regression4["level_degree"] != 4
        or regression4["completed_pair_count"] != 6
    ):
        raise AssertionError("degree-four regression failed")
    regression_obstruction = regression5["first_obstruction"]
    if (
        regression5["level_degree"] != 5
        or regression_obstruction is None
        or regression_obstruction["zero_levels"]
        != ["T^5+T^3+T^2+T+1"]
    ):
        raise AssertionError("known degree-five regression failed")

    return {
        "degree14_pairs": 171 + 14535,
        "degree13_zero_levels_each": 92,
    }


def audit() -> dict[str, object]:
    census_rows = read_csv(CENSUS)
    if len(census_rows) != 1352:
        raise AssertionError(f"expected 1352 census rows, found {len(census_rows)}")
    canonical: dict[tuple[int, str], dict[str, str]] = {}
    for row in census_rows:
        row_key = key(row)
        if row_key in canonical:
            raise AssertionError(f"duplicate census row {row_key}")
        canonical[row_key] = row

    family_counts: Counter[tuple[int, int]] = Counter()
    family_values: dict[tuple[int, int], Counter[int]] = defaultdict(Counter)
    for row in census_rows:
        family = int(row["q"]), int(row["deg_p"])
        family_counts[family] += 1
        family_values[family][int(row["E"])] += 1
    if len(family_counts) != 31:
        raise AssertionError(f"expected 31 families, found {len(family_counts)}")

    summary_rows = read_csv(FAMILY_SUMMARY)
    if len(summary_rows) != 31:
        raise AssertionError("family summary does not contain 31 rows")
    for row in summary_rows:
        family = int(row["q"]), int(row["deg_p"])
        if int(row["characteristics"]) != family_counts[family]:
            raise AssertionError(f"family size mismatch for {family}")
        expected_distribution = ", ".join(
            f"{value} ({count})"
            for value, count in sorted(family_values[family].items())
        )
        if row["E_distribution"] != expected_distribution:
            raise AssertionError(
                f"distribution mismatch for {family}: "
                f"{row['E_distribution']!r} != {expected_distribution!r}"
            )

    detailed_rows = read_csv(DETAILED)
    degree8_rows = read_csv(DEGREE_EIGHT)
    degree9_rows = read_csv(DEGREE_NINE)
    covered = set()
    covered |= compare_rows(canonical, detailed_rows, "detailed table")
    covered |= compare_rows(canonical, degree8_rows, "degree-eight table")
    covered |= compare_rows(canonical, degree9_rows, "degree-nine table")
    if covered != set(canonical):
        missing = sorted(set(canonical) - covered)
        raise AssertionError(f"canonical rows lack source coverage: {missing[:5]}")

    degree9_distribution = Counter(int(row["E"]) for row in degree9_rows)
    if degree9_distribution != Counter({16: 54, 14: 2}):
        raise AssertionError(f"unexpected degree-nine distribution {degree9_distribution}")

    independent = audit_independent_completion()
    return {
        "characteristics": len(census_rows),
        "families": len(family_counts),
        "binary_characteristics": sum(
            count for (q, _degree), count in family_counts.items() if q == 2
        ),
        "degree9_distribution": dict(sorted(degree9_distribution.items())),
        "degree9_manifest_files": audit_manifest(),
        "detailed_rows": len(detailed_rows),
        "degree8_rows": len(degree8_rows),
        "degree9_rows": len(degree9_rows),
        **independent,
    }


def main() -> None:
    print(json.dumps(audit(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
