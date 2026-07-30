#!/usr/bin/env python3
"""Audit result JSON files and create consolidated tables and a report."""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path

from sage.all import GF, PolynomialRing, sage_eval

from drinfeld_complete import build_supersingular_context


def divisors(n: int):
    return [d for d in range(1, n + 1) if n % d == 0]


def mobius(n: int) -> int:
    factors = 0
    candidate = 2
    while candidate * candidate <= n:
        if n % candidate == 0:
            n //= candidate
            factors += 1
            if n % candidate == 0:
                return 0
            while n % candidate == 0:
                n //= candidate
        candidate += 1
    if n > 1:
        factors += 1
    return -1 if factors % 2 else 1


def irreducible_count(q: int, degree: int) -> int:
    return sum(
        mobius(divisor) * q ** (degree // divisor)
        for divisor in divisors(degree)
    ) // degree


def theorem_C(q: int, degree: int) -> int:
    if degree < 3:
        return 0
    if degree % 2:
        return (q**degree - q) // (q - 1)
    return (q**degree - q**2) // (q**2 - 1)


def theorem_bound(q: int, degree: int) -> int:
    if degree <= 2:
        return 1
    C = theorem_C(q, degree)
    for exponent in range(1, 10 * degree + 20):
        Q = q**exponent
        if (Q + 1) ** 2 > 4 * C**2 * Q:
            return exponent
    raise AssertionError("cutoff search unexpectedly failed")


def parse_polynomial(A, text: str):
    T = A.gen()
    return A(
        sage_eval(
            text.replace("^", "**"),
            locals={"T": T, "a": A.base_ring().gen()},
        )
    )


def validate(record: dict, path: Path) -> None:
    q = int(record["q"])
    d = int(record["degree_p"])
    U = int(record["theorem_bound_U"])
    Fq = GF(q, name="a")
    A = PolynomialRing(Fq, "T")
    p = parse_polynomial(A, record["p"])
    assert p.is_monic() and p.is_irreducible() and int(p.degree()) == d, path
    assert U == theorem_bound(q, d), path

    # Rebuild the supersingular context instead of trusting the archived
    # construction data.
    context = build_supersingular_context(q, p)
    assert str(context.hasse_polynomial) == record["hasse_polynomial"], path
    assert (
        str(context.supersingular_polynomial)
        == record["supersingular_polynomial"]
    ), path
    assert [str(value) for value in context.j_invariants] == record["j_invariants"], path
    assert len(context.modules) == int(record["class_count"]), path
    assert [
        order // (q - 1) for order in context.automorphism_orders
    ] == record["weights"], path

    checks = record["construction_checks"]
    assert int(checks["hasse_degree"]) == int(checks["expected_hasse_degree"]), path
    assert int(checks["supersingular_degree"]) == int(
        checks["expected_class_number"]
    ), path
    for name in (
        "squarefree",
        "splits_over_q_2d",
        "recurrence_matches_hasse",
        "all_phi_p_pure_frobenius",
        "all_sage_supersingular",
    ):
        assert bool(checks[name]), (path, name)
    assert sorted(map(int, record["degree_summary"])) == list(range(1, U)), path

    by_degree = defaultdict(list)
    for level in record["levels"]:
        degree = int(level["degree"])
        ell = parse_polynomial(A, level["ell"])
        assert str(ell) == level["ell"], (path, level["ell"])
        assert ell.is_monic() and ell.is_irreducible(), (path, level["ell"])
        assert int(ell.degree()) == degree and ell != p, (path, level["ell"])
        by_degree[degree].append(level)
        matrix = level["brandt_matrix"]
        evidence = level["entry_evidence"]
        assert len(matrix) == len(context.modules), (path, level["ell"])
        assert len(evidence) == len(context.modules), (path, level["ell"])
        zeros = []
        for i, row in enumerate(matrix):
            assert len(row) == len(context.modules), (path, level["ell"])
            assert len(evidence[i]) == len(context.modules), (path, level["ell"])
            for j, value in enumerate(row):
                if value is None:
                    assert evidence[i][j] == "spectral_positive", (path, level["ell"])
                    M = Fraction(q**d - 1, q**2 - 1)
                    C2 = (
                        Fraction(record["weights"][i]) * M - 1
                    ) * (
                        Fraction(record["weights"][j]) * M - 1
                    )
                    Q = q**degree
                    assert Fraction((Q + 1) ** 2) > 4 * C2 * Q, (
                        path,
                        level["ell"],
                        i,
                        j,
                    )
                else:
                    assert evidence[i][j] == "exhaustive_hom", (path, level["ell"])
                    assert int(value) >= 0, (path, level["ell"])
                    if int(value) == 0:
                        zeros.append([i, j])
        assert zeros == level["zero_entries"], (path, level["ell"])
        assert bool(level["complete"]) == (not zeros), (path, level["ell"])
        full = all(value is not None for row in matrix for value in row)
        row_sums = full and all(sum(row) == q**degree + 1 for row in matrix)
        weighted_symmetry = full and all(
            record["weights"][j] * matrix[i][j]
            == record["weights"][i] * matrix[j][i]
            for i in range(len(matrix))
            for j in range(len(matrix))
        )
        archived_checks = level["checks"]
        assert bool(archived_checks["full"]) == full, (path, level["ell"])
        assert bool(archived_checks["row_sums"]) == row_sums, (path, level["ell"])
        assert (
            bool(archived_checks["weighted_symmetry"]) == weighted_symmetry
        ), (path, level["ell"])
        if full:
            assert row_sums and weighted_symmetry, (path, level["ell"])

    bad_degrees = []
    for degree in range(1, U):
        summary = record["degree_summary"][str(degree)]
        expected_levels = {
            str(ell)
            for ell in A.polynomials(of_degree=degree)
            if ell.is_monic() and ell.is_irreducible() and ell != p
        }
        actual_levels = {level["ell"] for level in by_degree[degree]}
        assert len(by_degree[degree]) == len(actual_levels), (path, degree)
        assert actual_levels == expected_levels, (path, degree)
        expected = len(expected_levels)
        assert expected == irreducible_count(q, degree) - int(degree == d)
        assert int(summary["level_count"]) == expected, (path, degree)
        incomplete = sorted(
            level["ell"] for level in by_degree[degree] if not level["complete"]
        )
        assert sorted(summary["incomplete_levels"]) == incomplete, (path, degree)
        assert int(summary["incomplete_count"]) == len(incomplete), (path, degree)
        assert int(summary["complete_count"]) + len(incomplete) == expected, (
            path,
            degree,
        )
        if incomplete:
            bad_degrees.append(degree)
    expected_E = 1 if d <= 2 else 1 + max(bad_degrees, default=0)
    assert int(record["E"]) == expected_E, path
    assert int(record["E"]) <= U, path


def write_csv(path: Path, fieldnames: list[str], rows: list[dict]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("result_dir", type=Path)
    args = parser.parse_args()
    files = sorted(args.result_dir.glob("p_*.json"))
    if not files:
        raise SystemExit(f"no p_*.json files in {args.result_dir}")

    records = []
    for path in files:
        record = json.loads(path.read_text())
        validate(record, path)
        records.append((path, record))
    records.sort(key=lambda item: (item[1]["degree_p"], item[1]["p"]))
    q_values = {int(record["q"]) for _, record in records}
    if len(q_values) != 1:
        raise AssertionError("one report directory must contain one q")
    q = q_values.pop()

    characteristic_rows = []
    level_rows = []
    incomplete_rows = []
    for path, record in records:
        bad = [level for level in record["levels"] if not level["complete"]]
        characteristic_rows.append(
            {
                "q": q,
                "p": record["p"],
                "deg_p": record["degree_p"],
                "class_count": record["class_count"],
                "U": record["theorem_bound_U"],
                "E": record["E"],
                "levels_scanned": len(record["levels"]),
                "incomplete_levels": len(bad),
                "largest_incomplete_degree": max(
                    (level["degree"] for level in bad), default=""
                ),
                "vectors_enumerated": record["total_vectors_enumerated"],
                "result_file": path.name,
            }
        )
        for degree in range(1, int(record["theorem_bound_U"])):
            summary = record["degree_summary"][str(degree)]
            level_rows.append(
                {
                    "q": q,
                    "p": record["p"],
                    "deg_p": record["degree_p"],
                    "deg_ell": degree,
                    "prime_levels": summary["level_count"],
                    "complete": summary["complete_count"],
                    "incomplete": summary["incomplete_count"],
                    "incomplete_levels": "; ".join(summary["incomplete_levels"]),
                }
            )
        for level in bad:
            incomplete_rows.append(
                {
                    "q": q,
                    "p": record["p"],
                    "deg_p": record["degree_p"],
                    "ell": level["ell"],
                    "deg_ell": level["degree"],
                    "zero_entry_count": len(level["zero_entries"]),
                    "zero_entries": json.dumps(level["zero_entries"]),
                }
            )

    characteristic_fields = [
        "q",
        "p",
        "deg_p",
        "class_count",
        "U",
        "E",
        "levels_scanned",
        "incomplete_levels",
        "largest_incomplete_degree",
        "vectors_enumerated",
        "result_file",
    ]
    level_fields = [
        "q",
        "p",
        "deg_p",
        "deg_ell",
        "prime_levels",
        "complete",
        "incomplete",
        "incomplete_levels",
    ]
    incomplete_fields = [
        "q",
        "p",
        "deg_p",
        "ell",
        "deg_ell",
        "zero_entry_count",
        "zero_entries",
    ]
    write_csv(
        args.result_dir / "characteristics.csv",
        characteristic_fields,
        characteristic_rows,
    )
    write_csv(
        args.result_dir / "degree_data.csv",
        level_fields,
        level_rows,
    )
    write_csv(
        args.result_dir / "incomplete_levels.csv",
        incomplete_fields,
        incomplete_rows,
    )

    degree_groups = defaultdict(list)
    for row in characteristic_rows:
        degree_groups[int(row["deg_p"])].append(row)

    lines = [
        f"# Exact Drinfeld completeness data over F_{q}[T]",
        "",
        "Every JSON file was audited by rebuilding its supersingular context, "
        "reconstructing the exact irreducible-level set, recomputing every "
        "spectral inequality, and directly checking zero-entry evidence and "
        "(where fully materialized) Brandt row sums and weighted symmetry.",
        "",
        "## Completeness numbers",
        "",
        "| deg p | characteristics completed | all characteristics | E(p) values | U(p) |",
        "|---:|---:|---:|---|---:|",
    ]
    for degree, rows in sorted(degree_groups.items()):
        expected = irreducible_count(q, degree)
        distribution = Counter(int(row["E"]) for row in rows)
        values = ", ".join(
            f"{value}" + (f" ({count}×)" if count > 1 else "")
            for value, count in sorted(distribution.items())
        )
        cutoffs = sorted({int(row["U"]) for row in rows})
        lines.append(
            f"| {degree} | {len(rows)} | {expected} | {values} | "
            f"{', '.join(map(str, cutoffs))} |"
        )

    lines.extend(
        [
            "",
            "## Per characteristic",
            "",
            "| p | deg p | h(p) | U(p) | E(p) | bad levels | largest bad degree |",
            "|---|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for row in characteristic_rows:
        largest = row["largest_incomplete_degree"] or "—"
        lines.append(
            f"| `{row['p']}` | {row['deg_p']} | {row['class_count']} | "
            f"{row['U']} | **{row['E']}** | {row['incomplete_levels']} | "
            f"{largest} |"
        )

    lines.extend(
        [
            "",
            "## Aggregate behavior by deg p and deg ell",
            "",
            "| deg p | deg ell | characteristics | prime levels tested | complete | incomplete |",
            "|---:|---:|---:|---:|---:|---:|",
        ]
    )
    aggregate = defaultdict(Counter)
    for row in level_rows:
        key = (int(row["deg_p"]), int(row["deg_ell"]))
        aggregate[key]["characteristics"] += 1
        aggregate[key]["prime_levels"] += int(row["prime_levels"])
        aggregate[key]["complete"] += int(row["complete"])
        aggregate[key]["incomplete"] += int(row["incomplete"])
    for (degree_p, degree_ell), counts in sorted(aggregate.items()):
        lines.append(
            f"| {degree_p} | {degree_ell} | {counts['characteristics']} | "
            f"{counts['prime_levels']} | {counts['complete']} | "
            f"{counts['incomplete']} |"
        )

    lines.extend(
        [
            "",
            "The scan is conclusive for every listed characteristic: all prime "
            "levels below U(p) were decided exactly, and Theorem 1.2 proves "
            "completeness for every level degree at least U(p).",
            "",
        ]
    )
    (args.result_dir / "REPORT.md").write_text("\n".join(lines))
    print(
        f"validated {len(records)} characteristics; wrote CSV tables and "
        f"{args.result_dir / 'REPORT.md'}"
    )


if __name__ == "__main__":
    main()
