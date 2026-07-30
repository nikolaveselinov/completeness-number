#!/usr/bin/env sage-python
"""Audit compact degree-seven archives and write family tables/report."""

from __future__ import annotations

import csv
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

REPOSITORY = Path(__file__).resolve().parents[1]
if str(REPOSITORY) not in sys.path:
    sys.path.insert(0, str(REPOSITORY))

from fast_q2.validate_compact_archive import validate_archive


def write_csv(path: Path, fields: list[str], rows: list[dict]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit(f"usage: {sys.argv[0]} RESULT_DIRECTORY")
    directory = Path(sys.argv[1])
    files = sorted(directory.glob("p_*.compact.json"))
    if not files:
        raise SystemExit(f"no compact archives in {directory}")

    records = []
    for path in files:
        record = json.loads(path.read_text())
        audit = validate_archive(record)
        records.append((path, record, audit))
        print(f"validated {path.name}: E={audit['E']}", flush=True)
    records.sort(key=lambda item: item[1]["p"])

    characteristics = []
    degree_rows = []
    incomplete_rows = []
    certificate_totals: Counter[str] = Counter()
    high_degree_exceptions = []
    for path, record, audit in records:
        certificate_totals.update(audit["certificate_counts"])
        characteristics.append(
            {
                "q": 2,
                "p": record["p"],
                "deg_p": 7,
                "class_count": audit["class_count"],
                "U": audit["theorem_bound_U"],
                "E": audit["E"],
                "levels_scanned": audit["level_count"],
                "incomplete_levels": audit["incomplete_level_count"],
                "largest_incomplete_degree": (
                    audit["largest_incomplete_degree"] or ""
                ),
                "unordered_pairs": audit["unordered_pair_count"],
                "pair_level_decisions": audit["pair_level_decisions"],
                "vectors_visited": audit["total_iterations"],
                "result_file": path.name,
            }
        )
        for degree in range(1, 16):
            summary = record["degree_summary"][str(degree)]
            degree_rows.append(
                {
                    "q": 2,
                    "p": record["p"],
                    "deg_p": 7,
                    "deg_ell": degree,
                    "prime_levels": summary["level_count"],
                    "complete": summary["complete_count"],
                    "incomplete": summary["incomplete_count"],
                    "incomplete_levels": "; ".join(
                        summary["incomplete_levels"]
                    ),
                }
            )
        for level in record["levels"]:
            if not level["complete"]:
                incomplete_rows.append(
                    {
                        "q": 2,
                        "p": record["p"],
                        "deg_p": 7,
                        "ell": level["ell"],
                        "deg_ell": level["degree"],
                        "zero_entry_count": len(level["zero_entries"]),
                        "zero_entries": json.dumps(level["zero_entries"]),
                    }
                )
                if int(level["degree"]) >= 12:
                    high_degree_exceptions.append(
                        {
                            "p": record["p"],
                            "ell": level["ell"],
                            "degree": int(level["degree"]),
                            "zero_entries": level["zero_entries"],
                        }
                    )

    write_csv(
        directory / "characteristics.csv",
        list(characteristics[0]),
        characteristics,
    )
    write_csv(directory / "degree_data.csv", list(degree_rows[0]), degree_rows)
    write_csv(
        directory / "incomplete_levels.csv",
        [
            "q",
            "p",
            "deg_p",
            "ell",
            "deg_ell",
            "zero_entry_count",
            "zero_entries",
        ],
        incomplete_rows,
    )

    distribution = Counter(row["E"] for row in characteristics)
    aggregate = defaultdict(Counter)
    for row in degree_rows:
        degree = int(row["deg_ell"])
        aggregate[degree]["levels"] += int(row["prime_levels"])
        aggregate[degree]["complete"] += int(row["complete"])
        aggregate[degree]["incomplete"] += int(row["incomplete"])
    lines = [
        "# Exact compact q=2, degree-seven completeness data",
        "",
        f"Validated characteristics: **{len(records)}/18**.",
        "",
        "Every archive rebuilds the supersingular context and exact 4,719-level "
        "universe. Its 946 unordered-pair certificates cover all 4,464,174 "
        "pair-level decisions by exhaustive Hom traversal, a successful "
        "all-target witness stop, or the exact pairwise spectral inequality. "
        "Weighted Brandt symmetry transports only zero versus positivity.",
        "",
        "## E(p) distribution",
        "",
        "| E(p) | characteristics |",
        "|---:|---:|",
    ]
    for value, count in sorted(distribution.items()):
        lines.append(f"| {value} | {count} |")
    if high_degree_exceptions:
        lines.extend(
            [
                "",
                "## High-degree exceptional levels",
                "",
                "| p | ell | deg ell | exact ordered zero entries |",
                "|---|---|---:|---|",
            ]
        )
        for item in high_degree_exceptions:
            lines.append(
                f"| `{item['p']}` | `{item['ell']}` | {item['degree']} | "
                f"`{json.dumps(item['zero_entries'])}` |"
            )
    lines.extend(
        [
            "",
            "## Per characteristic",
            "",
            "| p | U(p) | E(p) | incomplete levels | largest incomplete degree |",
            "|---|---:|---:|---:|---:|",
        ]
    )
    for row in characteristics:
        largest = row["largest_incomplete_degree"] or "—"
        lines.append(
            f"| `{row['p']}` | {row['U']} | **{row['E']}** | "
            f"{row['incomplete_levels']} | {largest} |"
        )
    lines.extend(
        [
            "",
            "## Aggregate by level degree",
            "",
            "| deg ell | prime levels tested | complete | incomplete |",
            "|---:|---:|---:|---:|",
        ]
    )
    for degree, counts in sorted(aggregate.items()):
        lines.append(
            f"| {degree} | {counts['levels']} | {counts['complete']} | "
            f"{counts['incomplete']} |"
        )
    lines.extend(
        [
            "",
            "## Certificate totals",
            "",
            "| certificate | pair-level decisions |",
            "|---|---:|",
        ]
    )
    for name, count in sorted(certificate_totals.items()):
        lines.append(f"| `{name}` | {count} |")
    lines.extend(
        [
            "",
            "The scan is conclusive below U(p)=16; Theorem 1.2 proves every "
            "prime level of degree at least 16 is complete.",
            "",
        ]
    )
    (directory / "REPORT.md").write_text("\n".join(lines))
    print(f"wrote family report and CSV tables in {directory}")


if __name__ == "__main__":
    main()
