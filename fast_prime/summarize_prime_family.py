#!/usr/bin/env sage-python
"""Audit an odd-prime compact family and create CSV/Markdown summaries."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

REPOSITORY = Path(__file__).resolve().parents[1]
if str(REPOSITORY) not in sys.path:
    sys.path.insert(0, str(REPOSITORY))

from fast_prime.validate_prime_archive import validate_archive


def _write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def summarize(directory: Path) -> None:
    files = sorted(directory.glob("p_*.compact.json"))
    if not files:
        raise SystemExit(f"no compact archives in {directory}")
    records = []
    audits = []
    for path in files:
        audit = validate_archive(path)
        record = json.loads(path.read_text())
        audits.append(audit)
        records.append(record)
        print(f"validated {path.name}: E={audit['E']}", flush=True)
    q_values = {audit["q"] for audit in audits}
    degrees = {audit["degree_p"] for audit in audits}
    if len(q_values) != 1 or len(degrees) != 1:
        raise ValueError("a family directory must have one q and one deg p")
    prime = next(iter(q_values))
    degree_p = next(iter(degrees))
    cutoff_values = {audit["theorem_bound_U"] for audit in audits}
    if len(cutoff_values) != 1:
        raise ValueError("family cutoff is not constant")
    cutoff = next(iter(cutoff_values))

    characteristic_rows = []
    degree_rows = []
    incomplete_rows = []
    certificate_totals = defaultdict(int)
    for path, record, audit in zip(files, records, audits):
        characteristic_rows.append(
            {
                "q": prime,
                "p": audit["p"],
                "deg_p": degree_p,
                "class_count": audit["class_count"],
                "U": cutoff,
                "E": audit["E"],
                "levels_scanned": audit["level_count"],
                "incomplete_levels": audit["incomplete_level_count"],
                "largest_incomplete_degree": (
                    audit["largest_incomplete_degree"] or ""
                ),
                "unordered_pairs": audit["unordered_pair_count"],
                "pair_level_decisions": audit["pair_level_decisions"],
                "vectors_visited": audit["total_vectors_visited"],
                "result_file": path.name,
            }
        )
        for degree in range(1, cutoff):
            summary = record["degree_summary"][str(degree)]
            degree_rows.append(
                {
                    "q": prime,
                    "p": audit["p"],
                    "deg_p": degree_p,
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
                        "q": prime,
                        "p": audit["p"],
                        "deg_p": degree_p,
                        "ell": level["ell"],
                        "deg_ell": level["degree"],
                        "zero_entry_count": len(level["zero_entries"]),
                        "zero_entries": json.dumps(
                            level["zero_entries"], separators=(",", ":")
                        ),
                    }
                )
        for certificate in record["pair_certificates"]:
            for degree in certificate["spectral_degrees"]:
                certificate_totals["spectral_positive"] += record[
                    "degree_summary"
                ][str(degree)]["level_count"]
            run = certificate["run"]
            if run is not None:
                for _, count in run["counts"]:
                    certificate_totals[
                        "exhaustive_positive" if count else "exhaustive_zero"
                    ] += 1

    _write_csv(directory / "characteristics.csv", characteristic_rows)
    _write_csv(directory / "degree_data.csv", degree_rows)
    _write_csv(directory / "incomplete_levels.csv", incomplete_rows)

    distribution = Counter(row["E"] for row in characteristic_rows)
    aggregate = defaultdict(lambda: {"levels": 0, "complete": 0, "incomplete": 0})
    for row in degree_rows:
        item = aggregate[int(row["deg_ell"])]
        item["levels"] += int(row["prime_levels"])
        item["complete"] += int(row["complete"])
        item["incomplete"] += int(row["incomplete"])
    lines = [
        f"# Exact compact q={prime}, degree-{degree_p} completeness data",
        "",
        f"Validated characteristics: **{len(records)}**.",
        "",
        "Every non-spectral zero follows complete bounded-Hom exhaustion; "
        "the remaining positive entries satisfy the exact pair-specific "
        "spectral inequality.",
        "",
        "## E(p) distribution",
        "",
        "| E(p) | characteristics |",
        "|---:|---:|",
    ]
    for value, count in sorted(distribution.items()):
        lines.append(f"| {value} | {count} |")
    lines.extend(
        [
            "",
            "## Per characteristic",
            "",
            "| p | U(p) | E(p) | incomplete levels | largest incomplete degree |",
            "|---|---:|---:|---:|---:|",
        ]
    )
    for row in characteristic_rows:
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
    total_vectors = sum(
        audit["total_vectors_visited"] for audit in audits
    )
    lines.extend(
        [
            "",
            f"Total exhaustive vector visits: **{total_vectors:,}**.",
            "",
            f"The scan is conclusive below U(p)={cutoff}; Theorem 1.2 proves "
            f"every prime level of degree at least {cutoff} is complete.",
            "",
        ]
    )
    (directory / "REPORT.md").write_text("\n".join(lines))
    print(f"wrote family report and CSV tables in {directory}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("result_directory", type=Path)
    args = parser.parse_args()
    summarize(args.result_directory)


if __name__ == "__main__":
    main()
