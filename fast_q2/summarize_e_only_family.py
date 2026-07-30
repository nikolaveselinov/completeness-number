#!/usr/bin/env sage-python
"""Audit a complete q=2 E-only family and write compact summary artifacts."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter
from pathlib import Path

REPOSITORY = Path(__file__).resolve().parents[1]
if str(REPOSITORY) not in sys.path:
    sys.path.insert(0, str(REPOSITORY))

from sage.all import GF, PolynomialRing

from compute_family import polynomial_slug
from drinfeld_complete import monic_irreducibles
from fast_q2.replay_e_only_archive import (
    receipt_path,
    sha256,
    validate_replay_receipt,
)
from fast_q2.validate_e_only_archive import validate_e_only_archive


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty CSV {path}")
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def write_sha256_manifest(
    path: Path,
    *,
    directory: Path,
    files: list[Path],
) -> None:
    """Atomically bind the immutable primary archives and replay receipts."""
    relative_files = sorted(
        (
            candidate.relative_to(directory)
            for candidate in files
        ),
        key=lambda candidate: candidate.as_posix(),
    )
    if not relative_files:
        raise ValueError("refusing to write an empty SHA-256 manifest")
    lines = [
        f"{sha256(directory / relative)}  {relative.as_posix()}"
        for relative in relative_files
    ]
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text("\n".join(lines) + "\n")
    temporary.replace(path)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--require-replay-receipts", action="store_true")
    args = parser.parse_args()
    directory = args.directory
    files = sorted(directory.glob("p_*.e-only.json"))
    if not files:
        raise SystemExit(f"no E-only archives in {directory}")

    records = []
    audits = []
    for path in files:
        record = json.loads(path.read_text())
        audit = validate_e_only_archive(record)
        records.append((path, record))
        audits.append(audit)
        print(f"validated {path.name}: E={audit['E']}", flush=True)
    q_values = {audit["q"] for audit in audits}
    degree_values = {audit["degree_p"] for audit in audits}
    if q_values != {2} or len(degree_values) != 1:
        raise ValueError("family must contain one q=2 characteristic degree")
    degree_p = next(iter(degree_values))
    A = PolynomialRing(GF(2), "T")
    expected_polynomials = monic_irreducibles(A, degree_p)
    expected_by_text = {str(p): p for p in expected_polynomials}
    expected = set(expected_by_text)
    actual_counts = Counter(audit["p"] for audit in audits)
    actual = set(actual_counts)
    duplicates = sorted(
        p_text for p_text, count in actual_counts.items() if count != 1
    )
    if actual != expected or duplicates or len(files) != len(expected):
        missing = sorted(expected - actual)
        extra = sorted(actual - expected)
        raise ValueError(
            "family is not a one-to-one characteristic census: "
            f"missing={missing}, unexpected={extra}, duplicates={duplicates}, "
            f"files={len(files)}, expected_files={len(expected)}"
        )
    replay_files = []
    for (path, _record), audit in zip(records, audits):
        expected_name = (
            f"p_{polynomial_slug(expected_by_text[audit['p']])}.e-only.json"
        )
        if path.name != expected_name:
            raise ValueError(
                f"noncanonical archive name {path.name}; "
                f"expected {expected_name}"
            )
        if args.require_replay_receipts:
            replay_path = receipt_path(path)
            if not replay_path.exists():
                raise ValueError(f"missing replay receipt {replay_path}")
            receipt = json.loads(replay_path.read_text())
            validate_replay_receipt(
                receipt,
                archive=path,
                structural_audit=audit,
            )
            replay_files.append(replay_path)
            print(f"verified replay receipt {replay_path.name}", flush=True)

    manifest_files = [path for path, _record in records]
    manifest_files.extend(replay_files)
    write_sha256_manifest(
        directory / "SHA256SUMS",
        directory=directory,
        files=manifest_files,
    )

    rows = []
    for (path, record), audit in sorted(
        zip(records, audits), key=lambda item: item[1]["p"]
    ):
        witness = record["witness"]
        rows.append(
            {
                "q": 2,
                "p": audit["p"],
                "deg_p": degree_p,
                "class_count": audit["class_count"],
                "U": audit["theorem_bound_U"],
                "E": audit["E"],
                "witness_degree": audit["witness_degree"],
                "witness_kind": audit["witness_kind"],
                "witness_pair": json.dumps(witness.get("pair", "")),
                "witness_ell": witness.get("ell", ""),
                "positive_degrees_certified": audit["positive_degrees"],
                "vectors_visited": audit["total_iterations"],
                "result_file": path.name,
            }
        )
    write_csv(directory / "characteristics.csv", rows)

    distribution = Counter(row["E"] for row in rows)
    lines = [
        f"# Exact descending q=2, degree-{degree_p} completeness data",
        "",
        f"Validated characteristics: **{len(rows)}/{len(expected)}**.",
        "",
        "Each archive proves every prime level in every degree strictly above "
        "its witness degree positive for every Brandt entry, using exact "
        "pair-specific spectral inequalities or deterministic bounded-Hom "
        "counter certificates.  It then proves one zero at degree E(p)-1 by "
        "complete Hom-space exhaustion (or by the exact row-sum obstruction).",
        "",
        (
            "Every retained counter run was independently recomputed and "
            "matched its hash-bound replay receipt."
            if args.require_replay_receipts
            else "This report records a structural archive audit; use "
            "`--require-replay-receipts` for release-grade independent replay."
        ),
        "",
        (
            "`SHA256SUMS` binds every immutable primary archive and replay "
            "receipt used by this release audit."
            if args.require_replay_receipts
            else "`SHA256SUMS` binds every immutable primary archive used by "
            "this structural audit."
        ),
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
            "| p | U(p) | E(p) | witness | pair |",
            "|---|---:|---:|---|---|",
        ]
    )
    for row in rows:
        pair = row["witness_pair"] or "—"
        witness = row["witness_ell"] or row["witness_kind"]
        lines.append(
            f"| `{row['p']}` | {row['U']} | **{row['E']}** | "
            f"`{witness}` | `{pair}` |"
        )
    lines.extend(
        [
            "",
            f"The scan is conclusive below each recorded U(p); Theorem 1.2 "
            f"proves every prime level of degree at least U(p) complete.",
            "",
        ]
    )
    (directory / "REPORT.md").write_text("\n".join(lines))
    print(f"wrote family report, CSV, and SHA-256 manifest in {directory}")


if __name__ == "__main__":
    main()
