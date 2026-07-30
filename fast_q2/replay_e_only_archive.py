#!/usr/bin/env sage-python
"""Independently replay one q=2 E-only archive and write a hash-bound receipt."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import uuid
from pathlib import Path
from typing import Any

REPOSITORY = Path(__file__).resolve().parents[1]
if str(REPOSITORY) not in sys.path:
    sys.path.insert(0, str(REPOSITORY))

from sage.all import GF, PolynomialRing

from compute_family import polynomial_slug
from drinfeld_complete import monic_irreducibles
from fast_q2.validate_e_only_archive import validate_e_only_archive


REPLAY_RECEIPT_KIND = "q2_e_only_replay_receipt_v1"
VALIDATOR_SOURCE = REPOSITORY / "fast_q2" / "validate_e_only_archive.py"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def receipt_path(archive: Path) -> Path:
    suffix = ".e-only.json"
    if not archive.name.endswith(suffix):
        raise ValueError(f"archive name must end in {suffix}")
    stem = archive.name[: -len(suffix)]
    return archive.parent / "replay_receipts" / f"{stem}.replay.json"


def validate_replay_receipt(
    receipt_value: Any,
    *,
    archive: Path,
    structural_audit: dict[str, Any],
) -> dict[str, Any]:
    if not isinstance(receipt_value, dict):
        raise ValueError("replay receipt must be a JSON object")
    expected_keys = {
        "receipt_kind",
        "schema_version",
        "archive_filename",
        "archive_sha256",
        "replay_driver_sha256",
        "validator_source_sha256",
        "audit",
    }
    if set(receipt_value) != expected_keys:
        raise ValueError(
            "replay receipt keys differ: "
            f"expected={sorted(expected_keys)}, actual={sorted(receipt_value)}"
        )
    if receipt_value["receipt_kind"] != REPLAY_RECEIPT_KIND:
        raise ValueError("unexpected replay receipt kind")
    if receipt_value["schema_version"] != 1:
        raise ValueError("unexpected replay receipt schema version")
    if receipt_value["archive_filename"] != archive.name:
        raise ValueError("replay receipt names the wrong archive")
    if receipt_value["archive_sha256"] != sha256(archive):
        raise ValueError("archive changed after the independent replay")
    replay_audit = receipt_value["audit"]
    if not isinstance(replay_audit, dict):
        raise ValueError("replay receipt audit must be a JSON object")
    if replay_audit.get("recompute") is not True:
        raise ValueError("receipt does not record an independent recomputation")
    if int(replay_audit.get("replayed_runs", 0)) < 1:
        raise ValueError("receipt did not replay any compiled-counter run")
    for key in (
        "q",
        "p",
        "degree_p",
        "theorem_bound_U",
        "E",
        "class_count",
        "positive_degrees",
        "witness_kind",
        "witness_degree",
        "total_iterations",
        "certificate_iterations",
        "discarded_witness_runs",
    ):
        if replay_audit.get(key) != structural_audit.get(key):
            raise ValueError(f"replay receipt disagrees on {key}")
    return receipt_value


def _atomic_write(path: Path, value: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(
        f"{path.name}.{os.getpid()}.{uuid.uuid4().hex}.tmp"
    )
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--degree", type=int, required=True)
    parser.add_argument(
        "--index",
        type=int,
        help="zero-based irreducible index (defaults to SLURM_ARRAY_TASK_ID)",
    )
    parser.add_argument("--result-dir", type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    index = args.index
    if index is None:
        if "SLURM_ARRAY_TASK_ID" not in os.environ:
            raise ValueError("--index or SLURM_ARRAY_TASK_ID is required")
        index = int(os.environ["SLURM_ARRAY_TASK_ID"])
    A = PolynomialRing(GF(2), "T")
    characteristics = monic_irreducibles(A, args.degree)
    if not 0 <= index < len(characteristics):
        raise IndexError(
            f"index {index} is outside 0..{len(characteristics) - 1}"
        )
    p = characteristics[index]
    directory = args.result_dir or Path(
        f"results/local/q2_degree{args.degree}_e_only"
    )
    archive = directory / f"p_{polynomial_slug(p)}.e-only.json"
    if not archive.exists():
        raise FileNotFoundError(f"missing archive {archive}")
    record = json.loads(archive.read_text())
    audit = validate_e_only_archive(record, recompute=True)
    receipt = {
        "receipt_kind": REPLAY_RECEIPT_KIND,
        "schema_version": 1,
        "archive_filename": archive.name,
        "archive_sha256": sha256(archive),
        "replay_driver_sha256": sha256(Path(__file__)),
        "validator_source_sha256": sha256(VALIDATOR_SOURCE),
        "audit": audit,
    }
    destination = receipt_path(archive)
    _atomic_write(destination, receipt)
    print(
        f"independently replayed index={index}, p={p}, "
        f"runs={audit['replayed_runs']} -> {destination}",
        flush=True,
    )


if __name__ == "__main__":
    main()
