#!/usr/bin/env sage-python
"""Validate one compact exact q=2, degree-eight archive."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPOSITORY = Path(__file__).resolve().parents[1]
if str(REPOSITORY) not in sys.path:
    sys.path.insert(0, str(REPOSITORY))

from fast_q2.validate_compact_archive import (
    ARCHIVE_KIND_DEGREE8,
    ArchiveValidationError,
    validate_archive,
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    parser.add_argument("--recompute", action="store_true")
    args = parser.parse_args()
    try:
        record = json.loads(args.archive.read_text())
        if record.get("archive_kind") != ARCHIVE_KIND_DEGREE8:
            raise ArchiveValidationError(
                f"archive_kind must be {ARCHIVE_KIND_DEGREE8!r}"
            )
        summary = validate_archive(record, recompute=args.recompute)
    except (OSError, json.JSONDecodeError, ArchiveValidationError) as error:
        raise SystemExit(f"validation failed: {error}") from error
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
