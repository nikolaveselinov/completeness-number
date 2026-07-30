#!/usr/bin/env sage-python
"""Compare standard ordered-matrix archives with compact unordered archives."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def compare(standard_directory: Path, compact_directory: Path) -> dict:
    standard_files = sorted(standard_directory.glob("p_*.json"))
    compact_files = sorted(compact_directory.glob("p_*.compact.json"))
    standard = {path.stem: path for path in standard_files}
    compact = {
        path.name.removesuffix(".compact.json"): path
        for path in compact_files
    }
    if set(standard) != set(compact):
        missing_standard = sorted(set(compact) - set(standard))
        missing_compact = sorted(set(standard) - set(compact))
        raise ValueError(
            f"family mismatch: missing standard={missing_standard}, "
            f"missing compact={missing_compact}"
        )

    level_comparisons = 0
    zero_entry_comparisons = 0
    for slug in sorted(standard):
        standard_path = standard[slug]
        compact_path = compact[slug]
        ordered = json.loads(standard_path.read_text())
        unordered = json.loads(compact_path.read_text())
        characteristic = ordered["p"]
        for field in (
            "q",
            "p",
            "degree_p",
            "theorem_bound_U",
            "E",
            "class_count",
            "weights",
            "j_invariants",
        ):
            if ordered[field] != unordered[field]:
                raise ValueError(
                    f"{characteristic}: {field} differs between "
                    f"{standard_path.name} and {compact_path.name}"
                )
        ordered_levels = {item["ell"]: item for item in ordered["levels"]}
        compact_levels = {item["ell"]: item for item in unordered["levels"]}
        if set(ordered_levels) != set(compact_levels):
            raise ValueError(f"{characteristic}: level universes differ")
        for level in ordered_levels:
            first = ordered_levels[level]
            second = compact_levels[level]
            if first["degree"] != second["degree"]:
                raise ValueError(
                    f"{characteristic}, {level}: degree differs"
                )
            if first["complete"] != second["complete"]:
                raise ValueError(
                    f"{characteristic}, {level}: completeness differs"
                )
            if first["zero_entries"] != second["zero_entries"]:
                raise ValueError(
                    f"{characteristic}, {level}: zero positions differ"
                )
            level_comparisons += 1
            zero_entry_comparisons += len(first["zero_entries"])
    return {
        "characteristics": len(standard),
        "levels_compared": level_comparisons,
        "zero_entry_lists_compared": level_comparisons,
        "zero_entries_compared": zero_entry_comparisons,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("standard_directory", type=Path)
    parser.add_argument("compact_directory", type=Path)
    args = parser.parse_args()
    print(
        json.dumps(
            compare(args.standard_directory, args.compact_directory),
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
