"""Python wrapper for the independent exact C++ Boolean–Walsh counter."""

from __future__ import annotations

import math
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence


HERE = Path(__file__).resolve().parent
SOURCE = HERE / "q2_walsh_counter.cpp"
BINARY = HERE / "q2_walsh_counter"


class CppWalshCounterError(RuntimeError):
    """Raised when the isolated C++ helper rejects a form or its output."""


@dataclass(frozen=True)
class CppWalshMultiplicityResult:
    """Complete exact value distribution returned by the C++ helper."""

    counts: tuple[int, ...]
    dimension: int
    max_degree: int
    character_sums_evaluated: int
    nonzero_character_sums: int
    seconds: float
    includes_zero_vector: bool

    def count(self, code: int) -> int:
        if code < 0 or code >= len(self.counts):
            raise ValueError(
                f"code {code} is outside 0 <= code < {len(self.counts)}"
            )
        return self.counts[code]


def compile_walsh_counter() -> None:
    """Build the isolated helper directly, without touching any Make target."""
    if BINARY.exists() and BINARY.stat().st_mtime >= SOURCE.stat().st_mtime:
        return
    subprocess.run(
        [
            "g++",
            "-O3",
            "-std=c++20",
            "-DNDEBUG",
            str(SOURCE),
            "-o",
            str(BINARY),
        ],
        check=True,
    )


def _words(value: int, count: int) -> list[str]:
    if int(value) < 0:
        raise ValueError("packed coefficients must be nonnegative")
    mask = (1 << 64) - 1
    return [
        hex((int(value) >> (64 * index)) & mask)
        for index in range(count)
    ]


def run_walsh_counter(
    diagonal: Sequence[int],
    cross: Sequence[Sequence[int]],
    *,
    max_degree: int,
    field_degree: int,
    include_zero_vector: bool = False,
) -> CppWalshMultiplicityResult:
    """Return all exact fibres of one validated binary-valued norm form."""
    compile_walsh_counter()
    dimension = len(diagonal)
    if len(cross) != dimension or any(
        len(row) != dimension for row in cross
    ):
        raise ValueError("cross must be a square dimension-by-dimension matrix")
    for index in range(dimension):
        if int(cross[index][index]) != 0:
            raise ValueError(f"cross[{index}][{index}] must be zero")
        for other in range(index + 1, dimension):
            if int(cross[index][other]) != int(cross[other][index]):
                raise ValueError(
                    f"cross[{index}][{other}] and "
                    f"cross[{other}][{index}] must agree"
                )
    if max_degree < 0:
        raise ValueError("max_degree must be nonnegative")
    if field_degree <= 0:
        raise ValueError("field_degree must be positive")

    word_count = math.ceil((max_degree + 1) * field_degree / 64)
    tokens = [
        str(dimension),
        str(max_degree),
        str(field_degree),
        str(word_count),
        str(int(include_zero_vector)),
    ]
    for value in diagonal:
        tokens.extend(_words(int(value), word_count))
    for left in range(dimension):
        for right in range(left + 1, dimension):
            tokens.extend(_words(int(cross[left][right]), word_count))

    completed = subprocess.run(
        [str(BINARY)],
        input=" ".join(tokens),
        text=True,
        capture_output=True,
    )
    if completed.returncode:
        message = completed.stderr.strip() or (
            f"helper exited with status {completed.returncode}"
        )
        raise CppWalshCounterError(message)

    lines = completed.stdout.splitlines()
    if not lines:
        raise CppWalshCounterError("the C++ Walsh counter returned no output")
    header = lines[0].split()
    if header[0] != "SUMMARY" or len(header) != 7:
        raise CppWalshCounterError(
            f"unexpected helper output: {lines[0]!r}"
        )
    reported_dimension = int(header[1])
    reported_max_degree = int(header[2])
    output_count = 1 << (max_degree + 1)
    if reported_dimension != dimension or reported_max_degree != max_degree:
        raise CppWalshCounterError("helper summary changed the input dimensions")
    if int(header[3]) != output_count:
        raise CppWalshCounterError("helper evaluated the wrong code space")
    if bool(int(header[5])) != bool(include_zero_vector):
        raise CppWalshCounterError("helper changed include_zero_vector")
    if len(lines) != output_count + 1:
        raise CppWalshCounterError(
            f"expected {output_count} COUNT lines, got {len(lines) - 1}"
        )

    counts = [0] * output_count
    for expected_code, line in enumerate(lines[1:]):
        fields = line.split()
        if (
            len(fields) != 3
            or fields[0] != "COUNT"
            or int(fields[1]) != expected_code
        ):
            raise CppWalshCounterError(
                f"unexpected helper output: {line!r}"
            )
        counts[expected_code] = int(fields[2])

    expected_total = (1 << dimension) - int(not include_zero_vector)
    if sum(counts) != expected_total:
        raise CppWalshCounterError("helper counts have the wrong total")
    return CppWalshMultiplicityResult(
        counts=tuple(counts),
        dimension=dimension,
        max_degree=max_degree,
        character_sums_evaluated=int(header[3]),
        nonzero_character_sums=int(header[4]),
        seconds=float(header[6]),
        includes_zero_vector=bool(include_zero_vector),
    )
