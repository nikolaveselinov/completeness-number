"""Exact checks for the isolated optimized C++ Boolean-Walsh counter."""

import pytest
from sage.all import GF, PolynomialRing

from drinfeld_complete import build_supersingular_context, monic_irreducibles
from fast_q2.compare_walsh_gray import run_read_only_gray_counter
from fast_q2.gray_counter import BINARY as GRAY_BINARY
from fast_q2.gray_counter import SOURCE as GRAY_SOURCE
from fast_q2.gray_counter import interpolate_pair
from fast_q2.q2_walsh_cpp import (
    BINARY as WALSH_BINARY,
    CppWalshCounterError,
    run_walsh_counter,
)
from fast_q2.walsh_norm_counter import walsh_multiplicities


def _gray_state() -> tuple[int, int]:
    if not GRAY_BINARY.exists():
        pytest.skip("read-only Gray comparison requires the existing binary")
    if GRAY_BINARY.stat().st_mtime < GRAY_SOURCE.stat().st_mtime:
        pytest.skip("existing Gray binary is older than its source")
    return GRAY_BINARY.stat().st_mtime_ns, GRAY_BINARY.stat().st_size


def _symmetric_cross(dimension: int) -> list[list[int]]:
    return [[0] * dimension for _ in range(dimension)]


def _check_all_three(
    diagonal: list[int],
    cross: list[list[int]],
    *,
    max_degree: int,
    field_degree: int,
) -> None:
    gray_before = _gray_state()
    cpp = run_walsh_counter(
        diagonal,
        cross,
        max_degree=max_degree,
        field_degree=field_degree,
    )
    python = walsh_multiplicities(
        diagonal,
        cross,
        max_degree=max_degree,
        field_degree=field_degree,
    )
    target_codes = list(range(1, 1 << (max_degree + 1)))
    gray = run_read_only_gray_counter(
        diagonal,
        cross,
        max_degree=max_degree,
        field_degree=field_degree,
        target_codes=target_codes,
    )

    assert GRAY_BINARY.exists()
    assert (
        GRAY_BINARY.stat().st_mtime_ns,
        GRAY_BINARY.stat().st_size,
    ) == gray_before
    assert cpp.counts == python.counts
    assert cpp.count(0) == 0
    assert gray.exhaustive
    assert gray.invalid_norms == 0
    assert {
        code: cpp.count(code) for code in target_codes
    } == gray.counts


def test_cpp_walsh_matches_python_and_gray_on_synthetic_form():
    dimension = 4
    max_degree = 4
    diagonal = [1 << index for index in range(dimension)]
    cross = _symmetric_cross(dimension)
    for left, right in [(0, 1), (0, 3), (1, 2), (2, 3)]:
        cross[left][right] = 1 << max_degree
        cross[right][left] = 1 << max_degree
    _check_all_three(
        diagonal,
        cross,
        max_degree=max_degree,
        field_degree=1,
    )


@pytest.mark.parametrize(
    ("degree_p", "max_degree"),
    [(3, 5), (7, 10)],
)
def test_cpp_walsh_matches_python_and_gray_on_real_form(
    degree_p, max_degree
):
    A = PolynomialRing(GF(2), "T")
    characteristic = monic_irreducibles(A, degree_p)[0]
    ctx = build_supersingular_context(2, characteristic)
    diagonal, cross = interpolate_pair(ctx, 0, 1, max_degree)
    _check_all_three(
        diagonal,
        cross,
        max_degree=max_degree,
        field_degree=2 * degree_p,
    )


def test_cpp_walsh_includes_zero_vector_on_request():
    result = run_walsh_counter(
        [1],
        [[0]],
        max_degree=0,
        field_degree=1,
        include_zero_vector=True,
    )
    assert result.counts == (1, 1)
    assert sum(result.counts) == 2


def test_cpp_walsh_rejects_nonbinary_extension_coefficient():
    with pytest.raises(
        CppWalshCounterError, match="not literally 0 or 1"
    ):
        run_walsh_counter(
            [2],
            [[0]],
            max_degree=0,
            field_degree=3,
        )


def test_cpp_build_does_not_change_gray_binary():
    before = _gray_state()
    assert WALSH_BINARY.exists()
    assert (
        GRAY_BINARY.stat().st_mtime_ns,
        GRAY_BINARY.stat().st_size,
    ) == before
