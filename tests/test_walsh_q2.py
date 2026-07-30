"""Exact checks for the experimental Boolean-Walsh q=2 counter."""

import pytest
from sage.all import GF, PolynomialRing

from drinfeld_complete import build_supersingular_context, monic_irreducibles
from fast_q2.compare_walsh_gray import run_read_only_gray_counter
from fast_q2.gray_counter import (
    BINARY,
    SOURCE,
    interpolate_pair,
)
from fast_q2.walsh_norm_counter import (
    NonBinaryNormFormError,
    quadratic_character_sum,
    walsh_multiplicities,
)


def _brute_character_sum(encoded: int, dimension: int) -> int:
    total = 0
    for vector in range(1 << dimension):
        value = (encoded & vector).bit_count() & 1
        position = dimension
        for left in range(dimension):
            for right in range(left + 1, dimension):
                value ^= (
                    ((encoded >> position) & 1)
                    & ((vector >> left) & 1)
                    & ((vector >> right) & 1)
                )
                position += 1
        total += -1 if value else 1
    return total


def _require_read_only_gray_binary() -> tuple[int, int]:
    """Require the live helper while preserving a cheap outer fingerprint."""
    if not BINARY.exists():
        pytest.skip("read-only Gray comparison requires the existing binary")
    if BINARY.stat().st_mtime < SOURCE.stat().st_mtime:
        pytest.skip("existing Gray binary is older than its source")
    return BINARY.stat().st_mtime_ns, BINARY.stat().st_size


def _assert_gray_binary_unchanged(before: tuple[int, int]) -> None:
    assert BINARY.exists()
    assert (BINARY.stat().st_mtime_ns, BINARY.stat().st_size) == before


def _symmetric_cross(dimension: int) -> list[list[int]]:
    return [[0] * dimension for _ in range(dimension)]


def test_quadratic_character_sum_exhaustive_through_dimension_four():
    for dimension in range(5):
        monomials = dimension + dimension * (dimension - 1) // 2
        for encoded in range(1 << monomials):
            assert quadratic_character_sum(
                encoded, dimension
            ) == _brute_character_sum(encoded, dimension)


def test_walsh_matches_gray_on_synthetic_vector_form():
    before = _require_read_only_gray_binary()
    dimension = 4
    max_degree = 4
    diagonal = [1 << index for index in range(dimension)]
    cross = _symmetric_cross(dimension)
    for left, right in [(0, 1), (0, 3), (1, 2), (2, 3)]:
        cross[left][right] = 1 << max_degree
        cross[right][left] = 1 << max_degree

    walsh = walsh_multiplicities(
        diagonal,
        cross,
        max_degree=max_degree,
        field_degree=1,
    )
    target_codes = list(range(1, 1 << (max_degree + 1)))
    gray = run_read_only_gray_counter(
        diagonal,
        cross,
        max_degree=max_degree,
        field_degree=1,
        target_codes=target_codes,
    )
    _assert_gray_binary_unchanged(before)

    assert gray.exhaustive
    assert gray.invalid_norms == 0
    assert walsh.count(0) == 0
    assert {
        code: walsh.count(code) for code in target_codes
    } == gray.counts


@pytest.mark.parametrize(
    ("degree_p", "max_degree"),
    [(3, 5), (7, 10)],
)
def test_walsh_matches_gray_on_real_norm_form(degree_p, max_degree):
    before = _require_read_only_gray_binary()
    A = PolynomialRing(GF(2), "T")
    characteristic = monic_irreducibles(A, degree_p)[0]
    ctx = build_supersingular_context(2, characteristic)
    diagonal, cross = interpolate_pair(ctx, 0, 1, max_degree)

    walsh = walsh_multiplicities(
        diagonal,
        cross,
        max_degree=max_degree,
        field_degree=2 * degree_p,
    )
    target_codes = list(range(1, 1 << (max_degree + 1)))
    gray = run_read_only_gray_counter(
        diagonal,
        cross,
        max_degree=max_degree,
        field_degree=2 * degree_p,
        target_codes=target_codes,
    )
    _assert_gray_binary_unchanged(before)

    assert gray.exhaustive
    assert gray.invalid_norms == 0
    assert walsh.count(0) == 0
    assert {
        code: walsh.count(code) for code in target_codes
    } == gray.counts


def test_walsh_rejects_nonbinary_extension_coefficient():
    with pytest.raises(NonBinaryNormFormError, match="not 0 or 1"):
        walsh_multiplicities(
            [2],
            [[0]],
            max_degree=0,
            field_degree=3,
        )
