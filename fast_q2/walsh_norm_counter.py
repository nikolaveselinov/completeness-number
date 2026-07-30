"""Exact Boolean–Walsh multiplicities for binary-valued norm forms.

This reference implementation applies when every extension-field coefficient
block in the interpolated raw norm form is zero or one. Under that
hypothesis, each output coefficient is a Boolean quadratic form, and its
complete value distribution follows from its Walsh spectrum. The module is
independent of Sage and of ``q2_norm_counter``.
"""

from __future__ import annotations

from dataclasses import dataclass
from time import perf_counter
from typing import Sequence


class NonBinaryNormFormError(ValueError):
    """Raised when a packed coefficient is not literally in F_2."""


@dataclass(frozen=True)
class WalshMultiplicityResult:
    """Complete exact value distribution of one vector-valued quadratic form."""

    counts: tuple[int, ...]
    dimension: int
    max_degree: int
    character_sums_evaluated: int
    nonzero_character_sums: int
    seconds: float
    includes_zero_vector: bool

    def count(self, code: int) -> int:
        """Return the multiplicity of one polynomial bit code."""
        if code < 0 or code >= len(self.counts):
            raise ValueError(
                f"code {code} is outside 0 <= code < {len(self.counts)}"
            )
        return self.counts[code]


def _monomial_count(dimension: int) -> int:
    return dimension + dimension * (dimension - 1) // 2


def quadratic_character_sum(encoded: int, dimension: int) -> int:
    r"""Evaluate an exact quadratic character sum over ``F_2^dimension``.

    ``encoded`` stores the linear coefficients first, followed by the
    square-free quadratic coefficients in lexicographic pair order
    ``(0,1), (0,2), ..., (dimension-2, dimension-1)``.  Thus this returns

    .. math::

       \sum_{x\in\mathbf F_2^r}(-1)^{q(x)}.

    Pair elimination computes the answer using only integer and bit
    operations; there is no floating-point approximation.
    """
    if dimension < 0:
        raise ValueError("dimension must be nonnegative")
    if encoded < 0:
        raise ValueError("encoded form must be nonnegative")
    monomial_count = _monomial_count(dimension)
    if encoded >> monomial_count:
        raise ValueError("encoded form has coefficients beyond its dimension")

    linear = encoded & ((1 << dimension) - 1)
    rows = [0] * dimension
    position = dimension
    for left in range(dimension):
        for right in range(left + 1, dimension):
            if encoded & (1 << position):
                rows[left] |= 1 << right
                rows[right] |= 1 << left
            position += 1

    active = (1 << dimension) - 1
    constant = 0
    eliminated_pairs = 0
    while True:
        left = -1
        right = -1
        remaining = active
        while remaining:
            low = remaining & -remaining
            candidate = low.bit_length() - 1
            neighbours = rows[candidate] & active
            if neighbours:
                left = candidate
                neighbour = neighbours & -neighbours
                right = neighbour.bit_length() - 1
                break
            remaining ^= low
        if left < 0:
            break

        left_bit = 1 << left
        right_bit = 1 << right
        left_constant = bool(linear & left_bit)
        right_constant = bool(linear & right_bit)
        left_terms = rows[left] & active & ~right_bit
        right_terms = rows[right] & active & ~left_bit

        active &= ~(left_bit | right_bit)
        linear &= active
        remaining = active
        while remaining:
            low = remaining & -remaining
            index = low.bit_length() - 1
            rows[index] &= active
            remaining ^= low

        # Eliminating x,y from
        #   xy + x*A(z) + y*B(z) + C(z)
        # contributes 2*(-1)^(A(z)B(z)+C(z)).
        constant ^= int(left_constant and right_constant)
        if left_constant:
            linear ^= right_terms
        if right_constant:
            linear ^= left_terms

        left_remaining = left_terms
        while left_remaining:
            left_low = left_remaining & -left_remaining
            left_index = left_low.bit_length() - 1
            right_remaining = right_terms
            while right_remaining:
                right_low = right_remaining & -right_remaining
                right_index = right_low.bit_length() - 1
                if left_index == right_index:
                    linear ^= left_low
                else:
                    rows[left_index] ^= right_low
                    rows[right_index] ^= left_low
                right_remaining ^= right_low
            left_remaining ^= left_low
        eliminated_pairs += 1

    if linear & active:
        return 0
    magnitude = 1 << (eliminated_pairs + active.bit_count())
    return -magnitude if constant else magnitude


def _binary_component_forms(
    diagonal: Sequence[int],
    cross: Sequence[Sequence[int]],
    *,
    max_degree: int,
    field_degree: int,
) -> tuple[int, ...]:
    """Transpose packed vector coefficients into scalar Boolean forms."""
    if max_degree < 0:
        raise ValueError("max_degree must be nonnegative")
    if field_degree <= 0:
        raise ValueError("field_degree must be positive")

    dimension = len(diagonal)
    if len(cross) != dimension or any(
        len(row) != dimension for row in cross
    ):
        raise ValueError("cross must be a square dimension-by-dimension matrix")

    packed_width = (max_degree + 1) * field_degree
    field_mask = (1 << field_degree) - 1
    components = [0] * (max_degree + 1)

    def add_packed(value: int, monomial_position: int, label: str) -> None:
        if value < 0:
            raise ValueError(f"{label} must be nonnegative")
        if value >> packed_width:
            raise ValueError(f"{label} exceeds the requested packed width")
        for degree in range(max_degree + 1):
            coefficient = (value >> (degree * field_degree)) & field_mask
            if coefficient not in (0, 1):
                raise NonBinaryNormFormError(
                    f"{label}, T^{degree} coefficient block is "
                    f"{coefficient}, not 0 or 1"
                )
            if coefficient:
                components[degree] |= 1 << monomial_position

    for index, value in enumerate(diagonal):
        add_packed(int(value), index, f"diagonal[{index}]")

    position = dimension
    for left in range(dimension):
        if int(cross[left][left]) != 0:
            raise ValueError(f"cross[{left}][{left}] must be zero")
        for right in range(left + 1, dimension):
            value = int(cross[left][right])
            if value != int(cross[right][left]):
                raise ValueError(
                    f"cross[{left}][{right}] and cross[{right}][{left}] "
                    "must agree"
                )
            add_packed(value, position, f"cross[{left}][{right}]")
            position += 1
    return tuple(components)


def _fwht(values: list[int]) -> None:
    """Apply the unnormalized Walsh-Hadamard transform in place."""
    width = 1
    while width < len(values):
        step = 2 * width
        for start in range(0, len(values), step):
            for offset in range(width):
                left = start + offset
                right = left + width
                left_value = values[left]
                right_value = values[right]
                values[left] = left_value + right_value
                values[right] = left_value - right_value
        width = step


def walsh_multiplicities(
    diagonal: Sequence[int],
    cross: Sequence[Sequence[int]],
    *,
    max_degree: int,
    field_degree: int,
    include_zero_vector: bool = False,
) -> WalshMultiplicityResult:
    r"""Return every exact output multiplicity via Fourier inversion.

    If ``Q: F_2^r -> F_2^(max_degree+1)`` is the binary-valued packed norm
    form, this computes all character sums

    ``sum_x (-1)^(s dot Q(x))``

    and applies one inverse Walsh-Hadamard transform.  By default the zero
    input vector is removed, matching the primary Gray counter.
    """
    started = perf_counter()
    dimension = len(diagonal)
    components = _binary_component_forms(
        diagonal,
        cross,
        max_degree=max_degree,
        field_degree=field_degree,
    )
    output_dimension = max_degree + 1
    output_count = 1 << output_dimension

    spectrum = [0] * output_count
    encoded_scalar_form = 0
    spectrum[0] = 1 << dimension
    for ordinal in range(1, output_count):
        flip = (ordinal & -ordinal).bit_length() - 1
        encoded_scalar_form ^= components[flip]
        character = quadratic_character_sum(encoded_scalar_form, dimension)
        gray_code = ordinal ^ (ordinal >> 1)
        spectrum[gray_code] = character

    nonzero_character_sums = sum(value != 0 for value in spectrum)
    _fwht(spectrum)
    for code, value in enumerate(spectrum):
        quotient, remainder = divmod(value, output_count)
        if remainder:
            raise ArithmeticError(
                f"Fourier inversion was nonintegral at code {code}"
            )
        if quotient < 0:
            raise ArithmeticError(
                f"Fourier inversion was negative at code {code}"
            )
        spectrum[code] = quotient

    if not include_zero_vector:
        if spectrum[0] < 1:
            raise ArithmeticError("zero input was missing from the zero fibre")
        spectrum[0] -= 1

    expected_total = (1 << dimension) - int(not include_zero_vector)
    if sum(spectrum) != expected_total:
        raise ArithmeticError("Walsh multiplicities have the wrong total")

    return WalshMultiplicityResult(
        counts=tuple(spectrum),
        dimension=dimension,
        max_degree=max_degree,
        character_sums_evaluated=output_count,
        nonzero_character_sums=nonzero_character_sums,
        seconds=perf_counter() - started,
        includes_zero_vector=include_zero_vector,
    )
