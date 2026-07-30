"""Cross-backend regression tests for the compiled odd-prime counter."""

from sage.all import GF, PolynomialRing

from drinfeld_complete.compact_prime import (
    interpolate_prime_form,
    polynomial_code,
    run_prime_counter,
)
from drinfeld_complete.core import (
    build_supersingular_context,
    enumerate_pair_norms_q3,
    monic_irreducibles,
)


def _q3_key_to_code(key: tuple[int, int], max_degree: int) -> int:
    ones, twos = key
    code = 0
    place = 1
    for degree in range(max_degree + 1):
        coefficient = 1 if (ones >> degree) & 1 else 0
        if (twos >> degree) & 1:
            coefficient = 2
        code += coefficient * place
        place *= 3
    return code


def test_compiled_prime_counter_matches_ternary_backend():
    field = GF(3)
    ring = PolynomialRing(field, "T")
    characteristic = ring("T^3 + 2*T + 1")
    context = build_supersingular_context(3, characteristic)
    max_degree = 5
    levels = [
        level
        for degree in range(1, max_degree + 1)
        for level in monic_irreducibles(ring, degree)
        if level != characteristic
    ]
    target_codes = [polynomial_code(level, 3) for level in levels]
    form = interpolate_prime_form(context, 0, 0, max_degree)
    compiled = run_prime_counter(3, form, max_degree, target_codes)
    ternary, _ = enumerate_pair_norms_q3(context, 0, 0, max_degree)
    expected = {
        _q3_key_to_code(key, max_degree): count
        for key, count in ternary.items()
    }
    assert compiled.iterations == 3 ** form["dimension"] - 1
    assert compiled.invalid_norms == 0
    assert compiled.counts == {
        code: expected.get(code, 0) for code in target_codes
    }
