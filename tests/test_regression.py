"""End-to-end regression checks from independently computed Brandt data."""

from collections import Counter
import itertools

from sage.all import GF, Hom, PolynomialRing

from drinfeld_complete import (
    build_supersingular_context,
    compute_completeness,
    monic_irreducibles,
)
from drinfeld_complete.core import (
    _poly_key_q4,
    _poly_key_q8,
    enumerate_pair_norms_q4,
    enumerate_pair_norms_q8,
)


def run_case(q, p_text, expected_E, expected_bad, *, full_matrices=True):
    Fq = GF(q, name="a")
    A = PolynomialRing(Fq, "T")
    T = A.gen()
    p = A(eval(p_text.replace("^", "**"), {"T": T}))
    ctx = build_supersingular_context(q, p)
    result = compute_completeness(ctx, full_matrices=full_matrices)
    bad = {
        level["ell"]
        for level in result["levels"]
        if not level["complete"]
    }
    assert result["E"] == expected_E
    assert bad == expected_bad
    return result


def test_degree_3_characteristic():
    result = run_case(
        2,
        "T^3 + T + 1",
        6,
        {
            "T",
            "T + 1",
            "T^2 + T + 1",
            "T^3 + T^2 + 1",
            "T^5 + T^3 + T^2 + T + 1",
        },
    )
    exceptional = next(
        level
        for level in result["levels"]
        if level["ell"] == "T^5 + T^3 + T^2 + T + 1"
    )
    # The two nonzero j-roots may be returned in either finite-field order.
    assert exceptional["brandt_matrix"] in (
        [
            [0, 15, 18],
            [5, 12, 16],
            [6, 16, 11],
        ],
        [
            [0, 18, 15],
            [6, 11, 16],
            [5, 16, 12],
        ],
    )


def test_degree_4_characteristic():
    run_case(
        2,
        "T^4 + T + 1",
        3,
        {"T", "T + 1", "T^2 + T + 1"},
    )


def test_q3_degree_3_spectral_backend():
    result = run_case(
        3,
        "T^3 + 2*T + 1",
        4,
        {
            "T",
            "T + 1",
            "T + 2",
            "T^2 + 1",
            "T^2 + T + 2",
            "T^2 + 2*T + 2",
            "T^3 + 2*T^2 + 1",
            "T^3 + 2*T^2 + T + 1",
            "T^3 + 2*T^2 + 2*T + 2",
        },
        full_matrices=False,
    )
    assert result["theorem_bound_U"] == 6
    assert result["total_vectors_enumerated"] == 99_128
    assert {
        item["packing_field_degree"] for item in result["enumerations"]
    } == {1}
    assert any(
        evidence == "spectral_positive"
        for level in result["levels"]
        for row in level["entry_evidence"]
        for evidence in row
    )


def test_q3_degree_4_full_brandt_backend():
    result = run_case(
        3,
        "T^4 + T + 2",
        5,
        {
            "T",
            "T + 1",
            "T + 2",
            "T^2 + 1",
            "T^2 + T + 2",
            "T^2 + 2*T + 2",
            "T^3 + 2*T + 1",
            "T^3 + 2*T + 2",
            "T^3 + T^2 + 2",
            "T^3 + T^2 + 2*T + 1",
            "T^3 + 2*T^2 + 1",
            "T^3 + 2*T^2 + 2*T + 2",
            "T^4 + T^2 + 2",
            "T^4 + T^2 + T + 1",
            "T^4 + 2*T^2 + 2",
        },
        full_matrices=False,
    )
    assert result["theorem_bound_U"] == 6
    assert result["total_vectors_enumerated"] == 1_968_200
    assert all(
        level["checks"]
        == {
            "full": True,
            "row_sums": True,
            "weighted_symmetry": True,
        }
        for level in result["levels"]
    )


def test_nonprime_base_field_splitting_constructor():
    Fq = GF(4, name="a")
    A = PolynomialRing(Fq, "T")
    p = monic_irreducibles(A, 3)[0]
    ctx = build_supersingular_context(4, p)
    assert len(ctx.modules) == 5
    assert int(ctx.splitting_field.cardinality()) == 4**6
    assert all(
        [
            ctx.checks["squarefree"],
            ctx.checks["splits_over_q_2d"],
            ctx.checks["recurrence_matches_hasse"],
            ctx.checks["all_phi_p_pure_frobenius"],
            ctx.checks["all_sage_supersingular"],
        ]
    )


def test_theorem_only_low_degree_archive_needs_no_exhaustive_backend():
    Fq = GF(5)
    A = PolynomialRing(Fq, "T")
    p = monic_irreducibles(A, 2)[0]
    ctx = build_supersingular_context(5, p)
    result = compute_completeness(ctx)
    assert result["theorem_bound_U"] == 1
    assert result["E"] == 1
    assert result["levels"] == []
    assert result["enumerations"] == []
    assert result["total_vectors_enumerated"] == 0


def test_q4_exhaustive_backend_and_direct_norm_regression():
    Fq = GF(4, name="a")
    a = Fq.gen()
    A = PolynomialRing(Fq, "T")
    T = A.gen()
    p = A(T**3 + a)
    ctx = build_supersingular_context(4, p)

    # Compare every vector in a small bounded Hom space against Sage's direct
    # norm, including all nonbinary F_4 coefficients and cross terms.
    source, target, max_degree = 1, 2, 2
    homset = Hom(ctx.modules[source], ctx.modules[target])
    basis = homset.basis(degree=max_degree)
    field_values = [Fq(0), Fq(1), a, a + 1]
    direct_counts = Counter()
    for coefficients in itertools.product(field_values, repeat=len(basis)):
        if all(not coefficient for coefficient in coefficients):
            continue
        ore = ctx.modules[source].ore_polring().zero()
        for coefficient, basis_element in zip(coefficients, basis):
            ore += coefficient * basis_element.ore_polynomial()
        key = _poly_key_q4(homset(ore).norm().gen())
        direct_counts[key] += 1
    packed_counts, metadata = enumerate_pair_norms_q4(
        ctx, source, target, max_degree, cross_check_samples=32
    )
    assert packed_counts == direct_counts
    assert metadata["vectors_enumerated"] == 4**len(basis) - 1 == 255

    # Exercise the complete Appendix A scan, including spectral certificates.
    result = compute_completeness(ctx)
    assert result["theorem_bound_U"] == 6
    assert result["E"] == 4
    assert result["total_vectors_enumerated"] == 1_085_415
    assert {
        degree: (
            summary["level_count"],
            summary["complete_count"],
            summary["incomplete_count"],
        )
        for degree, summary in result["degree_summary"].items()
    } == {
        "1": (4, 0, 4),
        "2": (6, 0, 6),
        "3": (19, 7, 12),
        "4": (60, 60, 0),
        "5": (204, 204, 0),
    }


def test_q8_exhaustive_backend_matches_every_direct_norm():
    Fq = GF(8, name="a")
    a = Fq.gen()
    A = PolynomialRing(Fq, "T")
    p = monic_irreducibles(A, 3)[0]
    ctx = build_supersingular_context(8, p)
    source, target, max_degree = 1, 2, 2
    homset = Hom(ctx.modules[source], ctx.modules[target])
    basis = homset.basis(degree=max_degree)
    field_values = []
    for code in range(8):
        value = Fq.zero()
        for exponent in range(3):
            if code & (1 << exponent):
                value += a**exponent
        field_values.append(value)

    direct_counts = Counter()
    for coefficients in itertools.product(field_values, repeat=len(basis)):
        if all(not coefficient for coefficient in coefficients):
            continue
        ore = ctx.modules[source].ore_polring().zero()
        for coefficient, basis_element in zip(coefficients, basis):
            ore += coefficient * basis_element.ore_polynomial()
        direct_counts[_poly_key_q8(homset(ore).norm().gen())] += 1

    packed_counts, metadata = enumerate_pair_norms_q8(
        ctx, source, target, max_degree, cross_check_samples=64
    )
    assert packed_counts == direct_counts
    assert metadata["vectors_enumerated"] == 8**len(basis) - 1 == 4095
    assert metadata["gray_radix"] == 8
