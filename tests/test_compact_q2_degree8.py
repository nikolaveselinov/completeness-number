"""Structural regression for the separate degree-eight compact family."""

from sage.all import GF, PolynomialRing

from drinfeld_complete import (
    build_supersingular_context,
    monic_irreducibles,
    theorem_bound,
)
from drinfeld_complete.core import pair_is_spectrally_positive
from fast_q2.validate_compact_archive import (
    ARCHIVE_KIND_DEGREE8,
    ARCHIVE_SPECS,
)


def test_degree8_compact_family_parameters():
    A = PolynomialRing(GF(2), "T")
    characteristics = monic_irreducibles(A, 8)
    assert len(characteristics) == 30
    assert theorem_bound(2, 8) == 15
    ctx = build_supersingular_context(2, characteristics[0])
    assert len(ctx.modules) == 85
    assert set(map(int, ctx.automorphism_orders)) == {1}
    assert not any(
        pair_is_spectrally_positive(ctx, 0, 0, degree)
        for degree in range(1, 15)
    )
    assert ARCHIVE_SPECS[ARCHIVE_KIND_DEGREE8] == {
        "degree_p": 8,
        "cutoff": 15,
        "level_count": 2537,
        "vertex_count": 85,
        "unordered_pair_count": 3655,
    }
