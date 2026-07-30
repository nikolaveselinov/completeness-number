"""Small invariants for the compact characteristic-two engine."""

from copy import deepcopy

import pytest

from drinfeld_complete.compact_q2 import _deterministic_sample_vectors
from drinfeld_complete.e_only_q2 import compute_e_only_q2, e_only_pair_order
from drinfeld_complete import (
    build_supersingular_context,
    monic_irreducibles,
    theorem_bound,
)
from fast_q2.gray_counter import run_counter
from fast_q2.validate_compact_archive import ArchiveValidationError
from fast_q2.validate_e_only_archive import validate_e_only_archive
from sage.all import GF, PolynomialRing


def test_direct_norm_samples_exercise_cross_terms():
    samples = _deterministic_sample_vectors(16, 8)
    assert len(samples) == 8
    assert len(samples) == len(set(samples))
    assert any(vector.bit_count() == 1 for vector in samples)
    assert any(vector.bit_count() > 1 for vector in samples)


def test_compiled_counter_accepts_degree_nine_packing():
    """Exercise a coherent coefficient beyond the former four-word limit."""
    max_degree = 19
    field_degree = 18
    target_code = 1 << max_degree
    result = run_counter(
        [1 << (max_degree * field_degree)],
        [[0]],
        max_degree=max_degree,
        field_degree=field_degree,
        target_codes=[target_code],
        stop_when_seen=False,
    )
    assert result.exhaustive
    assert result.iterations == 1
    assert result.seen == 1
    assert result.invalid_norms == 0
    assert result.counts == {target_code: 1}


def test_degree_nine_e_only_spectral_tiers():
    A = PolynomialRing(GF(2), "T")
    characteristics = monic_irreducibles(A, 9)
    assert len(characteristics) == 56
    assert theorem_bound(2, 9) == 20
    ctx = build_supersingular_context(2, characteristics[0])
    assert len(ctx.modules) == 171
    assert [
        (index, int(weight))
        for index, weight in enumerate(ctx.automorphism_orders)
        if int(weight) != 1
    ] == [(0, 3)]
    assert e_only_pair_order(ctx, 19) == [(0, 0)]
    assert len(e_only_pair_order(ctx, 18)) == 171
    assert len(e_only_pair_order(ctx, 17)) == 171
    assert len(e_only_pair_order(ctx, 16)) == 14706
    assert e_only_pair_order(ctx, 16)[:3] == [(0, 0), (1, 1), (2, 2)]


def test_e_only_work_counter_includes_discarded_witness_attempts(tmp_path):
    """The work total may exceed the retained mathematical certificate."""
    A = PolynomialRing(GF(2), "T")
    p = A("T^4 + T + 1")
    ctx = build_supersingular_context(2, p)
    record = compute_e_only_q2(
        ctx,
        checkpoint_path=tmp_path / "checkpoint.json",
    )
    audit = validate_e_only_archive(record)
    assert audit["E"] == 3
    assert audit["discarded_witness_runs"] == 1
    assert audit["total_iterations"] > audit["certificate_iterations"]

    tampered = deepcopy(record)
    tampered["total_vectors_visited"] = audit["certificate_iterations"]
    with pytest.raises(ArchiveValidationError, match="total_vectors_visited"):
        validate_e_only_archive(tampered)

    tampered = deepcopy(record)
    tampered["engine"] = {"fabricated": True}
    with pytest.raises(ArchiveValidationError, match="engine"):
        validate_e_only_archive(tampered)
