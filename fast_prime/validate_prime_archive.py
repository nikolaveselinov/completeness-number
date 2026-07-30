#!/usr/bin/env sage-python
"""Strict independent validator for odd-prime compact zero archives."""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

from sage.all import GF, PolynomialRing, is_prime, sage_eval

REPOSITORY = Path(__file__).resolve().parents[1]
if str(REPOSITORY) not in sys.path:
    sys.path.insert(0, str(REPOSITORY))

from drinfeld_complete.compact_prime import (
    COMPACT_KIND,
    COUNTER_BINARY,
    COUNTER_SOURCE,
    polynomial_code,
)
from drinfeld_complete.core import (
    build_supersingular_context,
    monic_irreducibles,
    pair_is_spectrally_positive,
    theorem_bound,
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def validate_archive(path: Path) -> dict:
    record = json.loads(path.read_text())
    _require(record["archive_kind"] == COMPACT_KIND, "wrong archive kind")
    _require(record["schema_version"] == 1, "wrong schema version")
    prime = int(record["q"])
    _require(bool(is_prime(prime)) and prime != 2, "q is not an odd prime")
    field = GF(prime)
    ring = PolynomialRing(field, "T")
    T = ring.gen()
    characteristic = ring(
        sage_eval(
            str(record["p"]).replace("^", "**"),
            locals={"T": T},
        )
    )
    _require(characteristic.is_monic(), "p is not monic")
    _require(characteristic.is_irreducible(), "p is not irreducible")
    degree_p = int(characteristic.degree())
    _require(record["degree_p"] == degree_p, "degree_p mismatch")
    cutoff = theorem_bound(prime, degree_p)
    _require(record["theorem_bound_U"] == cutoff, "cutoff mismatch")

    context = build_supersingular_context(prime, characteristic)
    class_count = len(context.modules)
    _require(record["class_count"] == class_count, "class count mismatch")
    expected_aut = [int(value) for value in context.automorphism_orders]
    expected_weights = [
        value // (prime - 1) for value in expected_aut
    ]
    _require(record["automorphism_orders"] == expected_aut, "Aut mismatch")
    _require(record["weights"] == expected_weights, "weight mismatch")
    _require(
        record["j_invariants"]
        == [str(value) for value in context.j_invariants],
        "j-invariant ordering mismatch",
    )

    expected_levels = []
    expected_by_code = {}
    codes_by_degree = {}
    for degree in range(1, cutoff):
        levels = [
            level
            for level in monic_irreducibles(ring, degree)
            if level != characteristic
        ]
        codes = sorted(polynomial_code(level, prime) for level in levels)
        codes_by_degree[degree] = codes
        for level in levels:
            code = polynomial_code(level, prime)
            item = {
                "ell": str(level),
                "code": code,
                "degree": degree,
            }
            expected_levels.append(item)
            expected_by_code[code] = item
    expected_levels.sort(key=lambda item: (item["degree"], item["code"]))

    pair_total = class_count * (class_count + 1) // 2
    _require(
        len(record["pair_certificates"]) == pair_total,
        "pair certificate count mismatch",
    )
    seen_pairs = set()
    zero_pairs_by_code = {code: [] for code in expected_by_code}
    certificate_decisions = 0
    vector_total = 0
    for certificate in record["pair_certificates"]:
        pair = tuple(int(value) for value in certificate["pair"])
        _require(len(pair) == 2, "malformed pair")
        source, target = pair
        _require(
            0 <= source <= target < class_count,
            f"invalid unordered pair {pair}",
        )
        _require(pair not in seen_pairs, f"duplicate pair {pair}")
        seen_pairs.add(pair)
        expected_spectral = [
            degree
            for degree in range(1, cutoff)
            if pair_is_spectrally_positive(
                context, source, target, degree
            )
        ]
        _require(
            certificate["spectral_degrees"] == expected_spectral,
            f"spectral degrees mismatch for {pair}",
        )
        non_spectral = [
            degree
            for degree in range(1, cutoff)
            if degree not in expected_spectral
        ]
        run = certificate["run"]
        if not non_spectral:
            _require(run is None, f"unneeded run for {pair}")
            certificate_decisions += sum(
                len(codes_by_degree[degree])
                for degree in expected_spectral
            )
            continue
        _require(isinstance(run, dict), f"missing run for {pair}")
        _require(
            run["target_degrees"] == non_spectral,
            f"target degree mismatch for {pair}",
        )
        max_degree = max(non_spectral)
        expected_dimension = (
            2 * (max_degree + 1) - (degree_p - 1)
        )
        _require(
            run["max_degree"] == max_degree,
            f"max degree mismatch for {pair}",
        )
        _require(
            run["dimension"] == expected_dimension,
            f"dimension mismatch for {pair}",
        )
        expected_iterations = prime**expected_dimension - 1
        _require(
            run["iterations"] == expected_iterations,
            f"iteration count mismatch for {pair}",
        )
        _require(run["exhaustive"] is True, f"nonexhaustive run for {pair}")
        _require(run["invalid_norms"] == 0, f"invalid norms for {pair}")
        _require(
            run["direct_sage_norm_cross_checks"] >= 1,
            f"missing direct checks for {pair}",
        )
        target_codes = [
            code
            for degree in non_spectral
            for code in codes_by_degree[degree]
        ]
        _require(
            run["target_count"] == len(target_codes),
            f"target count mismatch for {pair}",
        )
        counts = {
            int(code): int(count) for code, count in run["counts"]
        }
        _require(
            len(counts) == len(run["counts"]),
            f"duplicate count code for {pair}",
        )
        _require(
            set(counts) == set(target_codes),
            f"count universe mismatch for {pair}",
        )
        zeros = sorted(code for code, count in counts.items() if count == 0)
        _require(
            run["zero_codes"] == zeros,
            f"zero code mismatch for {pair}",
        )
        for code, count in counts.items():
            _require(count >= 0, f"negative count for {pair}, code {code}")
            _require(
                count % expected_aut[target] == 0,
                f"target Aut divisibility failed for {pair}, code {code}",
            )
            _require(
                count % expected_aut[source] == 0,
                f"source Aut divisibility failed for {pair}, code {code}",
            )
        for code in zeros:
            zero_pairs_by_code[code].append(pair)
        vector_total += expected_iterations
        certificate_decisions += sum(
            len(codes_by_degree[degree])
            for degree in expected_spectral
        )
        certificate_decisions += len(target_codes)

    _require(len(seen_pairs) == pair_total, "unordered pair coverage gap")
    _require(
        record["total_vectors_visited"] == vector_total,
        "total vector count mismatch",
    )
    _require(
        certificate_decisions == pair_total * len(expected_levels),
        "pair-level decision coverage mismatch",
    )

    archived_levels = sorted(
        record["levels"], key=lambda item: (item["degree"], item["code"])
    )
    _require(len(archived_levels) == len(expected_levels), "level count mismatch")
    archived_by_code = {
        int(item["code"]): item for item in archived_levels
    }
    _require(
        len(archived_by_code) == len(archived_levels),
        "duplicate archived level code",
    )
    _require(
        set(archived_by_code) == set(expected_by_code),
        "archived level universe mismatch",
    )
    largest_bad = 0
    incomplete_count = 0
    derived_degree_summary = {}
    for degree in range(1, cutoff):
        incomplete = []
        for code in codes_by_degree[degree]:
            expected = expected_by_code[code]
            item = archived_by_code[code]
            _require(item["ell"] == expected["ell"], f"ell mismatch at {code}")
            _require(item["degree"] == degree, f"degree mismatch at {code}")
            ordered = []
            for source, target in zero_pairs_by_code[code]:
                ordered.append([source, target])
                if source != target:
                    ordered.append([target, source])
            ordered.sort()
            _require(
                item["zero_entries"] == ordered,
                f"ordered zero mismatch at {code}",
            )
            _require(
                item["complete"] is (not ordered),
                f"complete flag mismatch at {code}",
            )
            if ordered:
                incomplete.append(expected["ell"])
                incomplete_count += 1
                largest_bad = max(largest_bad, degree)
        derived_degree_summary[str(degree)] = {
            "level_count": len(codes_by_degree[degree]),
            "complete_count": len(codes_by_degree[degree]) - len(incomplete),
            "incomplete_count": len(incomplete),
            "incomplete_levels": incomplete,
        }
    _require(
        record["degree_summary"] == derived_degree_summary,
        "degree summary mismatch",
    )
    expected_E = largest_bad + 1
    _require(record["E"] == expected_E, "E mismatch")
    _require(
        record["engine"]["counter_source_sha256"] == _sha256(COUNTER_SOURCE),
        "counter source hash mismatch",
    )
    _require(
        record["engine"]["counter_binary_sha256"] == _sha256(COUNTER_BINARY),
        "counter binary hash mismatch",
    )
    return {
        "q": prime,
        "p": str(characteristic),
        "degree_p": degree_p,
        "class_count": class_count,
        "theorem_bound_U": cutoff,
        "E": expected_E,
        "level_count": len(expected_levels),
        "incomplete_level_count": incomplete_count,
        "largest_incomplete_degree": largest_bad or None,
        "unordered_pair_count": pair_total,
        "pair_level_decisions": pair_total * len(expected_levels),
        "total_vectors_visited": vector_total,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("archive", type=Path)
    args = parser.parse_args()
    print(json.dumps(validate_archive(args.archive), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
