"""Exact completeness computations for supersingular rank-two Drinfeld modules."""

from .compact_prime import compute_compact_odd_prime
from .compact_q2 import compute_compact_q2_degree7
from .compact_q2_degree8 import compute_compact_q2_degree8
from .core import (
    build_supersingular_context,
    compute_completeness,
    monic_irreducibles,
    theorem_bound,
)
from .e_only_q2 import compute_e_only_q2

__version__ = "0.1.0"

__all__ = [
    "__version__",
    "build_supersingular_context",
    "compute_compact_odd_prime",
    "compute_compact_q2_degree7",
    "compute_compact_q2_degree8",
    "compute_completeness",
    "compute_e_only_q2",
    "monic_irreducibles",
    "theorem_bound",
]
