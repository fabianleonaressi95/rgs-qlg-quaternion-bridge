"""Randomized gauge-covariance tests on the Möbius subcomplex."""
import numpy as np

from src.mobius_complex import build_mobius_complex
from src.mobius_gauge_audit import (
    audit_mobius_gauge_covariance,
    build_random_mobius_lattice,
)
from src.holonomies import dagger


def test_randomized_mobius_gauge_audit_passes_with_fixed_seed():
    result = audit_mobius_gauge_covariance(seed=20261103, atol=1e-10)
    assert result["passed"] is True
    assert result["vertices"] == 12
    assert result["edges"] == 18
    assert result["faces_tested"] == 6
    for key in (
        "open_path_covariance_max",
        "closed_loop_conjugation_max",
        "wilson_trace_invariance_max",
        "inverse_link_consistency_max",
    ):
        assert result[key] <= result["tolerance"]


def test_randomized_audit_is_reproducible():
    first = audit_mobius_gauge_covariance(seed=17)
    second = audit_mobius_gauge_covariance(seed=17)
    assert first == second


def test_gauge_transformation_preserves_su2_links_and_reverse_edges():
    complex_ = build_mobius_complex()
    lattice, gauge = build_random_mobius_lattice(seed=31)
    transformed = lattice.gauge_transform(gauge)
    identity = np.eye(2, dtype=complex)
    for u, v in complex_.edges:
        for link in (transformed.link(u, v), transformed.link(v, u)):
            assert np.allclose(dagger(link) @ link, identity, atol=1e-12)
            assert np.isclose(np.linalg.det(link), 1.0, atol=1e-12)
        assert np.allclose(
            transformed.link(v, u), dagger(transformed.link(u, v)), atol=1e-12
        )
