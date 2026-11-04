import numpy as np
import pytest

from src.q6_complex import build_q6_complex
from src.q6_hodge_laplacian import (
    audit_q6_hodge_laplacians,
    build_q6_hodge_laplacians,
)


def test_unit_weights_match_incidence_hodge_formulas():
    q6 = build_q6_complex()
    hodge = build_q6_hodge_laplacians(complex_=q6)
    b1 = q6.boundary_1.astype(float)
    b2 = q6.boundary_2.astype(float)
    assert np.allclose(hodge.laplacian_0, b1 @ b1.T)
    assert np.allclose(hodge.laplacian_1, b1.T @ b1 + b2 @ b2.T)
    assert np.allclose(hodge.laplacian_2, b2.T @ b2)


def test_harmonic_dimensions_match_q6_cellular_betti_numbers():
    report = audit_q6_hodge_laplacians()
    assert report["cell_counts"] == {"C0": 64, "C1": 192, "C2": 240}
    assert report["boundary_ranks"] == {"rank_B1": 63, "rank_B2": 129}
    assert report["expected_betti_numbers"] == (1, 0, 111)
    assert report["observed_harmonic_dimensions"] == (1, 0, 111)
    assert report["boundary_squared_exact_zero"]
    assert report["pass"]


def test_positive_nonuniform_weights_preserve_harmonic_dimensions():
    rng = np.random.default_rng(20261104)
    hodge = build_q6_hodge_laplacians(
        weights_0=np.exp(rng.normal(scale=0.2, size=64)),
        weights_1=np.exp(rng.normal(scale=0.2, size=192)),
        weights_2=np.exp(rng.normal(scale=0.2, size=240)),
    )
    report = audit_q6_hodge_laplacians(hodge)
    assert report["observed_harmonic_dimensions"] == (1, 0, 111)
    assert min(report["minimum_eigenvalues"]) >= -1e-9
    assert report["pass"]


@pytest.mark.parametrize("kwargs", [
    {"weights_0": np.zeros(64)},
    {"weights_1": -np.ones(192)},
    {"weights_2": np.full(240, np.nan)},
    {"weights_0": np.ones(63)},
])
def test_invalid_weights_are_rejected(kwargs):
    with pytest.raises(ValueError):
        build_q6_hodge_laplacians(**kwargs)
