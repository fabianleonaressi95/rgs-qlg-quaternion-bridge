import numpy as np

from src.q6_complex import build_q6_complex
from src.q6_whitney_hodge import (
    audit_q6_whitney_hodge,
    build_q6_whitney_hodge,
    build_q6_whitney_metrics,
)


def test_whitney_gram_matrices_are_spd_and_have_expected_blocks():
    metrics = build_q6_whitney_metrics()
    assert metrics.cell_counts == (64, 192, 240)
    assert metrics.metric_0.shape == (64, 64)
    assert metrics.metric_1.shape == (192, 192)
    assert metrics.metric_2.shape == (240, 240)
    for h in (metrics.metric_0, metrics.metric_1, metrics.metric_2):
        assert np.allclose(h, h.T)
        assert np.linalg.eigvalsh(h)[0] > 0
    assert len(metrics.orientation_blocks[0]) == 1
    assert len(metrics.orientation_blocks[1]) == 6
    assert len(metrics.orientation_blocks[2]) == 15


def test_whitney_gram_entries_match_tensor_product_hat_integrals():
    metrics = build_q6_whitney_metrics()
    assert np.isclose(metrics.metric_0[0, 63], (1.0 / 6.0) ** 6)
    q6 = build_q6_complex()
    # The six edge orientations form six mutually orthogonal 32-edge blocks.
    assert len(q6.edges) == 192
    assert np.count_nonzero(metrics.metric_1) == 6 * 32 * 32
    # Fifteen face orientations form mutually orthogonal 16-face blocks.
    assert np.count_nonzero(metrics.metric_2) == 15 * 16 * 16


def test_whitney_hodge_is_symmetric_positive_and_matches_betti_numbers():
    report = audit_q6_whitney_hodge()
    assert report["cell_counts"] == {"C0": 64, "C1": 192, "C2": 240}
    assert report["cellular_betti_numbers"] == (1, 0, 111)
    assert report["harmonic_dimensions"] == (1, 0, 111)
    assert max(report["laplacian_symmetry_residuals"]) < 1e-9
    assert min(report["laplacian_min_eigenvalues"]) > -1e-9
    assert report["boundary_squared_exact_zero"]
    assert report["pass"]


def test_whitney_hodge_builder_accepts_explicit_metric_object():
    metrics = build_q6_whitney_metrics()
    hodge = build_q6_whitney_hodge(metrics)
    assert hodge.laplacian_0.shape == (64, 64)
    assert hodge.laplacian_1.shape == (192, 192)
    assert hodge.laplacian_2.shape == (240, 240)
