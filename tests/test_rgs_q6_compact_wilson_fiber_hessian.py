import numpy as np

from src.rgs_q6_compact_wilson_fiber_hessian import (
    compact_wilson_fiber_hessian,
    finite_fibre_log_partition,
)


def test_compact_wilson_fibre_hessian_matches_finite_difference():
    rng = np.random.default_rng(20261104)
    m = rng.normal(size=(17, 3))
    eta = rng.uniform(-np.pi, np.pi, size=(31, 17))
    logw = rng.normal(size=31)
    beta = 0.73
    report = compact_wilson_fiber_hessian(m, eta, logw, beta)
    h = np.asarray(report["effective_hessian"])
    eps = 2e-4
    f0 = finite_fibre_log_partition(np.zeros(3), m, eta, logw, beta)
    fd = np.zeros((3, 3))
    for i in range(3):
        ei = np.eye(3)[i] * eps
        for j in range(3):
            ej = np.eye(3)[j] * eps
            fd[i, j] = -(
                finite_fibre_log_partition(ei + ej, m, eta, logw, beta)
                - finite_fibre_log_partition(ei - ej, m, eta, logw, beta)
                - finite_fibre_log_partition(-ei + ej, m, eta, logw, beta)
                + finite_fibre_log_partition(-ei - ej, m, eta, logw, beta)
            ) / (4 * eps * eps)
    assert report["probability_sum_residual"] < 1e-12
    assert np.allclose(h, fd, atol=2e-6, rtol=2e-5)


def test_single_fibre_point_reduces_to_cosine_hessian():
    m = np.array([[1.0, 0.0], [0.0, 2.0], [1.0, 1.0]])
    eta = np.zeros((1, 3))
    report = compact_wilson_fiber_hessian(m, eta, np.zeros(1), 0.4)
    expected = 0.4 * (m.T @ m)
    assert np.allclose(report["effective_hessian"], expected, atol=1e-12)
    assert np.allclose(report["gradient_covariance"], np.zeros((2, 2)), atol=1e-12)
