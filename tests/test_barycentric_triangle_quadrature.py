import numpy as np
import pytest

from src.barycentric_triangle_quadrature import (
    quadratic_triangle_integral,
    triangle_barycentric_data,
    triangle_cubic_remainder_bound,
)

V = np.array([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 1.0, 0.0]])


def _degree_three_exact_integral(vertices, function):
    """Independent degree-three Dunavant rule (weights sum to one)."""
    area, _, _ = triangle_barycentric_data(vertices)
    barycentric = [
        (1 / 3, 1 / 3, 1 / 3),
        (0.6, 0.2, 0.2),
        (0.2, 0.6, 0.2),
        (0.2, 0.2, 0.6),
    ]
    weights = [-27 / 48, 25 / 48, 25 / 48, 25 / 48]
    points = [np.asarray(lam) @ vertices for lam in barycentric]
    return area * sum(w * function(x) for w, x in zip(weights, points))


def test_triangle_barycentre_is_vertex_average():
    area, centroid, _ = triangle_barycentric_data(V)
    assert area == pytest.approx(1.0)
    assert np.allclose(centroid, V.mean(axis=0))


def test_triangle_covariance_matches_barycentric_second_moments():
    _, centroid, covariance = triangle_barycentric_data(V)
    centred = V - centroid
    expected = sum(np.outer(y, y) for y in centred) / 12.0
    assert np.allclose(covariance, expected, atol=1e-15)
    assert np.allclose(covariance, covariance.T)
    assert np.linalg.eigvalsh(covariance).min() >= -1e-14


def test_constant_and_affine_triangle_integrals_are_barycentrically_exact():
    area, centroid, _ = triangle_barycentric_data(V)
    constant = 3.25
    linear = np.array([1.5, -2.0, 4.0])
    integral = quadratic_triangle_integral(V, constant, linear, np.zeros((3, 3)))
    assert integral == pytest.approx(area * (constant + linear @ centroid))


def test_quadratic_triangle_integral_is_exact():
    constant = -0.7
    linear = np.array([0.5, -1.25, 2.0])
    hessian = np.array([[2.0, 0.3, 0.0], [0.3, -1.0, 0.0], [0.0, 0.0, 0.8]])
    exact = quadratic_triangle_integral(V, constant, linear, hessian)
    independent = _degree_three_exact_integral(
        V, lambda x: constant + linear @ x + 0.5 * x @ hessian @ x
    )
    assert exact == pytest.approx(independent, rel=1e-13, abs=1e-13)


def test_cubic_remainder_respects_one_sixth_bound_and_scales_as_h5():
    ratios, remainders = [], []
    for scale in (1.0, 0.5, 0.25):
        vertices = scale * V
        area, centroid, covariance = triangle_barycentric_data(vertices)
        diameter = max(np.linalg.norm(vertices[i] - vertices[j])
                       for i in range(3) for j in range(i))
        g = lambda x: x[0] ** 3
        hessian = np.zeros((3, 3))
        hessian[0, 0] = 6.0 * centroid[0]
        quadrature = area * g(centroid)
        correction = 0.5 * area * np.sum(hessian * covariance)
        exact = _degree_three_exact_integral(vertices, g)
        remainder = exact - quadrature - correction
        bound = triangle_cubic_remainder_bound(area, diameter, 6.0)
        assert abs(remainder) <= bound + 1e-14
        ratios.append(abs(remainder) / bound)
        remainders.append(abs(remainder))
    assert max(ratios) <= 1.0 / 6.0 + 1e-14
    assert remainders[1] == pytest.approx(remainders[0] / 32.0)
    assert remainders[2] == pytest.approx(remainders[1] / 32.0)


def test_degenerate_triangle_and_invalid_bound_are_rejected():
    with pytest.raises(ValueError, match="degenerate"):
        triangle_barycentric_data(np.zeros((3, 3)))
    with pytest.raises(ValueError, match="finite and nonnegative"):
        triangle_cubic_remainder_bound(1.0, -1.0, 2.0)