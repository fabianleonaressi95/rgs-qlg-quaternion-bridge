"""Exact barycentric moments and quadrature bounds for planar triangles.

This module is deliberately local: it does not construct a Q6 surface or a
refinement family. It provides the cellwise algebra needed by Gate E.1.
"""
from __future__ import annotations

import numpy as np


def triangle_barycentric_data(vertices: np.ndarray) -> tuple[float, np.ndarray, np.ndarray]:
    """Return (area, areal centroid, ambient 3x3 covariance) of a triangle.

    Sigma = area**-1 * integral_F (x-c)(x-c)^T dA
          = (1/12) sum_i (v_i-c)(v_i-c)^T.
    """
    v = np.asarray(vertices, dtype=float)
    if v.shape != (3, 3):
        raise ValueError("vertices must have shape (3, 3)")
    if not np.isfinite(v).all():
        raise ValueError("vertices must be finite")
    e1, e2 = v[1] - v[0], v[2] - v[0]
    area = 0.5 * float(np.linalg.norm(np.cross(e1, e2)))
    if not np.isfinite(area) or area <= np.finfo(float).eps:
        raise ValueError("degenerate triangle")
    centroid = np.mean(v, axis=0)
    centred = v - centroid
    covariance = (centred.T @ centred) / 12.0
    return area, centroid, covariance


def quadratic_triangle_integral(vertices, constant, linear, hessian) -> float:
    """Integrate g(x)=constant+linear.x+1/2 x^T H x exactly on a triangle."""
    area, centroid, covariance = triangle_barycentric_data(vertices)
    b = np.asarray(linear, dtype=float)
    h = np.asarray(hessian, dtype=float)
    if b.shape != (3,) or h.shape != (3, 3):
        raise ValueError("linear and hessian must have shapes (3,) and (3, 3)")
    if not np.isfinite(b).all() or not np.isfinite(h).all() or not np.isfinite(constant):
        raise ValueError("polynomial coefficients must be finite")
    if not np.allclose(h, h.T, rtol=1e-12, atol=1e-12):
        raise ValueError("hessian must be symmetric")
    value_at_centroid = float(constant + b @ centroid + 0.5 * centroid @ h @ centroid)
    return float(area * (value_at_centroid + 0.5 * np.sum(h * covariance)))


def triangle_cubic_remainder_bound(area, diameter, third_derivative_supremum) -> float:
    """Return A*h^3*sup||D^3 g||/6, the local Taylor remainder bound."""
    values = np.asarray([area, diameter, third_derivative_supremum], dtype=float)
    if not np.isfinite(values).all() or np.any(values < 0.0):
        raise ValueError("area, diameter and derivative bound must be finite and nonnegative")
    return float(area * diameter**3 * third_derivative_supremum / 6.0)