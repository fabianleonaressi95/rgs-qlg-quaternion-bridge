"""Affine simplex barycentric geometry and P1 finite-element operators.

This module is geometric/numerical infrastructure. It does not identify the
simplex metric with a physical RGS spacetime metric.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
import numpy as np


@dataclass(frozen=True)
class SimplexGeometry:
    vertices: np.ndarray
    barycentric_gradients: np.ndarray
    volume: float
    metric: np.ndarray
    metric_inverse: np.ndarray

    @property
    def dimension(self) -> int:
        return int(self.vertices.shape[1])


def simplex_geometry(vertices: np.ndarray, rtol: float = 1e-12) -> SimplexGeometry:
    """Build exact-affine P1 geometry for d+1 vertices in R^d."""
    x = np.asarray(vertices, dtype=float)
    if x.ndim != 2 or x.shape[0] != x.shape[1] + 1:
        raise ValueError("vertices must have shape (d+1,d)")
    if not np.all(np.isfinite(x)):
        raise ValueError("vertices must be finite")
    d = x.shape[1]
    jac = (x[1:] - x[0]).T
    singular = np.linalg.svd(jac, compute_uv=False)
    if singular[-1] <= rtol * singular[0]:
        raise ValueError("degenerate simplex: affine Jacobian is rank deficient")
    metric = jac.T @ jac
    metric_inverse = np.linalg.inv(metric)
    augmented = np.vstack((np.ones(d + 1), x.T))
    inverse_augmented = np.linalg.inv(augmented)
    gradients = inverse_augmented[:, 1:]
    volume = abs(float(np.linalg.det(jac))) / math.factorial(d)
    return SimplexGeometry(x.copy(), gradients, volume, metric, metric_inverse)


def barycentric_coordinates(vertices: np.ndarray, points: np.ndarray) -> np.ndarray:
    """Return P1 barycentric coordinates for one point or an array of points."""
    geom = simplex_geometry(vertices)
    p = np.asarray(points, dtype=float)
    scalar = p.ndim == 1
    if scalar:
        p = p[None, :]
    if p.ndim != 2 or p.shape[1] != geom.dimension or not np.all(np.isfinite(p)):
        raise ValueError("points must have shape (n,d) and be finite")
    augmented = np.column_stack((np.ones(p.shape[0]), p))
    # lambda_i(x) = [1,x] dot row_i(inv([1;X^T]))
    inv_aug = np.linalg.inv(np.vstack((np.ones(len(geom.vertices)), geom.vertices.T)))
    result = augmented @ inv_aug.T
    return result[0] if scalar else result


def weighted_barycentre(vertices: np.ndarray, weights: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return centre, centred vertices, and weighted second-moment tensor."""
    x = np.asarray(vertices, dtype=float)
    w = np.asarray(weights, dtype=float)
    if x.ndim != 2 or w.shape != (x.shape[0],):
        raise ValueError("weights must have one entry per vertex")
    if not np.all(np.isfinite(x)) or not np.all(np.isfinite(w)) or np.any(w < 0):
        raise ValueError("vertices and nonnegative weights must be finite")
    total = float(np.sum(w))
    if total <= 0:
        raise ValueError("weights must have positive sum")
    w = w / total
    centre = w @ x
    centred = x - centre
    covariance = np.einsum("i,ij,ik->jk", w, centred, centred)
    return centre, centred, covariance


def p1_element_matrices(vertices: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return consistent P1 mass and Euclidean Dirichlet stiffness matrices."""
    geom = simplex_geometry(vertices)
    d = geom.dimension
    n = d + 1
    mass = geom.volume / ((d + 1) * (d + 2)) * (np.ones((n, n)) + np.eye(n))
    stiffness = geom.volume * (geom.barycentric_gradients @ geom.barycentric_gradients.T)
    return mass, stiffness


def audit_simplex(vertices: np.ndarray) -> dict[str, float | int | bool]:
    """Return reproducible affine and matrix-structure residuals."""
    geom = simplex_geometry(vertices)
    mass, stiffness = p1_element_matrices(vertices)
    grad_sum = np.sum(geom.barycentric_gradients, axis=0)
    ones = np.ones(len(vertices))
    return {
        "dimension": geom.dimension,
        "volume": geom.volume,
        "partition_of_unity_residual": float(np.max(np.abs(
            barycentric_coordinates(vertices, vertices) - np.eye(len(vertices))
        ))),
        "gradient_sum_residual": float(np.linalg.norm(grad_sum)),
        "mass_symmetry_residual": float(np.linalg.norm(mass - mass.T)),
        "stiffness_symmetry_residual": float(np.linalg.norm(stiffness - stiffness.T)),
        "stiffness_constant_mode_residual": float(np.linalg.norm(stiffness @ ones)),
        "mass_positive_definite": bool(np.min(np.linalg.eigvalsh(mass)) > 0.0),
        "stiffness_positive_semidefinite": bool(np.min(np.linalg.eigvalsh(stiffness)) > -1e-12),
    }
