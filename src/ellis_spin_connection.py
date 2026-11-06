"""Levi-Civita spin transport on the exact Ellis-RGS hypersurface.

The spatial metric is ds^2 = du^2 + (u^2+r0^2) dOmega^2, embedded
isometrically in R^4 by F(u,n)=(R(u)n, r0 asinh(u/r0)).
Parallel transport is integrated intrinsically through the embedded
hypersurface: for a tangent vector V and unit normal N,
V' = -<V,N'> N. The resulting SO(3) frame map is lifted to SU(2).

This is a geometric edge connection, not the polar coframe heuristic.
"""
from __future__ import annotations

from collections.abc import Iterable
import numpy as np

from src.discrete_dirac import _rotation_to_quaternion
from src.representations_2i import quaternion_to_fundamental_su2


def _point_data(point_r4: np.ndarray, r0: float):
    p = np.asarray(point_r4, dtype=float)
    if p.shape != (4,) or not np.all(np.isfinite(p)):
        raise ValueError("each Ellis point must be a finite 4-vector")
    z = float(p[3])
    u = r0 * np.sinh(z / r0)
    R = np.sqrt(u*u + r0*r0)
    n = p[:3] / R
    norm = float(np.linalg.norm(n))
    if abs(norm - 1.0) > 2e-7:
        raise ValueError("point does not lie on the declared Ellis hypersurface")
    n /= norm
    expected = np.r_[R*n, r0*np.arcsinh(u/r0)]
    if np.linalg.norm(expected-p) > 2e-7 * max(1.0, np.linalg.norm(p)):
        raise ValueError("point is inconsistent with the exact Ellis embedding")
    return u, n


def ellis_orthonormal_frame(u: float, n: np.ndarray, r0: float) -> np.ndarray:
    """Return an oriented 4x3 orthonormal tangent frame in the Ellis embedding."""
    u, r0 = float(u), float(r0)
    n = np.asarray(n, dtype=float)
    if n.shape != (3,) or not np.isfinite(u) or not np.isfinite(r0) or r0 <= 0:
        raise ValueError("expected finite u, positive r0 and n in R^3")
    norm = float(np.linalg.norm(n))
    if abs(norm-1.0) > 1e-9:
        raise ValueError("angular direction n must be unit")
    n = n/norm
    axis = np.eye(3)[int(np.argmin(np.abs(n)))]
    v1 = axis - float(np.dot(axis, n))*n
    v1 /= np.linalg.norm(v1)
    v2 = np.cross(n, v1)
    R = np.sqrt(u*u+r0*r0)
    radial = np.r_[(u/R)*n, r0/R]
    angular1 = np.r_[v1, 0.0]
    angular2 = np.r_[v2, 0.0]
    frame = np.column_stack((radial, angular1, angular2))
    if np.linalg.norm(frame.T@frame-np.eye(3)) > 2e-12:
        raise ArithmeticError("Ellis tangent frame lost orthonormality")
    return frame


def _normal(u: float, n: np.ndarray, r0: float) -> np.ndarray:
    R = np.sqrt(u*u+r0*r0)
    return np.r_[(r0/R)*n, -u/R]


def _normal_derivative(u: float, n: np.ndarray, ndot: np.ndarray,
                       udot: float, r0: float) -> np.ndarray:
    R = np.sqrt(u*u+r0*r0)
    return np.r_[
        (-r0*u*udot/R**3)*n + (r0/R)*ndot,
        -(r0*r0*udot/R**3),
    ]


def _great_circle(n0: np.ndarray, n1: np.ndarray):
    dot = float(np.clip(np.dot(n0, n1), -1.0, 1.0))
    theta = float(np.arccos(dot))
    if theta < 1e-13:
        return theta, lambda t: n0.copy(), lambda t: np.zeros(3)
    if np.pi-theta < 1e-8:
        raise ValueError("antipodal angular endpoints have no unique shortest path")
    tangent = (n1-dot*n0)/np.sin(theta)
    return (
        theta,
        lambda t: np.cos(theta*t)*n0 + np.sin(theta*t)*tangent,
        lambda t: theta*(-np.sin(theta*t)*n0 + np.cos(theta*t)*tangent),
    )


def ellis_edge_spin_link(point_i, point_j, r0: float, steps: int = 32):
    """Return U_ji in SU(2), parallel transport from vertex i to vertex j.

    The path is explicitly fixed: u is interpolated affinely and n follows
    the shorter great-circle arc. The ODE is integrated by classical RK4.
    """
    r0 = float(r0)
    steps = int(steps)
    if not np.isfinite(r0) or r0 <= 0 or steps < 1:
        raise ValueError("r0 must be positive and steps must be >= 1")
    ui, ni = _point_data(point_i, r0)
    uj, nj = _point_data(point_j, r0)
    theta, n_path, n_dot_path = _great_circle(ni, nj)
    udot = uj-ui
    Fi = ellis_orthonormal_frame(ui, ni, r0)
    Fj = ellis_orthonormal_frame(uj, nj, r0)
    V0 = Fi.copy()

    def rhs(t, V):
        u = ui+t*udot
        n = n_path(t)
        ndot = n_dot_path(t)
        N = _normal(u, n, r0)
        Ndot = _normal_derivative(u, n, ndot, udot, r0)
        return -np.outer(N, Ndot@V)

    h = 1.0/steps
    V = V0.copy()
    for k in range(steps):
        t = k*h
        k1 = rhs(t, V)
        k2 = rhs(t+h/2, V+h*k1/2)
        k3 = rhs(t+h/2, V+h*k2/2)
        k4 = rhs(t+h, V+h*k3)
        V = V + (h/6)*(k1+2*k2+2*k3+k4)

    N1 = _normal(uj, nj, r0)
    V = (np.eye(4)-np.outer(N1, N1))@V
    left, _, right = np.linalg.svd(V, full_matrices=False)
    V = left@right
    rotation = Fj.T@V
    # Numerical integration should preserve the oriented tangent frame.
    if np.linalg.det(rotation) < 0:
        raise ArithmeticError("parallel transport changed tangent-frame orientation")
    uquat = _rotation_to_quaternion(rotation)
    link = quaternion_to_fundamental_su2(uquat)
    orth = float(np.linalg.norm(rotation.T@rotation-np.eye(3), ord="fro"))
    det_res = float(abs(np.linalg.det(rotation)-1.0))
    return link, {
        "steps": steps,
        "angular_arc_length": theta,
        "so3_orthogonality_residual": orth,
        "so3_determinant_residual": det_res,
        "spin_unitarity_residual": float(np.linalg.norm(
            link.conj().T@link-np.eye(2), ord="fro")),
        "transport_equation": "V'=-<V,N'>N",
        "path": "affine u and shortest great-circle n",
    }


def ellis_edge_spin_links(vertices_r4, edges: Iterable[tuple[int, int]],
                          r0: float, steps: int = 32):
    """Build canonical i<j Levi-Civita SU(2) links on an Ellis mesh."""
    vertices = np.asarray(vertices_r4, dtype=float)
    if vertices.ndim != 2 or vertices.shape[1] != 4:
        raise ValueError("vertices_r4 must have shape (N,4)")
    links = {}
    max_residual = 0.0
    for raw_i, raw_j in edges:
        i, j = int(raw_i), int(raw_j)
        if not (0 <= i < j < len(vertices)):
            raise ValueError("edges must be canonical i<j vertex indices")
        link, audit = ellis_edge_spin_link(vertices[i], vertices[j], r0, steps)
        links[(i, j)] = link
        max_residual = max(max_residual, audit["so3_orthogonality_residual"],
                           audit["so3_determinant_residual"],
                           audit["spin_unitarity_residual"])
    return links, {
        "vertices": int(len(vertices)),
        "edges": int(len(links)),
        "r0": float(r0),
        "integration_steps_per_edge": int(steps),
        "max_transport_residual": max_residual,
        "connection": "Levi-Civita of exact Ellis-RGS hypersurface",
        "spin_lift": "SU(2) fundamental representation",
        "path_rule": "affine radial coordinate plus shortest spherical geodesic",
        "status": "GEOMETRIC_FINITE_MESH_SPIN_CONNECTION",
    }
