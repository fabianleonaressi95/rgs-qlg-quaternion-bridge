"""P1 barycentric finite elements on the geometric 1-skeleton of Q6.

The canonical embedding maps each six-bit vertex to its bit vector in R^6.
Each Q6 edge is treated as a one-dimensional affine simplex. The assembled
operators are therefore graph-skeleton FEM operators, not a six-dimensional
volume FEM and not a discretisation of curved RGS spacetime.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from src.barycentric_geometry import weighted_barycentre
from src.q6_complex import Q6Complex, build_q6_complex


@dataclass(frozen=True)
class Q6BarycentricFEM:
    """Consistent P1 mass and Dirichlet stiffness on all Q6 edges."""

    coordinates: np.ndarray
    edges: tuple[tuple[int, int], ...]
    edge_lengths: np.ndarray
    boundary_1: np.ndarray
    mass: np.ndarray
    stiffness: np.ndarray
    barycentre: np.ndarray
    covariance: np.ndarray


def q6_hypercube_coordinates() -> np.ndarray:
    """Return the canonical 64 x 6 bit-vector embedding (axis 0 is MSB)."""
    return np.asarray(
        [[(vertex >> (5 - axis)) & 1 for axis in range(6)] for vertex in range(64)],
        dtype=float,
    )


def build_q6_barycentric_fem(
    coordinates: np.ndarray | None = None,
    complex_: Q6Complex | None = None,
) -> Q6BarycentricFEM:
    """Assemble 1D P1 element matrices over the canonical Q6 edge set.

    For an edge of Euclidean length L, the local matrices are
        M_e = L/6 [[2,1],[1,2]],
        K_e = 1/L [[1,-1],[-1,1]].
    """
    q6 = build_q6_complex() if complex_ is None else complex_
    x = q6_hypercube_coordinates() if coordinates is None else np.asarray(coordinates, dtype=float)
    if x.ndim != 2 or x.shape[0] != len(q6.vertices) or x.shape[1] < 1:
        raise ValueError("coordinates must have shape (64,d) with d >= 1")
    if not np.all(np.isfinite(x)):
        raise ValueError("coordinates must be finite")

    n = len(q6.vertices)
    mass = np.zeros((n, n), dtype=float)
    stiffness = np.zeros((n, n), dtype=float)
    lengths = np.empty(len(q6.edges), dtype=float)
    for edge_id, (u, v) in enumerate(q6.edges):
        delta = x[v] - x[u]
        length = float(np.linalg.norm(delta))
        if not np.isfinite(length) or length <= 1e-12:
            raise ValueError(f"Q6 edge {(u, v)} is geometrically degenerate")
        lengths[edge_id] = length
        local_mass = (length / 6.0) * np.array([[2.0, 1.0], [1.0, 2.0]])
        local_stiffness = (1.0 / length) * np.array([[1.0, -1.0], [-1.0, 1.0]])
        ids = np.ix_([u, v], [u, v])
        mass[ids] += local_mass
        stiffness[ids] += local_stiffness

    weights = np.full(n, 1.0 / n)
    centre, _, covariance = weighted_barycentre(x, weights)
    return Q6BarycentricFEM(
        coordinates=x.copy(),
        edges=q6.edges,
        edge_lengths=lengths,
        boundary_1=np.asarray(q6.boundary_1, dtype=np.int8).copy(),
        mass=mass,
        stiffness=stiffness,
        barycentre=centre,
        covariance=covariance,
    )


def q6_edge_gradients(
    fem: Q6BarycentricFEM,
    nodal_values: np.ndarray,
) -> np.ndarray:
    """Return constant ambient gradients on edges for scalar or vector nodal data.

    Scalar input (64,) returns (192,d); vector-valued input (64,k) returns
    (192,k,d). Gradients are tangent to the corresponding embedded edge.
    """
    values = np.asarray(nodal_values, dtype=float)
    n, d = fem.coordinates.shape
    scalar = values.ndim == 1
    if scalar:
        values = values[:, None]
    if values.ndim != 2 or values.shape[0] != n or not np.all(np.isfinite(values)):
        raise ValueError("nodal_values must be finite with shape (64,) or (64,k)")
    gradients = np.empty((len(fem.edges), values.shape[1], d), dtype=float)
    for i, (u, v) in enumerate(fem.edges):
        tangent = fem.coordinates[v] - fem.coordinates[u]
        gradients[i] = (values[v] - values[u])[:, None] * tangent[None, :] / (
            fem.edge_lengths[i] ** 2
        )
    return gradients[:, 0, :] if scalar else gradients


def q6_dirichlet_energy(fem: Q6BarycentricFEM, nodal_values: np.ndarray) -> float:
    """Return 1/2 f^T K f, equal to half the sum of edgewise Dirichlet energies."""
    values = np.asarray(nodal_values, dtype=float)
    if values.ndim != 1 or values.shape != (len(fem.coordinates),) or not np.all(np.isfinite(values)):
        raise ValueError("nodal_values must be a finite vector of length 64")
    return 0.5 * float(values @ fem.stiffness @ values)


def audit_q6_barycentric_fem(
    fem: Q6BarycentricFEM | None = None,
    complex_: Q6Complex | None = None,
) -> dict[str, object]:
    """Report geometric, chain-complex and FEM invariants for the Q6 1-skeleton."""
    q6 = build_q6_complex() if complex_ is None else complex_
    model = build_q6_barycentric_fem(complex_=q6) if fem is None else fem
    b1 = model.boundary_1.astype(float)
    graph_stiffness = b1 @ np.diag(1.0 / model.edge_lengths) @ b1.T
    ones = np.ones(len(model.coordinates))
    mass_eigenvalues = np.linalg.eigvalsh(model.mass)
    stiffness_eigenvalues = np.linalg.eigvalsh(model.stiffness)
    chain_residual = (
        q6.boundary_1.astype(np.int64) @ q6.boundary_2.astype(np.int64)
    )
    return {
        "vertices": len(model.coordinates),
        "edges": len(model.edges),
        "embedding_dimension": int(model.coordinates.shape[1]),
        "edge_length_min": float(np.min(model.edge_lengths)),
        "edge_length_max": float(np.max(model.edge_lengths)),
        "barycentre": model.barycentre.tolist(),
        "covariance": model.covariance.tolist(),
        "mass_symmetry_residual": float(np.linalg.norm(model.mass - model.mass.T)),
        "stiffness_symmetry_residual": float(np.linalg.norm(model.stiffness - model.stiffness.T)),
        "mass_min_eigenvalue": float(mass_eigenvalues[0]),
        "mass_positive_definite": bool(mass_eigenvalues[0] > 0.0),
        "stiffness_min_eigenvalue": float(stiffness_eigenvalues[0]),
        "stiffness_semidefinite": bool(stiffness_eigenvalues[0] >= -1e-10),
        "stiffness_constant_mode_residual": float(np.linalg.norm(model.stiffness @ ones)),
        "stiffness_graph_factorization_residual": float(np.linalg.norm(model.stiffness - graph_stiffness)),
        "stiffness_kernel_dimension": int(np.count_nonzero(np.abs(stiffness_eigenvalues) < 1e-9)),
        "boundary_squared_exact_zero": bool(np.count_nonzero(chain_residual) == 0),
        "scope": "P1 FEM on Q6 geometric 1-skeleton only; no 2-cell Hodge star, continuum limit, or physical spacetime claim",
    }
