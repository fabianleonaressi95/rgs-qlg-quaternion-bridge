"""Common-L2 P1 reconstruction and exact finite-element inner product.

This implements G4 for the flat periodic Freudenthal benchmark only.
The map J_h sends nodal values to their continuous piecewise-linear
interpolant on the periodic simplicial mesh. Its L2 inner product is
assembled exactly, not replaced by the lumped G3 mass.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from src.flat_torus_refinement import FlatTorusMesh


def assemble_consistent_mass0(mesh: FlatTorusMesh) -> np.ndarray:
    """Assemble the exact P1 L2 mass matrix on the periodic mesh.

    On a d-simplex of volume V, the local matrix is
        V / ((d+1)(d+2)) * (I + 11^T).
    """
    n_vertices = int(mesh.vertices.shape[0])
    d = int(mesh.dimension)
    mass = np.zeros((n_vertices, n_vertices), dtype=float)
    factor_denominator = (d + 1) * (d + 2)
    for simplex in mesh.simplices:
        ids = np.asarray(simplex, dtype=int)
        origin = mesh.vertices[ids[0]]
        offsets = np.mod(mesh.vertices[ids] - origin, mesh.side_length)
        coords = origin + offsets
        A = (coords[1:] - coords[0]).T
        volume = abs(float(np.linalg.det(A))) / float(np.math.factorial(d))
        local = (volume / factor_denominator) * (
            np.eye(d + 1) + np.ones((d + 1, d + 1))
        )
        mass[np.ix_(ids, ids)] += local
    return mass


@dataclass(frozen=True)
class P1L2Embedding:
    """The explicit linear map J_h: R^V -> L2(T^d_L), represented by P1 data."""

    mesh: FlatTorusMesh
    nodal_values: np.ndarray

    def __post_init__(self) -> None:
        values = np.asarray(self.nodal_values)
        if values.ndim != 1 or values.shape[0] != self.mesh.vertices.shape[0]:
            raise ValueError("nodal_values must be a vector with one value per vertex")
        if not np.all(np.isfinite(values)):
            raise ValueError("nodal_values must be finite")
        object.__setattr__(self, "nodal_values", values.copy())

    def evaluate(self, points: np.ndarray, tolerance: float = 1e-10) -> np.ndarray:
        """Evaluate the periodic P1 interpolant at points in R^d.

        Coordinates are reduced modulo L. Element coordinates are locally
        unwrapped, so simplices crossing a periodic seam are evaluated using
        the corresponding translated point.
        """
        pts = np.asarray(points, dtype=float)
        scalar_input = pts.ndim == 1
        if scalar_input:
            pts = pts[None, :]
        if pts.ndim != 2 or pts.shape[1] != self.mesh.dimension:
            raise ValueError("points must have shape (n_points, dimension)")
        if not np.all(np.isfinite(pts)):
            raise ValueError("points must be finite")

        L = float(self.mesh.side_length)
        pts = np.mod(pts, L)
        result = np.empty(pts.shape[0], dtype=np.result_type(self.nodal_values, float))
        assigned = np.zeros(pts.shape[0], dtype=bool)
        d = self.mesh.dimension

        for simplex in self.mesh.simplices:
            ids = np.asarray(simplex, dtype=int)
            origin = self.mesh.vertices[ids[0]]
            offsets = np.mod(self.mesh.vertices[ids] - origin, L)
            coords = origin + offsets
            A = (coords[1:] - coords[0]).T
            try:
                inverse = np.linalg.inv(A)
            except np.linalg.LinAlgError as exc:
                raise RuntimeError("degenerate simplex in P1 reconstruction") from exc

            remaining = np.flatnonzero(~assigned)
            if remaining.size == 0:
                break
            candidates = pts[remaining]
            shifts = np.ceil((origin[None, :] - candidates) / L - 1e-13)
            unwrapped = candidates + L * shifts
            xi = (unwrapped - coords[0]) @ inverse.T
            bary = np.column_stack((1.0 - np.sum(xi, axis=1), xi))
            inside = np.min(bary, axis=1) >= -tolerance
            inside &= np.max(bary, axis=1) <= 1.0 + tolerance
            if np.any(inside):
                local_values = self.nodal_values[ids]
                result[remaining[inside]] = bary[inside] @ local_values
                assigned[remaining[inside]] = True

        if not np.all(assigned):
            missing = int(np.count_nonzero(~assigned))
            raise RuntimeError(f"could not locate {missing} point(s) in periodic mesh")
        return result[0] if scalar_input else result


def J_h(mesh: FlatTorusMesh, nodal_values: np.ndarray) -> P1L2Embedding:
    """Apply the common-Hilbert-space embedding by continuous P1 interpolation."""
    return P1L2Embedding(mesh=mesh, nodal_values=nodal_values)


def discrete_inner_product(
    mesh: FlatTorusMesh, left: np.ndarray, right: np.ndarray
) -> complex | float:
    """Exact L2 inner product of the two P1 reconstructions."""
    u = np.asarray(left)
    v = np.asarray(right)
    n = mesh.vertices.shape[0]
    if u.shape != (n,) or v.shape != (n,):
        raise ValueError("left and right must have one value per mesh vertex")
    mass = assemble_consistent_mass0(mesh)
    value = np.vdot(u, mass @ v)
    return float(np.real(value)) if np.isrealobj(u) and np.isrealobj(v) else complex(value)


def audit_g4_embedding(mesh: FlatTorusMesh) -> dict[str, float | int | bool | list[int]]:
    """Audit G4 map dimensions, injectivity proxy and exact L2 metric."""
    mass = assemble_consistent_mass0(mesh)
    eigenvalues = np.linalg.eigvalsh((mass + mass.T) / 2)
    constant = np.ones(mesh.vertices.shape[0], dtype=float)
    volume_from_mass = float(constant @ mass @ constant)
    expected_volume = float(mesh.side_length**mesh.dimension)
    symmetry = float(np.linalg.norm(mass - mass.T))
    return {
        "embedding_domain_dimension": int(mesh.vertices.shape[0]),
        "embedding_codomain": "L2(periodic flat torus)",
        "mass_matrix_shape": list(mass.shape),
        "mass_symmetry_residual": symmetry,
        "mass_min_eigenvalue": float(eigenvalues[0]),
        "mass_max_eigenvalue": float(eigenvalues[-1]),
        "mass_condition_number": float(eigenvalues[-1] / eigenvalues[0]),
        "mass_positive_definite": bool(eigenvalues[0] > 0.0),
        "constant_function_l2_norm_squared": volume_from_mass,
        "expected_torus_volume": expected_volume,
        "volume_relative_error": abs(volume_from_mass - expected_volume) / expected_volume,
    }
