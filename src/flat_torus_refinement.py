"""Periodic flat-torus simplicial refinement and P1 finite-element audit.

The construction is a controlled geometric benchmark for G0--G3.
It is not an RGS-curved refinement and does not claim a Mosco theorem.
"""
from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations, permutations
import math
import os
from typing import Any

import numpy as np


@dataclass(frozen=True)
class FlatTorusMesh:
    dimension: int
    cells_per_axis: int
    side_length: float
    vertices: np.ndarray
    simplices: tuple[tuple[int, ...], ...]
    edges: tuple[tuple[int, int], ...]
    triangles: tuple[tuple[int, int, int], ...]
    boundary_1: np.ndarray
    boundary_2: np.ndarray
    mass0: np.ndarray
    stiffness0: np.ndarray
    max_edge_length: float
    max_simplex_condition: float


def _vertex_id(index: tuple[int, ...], n: int) -> int:
    return int(np.ravel_multi_index(index, (n,) * len(index)))


def build_flat_torus_mesh(
    cells_per_axis: int,
    dimension: int = 2,
    side_length: float = 1.0,
) -> FlatTorusMesh:
    """Build a periodic Freudenthal triangulation of a flat d-torus.

    Each periodic cube is split into d! simplices, one for each axis
    permutation. Coordinates used for element geometry are unwrapped within
    the cell; vertex identifiers are periodic modulo cells_per_axis.
    """
    n = int(cells_per_axis)
    d = int(dimension)
    length = float(side_length)
    if d not in (2, 3):
        raise ValueError("dimension must be 2 or 3")
    if n < 3:
        raise ValueError("cells_per_axis must be >= 3 for this periodic mesh")
    if not np.isfinite(length) or length <= 0:
        raise ValueError("side_length must be finite and positive")

    h = length / n
    grid = np.indices((n,) * d).reshape(d, -1).T
    vertices = grid.astype(float) * h
    simplex_ids: list[tuple[int, ...]] = []
    stiffness = np.zeros((n**d, n**d), dtype=float)
    mass = np.zeros(n**d, dtype=float)
    edge_set: set[tuple[int, int]] = set()
    triangle_set: set[tuple[int, int, int]] = set()
    max_cond = 0.0
    max_edge = 0.0
    eye_steps = np.eye(d, dtype=int)

    for base in np.ndindex(*(n,) * d):
        base_arr = np.asarray(base, dtype=int)
        for perm in permutations(range(d)):
            offsets = [np.zeros(d, dtype=int)]
            current = np.zeros(d, dtype=int)
            for axis in perm:
                current = current + eye_steps[axis]
                offsets.append(current.copy())

            ids = tuple(
                _vertex_id(tuple((base_arr + off) % n), n)
                for off in offsets
            )
            if len(set(ids)) != d + 1:
                raise RuntimeError("degenerate periodic simplex: repeated vertex")
            simplex_ids.append(ids)

            coords = np.asarray(
                [(base_arr + off) * h for off in offsets], dtype=float
            )
            A = (coords[1:] - coords[0]).T
            det = float(np.linalg.det(A))
            volume = abs(det) / math.factorial(d)
            if volume <= 0:
                raise RuntimeError("degenerate simplex geometry")
            cond = float(np.linalg.cond(A))
            max_cond = max(max_cond, cond)
            invA = np.linalg.inv(A)
            grads = np.empty((d + 1, d), dtype=float)
            grads[1:, :] = invA
            grads[0, :] = -np.sum(grads[1:, :], axis=0)
            local_K = volume * (grads @ grads.T)
            for a, va in enumerate(ids):
                mass[va] += volume / (d + 1)
                for b, vb in enumerate(ids):
                    stiffness[va, vb] += local_K[a, b]

            for va, vb in combinations(ids, 2):
                edge = (min(va, vb), max(va, vb))
                edge_set.add(edge)
            for tri in combinations(ids, 3):
                triangle_set.add(tuple(sorted(tri)))

            # Physical edge lengths use the unwrapped element coordinates.
            for a, b in combinations(range(d + 1), 2):
                max_edge = max(max_edge, float(np.linalg.norm(coords[a] - coords[b])))

    edges = tuple(sorted(edge_set))
    triangles = tuple(sorted(triangle_set))
    edge_index = {edge: k for k, edge in enumerate(edges)}
    b1 = np.zeros((n**d, len(edges)), dtype=np.int8)
    for k, (u, v) in enumerate(edges):
        b1[u, k] = -1
        b1[v, k] = 1

    b2 = np.zeros((len(edges), len(triangles)), dtype=np.int8)
    for f, (a, b, c) in enumerate(triangles):
        # Boundary of [a,b,c] = [b,c] - [a,c] + [a,b].
        for edge, sign in (((b, c), 1), ((a, c), -1), ((a, b), 1)):
            b2[edge_index[edge], f] += sign

    return FlatTorusMesh(
        dimension=d,
        cells_per_axis=n,
        side_length=length,
        vertices=vertices,
        simplices=tuple(simplex_ids),
        edges=edges,
        triangles=triangles,
        boundary_1=b1,
        boundary_2=b2,
        mass0=mass,
        stiffness0=stiffness,
        max_edge_length=max_edge,
        max_simplex_condition=max_cond,
    )


def audit_flat_torus_mesh(mesh: FlatTorusMesh) -> dict[str, Any]:
    """Audit G0--G3 structural, geometric and scalar-Hodge data."""
    chain = (
        mesh.boundary_1.astype(np.int64)
        @ mesh.boundary_2.astype(np.int64)
    )
    K = mesh.stiffness0
    eig_min = float(np.linalg.eigvalsh((K + K.T) / 2)[0])
    return {
        "dimension": mesh.dimension,
        "cells_per_axis": mesh.cells_per_axis,
        "side_length": mesh.side_length,
        "spacing": mesh.side_length / mesh.cells_per_axis,
        "vertex_count": int(mesh.vertices.shape[0]),
        "edge_count": len(mesh.edges),
        "triangle_count": len(mesh.triangles),
        "simplex_count": len(mesh.simplices),
        "boundary_1_shape": list(mesh.boundary_1.shape),
        "boundary_2_shape": list(mesh.boundary_2.shape),
        "boundary_squared_max_abs": int(np.max(np.abs(chain))) if chain.size else 0,
        "boundary_squared_exact_zero": bool(np.count_nonzero(chain) == 0),
        "mass0_min": float(np.min(mesh.mass0)),
        "mass0_max": float(np.max(mesh.mass0)),
        "mass0_total": float(np.sum(mesh.mass0)),
        "stiffness_symmetry_residual": float(np.linalg.norm(K - K.T)),
        "stiffness_row_sum_residual": float(np.linalg.norm(np.sum(K, axis=1))),
        "stiffness_min_eigenvalue": eig_min,
        "max_edge_length": mesh.max_edge_length,
        "max_simplex_condition": mesh.max_simplex_condition,
        "hodge_star_0_positive": bool(np.all(mesh.mass0 > 0)),
        "stiffness_positive_semidefinite_with_tolerance": bool(eig_min > -1e-10),
    }


def _smooth_mode_energy_audit(mesh: FlatTorusMesh, mode: tuple[int, ...]) -> dict[str, float]:
    d = mesh.dimension
    if len(mode) != d:
        raise ValueError("mode must have one integer component per dimension")
    x = mesh.vertices
    k = 2.0 * np.pi * np.asarray(mode, dtype=float) / mesh.side_length
    u = np.sin(x @ k)
    discrete = float(u @ mesh.stiffness0 @ u)
    continuum = (
        0.5 * float(np.dot(k, k)) * mesh.side_length**d
        if any(mode) else 0.0
    )
    return {
        "discrete_energy": discrete,
        "continuum_energy": continuum,
        "relative_energy_error": (
            abs(discrete - continuum) / continuum if continuum > 0 else abs(discrete)
        ),
    }


def flat_torus_refinement_audit(
    levels: tuple[int, ...] = (4, 8, 16),
    dimension: int = 2,
    side_length: float = 1.0,
    mode: tuple[int, ...] | None = None,
) -> dict[str, Any]:
    """Build at least three mesh levels and emit a reproducible G0--G3 report."""
    ns = tuple(int(v) for v in levels)
    if len(ns) < 3 or any(b <= a for a, b in zip(ns, ns[1:])):
        raise ValueError("levels must contain at least three strictly increasing integers")
    if mode is None:
        mode = (1,) + (0,) * (dimension - 1)
    reports = []
    for n in ns:
        mesh = build_flat_torus_mesh(n, dimension, side_length)
        report = audit_flat_torus_mesh(mesh)
        report.update(_smooth_mode_energy_audit(mesh, tuple(mode)))
        reports.append(report)
    return {
        "audit": "FLAT_TORUS_FREUDENTHAL_G0_G3",
        "generator_version": "0.1.0",
        "provenance_commit": os.environ.get("GITHUB_SHA"),
        "geometry": "flat periodic torus (R/LZ)^d",
        "triangulation": "periodic Freudenthal simplices",
        "finite_element_space": "continuous piecewise-linear P1",
        "hodge_star_scope": "positive lumped 0-form Hodge star only",
        "levels": list(ns),
        "dimension": int(dimension),
        "side_length": float(side_length),
        "mode": list(mode),
        "refinements": reports,
        "scope_warning": (
            "G0-G3 flat benchmark only; not a curved RGS refinement and "
            "not a proof of Mosco convergence."
        ),
    }
