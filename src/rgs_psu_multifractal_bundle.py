"""Finite RGS vector-bundle connection and multifractal scaling diagnostics.

The base is the eight Z2^3 orbit labels.  The rank-three fibre is the
coordinate tangent/sign space with local frame D_b=diag(p_b), where p_b
is the positive RGS orbit representative.  This choice induces a
metric-compatible GL(3,R) connection U_{b<-a}=D_b^{-1}D_a.

Important: this canonical connection is pure gauge and therefore flat.
It is a baseline forced by the archived seed, not a discovered Yang-Mills
field.  Nonzero curvature requires independent dynamical link data.

The multifractal routine accepts a *nested sequence* of positive normalized
cell measures.  It estimates tau(q), D_q and D_1 from scale regressions.
It does not manufacture a physical measure from the finite 64-point seed.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import numpy as np


@dataclass(frozen=True)
class DiscreteBundleConnection:
    representatives: np.ndarray
    local_frames: np.ndarray
    fibre_metrics: np.ndarray
    base_edges: np.ndarray
    base_faces: np.ndarray
    transports: dict[tuple[int, int], np.ndarray]


def build_base_cube_2complex() -> tuple[np.ndarray, np.ndarray]:
    """Return the 12 edges and six square faces of the base Q3."""
    edges: list[tuple[int, int]] = []
    for v in range(8):
        for axis in range(3):
            w = v ^ (1 << (2 - axis))
            if v < w:
                edges.append((v, w))
    faces: list[tuple[int, int, int, int]] = []
    for a in range(3):
        ma = 1 << (2 - a)
        for b in range(a + 1, 3):
            mb = 1 << (2 - b)
            for v in range(8):
                if (v & ma) or (v & mb):
                    continue
                faces.append((v, v ^ ma, v ^ ma ^ mb, v ^ mb))
    return np.asarray(sorted(edges), dtype=np.int64), np.asarray(faces, dtype=np.int64)


def build_rgs_vector_bundle_connection(
    representatives: np.ndarray | None = None,
) -> DiscreteBundleConnection:
    """Build the seed-induced metric-compatible rank-three GL(3) connection."""
    if representatives is None:
        from .rgs_ellis_refinement import load_rgs_seed_data

        seed = load_rgs_seed_data()
        representatives = np.abs(seed.coordinates.reshape(8, 8, 3)[:, 0, :])
    p = np.asarray(representatives, dtype=np.float64)
    if p.shape != (8, 3) or not np.all(np.isfinite(p)) or np.any(p <= 0.0):
        raise ValueError("representatives must be a finite positive (8,3) array")

    frames = np.asarray([np.diag(row) for row in p])
    metrics = np.einsum("bij,bik->bjk", frames, frames)
    edges, faces = build_base_cube_2complex()
    transports: dict[tuple[int, int], np.ndarray] = {}
    for a, b in edges:
        # U_{b<-a} maps fibre coordinates at a into fibre coordinates at b.
        u_ba = np.linalg.solve(frames[b], frames[a])
        transports[(int(a), int(b))] = u_ba
        transports[(int(b), int(a))] = np.linalg.inv(u_ba)
    return DiscreteBundleConnection(p, frames, metrics, edges, faces, transports)


def face_holonomies(connection: DiscreteBundleConnection) -> np.ndarray:
    """Compute oriented Q3 face holonomies, based at each first vertex."""
    result = []
    for a, b, c, d in connection.base_faces:
        # Rightmost factor acts first: a -> b -> c -> d -> a.
        hol = (
            connection.transports[(int(a), int(d))]
            @ connection.transports[(int(d), int(c))]
            @ connection.transports[(int(c), int(b))]
            @ connection.transports[(int(b), int(a))]
        )
        result.append(hol)
    return np.asarray(result)


def gauge_transform_connection(
    connection: DiscreteBundleConnection, gauges: np.ndarray
) -> DiscreteBundleConnection:
    """Apply independent GL(3) changes of fibre frame at each base vertex."""
    g = np.asarray(gauges, dtype=np.float64)
    if g.shape != (8, 3, 3) or not np.all(np.isfinite(g)):
        raise ValueError("gauges must be finite with shape (8,3,3)")
    if any(np.linalg.matrix_rank(gb) != 3 for gb in g):
        raise ValueError("every vertex gauge matrix must be invertible")
    inv = np.linalg.inv(g)
    transformed = {
        (a, b): g[b] @ u @ inv[a]
        for (a, b), u in connection.transports.items()
    }
    metrics = np.asarray([
        inv[b].T @ connection.fibre_metrics[b] @ inv[b]
        for b in range(8)
    ])
    frames = np.asarray([connection.local_frames[b] @ inv[b] for b in range(8)])
    return DiscreteBundleConnection(
        connection.representatives.copy(), frames, metrics,
        connection.base_edges.copy(), connection.base_faces.copy(), transformed,
    )


def audit_rgs_vector_bundle_connection(
    connection: DiscreteBundleConnection | None = None,
    atol: float = 1e-11,
) -> dict:
    """Audit metric compatibility, reverse links, and canonical flatness."""
    c = build_rgs_vector_bundle_connection() if connection is None else connection
    metric_residuals = []
    inverse_residuals = []
    for a, b in c.base_edges:
        u = c.transports[(int(a), int(b))]
        metric_residuals.append(
            np.linalg.norm(u.T @ c.fibre_metrics[b] @ u - c.fibre_metrics[a])
        )
        inverse_residuals.append(
            np.linalg.norm(c.transports[(int(b), int(a))] @ u - np.eye(3))
        )
    hol = face_holonomies(c)
    curvature = hol - np.eye(3)[None, :, :]
    return {
        "schema": "RGS_PSU_RANK3_DISCRETE_BUNDLE_V1",
        "base_vertices": 8,
        "base_edges": int(len(c.base_edges)),
        "base_square_faces": int(len(c.base_faces)),
        "fibre_rank": 3,
        "structure_group": "GL(3,R), with metric-compatible transports",
        "max_metric_compatibility_residual": float(max(metric_residuals, default=0.0)),
        "max_reverse_link_residual": float(max(inverse_residuals, default=0.0)),
        "max_face_holonomy_residual": float(np.max(np.linalg.norm(curvature, axis=(1, 2)))),
        "face_holonomies": hol.tolist(),
        "metric_compatible": bool(max(metric_residuals, default=0.0) <= atol),
        "links_are_invertible": bool(max(inverse_residuals, default=0.0) <= atol),
        "canonical_connection_is_flat": bool(
            np.max(np.linalg.norm(curvature, axis=(1, 2))) <= atol
        ),
        "interpretation": (
            "The seed-induced connection is a pure-gauge GL(3) connection. "
            "It is a rigorously defined flat baseline; nonzero gauge curvature "
            "requires additional independent link/dynamical data."
        ),
    }


def multifractal_spectrum(
    scales: Iterable[float],
    level_measures: Iterable[np.ndarray],
    q_values: Iterable[float] = (-2.0, -1.0, 0.0, 0.5, 2.0, 3.0, 4.0),
) -> dict:
    """Estimate tau(q), D_q and D_1 from normalized measures on nested partitions.

    Every level measure must contain strictly positive masses summing to one.
    At least three decreasing scales are required.  The caller is responsible
    for ensuring that levels describe compatible nested partitions.
    """
    ell = np.asarray(tuple(scales), dtype=np.float64)
    measures = tuple(np.asarray(m, dtype=np.float64).reshape(-1) for m in level_measures)
    qs = np.asarray(tuple(q_values), dtype=np.float64)
    if ell.ndim != 1 or len(ell) < 3 or len(measures) != len(ell):
        raise ValueError("provide at least three scales and one measure per scale")
    if np.any(~np.isfinite(ell)) or np.any(ell <= 0.0) or np.any(np.diff(ell) >= 0.0):
        raise ValueError("scales must be finite, positive, and strictly decreasing")
    for mu in measures:
        if mu.size < 2 or np.any(~np.isfinite(mu)) or np.any(mu <= 0.0):
            raise ValueError("each level measure must contain at least two positive finite masses")
        if not np.isclose(mu.sum(), 1.0, rtol=1e-10, atol=1e-12):
            raise ValueError("each level measure must be normalized to one")
    if np.any(~np.isfinite(qs)) or len(np.unique(qs)) != len(qs):
        raise ValueError("q values must be finite and unique")

    log_ell = np.log(ell)
    log_inv_ell = -log_ell
    tau: dict[str, float] = {}
    dimensions: dict[str, float] = {}
    fits: dict[str, float] = {}
    entropy = np.asarray([-np.sum(mu * np.log(mu)) for mu in measures])
    d1, intercept1 = np.polyfit(log_inv_ell, entropy, 1)
    for q in qs:
        key = f"{q:g}"
        if np.isclose(q, 1.0):
            tau[key] = 0.0
            dimensions[key] = float(d1)
            fits[key] = float(intercept1)
            continue
        log_z = np.asarray([np.log(np.sum(mu ** q)) for mu in measures])
        slope, intercept = np.polyfit(log_ell, log_z, 1)
        tau[key] = float(slope)
        dimensions[key] = float(slope / (q - 1.0))
        fits[key] = float(intercept)
    return {
        "scales": ell.tolist(),
        "q_values": qs.tolist(),
        "tau_q": tau,
        "D_q": dimensions,
        "D_1_entropy_slope": float(d1),
        "fit_intercepts": fits,
        "level_cell_counts": [int(mu.size) for mu in measures],
        "nested_partition_compatibility": "CALLER_MUST_SUPPLY_AND_VERIFY_PARENT_CHILD_MAPS",
        "multifractality_diagnostic": (
            "NONLINEAR_DQ_OVER_Q_SUPPORTS_MULTIFRACTAL_SCALING"
            if np.ptp(list(dimensions.values())) > 1e-2
            else "DQ_APPROXIMATELY_CONSTANT_ON_TESTED_Q_RANGE"
        ),
    }
