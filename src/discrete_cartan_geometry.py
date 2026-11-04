"""Metric-compatible discrete connection and Cartan face holonomy.

Convention:
- G_v is a positive-definite metric matrix in the local frame at vertex v.
- U_(v<-u) maps vector components in frame u to frame v.
- Metric compatibility is U.T @ G_v @ U = G_u.
- Edge coframes theta_(u,v) are stored for canonical u<v in the source frame u.
- Face torsion is the signed sum of edge coframes parallel-transported to
  the first face vertex. Face curvature is the ordered connection holonomy.

This is a finite cochain/holonomy discretisation. It is not automatically
the unique Levi-Civita connection, and it does not imply a Pin lift.
"""
from __future__ import annotations

from typing import Hashable, Mapping
import numpy as np


def _spd_power(metric: np.ndarray, exponent: float) -> np.ndarray:
    g = np.asarray(metric, dtype=float)
    if g.ndim != 2 or g.shape[0] != g.shape[1] or not np.all(np.isfinite(g)):
        raise ValueError("Metric must be a finite square matrix.")
    if not np.allclose(g, g.T, atol=1e-12):
        raise ValueError("Metric must be symmetric.")
    values, vectors = np.linalg.eigh(g)
    if np.min(values) <= 0:
        raise ValueError("Metric must be positive definite.")
    return (vectors * (values ** exponent)[None, :]) @ vectors.T


def metric_compatible_edge_transport(
    metric_source: np.ndarray,
    metric_target: np.ndarray,
    orthogonal_factor: np.ndarray | None = None,
) -> np.ndarray:
    """Build U_(target<-source) with U.T G_target U = G_source.

    The orthogonal factor selects the frame-rotation freedom. Identity is the
    deterministic default; it is a prescription, not a uniqueness theorem.
    """
    gs = np.asarray(metric_source, dtype=float)
    gt = np.asarray(metric_target, dtype=float)
    if gs.shape != gt.shape or gs.ndim != 2 or gs.shape[0] != gs.shape[1]:
        raise ValueError("Source and target metrics must have equal square shapes.")
    dim = gs.shape[0]
    rotation = np.eye(dim) if orthogonal_factor is None else np.asarray(orthogonal_factor, dtype=float)
    if rotation.shape != (dim, dim):
        raise ValueError("Orthogonal factor has the wrong shape.")
    if not np.allclose(rotation.T @ rotation, np.eye(dim), atol=1e-10):
        raise ValueError("Orthogonal factor must be orthogonal.")
    if np.linalg.det(rotation) < 0:
        raise ValueError("Connection orthogonal factor must preserve orientation.")
    return _spd_power(gt, -0.5) @ rotation @ _spd_power(gs, 0.5)


def build_metric_compatible_connection(
    metrics: Mapping[Hashable, np.ndarray],
    edges: tuple[tuple[Hashable, Hashable], ...],
    orthogonal_factors: Mapping[tuple[Hashable, Hashable], np.ndarray] | None = None,
) -> dict[tuple[Hashable, Hashable], np.ndarray]:
    """Build both orientations of every metric-compatible edge transport."""
    factors = {} if orthogonal_factors is None else orthogonal_factors
    connection: dict[tuple[Hashable, Hashable], np.ndarray] = {}
    for source, target in edges:
        if source not in metrics or target not in metrics:
            raise KeyError("Every edge endpoint needs a local metric.")
        factor = factors.get((source, target))
        forward = metric_compatible_edge_transport(metrics[source], metrics[target], factor)
        reverse = np.linalg.inv(forward)
        connection[(source, target)] = forward
        connection[(target, source)] = reverse
    return connection


def audit_metric_compatibility(
    metrics: Mapping[Hashable, np.ndarray],
    connection: Mapping[tuple[Hashable, Hashable], np.ndarray],
    edges: tuple[tuple[Hashable, Hashable], ...],
) -> dict[str, object]:
    residuals = []
    inverse_residuals = []
    for source, target in edges:
        u = np.asarray(connection[(source, target)], dtype=float)
        v = np.asarray(connection[(target, source)], dtype=float)
        residuals.append(float(np.linalg.norm(
            u.T @ metrics[target] @ u - metrics[source], ord="fro"
        )))
        inverse_residuals.append(float(np.linalg.norm(v @ u - np.eye(u.shape[0]), ord="fro")))
    return {
        "edge_count": len(edges),
        "metric_compatibility_max_fro": max(residuals, default=0.0),
        "inverse_transport_max_fro": max(inverse_residuals, default=0.0),
    }


def transform_local_frame_data(
    metrics: Mapping[Hashable, np.ndarray],
    connection: Mapping[tuple[Hashable, Hashable], np.ndarray],
    edge_coframes: Mapping[tuple[Hashable, Hashable], np.ndarray],
    frame_changes: Mapping[Hashable, np.ndarray],
) -> tuple[dict, dict, dict]:
    """Apply arbitrary invertible local frame changes at vertices.

    G'_v=H_v^{-T}G_vH_v^{-1}; U'_(v<-u)=H_v U H_u^{-1};
    theta'_(u,v)=H_u theta_(u,v) for canonical oriented edges.
    """
    h_inv = {v: np.linalg.inv(h) for v, h in frame_changes.items()}
    new_metrics = {
        v: h_inv[v].T @ np.asarray(g) @ h_inv[v]
        for v, g in metrics.items()
    }
    new_connection = {
        (source, target): frame_changes[target] @ np.asarray(u) @ h_inv[source]
        for (source, target), u in connection.items()
    }
    new_edge_coframes = {
        edge: frame_changes[edge[0]] @ np.asarray(theta)
        for edge, theta in edge_coframes.items()
    }
    return new_metrics, new_connection, new_edge_coframes


def _path_transport(
    path: tuple[Hashable, ...],
    connection: Mapping[tuple[Hashable, Hashable], np.ndarray],
    dimension: int,
) -> np.ndarray:
    result = np.eye(dimension)
    for source, target in zip(path[:-1], path[1:]):
        result = np.asarray(connection[(source, target)]) @ result
    return result


def face_cartan_data(
    face: tuple[Hashable, Hashable, Hashable, Hashable],
    edge_coframes: Mapping[tuple[Hashable, Hashable], np.ndarray],
    connection: Mapping[tuple[Hashable, Hashable], np.ndarray],
) -> dict[str, np.ndarray | float]:
    """Return discrete torsion 2-cochain and curvature holonomy for one square."""
    if len(set(face)) != 4:
        raise ValueError("A square face must contain four distinct vertices.")
    sample = np.asarray(next(iter(edge_coframes.values())))
    theta_dim = sample.shape[0]
    connection_dim = next(iter(connection.values())).shape[0]
    torsion = np.zeros(theta_dim, dtype=float)
    holonomy = np.eye(connection_dim)
    for k in range(4):
        source, target = face[k], face[(k + 1) % 4]
        edge = (source, target) if source < target else (target, source)
        theta = np.asarray(edge_coframes[edge], dtype=float)
        if source < target:
            theta_at_source = theta
        else:
            theta_at_source = -np.asarray(connection[(edge[0], edge[1])]) @ theta
        path_to_base = (
            (source,) if k == 0
            else tuple([source] + [face[j] for j in range(k - 1, -1, -1)])
        )
        transport_to_base = _path_transport(path_to_base, connection, connection_dim)
        torsion += transport_to_base @ theta_at_source
        holonomy = np.asarray(connection[(source, target)]) @ holonomy
    curvature = holonomy - np.eye(connection_dim)
    return {
        "torsion": torsion,
        "curvature_holonomy": holonomy,
        "curvature_deviation_fro": float(np.linalg.norm(curvature, ord="fro")),
        "torsion_norm": float(np.linalg.norm(torsion)),
    }


def edge_coframes_from_embedding(
    points: np.ndarray,
    edges: tuple[tuple[int, int], ...],
    append_radial_difference: bool = False,
) -> dict[tuple[int, int], np.ndarray]:
    """Create canonical edge 1-coframes from an explicit vertex embedding."""
    x = np.asarray(points, dtype=float)
    if x.ndim != 2 or not np.all(np.isfinite(x)):
        raise ValueError("points must be a finite (N,d) array.")
    result = {}
    for source, target in edges:
        if source >= target:
            raise ValueError("Edges must be in canonical increasing order.")
        delta = x[target] - x[source]
        if append_radial_difference:
            delta = np.concatenate((
                delta,
                [np.linalg.norm(x[target]) - np.linalg.norm(x[source])],
            ))
        result[(source, target)] = delta
    return result


def audit_cartan_faces(
    faces: tuple[tuple[int, int, int, int], ...],
    edge_coframes: Mapping[tuple[int, int], np.ndarray],
    connection: Mapping[tuple[int, int], np.ndarray],
) -> dict[str, object]:
    data = [face_cartan_data(face, edge_coframes, connection) for face in faces]
    return {
        "audit": "DISCRETE_CARTAN_FACE_TORSION_CURVATURE",
        "faces": len(data),
        "max_torsion_norm": max((float(item["torsion_norm"]) for item in data), default=0.0),
        "max_curvature_deviation_fro": max(
            (float(item["curvature_deviation_fro"]) for item in data), default=0.0
        ),
        "definition": "face torsion is parallel-transported oriented edge-coframe sum; curvature is ordered face holonomy minus identity",
        "scope": "explicit finite Cartan discretisation; not automatically Levi-Civita or TEGR",
    }
