"""Tests for metric-compatible discrete connections and Cartan face data."""
import numpy as np
import pytest

from src.discrete_cartan_geometry import (
    audit_cartan_faces, audit_metric_compatibility,
    build_metric_compatible_connection, edge_coframes_from_embedding,
    face_cartan_data, metric_compatible_edge_transport, transform_local_frame_data,
)
from src.q6_complex import build_q6_complex


def test_metric_compatible_transport_and_reverse_are_exact_to_roundoff():
    rng = np.random.default_rng(41)
    a = rng.normal(size=(4, 4))
    b = rng.normal(size=(4, 4))
    gs = a.T @ a + np.eye(4)
    gt = b.T @ b + np.eye(4)
    u = metric_compatible_edge_transport(gs, gt)
    assert np.linalg.norm(u.T @ gt @ u - gs) < 2e-12
    assert np.linalg.norm(np.linalg.inv(u) @ u - np.eye(4)) < 2e-12


def test_q6_edge_connection_is_metric_compatible_on_all_192_edges():
    q6 = build_q6_complex()
    rng = np.random.default_rng(12)
    metrics = {}
    for vertex in q6.vertices:
        a = rng.normal(size=(4, 4))
        metrics[vertex] = a.T @ a + np.eye(4)
    connection = build_metric_compatible_connection(metrics, q6.edges)
    report = audit_metric_compatibility(metrics, connection, q6.edges)
    assert report["edge_count"] == 192
    assert report["metric_compatibility_max_fro"] < 2e-11
    assert report["inverse_transport_max_fro"] < 2e-11


def test_cartan_torsion_and_curvature_vanish_for_exact_flat_embedding():
    q6 = build_q6_complex()
    points = np.column_stack((np.arange(64, dtype=float), np.zeros(64),
                              np.zeros(64), np.zeros(64)))
    theta = edge_coframes_from_embedding(points, q6.edges)
    metrics = {v: np.eye(4) for v in q6.vertices}
    connection = build_metric_compatible_connection(metrics, q6.edges)
    report = audit_cartan_faces(q6.faces, theta, connection)
    assert report["faces"] == 240
    assert report["max_torsion_norm"] == 0.0
    assert report["max_curvature_deviation_fro"] == 0.0


def test_discrete_cartan_data_transform_covariantly_under_local_frame_changes():
    q6 = build_q6_complex()
    points = np.random.default_rng(3).normal(size=(64, 4))
    theta = edge_coframes_from_embedding(points, q6.edges)
    rng = np.random.default_rng(8)
    metrics, frames = {}, {}
    for vertex in q6.vertices:
        a = rng.normal(size=(4, 4))
        metrics[vertex] = a.T @ a + np.eye(4)
        q, _ = np.linalg.qr(rng.normal(size=(4, 4)))
        if np.linalg.det(q) < 0:
            q[:, 0] *= -1
        frames[vertex] = q
    connection = build_metric_compatible_connection(metrics, q6.edges)
    face = q6.faces[0]
    before = face_cartan_data(face, theta, connection)
    g2, u2, theta2 = transform_local_frame_data(metrics, connection, theta, frames)
    after = face_cartan_data(face, theta2, u2)
    expected_torsion = frames[face[0]] @ before["torsion"]
    expected_curvature = frames[face[0]] @ before["curvature_holonomy"] @ np.linalg.inv(frames[face[0]])
    assert np.linalg.norm(after["torsion"] - expected_torsion) < 2e-10
    assert np.linalg.norm(after["curvature_holonomy"] - expected_curvature) < 2e-10
    assert audit_metric_compatibility(g2, u2, q6.edges)["metric_compatibility_max_fro"] < 2e-10


def test_reflection_frame_changes_do_not_create_a_pin_lift():
    q6 = build_q6_complex()
    metrics = {v: np.eye(4) for v in q6.vertices}
    theta = edge_coframes_from_embedding(np.random.default_rng(2).normal(size=(64, 4)), q6.edges)
    connection = build_metric_compatible_connection(metrics, q6.edges)
    reflections = {v: np.diag([-1.0, 1.0, 1.0, 1.0]) for v in q6.vertices}
    transformed = transform_local_frame_data(metrics, connection, theta, reflections)
    assert all(np.linalg.det(h) == pytest.approx(-1.0) for h in reflections.values())
    assert audit_metric_compatibility(transformed[0], transformed[1], q6.edges)["metric_compatibility_max_fro"] < 1e-12
