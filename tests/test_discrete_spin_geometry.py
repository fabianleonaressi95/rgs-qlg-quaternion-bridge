"""Tests for Q6 frame alignment and Euclidean Spin(4)/Clifford transport."""
import numpy as np
import pytest
from scipy.linalg import expm

from src.discrete_spin_geometry import (
    build_q6_spin_transports,
    clifford_algebra_residual,
    edge_frame_spin_transport,
    euclidean_gamma_matrices,
    reverse_spin_transport,
)
from src.q6_complex import build_q6_complex
from src.unified_discrete_action import UnifiedQ6Fields


def _rotation(seed=3):
    rng = np.random.default_rng(seed)
    a = rng.normal(size=(4, 4))
    omega = 0.12 * (a - a.T)
    return expm(omega)


def test_euclidean_gamma_matrices_satisfy_clifford_algebra():
    assert clifford_algebra_residual() < 1e-13


def test_spin_lift_intertwines_a_nontrivial_so4_rotation():
    rng = np.random.default_rng(12)
    e0 = np.eye(4) + 0.03 * rng.normal(size=(4, 4))
    if np.linalg.det(e0) < 0:
        e0[0] *= -1
    rotation = _rotation()
    e1 = rotation @ e0
    result = edge_frame_spin_transport(e0, e1, 4, 9)
    assert np.allclose(result.vector_rotation, rotation, atol=2e-12)
    assert result.metric_mismatch < 2e-12
    assert result.coframe_alignment_residual < 2e-12
    assert result.clifford_residual < 3e-7
    assert np.allclose(
        result.spin_transport.conj().T @ result.spin_transport,
        np.eye(4), atol=2e-8,
    )


def test_reverse_transport_is_inverse():
    e0 = np.eye(4)
    e1 = _rotation(8)
    forward = edge_frame_spin_transport(e0, e1)
    backward = edge_frame_spin_transport(e1, e0)
    assert np.allclose(
        backward.spin_transport, reverse_spin_transport(forward), atol=2e-7
    )
    assert np.allclose(
        backward.vector_rotation, forward.vector_rotation.T, atol=2e-12
    )


def test_q6_builder_returns_all_192_edge_transports():
    q6 = build_q6_complex()
    e = np.repeat(np.eye(4)[None, :, :], 64, axis=0)
    fields = UnifiedQ6Fields(
        e,
        np.zeros((64, 2), dtype=complex),
        np.zeros((64, 2), dtype=complex),
    )
    transports = build_q6_spin_transports(fields, q6)
    assert len(transports) == 192
    assert all(item.clifford_residual < 3e-7 for item in transports.values())


def test_orientation_reversing_or_degenerate_coframes_are_rejected():
    with pytest.raises(ValueError, match="orientation"):
        edge_frame_spin_transport(np.eye(4), np.diag([-1.0, 1.0, 1.0, 1.0]))
    with pytest.raises(ValueError, match="non-degenerate"):
        edge_frame_spin_transport(np.eye(4), np.diag([1.0, 1.0, 1.0, 0.0]))


def test_metric_mismatch_is_reported_not_hidden():
    e0 = np.eye(4)
    e1 = np.diag([1.1, 1.0, 1.0, 1.0])
    result = edge_frame_spin_transport(e0, e1)
    assert result.metric_mismatch > 0.0
    assert result.coframe_alignment_residual > 0.0
