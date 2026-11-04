"""Tests for the coupled finite-dimensional Q6 action prototype."""
import numpy as np
import pytest

from src.holonomies import SU2Lattice
from src.unified_discrete_action import (
    ActionCouplings,
    UnifiedQ6Fields,
    action_components,
    gauge_transform_fields,
    identity_q6_links,
    unified_q6_action,
)
from src.q6_complex import build_q6_complex


def _su2(axis, angle):
    axis = np.asarray(axis, dtype=float)
    axis /= np.linalg.norm(axis)
    x, y, z = axis
    c, s = np.cos(angle / 2), np.sin(angle / 2)
    return np.array([
        [c - 1j*s*z, -s*y - 1j*s*x],
        [s*y - 1j*s*x, c + 1j*s*z],
    ], dtype=complex)


def _fields(seed=4):
    rng = np.random.default_rng(seed)
    e = np.repeat(np.eye(4)[None, :, :], 64, axis=0)
    e += 0.03 * rng.normal(size=e.shape)
    p = rng.normal(size=(64, 2)) + 1j*rng.normal(size=(64, 2))
    m = rng.normal(size=(64, 2)) + 1j*rng.normal(size=(64, 2))
    return UnifiedQ6Fields(e, p, m)


def _random_links(seed=9):
    q6 = build_q6_complex()
    rng = np.random.default_rng(seed)
    links = SU2Lattice(q6.vertices)
    for u, v in q6.edges:
        axis = rng.normal(size=3)
        angle = float(rng.uniform(-0.8, 0.8))
        links.add_link(u, v, _su2(axis, angle))
    return links


def test_identity_configuration_has_zero_action():
    e = np.repeat(np.eye(4)[None, :, :], 64, axis=0)
    zero = np.zeros((64, 2), dtype=complex)
    fields = UnifiedQ6Fields(e, zero, zero)
    parts = action_components(fields, identity_q6_links())
    assert all(abs(value) < 1e-12 for value in parts.values())


def test_action_is_nonnegative_for_nonnegative_couplings():
    parts = action_components(_fields(), _random_links())
    assert all(value >= -1e-10 for value in parts.values())


def test_local_su2_gauge_invariance_of_full_action():
    fields = _fields()
    links = _random_links()
    rng = np.random.default_rng(17)
    gauge = {
        v: _su2(rng.normal(size=3), float(rng.uniform(-1.0, 1.0)))
        for v in range(64)
    }
    transformed_links = links.gauge_transform(gauge)
    transformed_fields = gauge_transform_fields(fields, gauge)
    before = unified_q6_action(fields, links)
    after = unified_q6_action(transformed_fields, transformed_links)
    assert np.isclose(before, after, rtol=1e-11, atol=1e-10)


def test_action_decomposition_sums_to_total():
    parts = action_components(_fields(), _random_links())
    assert np.isclose(
        parts["total"],
        sum(parts[k] for k in ("coframe", "wilson", "hopping", "twin", "quartic")),
    )


def test_twin_term_vanishes_for_identical_sheets():
    fields = _fields()
    same = UnifiedQ6Fields(fields.coframe, fields.psi_plus, fields.psi_plus)
    parts = action_components(same, identity_q6_links())
    assert parts["twin"] == pytest.approx(0.0, abs=1e-12)


def test_degenerate_coframe_is_rejected():
    fields = _fields()
    e = fields.coframe.copy()
    e[0, 0] = e[0, 1]
    bad = UnifiedQ6Fields(e, fields.psi_plus, fields.psi_minus)
    with pytest.raises(ValueError, match="non-degenerate"):
        unified_q6_action(bad, identity_q6_links())


def test_wrong_link_graph_is_rejected():
    q6 = build_q6_complex()
    links = SU2Lattice(q6.vertices)
    with pytest.raises(ValueError, match="coincide exactly"):
        unified_q6_action(_fields(), links)


def test_negative_coupling_is_rejected():
    with pytest.raises(ValueError, match="non-negative"):
        ActionCouplings(twin_coupling=-1.0)
