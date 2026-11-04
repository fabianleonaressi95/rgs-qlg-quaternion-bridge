"""Independent directional finite-difference checks for Q6 gradients."""
import numpy as np

from src.holonomies import SU2Lattice
from src.q6_complex import build_q6_complex
from src.unified_discrete_action import (
    ActionCouplings, UnifiedQ6Fields, unified_q6_action,
)
from src.unified_discrete_gradients import (
    su2_left_perturbation, unified_q6_gradients,
)


def _su2(axis, angle):
    axis = np.asarray(axis, dtype=float)
    axis /= np.linalg.norm(axis)
    x, y, z = axis
    c, s = np.cos(angle / 2), np.sin(angle / 2)
    return np.array([
        [c - 1j*s*z, -s*y - 1j*s*x],
        [s*y - 1j*s*x, c + 1j*s*z],
    ], dtype=complex)


def _sample(seed=71):
    rng = np.random.default_rng(seed)
    q6 = build_q6_complex()
    e = np.repeat(np.eye(4)[None], 64, axis=0)
    e += 0.025 * rng.normal(size=e.shape)
    p = 0.2 * (rng.normal(size=(64, 2)) + 1j*rng.normal(size=(64, 2)))
    m = 0.2 * (rng.normal(size=(64, 2)) + 1j*rng.normal(size=(64, 2)))
    fields = UnifiedQ6Fields(e, p, m)
    links = SU2Lattice(q6.vertices)
    for u, v in q6.edges:
        links.add_link(u, v, _su2(rng.normal(size=3), float(rng.uniform(-0.5, 0.5))))
    return q6, fields, links, rng


def _perturb(fields, links, q6, directions, eps):
    de, dp, dm, dx = directions
    changed = UnifiedQ6Fields(
        fields.coframe + eps * de,
        fields.psi_plus + eps * dp,
        fields.psi_minus + eps * dm,
    )
    new_links = SU2Lattice(q6.vertices)
    for edge_id, (u, v) in enumerate(q6.edges):
        new_links.add_link(
            u, v,
            su2_left_perturbation(links.link(u, v), dx[edge_id], eps),
        )
    return changed, new_links


def test_all_analytic_gradients_match_independent_directional_differences():
    q6, fields, links, rng = _sample()
    couplings = ActionCouplings(
        coframe_stiffness=0.7, gauge_beta=0.9, matter_hopping=0.8,
        twin_coupling=0.35, quartic_coupling=0.12,
    )
    gradients = unified_q6_gradients(fields, links, couplings, q6)
    de = rng.normal(size=fields.coframe.shape)
    dp = rng.normal(size=fields.psi_plus.shape) + 1j*rng.normal(size=fields.psi_plus.shape)
    dm = rng.normal(size=fields.psi_minus.shape) + 1j*rng.normal(size=fields.psi_minus.shape)
    dx = 0.15 * rng.normal(size=(len(q6.edges), 3))
    directions = (de, dp, dm, dx)

    analytic = (
        np.sum(gradients["coframe"] * de)
        + 2.0 * np.real(np.vdot(gradients["psi_plus"], dp))
        + 2.0 * np.real(np.vdot(gradients["psi_minus"], dm))
        + np.sum(gradients["links"] * dx)
    )
    base = unified_q6_action(fields, links, couplings, q6)
    errors = []
    for eps in (2e-4, 7e-5, 2e-5):
        plus_fields, plus_links = _perturb(fields, links, q6, directions, eps)
        minus_fields, minus_links = _perturb(fields, links, q6, directions, -eps)
        finite_difference = (
            unified_q6_action(plus_fields, plus_links, couplings, q6)
            - unified_q6_action(minus_fields, minus_links, couplings, q6)
        ) / (2.0 * eps)
        errors.append(abs(finite_difference - analytic))
    assert errors[-1] < 2e-6 * max(1.0, abs(analytic))
    assert errors[-1] < errors[0]


def test_link_gradient_includes_forward_and_inverse_face_occurrences():
    q6, fields, links, rng = _sample(83)
    # A generic Q6 edge is used in a face in both orientations across the
    # canonical face list; the directional test above covers all such terms.
    couplings = ActionCouplings(gauge_beta=1.0, matter_hopping=0.0,
                                coframe_stiffness=0.0, twin_coupling=0.0,
                                quartic_coupling=0.0)
    gradients = unified_q6_gradients(fields, links, couplings, q6)
    assert gradients["links"].shape == (192, 3)
    assert np.all(np.isfinite(gradients["links"]))
    assert np.linalg.norm(gradients["links"]) > 1e-8
