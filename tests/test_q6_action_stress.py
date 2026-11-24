import numpy as np

from src.holonomies import SU2Lattice
from src.q6_complex import build_q6_complex
from src.unified_discrete_action import ActionCouplings, UnifiedQ6Fields
from src.q6_action_stress import (
    q6_action_shear_response,
    q6_action_shear_response_finite_difference,
    symmetric_shear_generator,
)


def sample_model():
    q6 = build_q6_complex()
    rng = np.random.default_rng(2406)
    coframe = np.broadcast_to(np.eye(4), (64, 4, 4)).copy()
    coframe += 0.025 * rng.normal(size=(64, 4, 4))
    psi_plus = rng.normal(size=(64, 2)) + 1j * rng.normal(size=(64, 2))
    psi_minus = rng.normal(size=(64, 2)) + 1j * rng.normal(size=(64, 2))
    links = SU2Lattice(q6.vertices)
    for u, v in q6.edges:
        links.add_link(u, v, np.eye(2, dtype=complex))
    return UnifiedQ6Fields(coframe, psi_plus, psi_minus), links


def test_symmetric_shear_is_traceless_and_default_axes_are_spatial():
    h = symmetric_shear_generator()
    assert np.allclose(h, h.T)
    assert np.isclose(np.trace(h), 0.0)
    assert h[0, 0] == h[0, 1] == h[0, 2] == h[0, 3] == 0.0
    assert h[1, 2] == h[2, 1] == 0.5


def test_analytic_action_shear_response_matches_central_difference():
    fields, links = sample_model()
    couplings = ActionCouplings(
        coframe_stiffness=0.7, gauge_beta=1.3, matter_hopping=0.8,
        twin_coupling=0.4, quartic_coupling=0.2,
    )
    exact = q6_action_shear_response(fields, links, 3.25, couplings=couplings)
    numeric = q6_action_shear_response_finite_difference(
        fields, links, 3.25, couplings=couplings, step=2e-6,
    )
    assert np.isclose(exact["Pi_action"], numeric, rtol=2e-6, atol=2e-7)


def test_metric_independent_wilson_and_volume_only_terms_have_zero_shear_part():
    fields, links = sample_model()
    couplings = ActionCouplings(
        coframe_stiffness=0.0, gauge_beta=1.0, matter_hopping=1.0,
        twin_coupling=1.0, quartic_coupling=1.0,
    )
    result = q6_action_shear_response(fields, links, 1.0, couplings=couplings)
    assert abs(result["Pi_wilson"]) < 1e-12
    assert abs(result["Pi_hopping"]) < 1e-12
    assert abs(result["Pi_twin"]) < 1e-12
    assert abs(result["Pi_quartic"]) < 1e-12
    assert abs(result["Pi_action"]) < 1e-12


def test_bad_reference_volume_is_rejected():
    fields, links = sample_model()
    for volume in (0.0, -1.0, np.nan):
        try:
            q6_action_shear_response(fields, links, volume)
        except ValueError:
            pass
        else:
            raise AssertionError("invalid reference volume must fail")


def test_inverse_square_edge_weight_derivative_matches_finite_difference():
    from src.q6_action_stress import inverse_square_edge_weight_shear_derivative
    rng = np.random.default_rng(192240)
    d = rng.normal(size=(17, 3))
    analytic = inverse_square_edge_weight_shear_derivative(d)
    h = symmetric_shear_generator(axis_a=0, axis_b=1, dimension=3)
    eps = 1e-6
    transform_plus = np.eye(3) + eps * h
    transform_minus = np.eye(3) - eps * h
    w_plus = 1.0 / np.sum((d @ transform_plus.T) ** 2, axis=1)
    w_minus = 1.0 / np.sum((d @ transform_minus.T) ** 2, axis=1)
    numeric = (w_plus - w_minus) / (2.0 * eps)
    assert np.allclose(analytic, numeric, rtol=2e-6, atol=2e-7)
