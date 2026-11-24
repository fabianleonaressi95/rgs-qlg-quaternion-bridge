import numpy as np

from src.holonomies import SU2Lattice
from src.q6_complex import build_q6_complex
from src.rgs_ellis_refinement import load_rgs_seed_data
from src.q6_metric_wilson_response import (
    face_area_log_shear_derivative,
    face_parallelogram_areas,
    metric_wilson_shear_finite_difference,
    metric_wilson_shear_response,
    shear_generator_3d,
)


def random_links(seed=240624):
    q6 = build_q6_complex()
    rng = np.random.default_rng(seed)
    lattice = SU2Lattice(q6.vertices)
    for u, v in q6.edges:
        z = rng.normal(size=4)
        z /= np.linalg.norm(z)
        a, b, c, d = z
        U = np.array([[a+1j*b, c+1j*d], [-c+1j*d, a-1j*b]])
        lattice.add_link(u, v, U)
    return lattice


def test_face_area_log_derivative_matches_finite_difference():
    x = np.asarray(load_rgs_seed_data().coordinates)
    h = shear_generator_3d(0, 1)
    exact = face_area_log_shear_derivative(x, h)
    eps = 1e-4
    ap = face_parallelogram_areas(((np.eye(3)+eps*h) @ x.T).T)
    am = face_parallelogram_areas(((np.eye(3)-eps*h) @ x.T).T)
    numeric = (np.log(ap)-np.log(am))/(2*eps)
    assert np.allclose(exact, numeric, rtol=2e-6, atol=2e-7)


def test_metric_wilson_stress_matches_action_finite_difference():
    x = np.asarray(load_rgs_seed_data().coordinates)
    links = random_links()
    exact = metric_wilson_shear_response(links, 2.75, beta=1.3, coordinates=x)
    numeric = metric_wilson_shear_finite_difference(
        links, 2.75, beta=1.3, coordinates=x, step=1e-4,
    )
    assert np.isclose(exact["Pi_xy"], numeric, rtol=3e-6, atol=3e-7)


def test_reference_configuration_has_unit_face_weights():
    x = np.asarray(load_rgs_seed_data().coordinates)
    from src.q6_metric_wilson_response import metric_weighted_plaquette_action
    result = metric_weighted_plaquette_action(x, random_links())
    assert np.allclose(result["weights"], 1.0)


def test_shear_generator_is_symmetric_traceless():
    h = shear_generator_3d()
    assert np.allclose(h, h.T)
    assert abs(np.trace(h)) < 1e-14


def test_metric_p1_matter_stress_matches_action_finite_difference():
    from src.unified_discrete_action import UnifiedQ6Fields
    from src.q6_metric_wilson_response import (
        metric_edge_matter_action,
        metric_edge_matter_shear_response,
        metric_edge_matter_shear_finite_difference,
    )
    rng = np.random.default_rng(91824)
    fields = UnifiedQ6Fields(
        np.broadcast_to(np.eye(4), (64,4,4)).copy(),
        rng.normal(size=(64,2)) + 1j*rng.normal(size=(64,2)),
        rng.normal(size=(64,2)) + 1j*rng.normal(size=(64,2)),
    )
    x = np.asarray(load_rgs_seed_data().coordinates)
    links = random_links(91825)
    action = metric_edge_matter_action(x, fields, links, coupling=0.73)
    assert action >= 0.0
    exact = metric_edge_matter_shear_response(
        fields, links, 2.75, coupling=0.73, coordinates=x,
    )
    numeric = metric_edge_matter_shear_finite_difference(
        fields, links, 2.75, coupling=0.73, coordinates=x, step=1e-6,
    )
    assert np.isclose(exact["Pi_xy"], numeric, rtol=3e-6, atol=3e-7)


def test_combined_candidate_stress_is_sum_of_gauge_and_matter_sectors():
    from src.unified_discrete_action import UnifiedQ6Fields
    from src.q6_metric_wilson_response import (
        q6_candidate_total_shear_response,
        metric_wilson_shear_response,
        metric_edge_matter_shear_response,
    )
    rng = np.random.default_rng(91826)
    fields = UnifiedQ6Fields(
        np.broadcast_to(np.eye(4), (64,4,4)).copy(),
        rng.normal(size=(64,2)) + 1j*rng.normal(size=(64,2)),
        rng.normal(size=(64,2)) + 1j*rng.normal(size=(64,2)),
    )
    x = np.asarray(load_rgs_seed_data().coordinates)
    links = random_links(91827)
    total = q6_candidate_total_shear_response(
        fields, links, 2.5, beta=0.9, matter_coupling=0.6, coordinates=x,
    )
    gauge = metric_wilson_shear_response(links, 2.5, beta=0.9, coordinates=x)
    matter = metric_edge_matter_shear_response(
        fields, links, 2.5, coupling=0.6, coordinates=x,
    )
    assert np.isclose(total["Pi_total"], gauge["Pi_xy"] + matter["Pi_xy"])
