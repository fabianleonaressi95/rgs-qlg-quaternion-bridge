"""Integration tests for Möbius/Q6, gauge, Clifford and multiscale audits."""
import numpy as np
import pytest

from src.discrete_differential import audit_discrete_differential_operators
from src.mobius_complex import build_mobius_complex
from src.mobius_gauge_audit import audit_mobius_gauge_covariance
from src.q6_complex import build_q6_complex
from src.q6_prism_map import audit_q6_prism_map, gray_code, inverse_gray_code
from src.rgs_multiscale import (
    equidistribution_moment_audit,
    finite_mosco_diagnostic,
    metric_correspondence_distortion,
    multifractal_partition_spectrum,
    sphere_moment_reference,
)
from src.binary_icosahedral import binary_icosahedral_elements
from src.clifford_qlg import clifford_grade_action_audit
from scripts.rgs_psu_integrated_audit import json_safe


def test_mobius_is_six_face_subcomplex_not_the_full_q6_240_face_complex():
    mobius, q6 = build_mobius_complex(), build_q6_complex()
    q6_faces = {frozenset(face) for face in q6.faces}
    assert (len(mobius.vertices), len(mobius.edges), len(mobius.faces)) == (12, 18, 6)
    assert (len(q6.vertices), len(q6.edges), len(q6.faces)) == (64, 192, 240)
    assert all(frozenset(face) in q6_faces for face in mobius.faces)
    assert len({frozenset(face) for face in mobius.faces}) == 6
    assert len(q6_faces) == 240


@pytest.mark.parametrize("seed", [7, 17, 31, 20261103])
def test_seeded_randomized_open_path_and_closed_loop_gauge_covariance(seed):
    report = audit_mobius_gauge_covariance(seed=seed, atol=1e-10)
    assert report["passed"] is True
    assert report["open_path_covariance_max"] <= 1e-10
    assert report["closed_loop_conjugation_max"] <= 1e-10
    assert report["wilson_trace_invariance_max"] <= 1e-10


def test_q6_to_prism_gray_map_is_injective_and_preserves_a_hamiltonian_cycle():
    assert all(inverse_gray_code(gray_code(i)) == i for i in range(64))
    report = audit_q6_prism_map()
    assert report["injective"] is True
    assert report["image_vertices"] == 64
    assert report["preserved_hamiltonian_cycle_edges"] == 64
    assert report["q6_edge_preservation_fraction"] == pytest.approx(1/3)
    assert report["preserves_all_q6_adjacencies"] is False


def test_discrete_de_rham_hodge_operators_on_mobius_and_q6():
    for complex_ in (build_mobius_complex(), build_q6_complex()):
        report = audit_discrete_differential_operators(
            complex_.boundary_1, complex_.boundary_2
        )
        assert report["passed"] is True
        assert report["d1_d0_nilpotency_fro"] == 0.0
        assert report["laplacian_min_eigenvalue"] >= -1e-12


def test_sphere_haar_moment_reference_values_in_four_dimensions():
    assert sphere_moment_reference((2, 0, 0, 0)) == pytest.approx(1/4)
    assert sphere_moment_reference((4, 0, 0, 0)) == pytest.approx(1/8)
    assert sphere_moment_reference((2, 2, 0, 0)) == pytest.approx(1/24)
    assert sphere_moment_reference((1, 1, 0, 0)) == 0.0


def test_finite_2i_moment_audit_is_labelled_as_finite_not_equdistributive_proof():
    result = equidistribution_moment_audit(
        binary_icosahedral_elements(), ((2, 0, 0, 0), (4, 0, 0, 0))
    )
    assert result["point_count"] == 120
    assert "not Weyl equidistribution" in result["scope"]


def test_multifractal_uniform_binary_partitions_have_dimension_one():
    scales = np.array([1/2, 1/4, 1/8, 1/16], dtype=float)
    partitions = [np.full(2**k, 1/(2**k)) for k in range(1, 5)]
    result = multifractal_partition_spectrum(scales, partitions, q_values=(0.0, 1.0, 2.0))
    assert result["spectrum"]["0.0"]["generalized_dimension"] == pytest.approx(1.0)
    assert result["spectrum"]["1.0"]["generalized_dimension"] == pytest.approx(1.0)
    assert result["spectrum"]["2.0"]["generalized_dimension"] == pytest.approx(1.0)


def test_finite_metric_correspondence_and_mosco_diagnostic_do_not_claim_theorems():
    d = np.array([[0.0, 1.0], [1.0, 0.0]])
    gh = metric_correspondence_distortion(d, d, np.array([[0, 0], [1, 1]]))
    assert gh["distortion"] == 0.0
    assert gh["half_distortion_bound"] == 0.0
    mosco = finite_mosco_diagnostic(np.array([1.0, 0.8]), np.array([1.0, 0.9]))
    assert mosco["finite_samples_nonnegative"] is True
    assert mosco["mosco_convergence_proved"] is False


def test_clifford_grade_action_preserves_grades_under_2i_spin_rotations():
    report = clifford_grade_action_audit()
    assert report["passed"] is True
    assert report["max_grade_leakage"] < 3e-9


def test_integrated_audit_json_safe_preserves_complex_components():
    payload = {"matrix": np.array([[1 + 2j, 3 - 4j]], dtype=complex),
               "nested": [np.complex128(5 + 6j), np.float64(2.5)]}
    safe = json_safe(payload)
    import json
    encoded = json.dumps(safe, allow_nan=False)
    decoded = json.loads(encoded)
    assert decoded["matrix"][0][0] == {"real": 1.0, "imag": 2.0}
    assert decoded["matrix"][0][1] == {"real": 3.0, "imag": -4.0}
    assert decoded["nested"][0] == {"real": 5.0, "imag": 6.0}
    assert decoded["nested"][1] == 2.5
