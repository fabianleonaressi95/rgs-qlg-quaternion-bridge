import numpy as np

from src.rgs_psu_cut_project import (
    PHI,
    PHI_STAR,
    audit_cut_project_frame,
    audit_model_set_patch,
    audit_torsion_coupled_internal_time_axis,
    build_cut_project_data,
    enumerate_model_set_patch,
    golden_cut_project_frame,
    icosahedral_axis_matrix,
)


def test_golden_conjugates_and_orthogonal_three_plus_three_split():
    assert np.isclose(PHI * PHI_STAR, -1.0, atol=1e-15)
    pp, pi, frame = golden_cut_project_frame()
    assert pp.shape == (3, 6)
    assert pi.shape == (3, 6)
    assert frame.shape == (6, 6)
    assert np.allclose(pp @ pp.T, np.eye(3), atol=1e-14)
    assert np.allclose(pi @ pi.T, np.eye(3), atol=1e-14)
    assert np.allclose(pp @ pi.T, 0.0, atol=1e-14)
    assert np.allclose(frame @ frame.T, np.eye(6), atol=1e-14)


def test_projection_matrix_is_the_recovered_icosahedral_convention():
    a = icosahedral_axis_matrix(PHI)
    expected = np.asarray([
        [0, 0, 1, -1, PHI, PHI],
        [1, -1, PHI, PHI, 0, 0],
        [PHI, PHI, 0, 0, 1, -1],
    ])
    assert np.array_equal(a, expected)


def test_projected_cube_is_the_rhombic_triacontahedral_window():
    data = build_cut_project_data()
    report = audit_cut_project_frame(data)
    assert report["pass"]
    assert report["acceptance_window"]["vertices"] == 32
    assert report["acceptance_window"]["facets"] == 30
    assert report["q6_seed"]["vertices"] == 64
    assert report["q6_seed"]["all_internal_points_inside_window"]


def test_finite_model_set_patch_is_inside_both_cuts():
    data = build_cut_project_data()
    patch = enumerate_model_set_patch(physical_radius=2.0, data=data)
    audit = audit_model_set_patch(patch)
    assert len(patch.lattice_points) > 0
    assert np.all(np.linalg.norm(patch.physical_points, axis=1) <= 2.0 + 1e-10)
    assert audit["rounded_projection_duplicates"] == 0
    assert audit["minimum_physical_separation"] > 0.0


def test_torsion_coupled_internal_time_axis_is_recoverable_and_stable():
    rng = np.random.default_rng(20261021)
    y = rng.normal(size=(64, 3))
    u = np.asarray([0.3, -0.4, 0.866025403784])
    u /= np.linalg.norm(u)
    latent = y @ u
    t = np.column_stack((
        latent + 0.05 * rng.normal(size=64),
        0.2 * rng.normal(size=64),
        -0.1 * latent + 0.05 * rng.normal(size=64),
        rng.normal(size=64),
    ))
    result = audit_torsion_coupled_internal_time_axis(
        y, t, n_permutations=99, n_bootstrap=100, seed=7
    )
    recovered = np.asarray(result["internal_time_axis_u_tau"])
    assert abs(float(recovered @ u)) > 0.98
    assert result["leading_coupling_strength"] > 0.0
    assert result["permutation_null"]["p_value"] <= 0.05
    assert result["bootstrap"]["n_valid"] > 80


def test_time_axis_rejects_missing_or_constant_torsion_data():
    y = np.arange(192, dtype=float).reshape(64, 3)
    with np.testing.assert_raises(ValueError):
        audit_torsion_coupled_internal_time_axis(y, np.ones((64, 4)))


def test_model_set_window_scaling_changes_the_finite_selection():
    data = build_cut_project_data()
    small = enumerate_model_set_patch(2.0, data, window_scale=0.8)
    large = enumerate_model_set_patch(2.0, data, window_scale=1.2)
    assert len(large.lattice_points) >= len(small.lattice_points)
