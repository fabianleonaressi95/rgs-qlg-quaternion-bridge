import numpy as np
import pytest

from src.rgs_ellis_refinement import (
    assemble_ellis_p1_operators,
    audit_ellis_rgs_refinement,
    build_ellis_rgs_refinement,
    load_rgs_seed_data,
)


def test_recovered_rgs_seed_and_q6_coframe_are_reconstructed():
    seed = load_rgs_seed_data()
    assert seed.coordinates.shape == (64, 3)
    assert seed.q6_edges.shape == (192, 3)
    assert seed.q6_raw_coframe.shape == (192, 4)
    assert seed.q6_unit_coframe.shape == (192, 4)
    assert seed.directions.shape == (56, 3)
    assert seed.vertex_to_direction.shape == (64,)
    assert np.allclose(np.linalg.norm(seed.directions, axis=1), 1.0)
    assert np.allclose(np.linalg.norm(seed.q6_unit_coframe, axis=1), 1.0)
    antipodal_distance = np.linalg.norm(
        seed.coordinates[:, None, :] + seed.coordinates[None, :, :], axis=2
    )
    antipodal = np.argmin(antipodal_distance, axis=1)
    assert np.array_equal(np.sort(antipodal), np.arange(64))
    assert np.max(antipodal_distance[np.arange(64), antipodal]) < 1e-14
    degree = (
        np.bincount(seed.q6_edges[:, 0], minlength=64)
        + np.bincount(seed.q6_edges[:, 1], minlength=64)
    )
    assert np.all(degree == 6)


def test_ellis_embedding_induces_the_target_warped_product_metric():
    r0 = 0.8292
    mesh = build_ellis_rgs_refinement(0, r0=r0, radial_base=8)
    u = mesh.radial_coordinates
    Rprime = u / np.sqrt(u * u + r0 * r0)
    zprime = r0 / np.sqrt(u * u + r0 * r0)
    assert np.max(np.abs(Rprime**2 + zprime**2 - 1.0)) < 2e-15
    assert mesh.vertices_r4.shape[1] == 4
    assert mesh.angular_directions.shape == (56, 3)
    assert mesh.angular_faces.shape == (108, 3)
    assert mesh.tetrahedra.shape == (2592, 4)
    assert np.min(mesh.element_gram[:, 0, 0]) > 0.0
    assert np.min(mesh.element_condition) > 1.0
    assert np.max(mesh.element_condition) < 100.0


def test_curved_mesh_chain_complex_and_p1_operators():
    mesh = build_ellis_rgs_refinement(
        0, r0=0.8292, radial_base=8, include_chain=True
    )
    assert mesh.chain is not None
    assert mesh.chain["chain_complex_exact"]
    assert mesh.chain["B1B2_nnz"] == 0
    assert mesh.chain["B2B3_nnz"] == 0

    mass, stiffness = assemble_ellis_p1_operators(mesh)
    n = mesh.vertices_r4.shape[0]
    assert mass.shape == (n, n)
    assert stiffness.shape == (n, n)
    assert np.linalg.norm((mass - mass.T).data) < 1e-12
    assert np.linalg.norm((stiffness - stiffness.T).data) < 1e-10
    assert np.min(mass.diagonal()) > 0.0
    one = np.ones(n)
    assert np.linalg.norm(stiffness @ one) < 1e-9
    assert np.isclose(one @ mass @ one, np.sum(mesh.element_volumes), rtol=1e-12)


def test_curved_rgs_refinement_converges_in_volume_and_diameter():
    report = audit_ellis_rgs_refinement(
        levels=(0, 1, 2), r0=0.8292, radial_base=8
    )
    rows = report["refinements"]
    errors = [row["volume_relative_error"] for row in rows]
    diameters = [row["h_max_chord"] for row in rows]
    conditions = [row["element_condition_max"] for row in rows]
    assert len(rows) == 3
    assert errors[0] > errors[1] > errors[2]
    assert errors[-1] < 0.01
    assert diameters[0] > diameters[1] > diameters[2]
    assert max(conditions) < 100.0
    assert report["seed"]["q6_coframe_reconstruction_residual"] == 0.0
    assert report["seed"]["q6_coframe_unit_norm_residual"] < 1e-14
    assert report["scope"]["Mosco_status"].startswith("OPEN")


def test_curved_refinement_rejects_invalid_parameters_and_levels():
    with pytest.raises(ValueError):
        build_ellis_rgs_refinement(-1)
    with pytest.raises(ValueError):
        build_ellis_rgs_refinement(0, r0=0.0)
    with pytest.raises(ValueError):
        audit_ellis_rgs_refinement(levels=(0, 2))
