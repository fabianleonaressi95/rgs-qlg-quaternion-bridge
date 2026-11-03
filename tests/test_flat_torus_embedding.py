import numpy as np
import pytest

from src.flat_torus_refinement import build_flat_torus_mesh
from src.flat_torus_embedding import (
    J_h,
    assemble_consistent_mass0,
    audit_g4_embedding,
    discrete_inner_product,
)


@pytest.mark.parametrize("dimension,n", [(2, 4), (2, 6), (3, 3)])
def test_g4_embedding_has_exact_positive_l2_mass(dimension, n):
    mesh = build_flat_torus_mesh(n, dimension=dimension, side_length=1.7)
    audit = audit_g4_embedding(mesh)
    mass = assemble_consistent_mass0(mesh)
    assert audit["mass_positive_definite"]
    assert audit["mass_symmetry_residual"] < 1e-12
    assert audit["volume_relative_error"] < 1e-12
    assert np.linalg.eigvalsh(mass)[0] > 0.0
    ones = np.ones(mesh.vertices.shape[0])
    assert np.isclose(ones @ mass @ ones, 1.7**dimension, rtol=1e-12, atol=1e-12)


def test_jh_reconstructs_nodal_values_and_is_periodic():
    mesh = build_flat_torus_mesh(5, dimension=2, side_length=2.0)
    x = mesh.vertices
    values = np.sin(2 * np.pi * x[:, 0] / 2.0) + 0.3 * np.cos(2 * np.pi * x[:, 1] / 2.0)
    reconstruction = J_h(mesh, values)

    assert np.allclose(reconstruction.evaluate(mesh.vertices), values, atol=1e-10)
    probes = np.array([[0.13, 0.27], [1.91, 0.04], [2.13, -1.73]])
    assert np.allclose(
        reconstruction.evaluate(probes),
        reconstruction.evaluate(probes + np.array([2.0, -4.0])),
        atol=1e-10,
    )


@pytest.mark.parametrize("dimension,n", [(2, 4), (3, 3)])
def test_jh_is_piecewise_affine_on_each_simplex(dimension, n):
    mesh = build_flat_torus_mesh(n, dimension=dimension)
    values = np.arange(mesh.vertices.shape[0], dtype=float) ** 2
    reconstruction = J_h(mesh, values)
    L = mesh.side_length

    for simplex in mesh.simplices:
        ids = np.asarray(simplex, dtype=int)
        origin = mesh.vertices[ids[0]]
        offsets = np.mod(mesh.vertices[ids] - origin, L)
        coords = origin + offsets
        barycenter = np.mod(np.mean(coords, axis=0), L)
        expected = float(np.mean(values[ids]))
        assert np.isclose(reconstruction.evaluate(barycenter), expected, atol=1e-9)


def test_discrete_inner_product_is_exact_p1_l2_product():
    mesh = build_flat_torus_mesh(4, dimension=2, side_length=1.3)
    u = np.sin(mesh.vertices[:, 0]) + mesh.vertices[:, 1]
    v = np.cos(mesh.vertices[:, 1]) - mesh.vertices[:, 0]
    mass = assemble_consistent_mass0(mesh)
    expected = float(u @ mass @ v)
    assert np.isclose(discrete_inner_product(mesh, u, v), expected, atol=1e-12)
    assert discrete_inner_product(mesh, u, u) > 0.0


def test_embedding_rejects_wrong_nodal_dimension_and_points():
    mesh = build_flat_torus_mesh(4, dimension=2)
    with pytest.raises(ValueError):
        J_h(mesh, np.zeros(3))
    reconstruction = J_h(mesh, np.zeros(mesh.vertices.shape[0]))
    with pytest.raises(ValueError):
        reconstruction.evaluate(np.zeros((2, 3)))
