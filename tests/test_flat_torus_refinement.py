import numpy as np
import pytest

from src.flat_torus_refinement import (
    audit_flat_torus_mesh,
    build_flat_torus_mesh,
    flat_torus_refinement_audit,
)


@pytest.mark.parametrize("dimension,n", [(2, 4), (2, 7), (3, 3), (3, 4)])
def test_periodic_freudenthal_mesh_has_valid_chain_complex(dimension, n):
    mesh = build_flat_torus_mesh(n, dimension=dimension)
    audit = audit_flat_torus_mesh(mesh)
    assert audit["boundary_squared_exact_zero"]
    assert audit["vertex_count"] == n**dimension
    assert audit["simplex_count"] == np.math.factorial(dimension) * n**dimension
    assert audit["hodge_star_0_positive"]
    assert audit["stiffness_positive_semidefinite_with_tolerance"]
    assert audit["stiffness_symmetry_residual"] < 1e-11
    assert audit["stiffness_row_sum_residual"] < 1e-11
    assert np.isclose(audit["mass0_total"], 1.0)


def test_2d_torus_expected_periodic_counts():
    mesh = build_flat_torus_mesh(5, dimension=2)
    assert len(mesh.vertices) == 25
    assert len(mesh.edges) == 75
    assert len(mesh.triangles) == 50
    assert len(mesh.simplices) == 50


def test_refinement_audit_runs_three_levels_and_energy_error_decreases():
    report = flat_torus_refinement_audit((4, 8, 16), dimension=2, mode=(1, 0))
    rows = report["refinements"]
    assert len(rows) == 3
    assert rows[-1]["spacing"] < rows[0]["spacing"]
    assert rows[-1]["relative_energy_error"] < rows[0]["relative_energy_error"]
    assert report["hodge_star_scope"] == "positive lumped 0-form Hodge star only"


def test_bad_mesh_parameters_are_rejected():
    with pytest.raises(ValueError):
        build_flat_torus_mesh(2, dimension=2)
    with pytest.raises(ValueError):
        build_flat_torus_mesh(4, dimension=4)
    with pytest.raises(ValueError):
        flat_torus_refinement_audit((4, 8), dimension=2)
