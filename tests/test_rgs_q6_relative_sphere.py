import numpy as np
import pytest

from src.rgs_q6_relative_sphere import (
    build_q6_a5_relative_sphere_map,
    q6_vertex_relative_direction,
)


def test_q6_vertices_map_equivariantly_to_finite_sphere_design():
    report = build_q6_a5_relative_sphere_map()
    assert report["linear_map_shape"] == [3, 6]
    assert report["linear_map_rank"] == 3
    assert report["linear_map_kernel_dimension"] == 3
    assert report["group_order"] == 60
    assert report["unique_direction_count"] == 32
    assert np.all(report["fibre_multiplicities"] == 2)
    assert report["antipodal_pair_count"] == 16
    assert report["all_checks_passed"]


def test_relative_direction_is_unit_and_vertex_indexed():
    directions = np.asarray([q6_vertex_relative_direction(v) for v in range(64)])
    assert directions.shape == (64, 3)
    assert np.allclose(np.linalg.norm(directions, axis=1), 1.0, atol=1e-13)
    assert np.allclose(directions[63], -directions[0], atol=1e-13)
    with pytest.raises(ValueError):
        q6_vertex_relative_direction(64)


def test_32_point_direction_set_is_spherical_five_design_but_not_six():
    report = build_q6_a5_relative_sphere_map()
    residuals = report["moment_residuals_to_uniform_S2"]
    assert report["design_strength_at_least"] == 5
    for degree in range(1, 6):
        assert residuals[str(degree)] < 1e-10
    assert residuals["6"] > 1e-5


def test_icosahedral_signed_axis_frame_is_tight_and_four_radial_shells():
    report = build_q6_a5_relative_sphere_map()
    assert report["axis_frame_residual"] < 1e-10
    assert len(report["radial_shells"]) == 4
    assert sum(report["radial_shells"].values()) == 64
