import numpy as np
import pytest

from src.kathara_square_similarity import (
    apply_stage, audit_kathara_square_similarity, kathara_centers_from_overlay,
    kathara_stage_transform,
)


def test_stage_map_is_square_diagonal_similarity():
    T = kathara_stage_transform()
    assert np.array_equal(T.T @ T, 2 * np.eye(2, dtype=int))
    assert int(round(np.linalg.det(T))) == 2
    assert np.array_equal(T @ T, np.array([[0, 2], [-2, 0]]))


def test_reported_center_overlays_determine_center9():
    a = audit_kathara_square_similarity(width=2.0)
    assert np.allclose(a.transform @ a.centers["9"], a.centers["12"])
    assert np.allclose(a.transform_squared @ a.centers["6"], a.centers["12"])
    assert np.allclose(a.centers["9"], a.transform @ a.centers["6"])
    assert a.residuals["central_column_ratio_residual"] < 1e-14
    assert a.residuals["stage_scale_ratio_6_to_9"] == pytest.approx(np.sqrt(2))
    assert a.residuals["stage_scale_ratio_9_to_12"] == pytest.approx(np.sqrt(2))


def test_two_stages_are_quarter_turn_and_double_scale():
    T = kathara_stage_transform()
    assert np.array_equal(np.linalg.matrix_power(T, 2), np.array([[0, 2], [-2, 0]]))
    assert np.allclose(apply_stage(np.array([-1.0, 0.0]), 2), [0.0, 2.0])


def test_lattice_image_is_checkerboard_index_two():
    T = kathara_stage_transform()
    for x in range(-4, 5):
        for y in range(-4, 5):
            u, v = T @ np.array([x, y])
            assert (u - v) % 2 == 0
    assert abs(round(np.linalg.det(T))) == 2


def test_invalid_width_rejected():
    with pytest.raises(ValueError):
        kathara_centers_from_overlay(0)
