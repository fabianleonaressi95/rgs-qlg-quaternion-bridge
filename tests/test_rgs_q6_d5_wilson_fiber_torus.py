import math
import sympy as sp
from src.rgs_q6_d5_wilson_fiber_torus import smith_torus_parameterization


def test_smith_parameterization_tracks_disconnected_annihilator():
    # L=< (2,0) > gives H={0,pi} x T^1, two components and dimension one.
    data = smith_torus_parameterization(sp.Matrix([[2, 0]]))
    assert data["torus_dimension"] == 1
    assert data["invariants"] == [2]
    assert abs(int(data["right_transform"].det())) == 1
    assert data["left_transform"]*sp.Matrix([[2,0]])*data["right_transform"] == data["smith_form"]


def test_smith_parameterization_for_primitive_lattice():
    data = smith_torus_parameterization(sp.Matrix([[1, 0, 0], [0, 1, 0]]))
    assert data["torus_dimension"] == 1
    assert data["invariants"] == [1, 1]
    assert data["smith_form"].shape == (2, 3)


def test_invalid_lattice_rank_rejected():
    import pytest
    with pytest.raises(ValueError):
        smith_torus_parameterization(sp.Matrix([[1,0],[2,0]]))
