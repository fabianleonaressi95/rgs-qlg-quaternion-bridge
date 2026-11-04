import math
import numpy as np
import sympy as sp
from scipy.special import i0
from src.rgs_q6_d5_wilson_fiber_torus import smith_torus_parameterization


def _integrate_toy_annihilator(data, beta, character):
    v = np.asarray(data["right_transform"].tolist(), dtype=float)
    u = np.arange(4096, dtype=float)/4096.0
    values = []
    for j in range(data["invariants"][0]):
        z = np.column_stack((
            np.full(u.size, 2.0*math.pi*j/data["invariants"][0]),
            2.0*math.pi*u,
        ))
        x = z @ v.T
        amplitude = np.exp(beta*np.sum(np.cos(x), axis=1))/(i0(beta)**2)
        phase = x @ np.asarray(character, dtype=float)
        values.append(np.mean(amplitude*np.exp(-1j*phase)))
    return np.mean(values)


def test_smith_parameterization_tracks_disconnected_annihilator():
    # L=< (2,0) > gives H={0,pi} x T^1, two components and dimension one.
    data = smith_torus_parameterization(sp.Matrix([[2, 0]]))
    assert data["torus_dimension"] == 1
    assert data["invariants"] == [2]
    assert data["left_transform"]*sp.Matrix([[2,0]])*data["right_transform"] == data["smith_form"]


def test_annihilator_integral_reproduces_exact_bessel_fiber_sums():
    beta = 1.25
    data = smith_torus_parameterization(sp.Matrix([[2, 0]]))
    zero_character = _integrate_toy_annihilator(data, beta, [0, 0])
    odd_character = _integrate_toy_annihilator(data, beta, [1, 0])
    assert abs(zero_character.imag) < 1e-12
    assert abs(odd_character.imag) < 1e-12
    assert math.isclose(zero_character.real, math.cosh(beta)/i0(beta), rel_tol=1e-10)
    assert math.isclose(odd_character.real, math.sinh(beta)/i0(beta), rel_tol=1e-10)


def test_smith_parameterization_for_primitive_lattice():
    data = smith_torus_parameterization(sp.Matrix([[1, 0, 0], [0, 1, 0]]))
    assert data["torus_dimension"] == 1
    assert data["invariants"] == [1, 1]
    assert data["smith_form"].shape == (2, 3)


def test_invalid_lattice_rank_rejected():
    import pytest
    with pytest.raises(ValueError):
        smith_torus_parameterization(sp.Matrix([[1,0],[2,0]]))
