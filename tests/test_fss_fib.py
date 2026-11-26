import numpy as np
import pytest
from src.fss_fib import (
    aubry_andre_hamiltonian, chi0_derivative, chi0_retarded,
    chemical_potential_for_filling, density_vertex, diagonalize_chain,
    fibonacci_approximant, fit_fss_power_law, nn_density_kernel,
    occupations_fermi_dirac, occupations_zero_temperature, rpa_pole_residue,
)

def test_fibonacci_approximants():
    assert fibonacci_approximant(10) == (55, 34, 34/55)
    assert fibonacci_approximant(12) == (144, 89, 89/144)
    with pytest.raises(ValueError):
        fibonacci_approximant(1)

def test_aa_hamiltonian_is_symmetric_and_critical_by_default():
    phi = (1 + np.sqrt(5))/2
    h = aubry_andre_hamiltonian(13, alpha=1/phi)
    assert np.allclose(h, h.T)
    assert np.allclose(np.diag(h), 2*np.cos(2*np.pi*np.arange(13)/phi))

def test_density_vertex_at_q_zero_is_identity():
    assert np.allclose(density_vertex(0.0, 7), np.eye(7))

def test_q_zero_response_vanishes_by_number_conservation():
    h = aubry_andre_hamiltonian(8, alpha=0.6180339887498948)
    e, u = diagonalize_chain(h)
    f = occupations_zero_temperature(e, 4)
    assert abs(chi0_retarded(e, u, 0.0, 0.3+0.1j, f)) < 1e-13

def test_response_derivative_matches_finite_difference():
    h = aubry_andre_hamiltonian(10, alpha=0.6180339887498948)
    e, u = diagonalize_chain(h)
    f = occupations_zero_temperature(e, 5)
    z, q, step = 0.7+0.2j, 0.7, 1e-6
    d = chi0_derivative(e, u, q, z, f)
    fd = (chi0_retarded(e,u,q,z+step,f)-chi0_retarded(e,u,q,z-step,f))/(2*step)
    assert d == pytest.approx(fd, rel=1e-6, abs=1e-8)

def test_finite_temperature_filling():
    e = np.array([-1.0, -0.2, 0.3, 1.2])
    mu = chemical_potential_for_filling(e, temperature=0.2, particle_number=2.0)
    f = occupations_fermi_dirac(e, temperature=0.2, chemical_potential=mu)
    assert np.sum(f) == pytest.approx(2.0, abs=1e-10)

def test_nn_kernel_and_pole_residue():
    assert nn_density_kernel(0.0, 1.5) == pytest.approx(3.0)
    assert rpa_pole_residue(0.5, -0.25, 2.0) == pytest.approx(1.0)

def test_fixed_exponent_fss_fit():
    sizes = [55, 144, 377, 987]
    values = [0.5 + 2.0/n**2 for n in sizes]
    fit = fit_fss_power_law(sizes, values, exponent=2.0)
    assert fit.limit == pytest.approx(0.5, abs=1e-12)
    assert fit.amplitude == pytest.approx(2.0, abs=1e-10)
    assert fit.residual_norm < 1e-12
