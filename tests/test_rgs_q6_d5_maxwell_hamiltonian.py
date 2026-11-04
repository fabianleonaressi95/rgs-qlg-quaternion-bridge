import math
import numpy as np
import pytest

from src.rgs_q6_d5_maxwell_hamiltonian import (
    build_q6_d5_maxwell_hamiltonian,
    maxwell_mode_parameters,
    phase_space_to_d5_maxwell_tangent,
)


def test_maxwell_legendre_coefficients_and_beta_independent_frequency():
    for beta in (0.125, 0.7, 1.0, 3.5, 12.0):
        p = maxwell_mode_parameters(beta, 12.0)
        assert p["electric_hessian"] == pytest.approx(1.0 / beta)
        assert p["magnetic_hessian"] == pytest.approx(12.0 * beta)
        assert p["oscillator_frequency"] == pytest.approx(math.sqrt(12.0))


def test_maxwell_normalized_d5_hamiltonian_audit_passes():
    for beta in (0.2, 1.0, 4.0):
        report = build_q6_d5_maxwell_hamiltonian(beta=beta)
        assert report["all_checks_passed"], report["checks"]
        assert report["liouville_jacobian"] == pytest.approx(1.0 / (32.0 * 10.0**5), rel=1e-12)


def test_d5_chart_energy_is_hamiltonian_for_multiple_beta():
    q = np.array([0.1, -0.2, 0.05, 0.03, -0.04])
    p = np.array([-0.03, 0.08, 0.11, -0.07, 0.02])
    for beta in (0.25, 1.0, 2.5):
        report = build_q6_d5_maxwell_hamiltonian(beta=beta)
        par = report["parameters"]
        h = 0.5 * (par["magnetic_hessian"] * np.dot(q, q) + par["electric_hessian"] * np.dot(p, p))
        z = phase_space_to_d5_maxwell_tangent(q, p, beta=beta)
        assert 10.0 * par["oscillator_frequency"] * np.vdot(z, z).real == pytest.approx(h)


@pytest.mark.parametrize("beta,eigenvalue", [(0.0, 12.0), (-1.0, 12.0), (1.0, 0.0)])
def test_invalid_maxwell_parameters_rejected(beta, eigenvalue):
    with pytest.raises(ValueError):
        maxwell_mode_parameters(beta, eigenvalue)
