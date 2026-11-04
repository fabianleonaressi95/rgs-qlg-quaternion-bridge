import numpy as np
import pytest

from src.rgs_ellis_spin_connection import ellis_coframe_jet, ellis_rgs_spin_connection


def test_ellis_analytic_coframe_reconstructs_metric_and_jet():
    u, theta, r0 = 0.37, 0.91, 0.8292
    e, de = ellis_coframe_jet(u, theta, r0)
    R = np.hypot(u, r0)
    expected = np.diag([-1.0, 1.0, R**2, R**2*np.sin(theta)**2])
    eta = np.diag([-1.0, 1.0, 1.0, 1.0])
    assert np.allclose(e.T @ eta @ e, expected, atol=1e-14)
    assert de[2, 2, 1] == pytest.approx(u/R)
    assert de[3, 3, 1] == pytest.approx(u*np.sin(theta)/R)
    assert de[3, 3, 2] == pytest.approx(R*np.cos(theta))


def test_ellis_spin_connection_is_torsion_free_and_metric_compatible():
    report = ellis_rgs_spin_connection(0.37, 0.91)
    assert report["passed"]
    assert report["antisymmetry_residual"] < 1e-12
    assert report["tetrad_postulate_residual"] < 1e-12
    assert report["torsion_free_residual"] < 1e-12
    omega = report["omega_lower"]
    assert np.linalg.norm(omega) > 0.1
    assert report["geometry"] == "ELLIS_RGS_ANALYTIC_CONTINUUM_CHART"


def test_ellis_chart_rejects_coordinate_axis_singularity():
    with pytest.raises(ValueError):
        ellis_coframe_jet(0.0, 0.0)
