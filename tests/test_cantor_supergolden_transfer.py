import numpy as np
import pytest
from src.cantor_supergolden_transfer import (
    PSI, audit_cantor_q6_supergolden, cylinder_probability,
    generalized_dimension, golden_branch_probabilities, perron_data,
    pressure, pressure_matrix, q6_laplacian, q6_spectrum_residual,
    supergolden_matrix,
)

def test_supergolden_identities():
    A=supergolden_matrix()
    assert abs(PSI**3-PSI**2-1.)<2e-14
    assert np.linalg.norm(np.linalg.matrix_power(A,3)-np.linalg.matrix_power(A,2)-np.eye(3))<1e-14

def test_parry_markov_data():
    d=perron_data()
    assert d["perron_residual"]<1e-12
    assert d["row_sum_residual"]<1e-12
    assert d["stationary_residual"]<1e-12
    assert np.all(d["P"]>=0.)

def test_golden_probabilities():
    p=golden_branch_probabilities()
    assert np.allclose(p["p0"]+p["p1"],1.,atol=1e-14)
    assert len(set(np.round(p["p0"],10)))==3

def test_pressure_special_values():
    assert np.isclose(pressure(1.),0.,atol=1e-12)
    assert np.isclose(max(abs(np.linalg.eigvals(pressure_matrix(0.)))),2.*PSI,atol=1e-12)

def test_generalized_dimensions():
    vals=[generalized_dimension(s) for s in (0.,.5,1.,2.,3.)]
    assert np.all(np.isfinite(vals))
    assert max(vals)-min(vals)>.05
    assert np.isclose(vals[0],np.log(2.*PSI)/np.log(3.*PSI))

def test_cylinder_probability():
    assert 0.<cylinder_probability([0,0,0],[0,1])<1.
    assert cylinder_probability([0,1],[0])==0.

def test_q6_spectrum():
    assert q6_laplacian().shape==(64,64)
    assert q6_spectrum_residual()<1e-12

def test_audit_keeps_physical_claim_open():
    r=audit_cantor_q6_supergolden()
    assert r["status"].startswith("symbolic model")
    assert "construct an intertwiner with the Q6 Laplacian" in r["open_obligations"]

@pytest.mark.parametrize("phi",[1.,0.,float("nan")])
def test_invalid_phi(phi):
    with pytest.raises(ValueError): golden_branch_probabilities(phi)
