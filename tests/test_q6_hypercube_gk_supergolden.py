import numpy as np
import pytest
from src.q6_hypercube_gk_supergolden import (
    green_kubo_ou, green_kubo_single_exponential, green_kubo_unitary_equivalence_residuals,
    q6_closed_two_cochain_projector, q6_hypercube_cosine_law,
    q6_projected_quadratic_flux_action, shear_stress_from_energy,
    supergolden_companion_matrix, supergolden_constant, supergolden_geometric_sum,
    supergolden_refinement_scales, supergolden_spectral_ratio_audit,
)
def test_hypercube_cosine_identity_and_taylor_expansion():
    r=q6_hypercube_cosine_law(np.linspace(-.02,.02,240),w=2.)
    assert r["identity_residual"]==pytest.approx(0.,abs=1e-13)
    assert r["F_Q6"]==pytest.approx(r["F_Q6_from_delta_mu4"])
    assert abs(r["taylor_remainder"])<1e-8 and r["plaquette_count"]==240
def test_closed_cochain_projector_is_orthogonal_rank_129():
    p=q6_closed_two_cochain_projector()
    assert np.allclose(p,p.T,atol=1e-10) and np.allclose(p@p,p,atol=1e-10)
    assert np.linalg.matrix_rank(p,tol=1e-8)==129
    assert q6_projected_quadratic_flux_action(np.ones(240))["rank_projector"]==129
def test_green_kubo_modal_sum():
    r=green_kubo_ou(np.diag([2.,4.]),np.eye(2),np.array([1.,2.]))
    assert r.convergent and r.correlation_integral==pytest.approx(1.5) and r.eta==pytest.approx(1.5)
    assert green_kubo_single_exponential(3.,2.)==pytest.approx(1.5)
def test_green_kubo_detects_nondecaying_zero_mode():
    r=green_kubo_ou(np.diag([0.,2.]),np.eye(2),np.array([1.,0.]))
    assert not r.convergent and np.isinf(r.eta)
def test_green_kubo_zero_mode_harmless_without_stress_projection():
    r=green_kubo_ou(np.diag([0.,2.]),np.eye(2),np.array([0.,1.]))
    assert r.convergent and r.eta==pytest.approx(.5)
def test_supergolden_polynomial_and_geometric_series():
    psi=supergolden_constant()
    assert psi**3-psi**2-1.==pytest.approx(0.,abs=2e-14)
    assert np.linalg.det(psi*np.eye(3)-supergolden_companion_matrix())==pytest.approx(0.,abs=1e-12)
    assert psi==pytest.approx(1.465571231876768,abs=1e-14)
    assert supergolden_geometric_sum(1.)==pytest.approx(psi/(psi-1.))
    assert supergolden_geometric_sum(1.,10000)==pytest.approx(psi/(psi-1.))
    s=supergolden_refinement_scales(1.,4)
    assert np.allclose(s[1:]/s[:-1],1./psi)
def test_spectral_ratio_audit_is_diagnostic():
    r=supergolden_spectral_ratio_audit(np.array([1.,2.,4.]))
    assert r["positive_eigenvalue_count"]==3 and "not a theorem" in r["interpretation"]
def test_shear_stress_is_negative_energy_derivative():
    assert shear_stress_from_energy(lambda g:3.+2.*g+5.*g*g)==pytest.approx(-2.,abs=1e-9)
def test_unitary_equivalence_transports_all_green_kubo_data():
    t=.37; u=np.array([[np.cos(t),-np.sin(t)],[np.sin(t),np.cos(t)]],complex)
    h=np.diag([1.,2.]).astype(complex); p=np.array([[0.,1.],[1.,0.]],complex)
    rho=np.diag([.7,.3]).astype(complex)
    r=green_kubo_unitary_equivalence_residuals(u,h,p,rho,u@h@u.conj().T,
        u@p@u.conj().T,u@rho@u.conj().T)
    assert max(r.values())<1e-12
