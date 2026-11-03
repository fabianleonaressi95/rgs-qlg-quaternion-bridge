import numpy as np
import pytest
from scipy.sparse.linalg import eigsh
from src.naressi_golden_bridge import (
    PHI, audit_golden_naressi_bridge, audit_qat_unitarity,
    golden_decay_to_q6_tau, golden_parameters, qat_hamiltonian,
    q6_entropy, tfim_golden_hamiltonian, tfim_ground_state_audit,
)

def test_capsule_golden_parameters():
    p=golden_parameters()
    assert np.isclose(p["q_phi"],0.42181947690765276,atol=1e-14)
    assert np.isclose(p["J_n"],p["q_phi"],atol=1e-14)
    assert np.isclose(p["chi_n"],-0.3110,atol=5e-5)
    assert np.isclose(p["omega_photon"],PHI**2)
    assert np.isclose(p["omega_phonon"],PHI)

def test_geometric_coupling_ratio_is_constant():
    j=np.array([golden_parameters(n=n)["J_n"] for n in range(1,7)])
    assert np.allclose(j[1:]/j[:-1],j[1]/j[0],atol=1e-14)

def test_qat_corrected_hamiltonian_is_hermitian_and_unitary():
    h=qat_hamiltonian(3)
    assert np.linalg.norm(h-h.conj().T)<1e-13
    report=audit_qat_unitarity()
    assert report["pass"]
    assert report["max_unitarity_residual"]<1e-10
    assert report["max_state_norm_residual"]<1e-10

def test_golden_decay_dictionary_matches_q6_spectral_ratio():
    tau=golden_decay_to_q6_tau()
    assert np.isclose(np.exp(-tau),golden_parameters()["q_phi"],atol=1e-14)
    assert np.isclose(q6_entropy(tau),2.5472014004557337,atol=1e-12)

def test_sparse_tfim_matches_dense_small_chain():
    h=tfim_golden_hamiltonian(4,g=0.7)
    assert np.linalg.norm((h-h.T).toarray())<1e-13
    dense=np.linalg.eigvalsh(h.toarray())[:2]
    sparse=np.sort(eigsh(h,k=2,which="SA",return_eigenvectors=False))
    assert np.allclose(sparse,dense,atol=1e-10)

def test_tfim_gap_and_entropy_are_finite():
    result=tfim_ground_state_audit(6,g=0.8)
    assert np.isfinite(result["ground_energy"])
    assert result["gap"]>0
    assert result["half_chain_entropy_nats"]>=0

def test_bridge_does_not_claim_geometric_derivation():
    report=audit_golden_naressi_bridge()
    assert report["q6_spectral_dictionary"]["status"].startswith("explicit spectral ansatz")
    assert report["evidence"]["tau_to_geometry"].startswith("OPEN")

@pytest.mark.parametrize("phi",[1.0,0.0,float("nan")])
def test_invalid_phi_is_rejected(phi):
    with pytest.raises(ValueError):
        golden_parameters(phi=phi)
