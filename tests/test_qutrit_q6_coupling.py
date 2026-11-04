import numpy as np
import pytest

from src.qutrit_q6_coupling import (
    build_q6_qutrit_chain,
    q6_gray_path,
    q6_qutrit_gauge_covariance_residual,
    q6_qutrit_spectrum,
)


def _su3_matrix(a, b, phase):
    r12 = np.array([
        [np.cos(a), np.sin(a), 0.0],
        [-np.sin(a), np.cos(a), 0.0],
        [0.0, 0.0, 1.0],
    ], dtype=complex)
    r23 = np.array([
        [1.0, 0.0, 0.0],
        [0.0, np.cos(b), np.sin(b)],
        [0.0, -np.sin(b), np.cos(b)],
    ], dtype=complex)
    diagonal = np.diag([
        np.exp(1j*phase), np.exp(-0.5j*phase), np.exp(-0.5j*phase)
    ])
    return r12 @ diagonal @ r23


def test_gray_path_is_embedded_in_q6_and_fibonacci_modulated():
    chain = build_q6_qutrit_chain(5, j_a=1.0, j_b=0.4)
    assert chain.site_vertices == (0, 1, 3, 2, 6)
    assert chain.word == "ABAA"
    assert chain.couplings == (1.0, 0.4, 1.0, 1.0)
    assert len(chain.links) == 4
    assert chain.hamiltonian.shape == (243, 243)
    assert np.allclose(chain.hamiltonian, chain.hamiltonian.conj().T)


def test_qutrit_q6_local_su3_gauge_covariance():
    chain = build_q6_qutrit_chain(4, j_a=1.0, j_b=0.37)
    gauge = {
        v: _su3_matrix(0.13*(i+1), -0.09*(i+1), 0.17*(i+1))
        for i, v in enumerate(chain.site_vertices)
    }
    links = {
        edge: _su3_matrix(0.08*(i+1), 0.05*(i+1), -0.12*(i+1))
        for i, edge in enumerate(chain.links)
    }
    chain = build_q6_qutrit_chain(4, j_a=1.0, j_b=0.37, links=links)
    residual = q6_qutrit_gauge_covariance_residual(chain, gauge)
    assert residual < 2e-10


def test_q6_embedded_qutrit_spectral_entanglement_audit():
    chain = build_q6_qutrit_chain(4, j_a=1.0, j_b=0.37)
    report = q6_qutrit_spectrum(chain)
    assert report["site_vertices"] == list(q6_gray_path(4))
    assert report["fibonacci_word"] == "ABA"
    assert report["ground_degeneracy"] >= 1
    assert report["gap_above_ground_manifold"] >= -1e-10
    assert report["schmidt_rank"] >= 1
    assert report["hermiticity_residual"] < 1e-12
    assert np.isfinite(report["entanglement_entropy_bits"])


def test_invalid_q6_embedding_and_link_support_rejected():
    with pytest.raises(ValueError):
        build_q6_qutrit_chain(8)
    with pytest.raises(ValueError):
        build_q6_qutrit_chain(3, links={(0, 63): np.eye(3)})
    with pytest.raises(ValueError):
        build_q6_qutrit_chain(3, links={(0, 1): 2*np.eye(3)})
