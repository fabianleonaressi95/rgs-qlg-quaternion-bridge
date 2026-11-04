import numpy as np

from src.q6_complex import build_q6_complex
from src.spinor_kirchhoff import (
    audit_spinor_kirchhoff, covariant_laplacian, dagger,
    edge_currents, gauge_transform_links, gauge_transform_spinors,
    random_q6_spinor_data,
)


def test_seeded_spinor_kirchhoff_audit_passes():
    report = audit_spinor_kirchhoff(seed=20261104, gauge_seed=20261105, atol=1e-10)
    assert report["passed"], report
    assert report["dimensions"]["vertices"] == 64
    assert report["dimensions"]["edges"] == 192
    assert report["dimensions"]["faces"] == 240


def test_edge_currents_are_local_gauge_invariant():
    q6, edges, links, psi, weights, gauge = random_q6_spinor_data(7123)
    links_g = gauge_transform_links(edges, links, gauge)
    psi_g = gauge_transform_spinors(psi, gauge)
    J = edge_currents(edges, links, psi, weights)
    Jg = edge_currents(edges, links_g, psi_g, weights)
    assert np.max(np.abs(J-Jg)) < 1e-11


def test_covariant_laplacian_is_hermitian_and_gauge_intertwines():
    q6, edges, links, psi, weights, gauge = random_q6_spinor_data(8123)
    links_g = gauge_transform_links(edges, links, gauge)
    H = covariant_laplacian(edges, links, weights)
    Hg = covariant_laplacian(edges, links_g, weights)
    G = np.zeros((128, 128), dtype=complex)
    for v in range(64):
        G[2*v:2*v+2, 2*v:2*v+2] = gauge[v]
    assert np.linalg.norm(H-H.conj().T) < 1e-11
    assert np.linalg.norm(Hg-G@H@G.conj().T) < 1e-10
