"""Exact cubical-product spectrum behind the Q6 five-mode Maxwell sector."""
import numpy as np

from src.rgs_q6_equivariant_cochains import build_q6_cochain_complex


def test_q6_coexact_hodge_spectrum_matches_product_complex_theorem():
    complex_data = build_q6_cochain_complex()
    d0 = np.asarray(complex_data["D0"], dtype=float)
    d1 = np.asarray(complex_data["D1"], dtype=float)

    assert d0.shape == (192, 64)
    assert d1.shape == (240, 192)
    assert np.linalg.matrix_rank(d0, tol=1e-10) == 63
    assert np.linalg.matrix_rank(d1, tol=1e-10) == 129
    assert np.max(np.abs(d1 @ d0)) == 0.0

    # Orthonormal basis of ker(D0^T), the 129-dimensional Coulomb tangent.
    _u, singular_values, vh = np.linalg.svd(d0.T, full_matrices=True)
    rank_d0 = int(np.sum(singular_values > 1e-10))
    coulomb_basis = vh[rank_d0:].T
    assert coulomb_basis.shape == (192, 129)
    assert np.linalg.norm(d0.T @ coulomb_basis) < 1e-10

    hodge = d1.T @ d1
    reduced = coulomb_basis.T @ hodge @ coulomb_basis
    reduced = 0.5 * (reduced + reduced.T)
    spectrum = np.linalg.eigvalsh(reduced)

    # Tensor-product proof: for k=0,...,5, lambda=2+2k and
    # multiplicity 6*C(5,k)-C(6,k+1) after removing exact 1-forms.
    expected = {4.0: 15, 6.0: 40, 8.0: 45, 10.0: 24, 12.0: 5}
    observed = {
        value: int(np.count_nonzero(np.isclose(spectrum, value, atol=1e-8, rtol=0.0)))
        for value in expected
    }
    assert observed == expected
    assert sum(observed.values()) == 129

    top = np.flatnonzero(np.isclose(spectrum, 12.0, atol=1e-8, rtol=0.0))
    assert top.size == 5
