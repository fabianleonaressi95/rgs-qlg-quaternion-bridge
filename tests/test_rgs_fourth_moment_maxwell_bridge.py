import numpy as np
import pytest

from src.rgs_fourth_moment_maxwell_bridge import (
    direct_antisymmetric_contraction,
    fourth_moment_anisotropy,
    maxwell_constitutive_response,
    maxwell_density,
)


def _symmetric_rank4():
    rng = np.random.default_rng(20261104)
    points = rng.normal(size=(40, 3))
    return np.einsum("ni,nj,nk,nl->ijkl", points, points, points, points)


def test_direct_contraction_with_antisymmetric_field_is_zero():
    r = _symmetric_rank4()
    f = np.array([[0.0, 0.3, -0.7], [-0.3, 0.0, 0.2], [0.7, -0.2, 0.0]])
    assert abs(direct_antisymmetric_contraction(r, f)) < 1e-9


def test_induced_response_is_stf_and_positive():
    r = _symmetric_rank4()
    q = np.diag([1.0, -1.0, 0.0]) / np.sqrt(2.0)
    raw = np.einsum("ijkl,kl->ij", r, q)
    raw = 0.5 * (raw + raw.T)
    raw -= np.eye(3) * np.trace(raw) / 3.0
    coupling = 0.2 / max(np.linalg.norm(raw, 2), 1.0)
    out = maxwell_constitutive_response(
        r, q, electric_coupling=coupling, magnetic_coupling=-coupling
    )
    assert abs(np.trace(out["anisotropy_operator"])) < 1e-10
    assert out["electric_eigenvalues"][0] > 0
    assert out["magnetic_eigenvalues"][0] > 0
    assert np.isfinite(maxwell_density(
        np.array([1.0, 0.2, -0.1]), np.array([0.3, -0.4, 0.5]),
        out["electric_matrix"], out["magnetic_matrix"]
    ))


def test_nonpositive_constitutive_matrix_is_rejected():
    r = _symmetric_rank4()
    q = np.diag([1.0, -1.0, 0.0]) / np.sqrt(2.0)
    with pytest.raises(ValueError, match="not positive definite"):
        maxwell_constitutive_response(
            r, q, electric_coupling=1e6, magnetic_coupling=0.0
        )


def test_tensor_covariance_under_orthogonal_frame_change():
    r = _symmetric_rank4()
    q = np.diag([1.0, -1.0, 0.0]) / np.sqrt(2.0)
    theta = 0.37
    g = np.array([[np.cos(theta), -np.sin(theta), 0.0],
                  [np.sin(theta), np.cos(theta), 0.0],
                  [0.0, 0.0, 1.0]])
    rg = np.einsum("ia,jb,kc,ld,abcd->ijkl", g, g, g, g, r)
    qg = g @ q @ g.T
    a = np.einsum("ijkl,kl->ij", fourth_moment_anisotropy(r), q)
    ag = np.einsum("ijkl,kl->ij", fourth_moment_anisotropy(rg), qg)
    assert np.linalg.norm(ag - g @ a @ g.T) < 1e-8
