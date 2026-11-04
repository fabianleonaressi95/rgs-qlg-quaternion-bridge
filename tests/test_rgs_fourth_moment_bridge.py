import numpy as np

from src.rgs_ellis_refinement import load_rgs_seed_data
from src.rgs_fourth_moment_bridge import (
    anisotropic_response_operator,
    audit_rgs_fourth_moment,
    fourth_moment_decomposition,
    isotropic_rank4_tensor,
)


def test_recovered_seed_has_second_moment_80_identity():
    x = load_rgs_seed_data().coordinates
    m2 = np.einsum("ni,nj->ij", x, x)
    assert np.allclose(m2, 80.0 * np.eye(3), atol=1e-10, rtol=0.0)


def test_fourth_moment_projection_is_trace_free_and_orthogonal():
    x = load_rgs_seed_data().coordinates
    data = fourth_moment_decomposition(x)
    r = data["M4_anisotropic"]
    assert np.linalg.norm(np.einsum("iikl->kl", r)) < 1e-9
    assert abs(data["isotropic_orthogonality"]) < 1e-8
    assert np.allclose(
        data["M4"], data["M4_isotropic"] + r, atol=1e-12, rtol=0.0
    )


def test_fourth_moment_has_nonzero_cubic_harmonic():
    x = load_rgs_seed_data().coordinates
    data = fourth_moment_decomposition(x)
    assert data["anisotropic_frobenius_norm"] > 1.0
    assert np.linalg.norm(data["M4_anisotropic"]) > 1.0


def test_anisotropic_response_splits_Eg_and_T2g_without_mixing():
    x = load_rgs_seed_data().coordinates
    data = fourth_moment_decomposition(x)
    response = anisotropic_response_operator(data["M4_anisotropic"])
    assert response["Eg_T2g_mixing_norm"] < 1e-9
    assert np.allclose(response["Eg_eigenvalues"], response["Eg_eigenvalues"][0], atol=1e-9)
    assert np.allclose(response["T2g_eigenvalues"], response["T2g_eigenvalues"][0], atol=1e-9)
    assert not np.isclose(response["Eg_eigenvalues"][0], response["T2g_eigenvalues"][0])


def test_full_fourth_moment_bridge_audit_passes():
    report = audit_rgs_fourth_moment()
    assert report["pass"]
    assert report["O_h"]["group_order"] == 48
    assert report["status"]["torsion_derivation"] == "OPEN"
