import numpy as np

from scripts.audit_cartan_refinement_consistency import (
    A_CURV,
    B_CURV,
    TORSION_VECTOR,
    convergence_table,
    edge_transport,
    face_observables,
)


def test_manufactured_cartan_pair_has_nonzero_curvature_and_torsion():
    assert np.linalg.norm(A_CURV - B_CURV) > 0.0
    assert np.linalg.norm(TORSION_VECTOR) > 0.0


def test_edge_transport_is_affine_and_has_homogeneous_last_row():
    G = edge_transport((0.0, 0.0), (0.1, 0.0))
    assert np.allclose(G[4], np.array([0.0, 0.0, 0.0, 0.0, 1.0]), atol=1e-13)


def test_face_observables_are_small_and_finite():
    log_H, tau = face_observables(0.05)
    assert np.isfinite(log_H).all()
    assert np.isfinite(tau).all()
    assert np.linalg.norm(log_H) > 0.0
    assert np.linalg.norm(tau) > 0.0


def test_face_consistency_has_expected_asymptotic_orders():
    rows = convergence_table()
    curvature_orders = [r["curvature_observed_order"] for r in rows[1:]]
    torsion_orders = [r["torsion_observed_order"] for r in rows[1:]]

    # This manufactured field has a cancellation that raises the curvature
    # face-error order to four; the generic local bound required is >= three.
    assert min(curvature_orders[-2:]) > 3.8
    assert min(torsion_orders[-2:]) > 2.8

    # Dividing by the triangular area gives density errors of orders two and one.
    assert rows[-1]["curvature_density_error"] < rows[-2]["curvature_density_error"] / 3.5
    assert rows[-1]["torsion_density_error"] < rows[-2]["torsion_density_error"] / 1.8
