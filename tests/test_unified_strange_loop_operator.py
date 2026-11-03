import numpy as np
import pytest

from src.unified_strange_loop_operator import (
    audit_unified_operator,
    barycentric_gradient_moire_audit,
    build_doppler_moire_operator,
    build_strange_loop_operator,
    cantor_parry_q6_laplacian,
    doppler_moire_phase_audit,
    escher_lifted_torus_laplacian,
    moving_barycenter_audit,
    torus_q6_embedding,
)


def test_operator_is_symmetric_psd_and_constant_mode_kernel():
    H, _ = build_strange_loop_operator()
    assert H.shape == (384, 384)
    assert np.linalg.norm(H-H.T) < 1e-12
    assert np.linalg.norm(H @ np.ones(384)) < 1e-11
    assert np.linalg.eigvalsh(H)[0] > -1e-10


def test_torus_gray_embedding_and_lift():
    ids = torus_q6_embedding()
    assert len(ids) == 16 and len(set(ids.tolist())) == 16
    L = escher_lifted_torus_laplacian()
    assert np.linalg.norm(L-L.T) < 1e-12
    assert np.linalg.norm(L @ np.ones(64)) < 1e-12


def test_unified_audit_and_prism_constant_intertwiner():
    r = audit_unified_operator()
    assert r["dimension"] == 384
    assert r["constant_mode_residual"] < 1e-10
    assert r["prism_intertwiner"]["rank"] == 1
    assert r["prism_intertwiner"]["intertwining_residual"] < 1e-10
    assert "only" in r["prism_intertwiner"]["interpretation"]


def test_invalid_weights():
    with pytest.raises(ValueError):
        build_strange_loop_operator(kappa=-1)


def test_cantor_emissions_are_coupled_to_parry_q6_and_stochastic():
    L, m = cantor_parry_q6_laplacian()
    assert L.shape == (192, 192)
    assert np.linalg.norm(L-L.T) < 1e-12
    assert np.linalg.norm(L @ np.ones(192)) < 1e-11
    assert np.all(m["cantor_p0"] > 0) and np.all(m["cantor_p1"] > 0)
    assert np.max(np.abs(m["cantor_p0"]+m["cantor_p1"]-1)) < 1e-12
    assert np.max(np.abs(m["parry_transition"].sum(axis=1)-1)) < 1e-12


def test_doppler_moire_phase_separation_and_relativistic_factor():
    r = doppler_moire_phase_audit(
        omega0=2.0, beta=0.3, k1=(1., 0.), k2=(1.1, 0.1), t=2.
    )
    expected = (1./np.sqrt(1.-.3**2))*(1.-.3)
    assert abs(r["doppler_factor"]-expected) < 1e-12
    assert abs(r["moire_length"]-2.*np.pi/np.sqrt(.02)) < 1e-10
    assert abs(r["delta_omega_observed"]+expected*.1) < 1e-12
    assert r["beat_period"] is not None
    assert len(r["moire_phase"]) == 64
    assert len(r["spacetime_phase"]) == 64
    assert "temporal" in r["phase_identity"] and "spatial" in r["phase_identity"]


def test_doppler_moire_extended_operator_is_hermitian():
    H, omega, _ = build_doppler_moire_operator(moire_strength=.25, beta=.2)
    assert H.shape == (384, 384)
    assert np.linalg.norm(H-H.T) < 1e-11
    assert np.linalg.eigvalsh(H)[0] > -1e-10
    assert omega.shape == (384,)
    assert np.allclose(omega, omega[0])


def test_barycenter_first_moment_covariance_and_exact_moire_identity():
    X = np.array([[0., 0.], [1., 0.], [0., 2.], [1., 2.]])
    w = np.array([1., 2., 3., 4.])
    dk = np.array([.12, -.08])
    g = np.array([[2., .2], [.2, 1.5]])
    r = barycentric_gradient_moire_audit(
        X, w, dk, metric=g, delta_omega=.3, phase0=.2, t=.7, eta=.4
    )
    assert r["first_moment_residual"] < 1e-13
    assert np.linalg.eigvalsh(np.asarray(r["covariance"]))[0] >= -1e-13
    assert abs(r["moire_average_exact"]-r["moire_average_characteristic"]) < 1e-12
    assert np.allclose(r["metric_phase_gradient"], np.linalg.solve(g, -dk))
    assert np.isfinite(r["moire_average_quadratic"])


def test_barycentric_quadratic_expansion_converges_for_small_wavevector_gap():
    X = np.array([[0., 0.], [1., 0.], [0., 1.], [1., 1.]])
    w = np.array([.1, .2, .3, .4])
    errors = []
    for eps in (1e-2, 5e-3):
        r = barycentric_gradient_moire_audit(X, w, [eps, .7*eps], phase0=.4)
        errors.append(abs(r["moire_average_exact"]-r["moire_average_quadratic"]))
    assert errors[1] < errors[0]


def test_moving_barycenter_includes_weight_transport():
    X = np.array([[0., 0.], [2., 1.], [1., 3.]])
    w = np.array([.2, .3, .5])
    V = np.array([[1., 0.], [0., 1.], [-1., 0.]])
    wd = np.array([.1, -.04, -.06])
    r = moving_barycenter_audit(X, w, V, wd)
    assert abs(r["normalized_weight_rate_sum"]) < 1e-14
    assert r["derivative_identity_residual"] < 1e-12


def test_barycenter_rejects_indefinite_metric_and_negative_weights():
    with pytest.raises(ValueError):
        barycentric_gradient_moire_audit([[0., 0.], [1., 1.]], [1., 1.], [1., 0.],
                                          metric=[[1., 2.], [2., 1.]])
    with pytest.raises(ValueError):
        barycentric_gradient_moire_audit([[0., 0.], [1., 1.]], [1., -1.], [1., 0.])



def test_q6_projection_recovers_affine_rgs_seed_exactly():
    from src.unified_strange_loop_operator import derive_q6_projection_from_rgs_seed
    B=np.array([[(g>>a)&1 for a in range(6)] for g in range(64)],dtype=float)
    P=np.array([[.3,-.2,.1,.4,.0,-.1],
                [.1,.2,-.3,.0,.5,.2],
                [-.2,.1,.0,.3,-.1,.4]])
    offset=np.array([1.2,-.7,.4])
    X=B@P.T+offset
    r=derive_q6_projection_from_rgs_seed(X)
    assert r["rank"]==6
    assert r["weighted_relative_residual"]<1e-12
    assert np.allclose(r["projection_3x6"],P,atol=1e-12)
    assert np.allclose(r["intercept_3"],offset,atol=1e-12)


def test_q6_projection_reports_non_affine_geometry_instead_of_hiding_it():
    from src.unified_strange_loop_operator import derive_q6_projection_from_rgs_seed
    B=np.array([[(g>>a)&1 for a in range(6)] for g in range(64)],dtype=float)
    X=np.column_stack((B[:,0],B[:,1],B[:,2]))
    X[:,2]+=0.07*((B[:,3]+B[:,4]+B[:,5])%2)
    r=derive_q6_projection_from_rgs_seed(X)
    assert r["rank"]==6
    assert r["weighted_relative_residual"]>1e-3
    assert r["max_point_residual"]>0


def test_q6_projection_rejects_wrong_seed_shape_and_bad_weights():
    from src.unified_strange_loop_operator import derive_q6_projection_from_rgs_seed
    with pytest.raises(ValueError):
        derive_q6_projection_from_rgs_seed(np.zeros((63,3)))
    with pytest.raises(ValueError):
        derive_q6_projection_from_rgs_seed(np.zeros((64,3)),weights=np.ones(63))
    L,m=cantor_parry_q6_laplacian()
    assert L.shape==(192,192)
    assert np.linalg.norm(L-L.T)<1e-12
    assert np.linalg.norm(L@np.ones(192))<1e-11
    assert np.all(m["cantor_p0"]>0) and np.all(m["cantor_p1"]>0)
    assert np.max(np.abs(m["cantor_p0"]+m["cantor_p1"]-1))<1e-12
    assert np.max(np.abs(m["parry_transition"].sum(axis=1)-1))<1e-12

def test_doppler_moire_phase_separation_and_relativistic_factor():
    r=doppler_moire_phase_audit(omega0=2.0,beta=0.3,k1=(1.,0.),k2=(1.1,0.1),t=2.)
    expected=(1./np.sqrt(1.-.3**2))*(1.-.3)
    assert abs(r["doppler_factor"]-expected)<1e-12
    assert abs(r["moire_length"]-2.*np.pi/np.sqrt(.02))<1e-10
    # Lorentz-transform the frequency difference: Delta omega' = gamma(Delta omega - beta Delta k_x).
    gamma=1./np.sqrt(1.-.3**2)
    assert abs(r["delta_omega_observed"]+gamma*.3*.1)<1e-12
    assert r["beat_period"] is not None
    assert len(r["moire_phase"])==64
    assert len(r["spacetime_phase"])==64
    assert "temporal" in r["phase_identity"] and "spatial" in r["phase_identity"]

def test_doppler_moire_extended_operator_is_hermitian():
    H,omega,_=build_doppler_moire_operator(moire_strength=.25,beta=.2)
    assert H.shape==(384,384)
    assert np.linalg.norm(H-H.T)<1e-11
    assert np.linalg.eigvalsh(H)[0]>-1e-10
    assert omega.shape==(384,)
    assert np.allclose(omega,omega[0])
