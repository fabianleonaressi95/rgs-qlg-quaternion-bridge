import numpy as np
import pytest
from src.rgs_q6_euclidean_time_wilson import (
    audit_q6_euclidean_time_wilson, build_q6_euclidean_time_data,
    q6_euclidean_gauge_transform, q6_euclidean_plaquette_fluxes,
    q6_euclidean_wilson_action, q6_lattice_gauss_law,
)


def test_periodic_time_cell_counts_and_boundary_identity():
    d = build_q6_euclidean_time_data(4)
    assert (d["n_spatial_links"], d["n_temporal_links"]) == (192, 64)
    assert (d["n_spatial_plaquettes"], d["n_temporal_plaquettes"]) == (240, 192)
    assert np.max(np.abs(d["D1"] @ d["D0"])) == 0


def test_plaquettes_and_action_are_gauge_invariant():
    nt = 3
    rng = np.random.default_rng(17)
    s, t = rng.uniform(-1, 1, (nt,192)), rng.uniform(-1, 1, (nt,64))
    lam = rng.uniform(-np.pi, np.pi, (nt,64))
    d = build_q6_euclidean_time_data(nt)
    sm, tm = q6_euclidean_gauge_transform(s,t,lam,d)
    f, fm = q6_euclidean_plaquette_fluxes(s,t,d), q6_euclidean_plaquette_fluxes(sm,tm,d)
    assert np.allclose(np.exp(1j*f["spatial"]),np.exp(1j*fm["spatial"]),atol=1e-10)
    assert np.allclose(np.exp(1j*f["temporal"]),np.exp(1j*fm["temporal"]),atol=1e-10)
    assert q6_euclidean_wilson_action(s,t,1.3,.8,d) == pytest.approx(
        q6_euclidean_wilson_action(sm,tm,1.3,.8,d),abs=1e-10)


def test_temporal_link_variation_is_discrete_gauss_law():
    nt=3; rng=np.random.default_rng(19)
    s,t=rng.uniform(-.4,.4,(nt,192)),rng.uniform(-.4,.4,(nt,64))
    d=build_q6_euclidean_time_data(nt); gauss=q6_lattice_gauss_law(s,t,d)
    eps,ti,vi,beta_t=1e-6,1,9,1.7
    tp,tm=t.copy(),t.copy(); tp[ti,vi]+=eps; tm[ti,vi]-=eps
    fd=(q6_euclidean_wilson_action(s,tp,1.2,beta_t,d)
        -q6_euclidean_wilson_action(s,tm,1.2,beta_t,d))/(2*eps)
    assert fd == pytest.approx(beta_t*gauss[ti,vi],abs=2e-7)


def test_full_euclidean_time_gauss_audit():
    report=audit_q6_euclidean_time_wilson(nt=3)
    assert report["all_checks_passed"], report["checks"]
    assert report["spacetime_interpretation"].startswith("6+1-dimensional")


@pytest.mark.parametrize("nt",[0,1,1.5])
def test_invalid_time_extent_rejected(nt):
    with pytest.raises(ValueError):
        build_q6_euclidean_time_data(nt)
