import math
import numpy as np
import pytest
from src.rgs_q6_d5_liouville_measure import d5_bergman_metric, d5_generic_norm, q6_d5_bergman_liouville_report, q6_maxwell_to_d5_chart, q6_d5_liouville_flow_report
from src.rgs_q6_d5_liouville_measure import (
    d5_bergman_metric, d5_generic_norm, d5_bergman_cutoff_volume,
    q6_d5_bergman_liouville_report, q6_maxwell_to_d5_chart, q6_d5_liouville_flow_report,
)


def test_origin_metric_and_volume_match_liouville():
    r = q6_d5_bergman_liouville_report(q=np.zeros(5), p=np.zeros(5))
    assert r["all_checks_passed"], r["checks"]
    assert np.allclose(d5_bergman_metric(np.zeros(5, dtype=complex)), 10*np.eye(5))
    assert r["bergman_volume_pullback_density"] == pytest.approx(1.)


@pytest.mark.parametrize("beta,eigenvalue", [(.2,12.), (1.,12.), (3.7,8.), (11.,4.)])
def test_exact_density_and_maxwell_hamiltonian_law(beta, eigenvalue):
    q = np.array([.07,-.03,.02,.01,-.04])
    p = np.array([-.02,.05,.01,-.06,.03])
    r = q6_d5_bergman_liouville_report(q, p, beta, eigenvalue)
    assert r["all_checks_passed"], r["checks"]
    assert r["bergman_volume_pullback_density"] == pytest.approx(r["generic_norm"]**-5, rel=2e-9)
    assert r["Maxwell_Hamiltonian"] == pytest.approx(r["chart_Hamiltonian"], rel=2e-14)
    assert r["frequency"] == pytest.approx(np.sqrt(eigenvalue))


def test_type_iv_determinant_formula():
    z = np.array([.05+.02j,.1-.03j,0,0,0], complex)
    n, _ = d5_generic_norm(z)
    assert np.linalg.det(d5_bergman_metric(z)).real == pytest.approx(10**5/n**5, rel=2e-9)


def test_exact_hua_cutoff_volume_in_liouville_units():
    epsilon = .05
    report = d5_bergman_cutoff_volume(epsilon)
    expected = 20.**5 * (math.pi**5 / 1920.) * report["normalized_Euclidean_Hua_factor"]
    assert report["Liouville_Bergman_volume"] == pytest.approx(expected, rel=1e-14)
    assert report["asymptotic_leading_coefficient"] == pytest.approx(25000*math.pi**5/9, rel=1e-14)
    assert report["Liouville_Bergman_volume"] > 0
    with pytest.raises(ValueError):
        d5_bergman_cutoff_volume(0.)


def test_invalid_inputs_rejected():
    with pytest.raises(ValueError): q6_maxwell_to_d5_chart(np.zeros(5),np.zeros(5),beta=0)
    with pytest.raises(ValueError): q6_maxwell_to_d5_chart(np.zeros(4),np.zeros(5))
    with pytest.raises(ValueError): q6_d5_bergman_liouville_report(np.ones(5)*2,np.zeros(5))


@pytest.mark.parametrize("beta,eigenvalue",[(.2,12.),(1.,12.),(3.7,8.),(11.,4.)])
    with pytest.raises(ValueError):
        q6_maxwell_to_d5_chart(np.zeros(5), np.zeros(5), beta=0)
    with pytest.raises(ValueError):
        q6_maxwell_to_d5_chart(np.zeros(4), np.zeros(5))
    with pytest.raises(ValueError):
        q6_d5_bergman_liouville_report(np.ones(5)*2, np.zeros(5))


@pytest.mark.parametrize("beta,eigenvalue", [(.2,12.), (1.,12.), (3.7,8.), (11.,4.)])
def test_maxwell_flow_is_common_phase_and_preserves_bergman_liouville_weight(beta,eigenvalue):
    r=q6_d5_liouville_flow_report(beta=beta,eigenvalue=eigenvalue,time=.231)
    assert r["all_checks_passed"], r["checks"]
    assert r["generic_norm_initial"] == pytest.approx(r["generic_norm_final"],rel=2e-12)
    assert r["bergman_density_initial"] == pytest.approx(r["bergman_density_final"],rel=2e-11)


def test_chart_uses_maxwell_canonical_reciprocal_scale():
    beta,lam=3.7,8.
    q=np.array([.07,-.03,.02,.01,-.04])
    p=np.array([-.02,.05,.01,-.06,.03])
    z=q6_maxwell_to_d5_chart(q,p,beta,lam)
    kappa=beta*np.sqrt(lam)
    expected=(np.sqrt(kappa)*q+1j*p/np.sqrt(kappa))/np.sqrt(20.)
    assert np.allclose(z,expected,rtol=0,atol=1e-14)
