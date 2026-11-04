import numpy as np
import pytest
from src.rgs_q6_d5_liouville_measure import d5_bergman_metric, d5_generic_norm, q6_d5_bergman_liouville_report, q6_maxwell_to_d5_chart

def test_origin_metric_and_volume_match_liouville():
    r=q6_d5_bergman_liouville_report(q=np.zeros(5),p=np.zeros(5))
    assert r["all_checks_passed"]
    assert np.allclose(d5_bergman_metric(np.zeros(5,dtype=complex)),10*np.eye(5))
    assert r["bergman_volume_pullback_density"] == pytest.approx(1.)

@pytest.mark.parametrize("beta,eigenvalue",[(.2,12.),(1.,12.),(3.7,8.),(11.,4.)])
def test_exact_density_law(beta,eigenvalue):
    q=np.array([.07,-.03,.02,.01,-.04]); p=np.array([-.02,.05,.01,-.06,.03])
    r=q6_d5_bergman_liouville_report(q,p,beta,eigenvalue)
    assert r["all_checks_passed"]
    assert r["bergman_volume_pullback_density"] == pytest.approx(r["generic_norm"]**-5,rel=2e-9)

def test_type_iv_determinant_formula():
    z=np.array([.05+.02j,.1-.03j,0,0,0],complex)
    n,_=d5_generic_norm(z)
    assert np.linalg.det(d5_bergman_metric(z)).real == pytest.approx(10**5/n**5,rel=2e-9)

def test_invalid_inputs_rejected():
    with pytest.raises(ValueError): q6_maxwell_to_d5_chart(np.zeros(5),np.zeros(5),beta=0)
    with pytest.raises(ValueError): q6_maxwell_to_d5_chart(np.zeros(4),np.zeros(5))
    with pytest.raises(ValueError): q6_d5_bergman_liouville_report(np.ones(5)*2,np.zeros(5))
