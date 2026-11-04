import numpy as np
import pytest
from src.rgs_q6_d5_a5_degree_six import (_invariant_harmonic_values,
    _uniform_sphere_harmonic_rms,build_q6_d5_degree_six_coefficient_audit)
from src.rgs_icosahedral_q6_intertwiner import icosahedral_signed_axis_representation

def test_reynolds_degree_six_harmonic_is_a5_invariant():
    r=np.asarray(icosahedral_signed_axis_representation()["rotations"],dtype=float)
    n=np.asarray([[.2,.3,.9]]); n/=np.linalg.norm(n,axis=1)[:,None]
    value=_invariant_harmonic_values(n,r)[0]
    assert _invariant_harmonic_values((r[17]@n[0])[None,:],r)[0]==pytest.approx(value,abs=2e-12)

def test_uniform_sphere_harmonic_rms_is_positive():
    r=np.asarray(icosahedral_signed_axis_representation()["rotations"],dtype=float)
    assert _uniform_sphere_harmonic_rms(r)>1e-8

def test_degree_six_coefficient_audit_is_reproducible_and_finite():
    args=dict(samples_per_batch=512,batches=4,seed=12345)
    a=build_q6_d5_degree_six_coefficient_audit(**args)
    b=build_q6_d5_degree_six_coefficient_audit(**args)
    assert a["all_checks_passed"]
    assert a["finite_q6_design"]["direction_count"]==32
    assert np.isfinite(a["continuous_s9_pushforward"]["normalized_degree_six_coefficient"])
    assert a["continuous_s9_pushforward"]["standard_error"]>0
    assert a["difference_continuous_minus_finite"]==b["difference_continuous_minus_finite"]
