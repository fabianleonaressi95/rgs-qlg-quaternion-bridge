import numpy as np
from src.rgs_ellis_refinement import load_rgs_seed_data
from src.rgs_fourth_moment_bridge import fourth_moment_decomposition
from src.rgs_fourth_moment_connection import (
    positive_constitutive_metric, covariant_derivative_q,
    gradient_energy_density, spin_current, audit_connection_current,
)

def fixture():
    x=load_rgs_seed_data().coordinates
    r=fourth_moment_decomposition(x)["M4_anisotropic"]
    q=np.diag([.4,-.1,-.3])
    p=np.array([[.2,.1,0.],[.1,-.05,.03],[0.,.03,-.15]])
    w=np.array([[0.,.12,-.04],[-.12,0.,.07],[.04,-.07,0.]])
    return r,q,p,w

def test_shifted_fourth_moment_metric_is_positive():
    r,_,_,_=fixture(); out=positive_constitutive_metric(r,beta=1.0)
    assert out["minimum_eigenvalue"]>0
    assert np.all(np.linalg.eigvalsh(out["operator"])>0)

def test_covariant_derivative_remains_symmetric_tracefree():
    _,q,p,w=fixture(); d=covariant_derivative_q(p,w,q)
    assert np.allclose(d,d.T,atol=1e-12)
    assert abs(np.trace(d))<1e-12

def test_gradient_energy_is_nonnegative():
    r,q,p,w=fixture(); h=positive_constitutive_metric(r)["operator"]
    e=gradient_energy_density(np.stack([p,2*p]),np.stack([w,w]),q,h,.8)
    assert e>=0

def test_spin_current_is_antisymmetric_and_matches_connection_variation():
    r,q,p,w=fixture()
    out=audit_connection_current(r,q,p,w)
    assert out["minimum_positive_eigenvalue"]>0
    assert out["current_antisymmetry_residual"]<1e-12
    assert out["connection_fd_max_abs"]<2e-7
    assert out["status"]["teleparallel_flatness_constraint"]=="OPEN"

def test_zero_gradient_has_zero_spin_current():
    r,q,_,w=fixture(); h=positive_constitutive_metric(r)["operator"]
    j=spin_current(q,np.zeros((3,3)),h)
    assert np.linalg.norm(j)==0.0
