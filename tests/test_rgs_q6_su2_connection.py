import numpy as np
import pytest
from src.rgs_q6_su2_connection import (
    audit_q6_su2_connection, random_su2_matrices, validate_su2_field,
    covariant_coboundary, dirac_kaehler_operator, gauge_transform_links,
    gauge_covariance_matrices, plaquette_holonomies,
)

def test_random_links_are_su2_and_invalid_links_rejected():
    u=random_su2_matrices(192,seed=10)
    r=validate_su2_field(u,192)
    assert r["unitarity_residual"]<1e-12 and r["determinant_residual"]<1e-12
    with pytest.raises(ValueError): validate_su2_field(np.zeros((192,2,2)),192)

def test_covariant_incidence_and_dirac_kaehler_hermiticity():
    u=random_su2_matrices(192,seed=11)
    d=covariant_coboundary(u)
    k=dirac_kaehler_operator(u)
    assert d.shape==(384,128) and k.shape==(512,512)
    assert np.linalg.norm(k-k.conj().T)<1e-12
    assert np.linalg.eigvalsh(k)[0]>-1e-10

def test_gauge_covariance_and_wilson_loop_invariance():
    u=random_su2_matrices(192,seed=12)
    g=random_su2_matrices(64,seed=13)
    ug=gauge_transform_links(u,g)
    d=covariant_coboundary(u); dg=covariant_coboundary(ug)
    p0,p1=gauge_covariance_matrices(g)
    assert np.linalg.norm(dg-p1@d@p0.conj().T)<1e-10
    h=plaquette_holonomies(u); hg=plaquette_holonomies(ug)
    assert np.max(np.abs(np.trace(h,axis1=1,axis2=2)-np.trace(hg,axis1=1,axis2=2)))<1e-10

def test_full_su2_connection_audit():
    report=audit_q6_su2_connection()
    assert report["all_checks_passed"]
    assert report["plaquettes"]==240
