import math
import numpy as np
import pytest
from src.rgs_q6_u1_quotient import (audit_q6_u1_quotient,gauge_transform_u1,
 q6_u1_quotient_data,sample_quotient_haar,tree_gauge_fix,u1_quadratic_action,u1_wilson_action,lift_chord_angles)


def test_q6_quotient_dimensions_and_rank():
 q=q6_u1_quotient_data()
 assert len(q["tree_edges"])==63 and len(q["chord_edges"])==129
 assert q["plaquette_map_rank"]==129
 assert np.max(np.abs(q["D1"]@q["D0"]))==0


def test_tree_gauge_fix_is_orbit_invariant():
 q=q6_u1_quotient_data(); rng=np.random.default_rng(11)
 x=rng.uniform(-1,1,192); lam=rng.uniform(-math.pi,math.pi,64)
 a,_=tree_gauge_fix(x,q); b,_=tree_gauge_fix(gauge_transform_u1(x,lam,q),q)
 assert np.allclose(np.exp(1j*a[list(q["chord_edges"])]),np.exp(1j*b[list(q["chord_edges"])]),atol=1e-10)


def test_gauge_and_icosahedral_audits():
 r=audit_q6_u1_quotient()
 assert r["all_checks_passed"]
 assert r["gauge_residual"]<1e-10 and r["icosahedral_residual"]<1e-10


def test_wilson_action_has_quadratic_small_field_limit():
 q=q6_u1_quotient_data(); x=1e-4*np.random.default_rng(12).normal(size=192)
 a=u1_wilson_action(x,1.7,q); b=u1_quadratic_action(x,1.7,q)
 assert b>0 and abs(a-b)/b<1e-7


def test_normalized_haar_coordinates_and_lift():
 x=sample_quotient_haar(7,13); q=q6_u1_quotient_data(); y=lift_chord_angles(x,q)
 assert x.shape==(7,129) and y.shape==(7,192)
 assert np.all(x>=-math.pi) and np.all(x<math.pi)
 assert np.all(y[...,list(q["tree_edges"])]==0)


@pytest.mark.parametrize("beta",[0.,-1.,math.inf,math.nan])
def test_invalid_beta_rejected(beta):
 with pytest.raises(ValueError): u1_wilson_action(np.zeros(192),beta)
