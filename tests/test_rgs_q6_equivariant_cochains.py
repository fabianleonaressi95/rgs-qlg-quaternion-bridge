import numpy as np
from src.rgs_q6_equivariant_cochains import build_q6_cochain_complex, induced_cochain_actions, audit_q6_equivariant_cochains
from src.rgs_icosahedral_q6_intertwiner import icosahedral_signed_axis_representation

def test_q6_cubical_cochain_dimensions_and_boundary_squared_zero():
    c=build_q6_cochain_complex()
    assert c["D0"].shape==(192,64) and c["D1"].shape==(240,192)
    assert np.array_equal(c["D1"]@c["D0"],np.zeros((240,64),dtype=np.int8))

def test_icosahedral_action_lifts_to_orthogonal_cochain_representations():
    c=build_q6_cochain_complex()
    for s in icosahedral_signed_axis_representation()["signed_permutations"]:
        p0,p1,p2=induced_cochain_actions(s,c)
        assert np.array_equal(p0.T@p0,np.eye(64,dtype=np.int8))
        assert np.array_equal(p1.T@p1,np.eye(192,dtype=np.int8))
        assert np.array_equal(p2.T@p2,np.eye(240,dtype=np.int8))

def test_full_equivariant_cochain_audit():
    r=audit_q6_equivariant_cochains()
    assert r["rank_D0"]==63 and r["rank_D1"]==129
    assert r["d1_d0_max_abs"]==0
    assert r["D0_equivariance_max_abs"]==0 and r["D1_equivariance_max_abs"]==0
    assert r["cochain_orthogonality_failures"]==0 and r["all_checks_passed"]
