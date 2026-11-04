import numpy as np
from src.rgs_fourth_moment_bridge import fourth_moment_decomposition
from src.rgs_ellis_refinement import load_rgs_seed_data
from src.rgs_fourth_moment_variation import (
    anisotropic_potential, coframe_action_density, coframe_euler_derivative,
    lapse_euler_derivative, q_euler_derivative, rotate_rank2, rotate_rank4,
    audit_variational_bridge,
)

def fixture():
    r=fourth_moment_decomposition(load_rgs_seed_data().coordinates)["M4_anisotropic"]
    return r, np.diag([.4,-.1,-.3])

def test_coframe_variation_matches_central_difference():
    r,q=fixture(); e=np.array([[1.1,.1,0.],[0.,.9,.05],[.02,0.,1.2]])
    h=1e-6; analytic=coframe_euler_derivative(e,1.2,r,q,.7); numeric=np.zeros((3,3))
    for a in range(3):
        for i in range(3):
            d=np.zeros((3,3)); d[a,i]=h
            numeric[a,i]=(coframe_action_density(e+d,1.2,r,q,.7)-coframe_action_density(e-d,1.2,r,q,.7))/(2*h)
    assert np.allclose(analytic,numeric,atol=2e-7,rtol=2e-7)

def test_lapse_and_q_variations():
    r,q=fixture(); e=np.diag([1.1,.9,1.2]); h=1e-6
    fd=(coframe_action_density(e,1.2+h,r,q)-coframe_action_density(e,1.2-h,r,q))/(2*h)
    assert abs(fd-lapse_euler_derivative(e,r,q))<1e-7
    b=np.diag([1.,-1.,0.])/np.sqrt(2); g=q_euler_derivative(r,q,.7)
    fdq=(anisotropic_potential(r,q+h*b,.7)-anisotropic_potential(r,q-h*b,.7))/(2*h)
    assert abs(fdq-np.einsum("ab,ab",g,b))<1e-7

def test_simultaneous_internal_tensor_rotation_invariance():
    r,q=fixture(); t=.37
    u=np.array([[np.cos(t),-np.sin(t),0.],[np.sin(t),np.cos(t),0.],[0.,0.,1.]])
    assert np.isclose(anisotropic_potential(r,q),anisotropic_potential(rotate_rank4(r,u),rotate_rank2(q,u)),atol=1e-9)

def test_complete_variational_audit():
    r,q=fixture(); out=audit_variational_bridge(r,q)
    assert out["coframe_gradient_fd_max_abs"]<2e-7
    assert out["lapse_gradient_fd_abs"]<1e-7
    assert out["q_gradient_sym0_fd_max_abs"]<1e-7
    assert out["status"]["teleparallel_torsion_generation"]=="NOT_GENERATED_BY_THIS_ALGEBRAIC_TERM"
