import math
import numpy as np
from src.rgs_smilga_two_particle_audit import (
    audit, boost, cm_pair, dot, from_total_relative, integrated_phase_space,
    jacobian_abs_8d, phase_space_density_dOmega, residuals,
    to_total_relative, wyler_alpha_inverse,
)


def test_transform_inverse_and_jacobian():
    rng = np.random.default_rng(7)
    a, b = rng.normal(size=(12,4)), rng.normal(size=(12,4))
    P, q = to_total_relative(a,b)
    aa, bb = from_total_relative(P,q)
    assert np.allclose(aa,a,atol=1e-14)
    assert np.allclose(bb,b,atol=1e-14)
    assert math.isclose(jacobian_abs_8d(),1.,abs_tol=1e-14)


def test_mass_shell_and_casimir_are_lorentz_covariant():
    dirs = np.array([[1.,0.,0.],[0.,1.,0.],[1.,2.,3.]])
    a,b = cm_pair(3.,1.,dirs)
    a,b = boost(a,np.array([.2,-.1,.15])), boost(b,np.array([.2,-.1,.15]))
    for x in residuals(a,b,1.).values():
        assert np.max(np.abs(x)) < 2e-13
    P,q = to_total_relative(a,b)
    assert np.allclose(dot(P,P),9.,atol=2e-13)
    assert np.allclose(dot(q,q),-1.25,atol=2e-13)


def test_two_body_phase_space_normalization():
    M,m=3.,1.
    assert math.isclose(4*math.pi*phase_space_density_dOmega(M,m),
                        integrated_phase_space(M,m),rel_tol=1e-14)
    assert math.isclose(integrated_phase_space(M,m),
                        math.sqrt(M*M/4-m*m)/(4*math.pi*M))
    assert integrated_phase_space(2.,1.) == 0.


def test_wyler_numeric_is_separate_from_phase_space_claim():
    assert math.isclose(wyler_alpha_inverse(),137.0360824481643,rel_tol=2e-14)
    assert "open" in audit()
