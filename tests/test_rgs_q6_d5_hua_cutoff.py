import math
import numpy as np
import pytest
from src.rgs_q6_d5_hua_cutoff import (
    d5_hua_generic_norm_density, d5_hua_cutoff_normalizer,
    d5_hua_cutoff_probability_density, d5_hua_cutoff_report,
)

def test_exact_factorized_hua_density_and_endpoints():
    t=np.linspace(0.,1.,1001)
    f=d5_hua_generic_norm_density(t)
    assert np.min(f)>=0.
    assert d5_hua_generic_norm_density(0.)==pytest.approx(20./3.)
    assert d5_hua_generic_norm_density(1.)==0.
    with pytest.raises(ValueError):
        d5_hua_generic_norm_density(-0.1)

@pytest.mark.parametrize("s",[-0.9,-0.5,0.,1.,2.,5.])
def test_density_mellin_transform_is_hua_ratio(s):
    mellin=(20./3.)/(s+1.)-60./(s+2.)-60./(s+3.)+(320./3.)/(s+2.5)+(20./3.)/(s+4.)
    chi=math.prod(s+k for k in (1.,2.,3.,4.))*(s+2.5)
    assert mellin==pytest.approx(60./chi,rel=2.e-13)

@pytest.mark.parametrize("epsilon",[0.5,0.2,0.1,0.05,0.01])
def test_exact_cutoff_normalizer_matches_quadrature(epsilon):
    nodes,weights=np.polynomial.legendre.leggauss(256)
    t=epsilon+(nodes+1.)*(1.-epsilon)/2.
    w=weights*(1.-epsilon)/2.
    numerical=float(np.dot(w,d5_hua_generic_norm_density(t)*t**-5))
    assert d5_hua_cutoff_normalizer(epsilon)==pytest.approx(numerical,rel=2.e-9)

def test_cutoff_measure_is_finite_and_normalized():
    epsilon=0.1
    nodes,weights=np.polynomial.legendre.leggauss(256)
    t=epsilon+(nodes+1.)*(1.-epsilon)/2.
    w=weights*(1.-epsilon)/2.
    density=d5_hua_cutoff_probability_density(t,epsilon)
    assert np.all(np.isfinite(density))
    assert float(np.dot(w,density))==pytest.approx(1.,rel=2.e-9)
    assert d5_hua_cutoff_probability_density(epsilon/2.,epsilon)==0.

def test_cutoff_asymptotic_and_report():
    epsilon=1.e-5
    assert d5_hua_cutoff_normalizer(epsilon)*epsilon**4==pytest.approx(5./3.,rel=2.e-4)
    assert d5_hua_cutoff_report()["all_checks_passed"]

@pytest.mark.parametrize("epsilon",[0.,1.,-0.1,float("nan"),float("inf")])
def test_invalid_cutoff_rejected(epsilon):
    with pytest.raises(ValueError):
        d5_hua_cutoff_normalizer(epsilon)
