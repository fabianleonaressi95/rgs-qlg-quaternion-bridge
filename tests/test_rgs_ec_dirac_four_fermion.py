import numpy as np
import pytest
from src.rgs_ec_algebraic_torsion import axial_torsion_from_dirac_spinor, dirac_axial_current
from src.rgs_ec_dirac_four_fermion import (
    audit_four_fermion_reduction, effective_eom_source,
    effective_four_fermion_coefficient, effective_four_fermion_density,
    multi_flavour_effective_density, raw_to_real_axial_current,
    real_to_raw_axial_current,
)


def test_raw_and_real_current_conventions_have_same_density():
    j=np.array([0.2,-0.1,0.3,0.4])
    a=audit_four_fermion_reduction(j,8.0)
    assert a["coefficient_for_real_J5_squared"] == 1.5
    assert a["coefficient_for_raw_A_squared"] == -1.5
    assert a["convention_conversion_residual"] < 1e-14


def test_spinor_current_to_torsion_and_four_fermion_pipeline():
    psi=np.array([0.2+0.1j,0.5-0.3j,-0.4+0.2j,0.7+0.6j])
    j=dirac_axial_current(psi)
    t=axial_torsion_from_dirac_spinor(psi,8.0)
    assert np.allclose(t["spinor_axial_current"],j)
    assert np.isfinite(effective_four_fermion_density(j,8.0))


def test_effective_source_is_gradient_of_density():
    j=np.array([0.2,-0.1,0.3,0.4]); k=2.7; eps=1e-6
    analytic=effective_eom_source(j,k); numerical=np.empty(4)
    for i in range(4):
        d=np.zeros(4); d[i]=eps
        numerical[i]=(effective_four_fermion_density(j+d,k)-effective_four_fermion_density(j-d,k))/(2*eps)
    assert np.linalg.norm(analytic-numerical)<1e-9


def test_multiflavour_current_is_summed_before_squaring():
    r=multi_flavour_effective_density([[1.,0,0,0],[.5,0,0,0]],1.0)
    assert np.allclose(r["total_axial_current"],[1.5,0,0,0])
    assert r["includes_cross_flavour_terms"]
    assert np.isclose(r["density"],-3/16*1.5**2)


def test_raw_real_round_trip_and_bad_kappa():
    j=np.array([.1,.2,-.3,.4])
    assert np.allclose(raw_to_real_axial_current(real_to_raw_axial_current(j)),j)
    with pytest.raises(ValueError): effective_four_fermion_coefficient(0)
