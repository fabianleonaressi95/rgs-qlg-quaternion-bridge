import numpy as np
import pytest

from src.rgs_ec_algebraic_torsion import (
    DIRAC_EC_CONVENTION_ID, axial_torsion_from_current,
)
from src.rgs_ec_spin_ensemble import (
    audit_spin_ensemble, build_cosmological_spin_state,
    cartan_normalized_spin_source, spatial_isotropy_residual,
    torsion_from_cartan_source,
)


def test_cartan_normalized_source_reproduces_action_derived_torsion():
    j=np.array([0.3,-0.2,0.1,0.4]); kappa=2.3
    expected=axial_torsion_from_current(j,kappa)["torsion"]
    assert np.allclose(torsion_from_cartan_source(j,kappa),expected)
    assert np.allclose(kappa*cartan_normalized_spin_source(j),expected)


def test_unpolarized_pair_has_zero_mean_but_nonzero_spin_fluctuation_contact():
    j=np.array([0.2,0.1,-0.3,0.4])
    result=audit_spin_ensemble([j,-j],1.7)
    assert np.allclose(result["mean_axial_current"],0)
    assert np.linalg.norm(result["mean_torsion"]) == 0
    assert result["mean_of_lorentz_norm_squared"] == pytest.approx(
        -j[0]**2+j[1]**2+j[2]**2+j[3]**2)
    assert result["ensemble_mean_four_fermion_density"] != pytest.approx(0)
    assert result["fluctuation_contact_contribution"] == pytest.approx(
        result["ensemble_mean_four_fermion_density"])
    assert result["moment_decomposition_residual"] < 1e-14


def test_isotropic_spatial_covariance_diagnostic():
    c=np.diag([0.2,0.5,0.5,0.5])
    r=spatial_isotropy_residual(c)
    assert r["spatial_isotropy_residual"] == pytest.approx(0)
    assert r["time_space_covariance_norm"] == 0


def test_weighted_moments_and_flavour_averaging():
    x=np.array([[1.,0,0,0],[-1.,0,0,0],[3.,0,0,0]])
    r=audit_spin_ensemble(x,1.0,weights=[1,1,2])
    assert np.allclose(r["mean_axial_current"],[1.5,0,0,0])
    assert r["sample_count"] == 3


def test_adapter_refuses_implicit_units_and_does_not_invent_beta():
    with pytest.raises(TypeError):
        build_cosmological_spin_state([[1.,0,0,0]],1.0,scale_factor=1,
            hubble_rate=1,energy_density=1,pressure=0)
    s=build_cosmological_spin_state([[1.,0,0,0]],1.0,scale_factor=1,
        hubble_rate=.1,energy_density=2,pressure=.3,
        units_id="natural-units",derivation_id="declared-test")
    assert s["torsion_lambda_beta"] is None
    assert s["effective_energy_density"] is None
    assert s["status"] == "SPIN_SOURCE_READY_COSMOLOGICAL_REDUCTION_OPEN"


def test_rejects_bad_weights_and_nonpositive_scale_factor():
    with pytest.raises(ValueError): audit_spin_ensemble([[1.,0,0,0]],1,weights=[-1])
    with pytest.raises(ValueError):
        build_cosmological_spin_state([[1.,0,0,0]],1,scale_factor=0,
            hubble_rate=0,energy_density=0,pressure=0,
            units_id="u",derivation_id="d")
