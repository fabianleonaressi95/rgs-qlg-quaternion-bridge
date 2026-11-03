import numpy as np
import pytest
from src.thermodynamics import compare_boundary_spectra, fermion_partition, modular_partition

def test_fermion_partition_closed_form_and_thermodynamic_identity():
    t=.7; e=np.array([-.4,.4]); r=fermion_partition(e,t); beta=1/t
    expected=2*np.log1p(np.exp(-.4*beta))-.4*beta
    assert np.isclose(r["logZ_vacuum_normalized"],expected)
    assert np.isclose(r["entropy_thermodynamic"],r["logZ_grand_canonical_raw"]+beta*r["energy_raw"])
    assert r["heat_capacity"]>0

def test_repeated_modes_equal_explicit_degeneracy():
    a=fermion_partition([-.2,-.2,.2,.2],.5)
    b=fermion_partition([-.2,.2],.5,[2.,2.])
    for k in ("logZ_vacuum_normalized","energy_above_filled_vacuum","entropy_thermodynamic","heat_capacity"):
        assert np.isclose(a[k],b[k],atol=1e-13)

def test_modular_endpoints_and_roundoff_are_handled():
    r=modular_partition([-2e-16,.25,.75,1+2e-16])
    expected=-2*(.25*np.log(.25)+.75*np.log(.75))
    assert r["projected_zero_modes"]==1 and r["projected_unit_modes"]==1
    assert np.isclose(r["entropy_entanglement"],expected)
    assert r["logZ_modular"] is None and r["modular_energies"][0] is None

def test_modular_rejects_non_roundoff_out_of_range():
    with pytest.raises(ValueError): modular_partition([-1e-4,.5])

def test_reference_json_reproduces_entropies_and_partition_functions():
    import json
    from pathlib import Path
    data=json.loads(Path("results/json/mobius_hawking_spectrum_reference.json").read_text())
    r=compare_boundary_spectra(data)
    assert r["naive"]["thermodynamics"]["mode_count"]==24
    assert r["wilson"]["thermodynamics"]["mode_count"]==24
    assert r["wilson"]["modular"]["projected_zero_modes"]==8
    assert r["entanglement_entropy_consistency"]["naive"]["absolute_residual"]<1e-10
    assert r["entanglement_entropy_consistency"]["wilson"]["absolute_residual"]<1e-8
    assert np.isclose(r["naive"]["thermodynamics"]["Z_vacuum_normalized"],1.1072369002993039,rtol=1e-11)
    assert np.isclose(r["wilson"]["thermodynamics"]["Z_vacuum_normalized"],1.2426750246234262,rtol=1e-11)

def test_invalid_fermion_inputs_are_rejected():
    with pytest.raises(ValueError): fermion_partition([],1.)
    with pytest.raises(ValueError): fermion_partition([1.],0.)
    with pytest.raises(ValueError): fermion_partition([1.],1.,[0.])
