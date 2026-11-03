import numpy as np
import pytest
from src.mosco_torus import euclidean_wilson_form_symbol, torus_mosco_sampling_audit


def test_flat_torus_wilson_form_converges_for_fixed_physical_mode():
    audit=torus_mosco_sampling_audit((16,32,64,128,256),dimension=3,mode=(1,0,0),r_exponent=0.0)
    rows=audit["refinements"]
    assert rows[-1]["relative_form_error"]<rows[0]["relative_form_error"]
    assert rows[-1]["relative_form_error"]<1e-3
    assert audit["fixed_mode_form_converges"]
    assert audit["doubler_mass_diverges"]


def test_flat_torus_doubler_mass_diverges_for_alpha_inside_window():
    audit=torus_mosco_sampling_audit((16,32,64,128),dimension=3,mode=(1,0,0),r_exponent=0.5)
    rows=audit["refinements"]
    assert rows[-1]["minimum_nonphysical_corner_mass"]>rows[0]["minimum_nonphysical_corner_mass"]


def test_flat_torus_scaling_endpoints_do_not_decouple_doublers_or_low_modes():
    low=torus_mosco_sampling_audit((16,32),dimension=3,r_exponent=-1.0)
    high=torus_mosco_sampling_audit((16,32),dimension=3,r_exponent=1.0)
    assert not low["fixed_mode_form_converges"]
    assert high["fixed_mode_form_converges"]
    assert not high["doubler_mass_diverges"]


def test_euclidean_wilson_form_has_only_the_physical_corner_zero():
    a=0.125
    assert euclidean_wilson_form_symbol([0.0],a,1.0)==0.0
    assert np.isclose(euclidean_wilson_form_symbol([np.pi/a],a,1.0),(2/a)**2)
    with pytest.raises(ValueError):
        torus_mosco_sampling_audit((8,16),dimension=3,mode=(1,0))
