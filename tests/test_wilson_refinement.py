import numpy as np
import pytest

from src.wilson_refinement import (
    naive_symbol_energy, refinement_audit, wilson_brillouin_corners,
    wilson_symbol,
)


def test_fixed_physical_momentum_has_second_order_naive_dispersion_error():
    audit=refinement_audit((16,32,64,128),dimension=3,mode=(1,0,0),r_exponent=0.0)
    rows=audit["refinements"]
    assert rows[-1]["fixed_mode_naive_relative_error"] < rows[0]["fixed_mode_naive_relative_error"]
    assert rows[-1]["fixed_mode_wilson_correction"] < rows[0]["fixed_mode_wilson_correction"]
    assert rows[-1]["fixed_mode_wilson_correction"] < 0.1
    assert audit["simultaneous_scaling_window"]


def test_three_dimensional_wilson_corners_lift_all_seven_doublers():
    a=0.125; r=1.0
    corners=wilson_brillouin_corners(3,a,r)
    assert len(corners)==8
    assert sum(c["pi_directions"]==0 for c in corners)==1
    assert sum(c["pi_directions"]>0 for c in corners)==7
    for c in corners:
        assert c["naive_energy"] < 1e-12
        assert np.isclose(c["wilson_mass"],2*r*c["pi_directions"]/a,atol=1e-11)


def test_wilson_parameter_scaling_window_is_minus_one_to_one():
    good=refinement_audit((16,32),dimension=3,r_exponent=0.5)
    bad_low=refinement_audit((16,32),dimension=3,r_exponent=-1.0)
    bad_doubler=refinement_audit((16,32),dimension=3,r_exponent=1.0)
    assert good["simultaneous_scaling_window"]
    assert not bad_low["low_momentum_wilson_term_vanishes"]
    assert not bad_doubler["doubler_mass_diverges"]


def test_wilson_symbol_has_expected_small_momentum_and_corner_limits():
    p=np.array([0.2,0.0,0.0]); a=1e-4; r=1.0
    w=wilson_symbol(p,a,r)
    assert np.isclose(w,0.5*r*a*np.dot(p,p),rtol=1e-6)
    assert np.isclose(naive_symbol_energy([np.pi/a],a),0.0,atol=1e-2)


def test_invalid_refinement_inputs_are_rejected():
    with pytest.raises(ValueError):
        refinement_audit((8,16),dimension=3,mode=(1,0))
    with pytest.raises(ValueError):
        wilson_symbol([1.0],0.0,1.0)
