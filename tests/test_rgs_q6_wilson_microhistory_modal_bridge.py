from fractions import Fraction

import pytest

from src.rgs_q6_wilson_microhistory_modal_bridge import (
    OrbitKripkeFrame,
    audit_q6_wilson_microhistory_modal_bridge,
    enumerate_orbit,
    flux_l1_order,
    labelled_microhistory,
    microhistory_orbit_stabilizer,
    symmetrized_microhistory_weight,
    wilson_factorial_majorant,
)


def test_exact_orbit_stabilizer_and_factorial_weight_identity():
    flux = (3, -2, 0, 1)
    data = microhistory_orbit_stabilizer(flux, edge_count=4)
    assert data["order"] == 6
    assert data["orbit_order"] == 60
    assert data["stabilizer_order"] == 12
    assert data["orbit_order"] * data["stabilizer_order"] == data["s_factorial"]
    beta = Fraction(7, 3)
    assert symmetrized_microhistory_weight(flux, beta, 4) == \
        wilson_factorial_majorant(flux, beta, 4)


def test_signed_tokens_keep_orientation_in_the_quotient():
    history = labelled_microhistory((2, -1, 0), edge_count=3)
    assert history == (1, 1, -2)
    assert len(enumerate_orbit(history)) == 3


def test_orbit_frame_is_s5_and_box_means_invariance_on_orbit():
    history = (1, 1, -2)
    frame = OrbitKripkeFrame.from_history(history)
    assert frame.is_equivalence_frame()
    invariant = lambda w: sorted(w) == sorted(history)
    noninvariant = lambda w: w[0] == 1
    assert frame.satisfies_box(history, invariant)
    assert not frame.satisfies_box(history, noninvariant)


def test_audit_and_validation():
    report = audit_q6_wilson_microhistory_modal_bridge()
    assert report["checks_passed"], report
    assert flux_l1_order((1, -2, 0), 3) == 3
    with pytest.raises(ValueError):
        flux_l1_order((1, 2), 3)
    with pytest.raises(TypeError):
        flux_l1_order((1, 2.0, 0), 3)
