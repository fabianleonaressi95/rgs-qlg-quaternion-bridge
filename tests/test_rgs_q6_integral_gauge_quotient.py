import sympy as sp

from src.rgs_q6_integral_gauge_quotient import (
    build_q6_integral_gauge_quotient_splitting,
)


def test_q6_integral_gauge_quotient_has_unimodular_five_mode_complement():
    out = build_q6_integral_gauge_quotient_splitting()

    assert out["H_cocharacter_basis"].shape == (192, 68)
    assert out["gauge_coordinates_in_H"].shape == (68, 64)
    assert out["gauge_adapted_edge_basis"].shape == (192, 68)
    assert out["gauge_edge_basis"].shape == (192, 63)
    assert out["quotient_cocharacters"].shape == (192, 5)
    assert out["quotient_characters"].shape == (192, 5)

    assert out["annihilator_lattice"] * out["gauge_adapted_edge_basis"] == sp.zeros(124, 68)
    assert out["H_cocharacter_basis"] * out["gauge_coordinates_in_H"] == out["incidence"]
    assert abs(int(out["quotient_pairing"].det())) == 1
    assert all(out["checks"].values())
