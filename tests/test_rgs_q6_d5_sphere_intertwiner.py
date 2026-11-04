import numpy as np

from src.rgs_q6_d5_sphere_intertwiner import (
    build_q6_d5_relative_direction_intertwiner,
    relative_direction_from_cotangent,
)


def test_exterior_square_contains_geometric_three_dimensional_intertwiner():
    report = build_q6_d5_relative_direction_intertwiner()
    assert report["source_irrep_dimension"] == 5
    assert report["exterior_square_dimension"] == 10
    assert report["target_geometric_irrep_dimension"] == 3
    assert report["intertwiner_nullity"] == 1
    assert report["intertwiner_rank"] == 3
    assert report["all_checks_passed"]


def test_cotangent_wedge_direction_is_A5_equivariant():
    report = build_q6_d5_relative_direction_intertwiner()
    rng = np.random.default_rng(20261104)
    q = rng.normal(size=5)
    p = rng.normal(size=5)
    base = relative_direction_from_cotangent(q, p, report["intertwiner"])
    assert base["direction"] is not None
    for index in (1, 7, 19, 43, 59):
        rho = report["A5_mode_representation"][index]
        rotation = report["geometric_rotations"][index]
        moved = relative_direction_from_cotangent(
            rho @ q, rho @ p, report["intertwiner"]
        )
        assert np.allclose(moved["vector"], rotation @ base["vector"], atol=2e-8)
        assert np.allclose(moved["direction"], rotation @ base["direction"], atol=2e-8)


def test_cotangent_direction_is_invariant_under_oscillator_SO2_phase():
    report = build_q6_d5_relative_direction_intertwiner()
    rng = np.random.default_rng(37)
    q, p = rng.normal(size=5), rng.normal(size=5)
    angle = 0.731
    q2 = np.cos(angle)*q - np.sin(angle)*p
    p2 = np.sin(angle)*q + np.cos(angle)*p
    a = relative_direction_from_cotangent(q, p, report["intertwiner"])
    b = relative_direction_from_cotangent(q2, p2, report["intertwiner"])
    assert np.allclose(a["wedge_coordinates"], b["wedge_coordinates"], atol=1e-12)
    assert np.allclose(a["vector"], b["vector"], atol=1e-10)
    assert np.allclose(a["direction"], b["direction"], atol=1e-10)


def test_collinear_cotangent_pair_has_zero_projected_direction_or_is_rejected_cleanly():
    report = build_q6_d5_relative_direction_intertwiner()
    q = np.arange(1., 6.)
    p = 2.0*q
    result = relative_direction_from_cotangent(q, p, report["intertwiner"])
    assert result["norm"] < 1e-12
    assert result["direction"] is None
