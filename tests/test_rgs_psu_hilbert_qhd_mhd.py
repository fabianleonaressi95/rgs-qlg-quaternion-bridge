import numpy as np

from src.rgs_psu_hilbert_qhd_mhd import (
    audit_hilbert_cantor_qhd_gauge_mhd,
    cantor_cylinder_mass,
    hilbert_d2xy,
    monopole_chern_number,
    q6_hilbert_cantor_table,
    u1_covariant_difference,
    u1_gauge_transform,
)


def test_q6_hilbert_curve_is_bijective():
    pts = [hilbert_d2xy(3, d) for d in range(64)]
    assert len(set(pts)) == 64
    assert all(0 <= x < 8 and 0 <= y < 8 for x, y in pts)
    assert q6_hilbert_cantor_table()["hilbert_bijection"]


def test_cantor_cylinder_measure_is_normalized_at_each_level():
    for n in range(1, 8):
        masses = [
            cantor_cylinder_mass([(k >> j) & 1 for j in range(n)], (0.37, 0.63))
            for k in range(2**n)
        ]
        assert np.isclose(sum(masses), 1.0)


def test_u1_covariant_difference_is_gauge_covariant():
    rng = np.random.default_rng(1204)
    psi = rng.normal(size=16) + 1j * rng.normal(size=16)
    links = np.exp(1j * rng.normal(size=16))
    phases = rng.normal(size=16)
    d0 = u1_covariant_difference(psi, links)
    psi_g, links_g = u1_gauge_transform(psi, links, phases)
    d1 = u1_covariant_difference(psi_g, links_g)
    assert np.allclose(d1, np.exp(1j * phases) * d0, atol=1e-12)


def test_monopole_chern_flux_quantization():
    report = monopole_chern_number(np.full(8, 2.0 * np.pi / 4.0))
    assert report["quantized_within_tolerance"]
    assert report["nearest_integer_charge"] == 2


def test_integrated_periodic_benchmark_residuals():
    report = audit_hilbert_cantor_qhd_gauge_mhd()
    assert report["u1_gauge"]["covariance_residual"] < 1e-11
    assert report["monopole"]["quantized_within_tolerance"]
    assert report["qhd"]["continuity_mass_residual"] < 1e-10
    assert report["mhd"]["induction_mass_residual"] < 1e-10
