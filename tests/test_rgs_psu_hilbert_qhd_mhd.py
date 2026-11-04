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


def test_fluxon_winding_is_quantized_and_gauge_invariant():
    from src.rgs_psu_hilbert_qhd_mhd import (
        fluxon_topological_charge,
        u1_gauge_transform,
    )
    n = 64
    theta = 2.0 * np.pi * np.arange(n) / n
    links = np.zeros(n)
    report = fluxon_topological_charge(theta, links)
    assert report["quantized"]
    assert report["nearest_integer_fluxon_number"] == 1

    rng = np.random.default_rng(126528)
    chi = rng.normal(size=n)
    psi = np.exp(1j * theta)
    psi_g, links_g = u1_gauge_transform(psi, np.exp(1j * links), chi)
    transformed = fluxon_topological_charge(
        np.angle(psi_g), np.angle(links_g)
    )
    assert transformed["quantized"]
    assert transformed["nearest_integer_fluxon_number"] == 1
    assert np.allclose(
        report["gauge_invariant_phase_drops"],
        transformed["gauge_invariant_phase_drops"],
        atol=1e-12,
    )


def test_sine_gordon_fluxon_has_one_flux_quantum_and_small_equation_residual():
    from src.rgs_psu_hilbert_qhd_mhd import (
        sine_gordon_fluxon_profile,
        sine_gordon_static_residual,
    )
    x = np.linspace(-20.0, 20.0, 4001)
    lam = 1.3
    kink = sine_gordon_fluxon_profile(x, center=0.2, josephson_length=lam)
    assert abs(kink["phase_jump"] - 2.0 * np.pi) < 1e-8
    # Integrating B/Phi0 = (1/2pi) dphi/dx gives one flux quantum.
    flux_quanta = np.sum(0.5 * (kink["normalized_magnetic_field"][:-1]\n                               + kink["normalized_magnetic_field"][1:]) * np.diff(x))
    assert abs(flux_quanta - 1.0) < 2e-5
    residual = sine_gordon_static_residual(kink["phase"], x[1] - x[0], lam)
    assert np.max(np.abs(residual)) < 2e-5


def test_fluxon_charge_is_not_declared_quantized_with_nontrivial_background_holonomy():
    from src.rgs_psu_hilbert_qhd_mhd import fluxon_topological_charge
    theta = np.linspace(0.0, 2.0 * np.pi, 32, endpoint=False)
    links = np.full(32, 0.03)
    report = fluxon_topological_charge(theta, links)
    assert not report["quantized"]
