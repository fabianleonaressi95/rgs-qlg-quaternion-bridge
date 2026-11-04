"""Compact U(1) Wilson theory on the periodic Euclidean-time product Q6 x C_Nt.

Q6 is the six-dimensional cubical hypercube graph. Adding one periodic
coordinate therefore gives a 6+1-dimensional Euclidean lattice, NOT a
physical 3+1-dimensional spacetime. This module is an exact gauge-invariant
temporal extension and a laboratory for the lattice Gauss constraint.
"""
from __future__ import annotations
import math
import numpy as np
from src.rgs_q6_equivariant_cochains import build_q6_cochain_complex


def _wrap_angle(x):
    return (np.asarray(x, dtype=float) + math.pi) % (2.0 * math.pi) - math.pi


def build_q6_euclidean_time_data(nt, complex_data=None):
    if not isinstance(nt, (int, np.integer)) or nt < 2:
        raise ValueError("nt must be an integer >= 2")
    c = build_q6_cochain_complex() if complex_data is None else complex_data
    return {"nt": int(nt), "complex": c, "D0": np.asarray(c["D0"], float),
            "D1": np.asarray(c["D1"], float), "n_spatial_links": 192,
            "n_temporal_links": 64, "n_spatial_plaquettes": 240,
            "n_temporal_plaquettes": 192}


def validate_q6_euclidean_links(spatial_links, temporal_links, nt):
    s, t = np.asarray(spatial_links, float), np.asarray(temporal_links, float)
    if s.shape != (nt, 192) or t.shape != (nt, 64):
        raise ValueError("expected spatial (nt,192) and temporal (nt,64) link arrays")
    if not np.all(np.isfinite(s)) or not np.all(np.isfinite(t)):
        raise ValueError("link angles must be finite")
    return s, t


def q6_euclidean_gauge_transform(spatial_links, temporal_links, gauge_angles, data=None):
    nt = np.asarray(spatial_links).shape[0]
    d = build_q6_euclidean_time_data(nt) if data is None else data
    s, t = validate_q6_euclidean_links(spatial_links, temporal_links, d["nt"])
    lam = np.asarray(gauge_angles, float)
    if lam.shape != (d["nt"], 64) or not np.all(np.isfinite(lam)):
        raise ValueError("gauge_angles must be finite with shape (nt,64)")
    moved_s = s + lam @ d["D0"].T
    moved_t = t + np.roll(lam, -1, axis=0) - lam
    return _wrap_angle(moved_s), _wrap_angle(moved_t)


def q6_euclidean_plaquette_fluxes(spatial_links, temporal_links, data=None):
    nt = np.asarray(spatial_links).shape[0]
    d = build_q6_euclidean_time_data(nt) if data is None else data
    s, t = validate_q6_euclidean_links(spatial_links, temporal_links, d["nt"])
    spatial_flux = s @ d["D1"].T
    temporal_flux = np.empty((nt, 192), dtype=float)
    edges = d["complex"]["edges"]
    for k in range(nt):
        kp = (k + 1) % nt
        for e, (a, b, _axis) in enumerate(edges):
            temporal_flux[k, e] = s[k, e] + t[k, b] - s[kp, e] - t[k, a]
    return {"spatial": spatial_flux, "temporal": temporal_flux}


def q6_euclidean_wilson_action(spatial_links, temporal_links,
                                beta_space=1.0, beta_time=1.0, data=None):
    if (not np.isfinite(beta_space) or beta_space <= 0
            or not np.isfinite(beta_time) or beta_time <= 0):
        raise ValueError("both Wilson couplings must be finite and positive")
    flux = q6_euclidean_plaquette_fluxes(spatial_links, temporal_links, data)
    return float(beta_space * np.sum(1.0 - np.cos(flux["spatial"]))
                 + beta_time * np.sum(1.0 - np.cos(flux["temporal"])))


def q6_lattice_gauss_law(spatial_links, temporal_links, data=None):
    """Euler-Lagrange derivative with respect to A_0(v,t): D0^T sin(F_0i)=0."""
    flux = q6_euclidean_plaquette_fluxes(spatial_links, temporal_links, data)
    d = build_q6_euclidean_time_data(flux["temporal"].shape[0]) if data is None else data
    return np.sin(flux["temporal"]) @ d["D0"]


def audit_q6_euclidean_time_wilson(nt=3, seed=20261104):
    d = build_q6_euclidean_time_data(nt)
    rng = np.random.default_rng(seed)
    s = rng.uniform(-0.5, 0.5, (nt, 192))
    t = rng.uniform(-0.5, 0.5, (nt, 64))
    lam = rng.uniform(-math.pi, math.pi, (nt, 64))
    bs, bt = 1.7, 0.9
    moved_s, moved_t = q6_euclidean_gauge_transform(s, t, lam, d)
    action = q6_euclidean_wilson_action(s, t, bs, bt, d)
    moved_action = q6_euclidean_wilson_action(moved_s, moved_t, bs, bt, d)
    flux = q6_euclidean_plaquette_fluxes(s, t, d)
    moved_flux = q6_euclidean_plaquette_fluxes(moved_s, moved_t, d)
    spatial_gauge = float(np.max(np.abs(
        np.exp(1j * flux["spatial"]) - np.exp(1j * moved_flux["spatial"]))))
    temporal_gauge = float(np.max(np.abs(
        np.exp(1j * flux["temporal"]) - np.exp(1j * moved_flux["temporal"]))))
    gauss = q6_lattice_gauss_law(s, t, d)
    # A0 is a Lagrange multiplier: dS/dA0(v,t) = beta_time * Gauss(v,t).
    eps = 1e-6
    ti, vi = 1 % nt, 7
    tp, tm = t.copy(), t.copy()
    tp[ti, vi] += eps
    tm[ti, vi] -= eps
    fd = (q6_euclidean_wilson_action(s, tp, bs, bt, d)
          - q6_euclidean_wilson_action(s, tm, bs, bt, d)) / (2.0 * eps)
    analytic = bt * gauss[ti, vi]
    derivative_res = abs(fd - analytic)
    checks = {
        "cell_counts": d["n_spatial_links"] == 192 and d["n_temporal_links"] == 64
            and d["n_spatial_plaquettes"] == 240 and d["n_temporal_plaquettes"] == 192,
        "D1D0_exact": bool(np.max(np.abs(d["D1"] @ d["D0"])) == 0),
        "spatial_flux_gauge_invariant": spatial_gauge < 1e-10,
        "temporal_flux_gauge_invariant": temporal_gauge < 1e-10,
        "wilson_action_gauge_invariant": abs(action - moved_action) < 1e-10,
        "gauss_law_is_A0_variation": derivative_res < 2e-7,
        "gauss_law_finite": bool(np.all(np.isfinite(gauss))),
    }
    return {
        "audit": "Q6_X_PERIODIC_TIME_U1_WILSON_GAUSS",
        "spacetime_interpretation": "6+1-dimensional Euclidean product lattice; not 3+1",
        "nt": nt, "spatial_links": nt * 192, "temporal_links": nt * 64,
        "spatial_plaquettes": nt * 240, "temporal_plaquettes": nt * 192,
        "wilson_action": "beta_s sum(1-cos F_ij) + beta_t sum(1-cos F_0i)",
        "gauge_action_residual": abs(action - moved_action),
        "spatial_flux_gauge_residual": spatial_gauge,
        "temporal_flux_gauge_residual": temporal_gauge,
        "gauss_law_shape": list(gauss.shape),
        "gauss_law_max_abs": float(np.max(np.abs(gauss))),
        "A0_variation_finite_difference": fd,
        "beta_time_times_gauss": analytic,
        "A0_variation_residual": derivative_res,
        "checks": checks, "all_checks_passed": bool(all(checks.values())),
        "interpretation": (
            "The temporal link is a compact lattice Lagrange multiplier; "
            "its Euler-Lagrange equation is the nonlinear lattice Gauss law. "
            "This is a finite 6+1D Euclidean gauge model, not yet the desired "
            "3+1D PSU-anti-twin dynamics."
        ),
        "open": [
            "derive a physical three-dimensional spatial subcomplex or justified dimensional reduction",
            "derive beta_s and beta_t and their anisotropic continuum scaling from the RGS-PSU action",
            "construct the Hamiltonian transfer matrix and identify the electric momentum spectrum",
            "couple the reduced five-mode sector to the PSU anti-twin vierbein",
        ],
    }
