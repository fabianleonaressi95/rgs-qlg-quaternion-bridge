"""Wilson transfer operator in the integer electric-flux representation on Q6.

The exact temporal kernel is diagonal in divergence-free integer electric
fluxes. A finite list of flux states gives a controlled Galerkin truncation
of that kernel; it is not the full infinite-dimensional transfer operator.
The spatial Wilson potential is kept as a multiplication operator in link
angle representation.

Q6 has six spatial dimensions. Continuum matching is therefore parameterized
by spatial_dim=6 and a dimensionful g^2, and must not be confused with 3+1D.
"""
from __future__ import annotations

import math
import numpy as np
from scipy.special import ive

from src.rgs_q6_euclidean_time_wilson import build_q6_euclidean_time_data
from src.rgs_q6_u1_character_electric_flux import (
    build_q6_integer_cycle_basis,
    q6_u1_integer_gauss_law,
)


def _positive_finite(name, value):
    x = float(value)
    if not math.isfinite(x) or x <= 0.0:
        raise ValueError(f"{name} must be finite and positive")
    return x


def validate_q6_integer_flux_states(flux_states, data=None):
    d = build_q6_euclidean_time_data(2) if data is None else data
    states = np.asarray(flux_states)
    if states.ndim != 2 or states.shape[1] != 192 or states.shape[0] < 1:
        raise ValueError("flux_states must have shape (n_states,192)")
    if not np.all(np.isfinite(states)) or not np.all(states == np.rint(states)):
        raise ValueError("flux_states must be finite integer vectors")
    states = states.astype(np.int64)
    div = np.asarray(d["D0"], dtype=np.int64).T @ states.T
    if np.any(div):
        raise ValueError("every flux state must satisfy the exact integer Gauss law")
    return states


def q6_u1_electric_transfer_eigenvalues(flux_states, beta_time, data=None):
    """Normalized temporal transfer eigenvalues prod_e I_|n_e|(beta)/I_0(beta)."""
    beta = _positive_finite("beta_time", beta_time)
    states = validate_q6_integer_flux_states(flux_states, data)
    orders = np.abs(states)
    # The common exp(-beta) scaling cancels in I_n/I_0.
    ratios = ive(orders, beta) / ive(0, beta)
    if np.any(ratios <= 0.0) or not np.all(np.isfinite(ratios)):
        raise FloatingPointError("scaled Bessel ratios are not finite and positive")
    return np.exp(np.sum(np.log(ratios), axis=1))


def q6_u1_electric_hamiltonian_eigenvalues(
    flux_states, beta_time, temporal_spacing, data=None
):
    """Exact finite-beta energies -log(lambda_n)/a_t, with vacuum energy removed."""
    a_t = _positive_finite("temporal_spacing", temporal_spacing)
    lam = q6_u1_electric_transfer_eigenvalues(flux_states, beta_time, data)
    return -np.log(lam) / a_t


def q6_u1_spatial_wilson_potential(spatial_links, beta_space, data=None):
    """V(A)=beta_s sum_f (1-cos((D1 A)_f)); gauge invariant on Q6."""
    beta = _positive_finite("beta_space", beta_space)
    d = build_q6_euclidean_time_data(2) if data is None else data
    a = np.asarray(spatial_links, dtype=float)
    if a.shape != (192,) or not np.all(np.isfinite(a)):
        raise ValueError("spatial_links must be a finite vector of shape (192,)")
    flux = np.asarray(d["D1"], dtype=float) @ a
    return float(beta * np.sum(1.0 - np.cos(flux)))


def q6_u1_symmetric_transfer_kernel(
    spatial_links_next,
    spatial_links_previous,
    flux_states,
    beta_time,
    beta_space,
    data=None,
):
    """Finite-character symmetric kernel <A'|T|A>.

    T = exp(-V/2) K_t exp(-V/2), with K_t represented by the supplied
    divergence-free integer flux states. The exact kernel is recovered only
    when the complete infinite cycle lattice is summed.
    """
    d = build_q6_euclidean_time_data(2) if data is None else data
    ap = np.asarray(spatial_links_next, dtype=float)
    a = np.asarray(spatial_links_previous, dtype=float)
    if ap.shape != (192,) or a.shape != (192,):
        raise ValueError("both link configurations must have shape (192,)")
    if not np.all(np.isfinite(ap)) or not np.all(np.isfinite(a)):
        raise ValueError("link configurations must be finite")
    states = validate_q6_integer_flux_states(flux_states, d)
    lam = q6_u1_electric_transfer_eigenvalues(states, beta_time, d)
    phase = np.exp(1j * (states @ (ap - a)))
    temporal = np.dot(lam, phase)
    vp = q6_u1_spatial_wilson_potential(ap, beta_space, d)
    va = q6_u1_spatial_wilson_potential(a, beta_space, d)
    return complex(np.exp(-0.5 * (vp + va)) * temporal)


def q6_u1_anisotropic_wilson_couplings(
    spatial_spacing, temporal_spacing, g_squared, spatial_dim=6
):
    """Continuum-matching Wilson couplings for d spatial dimensions.

    For S=(1/4g^2) integral F_mu_nu F_mu_nu, the plaquette coefficients are
      beta_t = a_s^(d-2)/(g^2 a_t),
      beta_s = a_t a_s^(d-4)/g^2.
    In d=3 this reduces to the familiar beta_t=a_s/(g^2 a_t),
    beta_s=a_t/(g^2 a_s). For Q6, d=6 and [g^2]=length^3.
    """
    a_s = _positive_finite("spatial_spacing", spatial_spacing)
    a_t = _positive_finite("temporal_spacing", temporal_spacing)
    g2 = _positive_finite("g_squared", g_squared)
    if not isinstance(spatial_dim, (int, np.integer)) or spatial_dim < 2:
        raise ValueError("spatial_dim must be an integer >= 2")
    beta_t = a_s ** (spatial_dim - 2) / (g2 * a_t)
    beta_s = a_t * a_s ** (spatial_dim - 4) / g2
    return {
        "spatial_dim": int(spatial_dim),
        "spatial_spacing": a_s,
        "temporal_spacing": a_t,
        "g_squared": g2,
        "beta_time": float(beta_t),
        "beta_space": float(beta_s),
        "electric_hamiltonian_coefficient": float(1.0 / (2.0 * beta_t * a_t)),
        "magnetic_potential_coefficient": float(beta_s / a_t),
        "electric_continuum_coefficient": float(g2 / (2.0 * a_s ** (spatial_dim - 2))),
        "magnetic_continuum_coefficient": float(a_s ** (spatial_dim - 4) / g2),
        "convention": "S=(1/(4 g^2)) integral d^(d+1)x F_mu_nu F_mu_nu",
    }


def audit_q6_u1_transfer_operator(beta_time=80.0, beta_space=0.7, a_s=0.2, a_t=0.05,
                                  g_squared=0.003, seed=20261104):
    d = build_q6_euclidean_time_data(2)
    cycles = build_q6_integer_cycle_basis(d)
    # A finite symmetric Galerkin set: vacuum plus +/- fundamental cycles.
    states = np.vstack([np.zeros((1, 192), dtype=np.int64), cycles, -cycles])
    rng = np.random.default_rng(seed)
    ap = rng.uniform(-0.3, 0.3, 192)
    a = rng.uniform(-0.3, 0.3, 192)
    lam = q6_u1_electric_transfer_eigenvalues(states, beta_time, d)
    kernel_ab = q6_u1_symmetric_transfer_kernel(
        ap, a, states, beta_time, beta_space, d
    )
    kernel_ba = q6_u1_symmetric_transfer_kernel(
        a, ap, states, beta_time, beta_space, d
    )
    couplings = q6_u1_anisotropic_wilson_couplings(
        a_s, a_t, g_squared, spatial_dim=6
    )
    one_cycle = cycles[:1]
    exact_energy = q6_u1_electric_hamiltonian_eigenvalues(
        one_cycle, beta_time, a_t, d
    )[0]
    leading_energy = float(np.sum(one_cycle[0] ** 2) / (2.0 * beta_time * a_t))
    checks = {
        "all_states_integer_gauss": bool(np.max(np.abs(
            q6_u1_integer_gauss_law(states[1], d)
        )) == 0),
        "transfer_weights_positive": bool(np.all(lam > 0.0)),
        "vacuum_transfer_eigenvalue_one": math.isclose(float(lam[0]), 1.0, rel_tol=1e-13),
        "transfer_eigenvalues_at_most_one": bool(np.max(lam) <= 1.0 + 1e-12),
        "symmetric_kernel_hermitian": abs(kernel_ab - np.conjugate(kernel_ba)) < 1e-11,
        "anisotropic_scaling_positive": (
            couplings["beta_time"] > 0.0 and couplings["beta_space"] > 0.0
        ),
        "electric_weak_coupling_asymptotic": abs(exact_energy - leading_energy)
            / max(leading_energy, 1e-15) < 0.03,
        "d6_scaling_exponents": (
            math.isclose(couplings["beta_time"],
                         a_s ** 4 / (g_squared * a_t), rel_tol=1e-14)
            and math.isclose(couplings["beta_space"],
                             a_t * a_s ** 2 / g_squared, rel_tol=1e-14)
        ),
    }
    return {
        "audit": "Q6_U1_INTEGER_FLUX_TRANSFER_OPERATOR",
        "n_flux_states": int(states.shape[0]),
        "flux_state_dimension": 192,
        "integer_cycle_rank": 129,
        "beta_time": float(beta_time),
        "beta_space": float(beta_space),
        "electric_energy_exact": float(exact_energy),
        "electric_energy_leading_n2": leading_energy,
        "electric_asymptotic_relative_error": float(
            abs(exact_energy - leading_energy) / max(leading_energy, 1e-15)
        ),
        "kernel_hermiticity_residual": float(
            abs(kernel_ab - np.conjugate(kernel_ba))
        ),
        "anisotropic_couplings": couplings,
        "checks": checks,
        "all_checks_passed": bool(all(checks.values())),
        "scope": (
            "Finite Galerkin evaluation of the exact character-diagonal temporal "
            "kernel and symmetric Wilson transfer factorization on the 6D spatial Q6. "
            "The supplied flux list is truncated; no physical 3+1D reduction is claimed."
        ),
        "open": [
            "control convergence as the integer cycle-lattice basis is enlarged",
            "derive the RGS-PSU microscopic action that fixes g_squared",
            "construct a physically justified 3D spatial subcomplex before 3+1D claims",
            "analyze the continuum limit and renormalization of the 6+1D gauge theory",
        ],
    }
