"""Coefficient-level matching of the Q6 Wilson Hamiltonian to the D5 Maxwell sector.

The matching is exact for the real five-dimensional coexact spectral
subspace after the fixed-flux, large-beta_t Hamiltonian limit. It does not
identify that real subspace with an integer sublattice of Z^192.
"""
from __future__ import annotations

import math
import numpy as np

from src.rgs_q6_d5_phase_space import build_q6_d5_phase_space_sector
from src.rgs_q6_u1_transfer_operator import q6_u1_anisotropic_wilson_couplings
from src.rgs_q6_d5_maxwell_hamiltonian import maxwell_mode_parameters


def q6_u1_maxwell_matching_report(
    spatial_spacing=1.0,
    temporal_spacing=0.05,
    g_squared=1.0,
    spatial_dim=6,
    hodge_eigenvalue=12.0,
):
    """Derive and audit the projected five-mode Maxwell coefficient identity.

    With beta_eff=a_s^(d-2)/g^2 and lambda_eff=lambda/a_s^2,
      H_5 = 1/(2 beta_eff) |p|^2
            + beta_eff lambda_eff/2 |q|^2.
    This is precisely the Maxwell oscillator convention used by the D5
    five-mode Hamiltonian. The equality is a projected continuum-Hamiltonian
    statement, not an exact equality of compact integer-flux Hilbert spaces.
    """
    vals = (spatial_spacing, temporal_spacing, g_squared, hodge_eigenvalue)
    if not all(math.isfinite(float(x)) and float(x) > 0 for x in vals):
        raise ValueError("all scales and the Hodge eigenvalue must be finite and positive")
    if not isinstance(spatial_dim, int) or spatial_dim < 2:
        raise ValueError("spatial_dim must be an integer >= 2")

    a_s, a_t, g2, lam = map(float, vals)
    sector = build_q6_d5_phase_space_sector(beta=1.0, eigenvalue=lam)
    basis = np.asarray(sector["selected_mode_basis"], dtype=float)
    complex_data = __import__(
        "src.rgs_q6_equivariant_cochains", fromlist=["build_q6_cochain_complex"]
    ).build_q6_cochain_complex()
    d0 = np.asarray(complex_data["D0"], dtype=float)
    d1 = np.asarray(complex_data["D1"], dtype=float)
    hodge = d1.T @ d1

    couplings = q6_u1_anisotropic_wilson_couplings(
        a_s, a_t, g2, spatial_dim=spatial_dim
    )
    beta_eff = a_s ** (spatial_dim - 2) / g2
    lambda_eff = lam / (a_s ** 2)
    maxwell = maxwell_mode_parameters(beta=beta_eff, eigenvalue=lambda_eff)

    # Transfer-Hamiltonian quadratic coefficients on the five normalized modes.
    electric_from_transfer = g2 / (a_s ** (spatial_dim - 2))
    magnetic_from_transfer = (a_s ** (spatial_dim - 4) / g2) * lam
    electric_from_d5 = maxwell["electric_hessian"]
    magnetic_from_d5 = maxwell["magnetic_hessian"]
    frequency_from_transfer = math.sqrt(
        electric_from_transfer * magnetic_from_transfer
    )
    frequency_from_d5 = maxwell["oscillator_frequency"]

    projected_hodge = basis.T @ hodge @ basis
    checks = {
        "five_mode_basis_shape": basis.shape == (192, 5),
        "basis_orthonormal": np.linalg.norm(basis.T @ basis - np.eye(5)) < 1e-10,
        "basis_satisfies_gauss_tangent": np.linalg.norm(d0.T @ basis) < 1e-10,
        "basis_is_hodge_eigenspace": np.linalg.norm(
            hodge @ basis - lam * basis
        ) < 1e-8,
        "projected_hodge_is_lambda_identity": np.linalg.norm(
            projected_hodge - lam * np.eye(5)
        ) < 1e-8,
        "electric_coefficient_matches_D5": math.isclose(
            electric_from_transfer, electric_from_d5, rel_tol=1e-12
        ),
        "magnetic_coefficient_matches_D5": math.isclose(
            magnetic_from_transfer, magnetic_from_d5, rel_tol=1e-12
        ),
        "frequency_matches": math.isclose(
            frequency_from_transfer, frequency_from_d5, rel_tol=1e-12
        ),
        "temporal_spacing_cancels_from_Hamiltonian_coefficients": (
            math.isclose(couplings["electric_hamiltonian_coefficient"],
                         electric_from_transfer / 2.0, rel_tol=1e-12)
            and math.isclose(couplings["magnetic_potential_coefficient"],
                             a_s ** (spatial_dim - 4) / g2, rel_tol=1e-12)
        ),
    }
    return {
        "audit": "Q6_WILSON_TO_D5_MAXWELL_COEFFICIENT_MATCHING",
        "spatial_dimension": spatial_dim,
        "mode_multiplicity": 5,
        "hodge_eigenvalue_lattice": lam,
        "spatial_spacing": a_s,
        "temporal_spacing": a_t,
        "g_squared": g2,
        "effective_maxwell_beta": beta_eff,
        "effective_maxwell_hodge_eigenvalue": lambda_eff,
        "electric_hessian_transfer": electric_from_transfer,
        "electric_hessian_D5": electric_from_d5,
        "magnetic_hessian_transfer": magnetic_from_transfer,
        "magnetic_hessian_D5": magnetic_from_d5,
        "frequency_transfer": frequency_from_transfer,
        "frequency_D5": frequency_from_d5,
        "projected_hodge_residual": float(np.linalg.norm(
            projected_hodge - lam * np.eye(5)
        )),
        "continuum_couplings": couplings,
        "checks": checks,
        "all_checks_passed": bool(all(checks.values())),
        "proved_within_assumptions": [
            "the normalized five-dimensional coexact eigenspace is a Gauss-tangent subspace",
            "the Wilson quadratic magnetic form restricts to lambda times the Euclidean norm",
            "the fixed-flux large-beta_t electric Hamiltonian coefficient matches the Maxwell Legendre coefficient",
            "the projected five-mode Hamiltonian is exactly the D5 Maxwell oscillator after beta and eigenvalue rescaling",
            "the oscillator frequency is sqrt(lambda)/a_s and is independent of g_squared",
        ],
        "not_proved": [
            "the compact integer cycle lattice contains this five-dimensional real eigenspace as an exact rank-five integer sublattice",
            "the finite-beta_t Bessel electric spectrum is exactly quadratic",
            "the global compact Q6 gauge theory equals the D5 bounded symmetric domain",
            "a physical 3+1D reduction or a value for the fine-structure constant",
            "the microscopic RGS-PSU action fixes g_squared or the electric charge normalization",
        ],
        "scope": (
            "Coefficient-level theorem for the real five-mode spectral projection "
            "of the continuum-matched 6+1D Q6 U(1) Hamiltonian. The passage from "
            "compact integer fluxes to this real Gaussian mode sector uses the "
            "fixed-flux weak-coupling/large-beta_t limit."
        ),
    }
