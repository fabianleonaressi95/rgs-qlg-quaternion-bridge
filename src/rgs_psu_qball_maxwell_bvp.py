"""Flat-space stationary gauged Q-ball coupled self-consistently to Maxwell A_0.

Conventions: metric (-,+,+,+), D_mu = partial_mu - i q A_mu,
Phi=f(r) exp(i omega t), A=A_0(r) dt, Omega=omega-q A_0.
The BVP is electrostatic (no magnetic field) and has a finite outer cutoff.
It does not include gravitational/vierbein backreaction or stability analysis.
"""
from __future__ import annotations

import numpy as np
from scipy.integrate import solve_bvp

from src.rgs_psu_qball_bvp import qball_stationary_bvp


def qball_maxwell_stationary_bvp(
    *,
    gauge_charge: float = 0.05,
    central_amplitude: float = 1.0,
    radius_max: float = 20.0,
    points: int = 180,
    mass_scale: float = 1.0,
    quartic: float = 1.0,
    sextic: float = 0.5,
    tolerance: float = 1.0e-7,
    max_nodes: int = 16000,
) -> dict:
    """Solve the coupled radial scalar-Maxwell BVP with omega as eigenparameter.

    Scalar:
      f''+2 f'/r = (m^2-Omega^2) f - 2 lambda f^3 + 3 g f^5
    Maxwell:
      A0''+2 A0'/r = +2 q Omega f^2
    Omega = omega-q A0.

    Conditions:
      f(0)=f0, f'(0)=0, f(R)=0, A0'(0)=0, A0(R)=0.
    The last condition fixes the residual constant-gauge freedom by setting
    the electrostatic potential to zero at the finite outer boundary.
    """
    q = float(gauge_charge)
    m, lam, g = float(mass_scale), float(quartic), float(sextic)
    f0, rmax, tol = float(central_amplitude), float(radius_max), float(tolerance)
    if not np.all(np.isfinite([q, m, lam, g, f0, rmax, tol])):
        raise ValueError("all parameters must be finite")
    if m <= 0.0 or lam <= 0.0 or g <= 0.0:
        raise ValueError("mass_scale, quartic and sextic must be positive")
    if f0 <= 0.0 or rmax <= 0.0 or tol <= 0.0:
        raise ValueError("central_amplitude, radius_max and tolerance must be positive")
    if points < 40 or max_nodes < points:
        raise ValueError("require points>=40 and max_nodes>=points")

    lower_sq = m * m - lam * lam / (4.0 * g)
    if lower_sq <= 0.0:
        raise ValueError("potential has no positive Q-ball frequency lower bound")
    omega_min = float(np.sqrt(lower_sq))

    # Use the already-tested global Q-ball as a continuation seed at q=0.
    seed = qball_stationary_bvp(
        central_amplitude=f0, radius_max=rmax, points=points,
        mass_scale=m, quartic=lam, sextic=g,
        tolerance=tol, max_nodes=max_nodes,
    )
    radius = np.linspace(0.0, rmax, int(points))
    f_seed = np.interp(radius, seed["radius"], seed["amplitude"])
    fp_seed = np.gradient(f_seed, radius, edge_order=2)
    # A smooth, boundary-compatible electrostatic initial iterate. Its scale
    # is O(q); the coupled BVP, not this seed, determines the final potential.
    shape = np.exp(-radius / max(1.0, rmax / 8.0)) * (1.0 - radius / rmax) ** 2
    a_seed = q * f0**2 * shape
    ap_seed = np.gradient(a_seed, radius, edge_order=2)
    initial = np.vstack((f_seed, fp_seed, a_seed, ap_seed))

    def equations(r, y, parameter):
        omega = parameter[0]
        f, fp, a0, ap = y
        local_frequency = omega - q * a0
        fpp = ((m * m - local_frequency**2) * f
               - 2.0 * lam * f**3 + 3.0 * g * f**5)
        app = 2.0 * q * local_frequency * f**2
        return np.vstack((fp, fpp, ap, app))

    def boundary(ya, yb, parameter):
        return np.asarray((
            ya[1], ya[0] - f0, yb[0], ya[3], yb[2]
        ), dtype=float)

    singular = np.zeros((4, 4), dtype=float)
    singular[1, 1] = -2.0
    singular[3, 3] = -2.0
    solution = solve_bvp(
        equations, boundary, radius, initial,
        p=np.asarray((seed["omega"],)), S=singular,
        tol=tol, max_nodes=int(max_nodes), verbose=0,
    )
    if solution.status != 0:
        raise RuntimeError(
            f"gauged Q-ball/Maxwell BVP did not converge: "
            f"status={solution.status}; {solution.message}"
        )

    omega = float(solution.p[0])
    if not omega_min < omega < m:
        raise RuntimeError(
            f"asymptotic frequency {omega:.12g} outside "
            f"necessary window ({omega_min:.12g}, {m:.12g})"
        )

    audit_radius = np.linspace(0.0, rmax, max(4 * int(points), 801))
    values = solution.sol(audit_radius)
    first = solution.sol(audit_radius, 1)
    second = solution.sol(audit_radius, 2)
    f, fp, a0, ap = values
    fpp, app = second[0], second[2]
    local_frequency = omega - q * a0

    scalar_rhs = ((m * m - local_frequency**2) * f
                  - 2.0 * lam * f**3 + 3.0 * g * f**5)
    scalar_laplacian = fpp.copy()
    scalar_laplacian[1:] += 2.0 * fp[1:] / audit_radius[1:]
    scalar_laplacian[0] = 3.0 * fpp[0]
    scalar_residual = scalar_laplacian - scalar_rhs

    maxwell_rhs = -2.0 * q * local_frequency * f**2
    maxwell_laplacian = app.copy()
    maxwell_laplacian[1:] += 2.0 * ap[1:] / audit_radius[1:]
    maxwell_laplacian[0] = 3.0 * app[0]
    maxwell_residual = maxwell_laplacian - maxwell_rhs

    potential = m * m * f**2 - lam * f**4 + g * f**6
    measure = 4.0 * np.pi * audit_radius**2
    scalar_energy_density = local_frequency**2 * f**2 + fp**2 + potential
    em_energy_density = 0.5 * ap**2
    energy_scalar = float(np.trapezoid(measure * scalar_energy_density, audit_radius))
    energy_maxwell = float(np.trapezoid(measure * em_energy_density, audit_radius))
    energy = energy_scalar + energy_maxwell

    # Q_particle is the phase Noether charge; Q_e=q Q_particle is electric charge.
    charge_density = 2.0 * local_frequency * f**2
    particle_charge = float(np.trapezoid(measure * charge_density, audit_radius))
    electric_charge = -q * particle_charge
    energy_per_particle = energy / particle_charge if particle_charge > 0.0 else float("inf")

    gauss_flux_charge = float(-4.0 * np.pi * rmax**2 * solution.sol(rmax)[3])
    gauss_residual = abs(gauss_flux_charge - electric_charge)
    return {
        "schema": "RGS_PSU_GAUGED_U1_QBALL_MAXWELL_BVP_V1",
        "converged": True,
        "solver_status": int(solution.status),
        "solver_message": str(solution.message),
        "gauge_charge": q,
        "omega": omega,
        "frequency_window": {
            "lower_open": omega_min, "upper_open": m, "satisfied": True,
            "interpretation": "necessary asymptotic window with A0(Rmax)=0",
        },
        "central_amplitude": f0,
        "radius_max": rmax,
        "mesh_nodes": int(solution.x.size),
        "audit_points": int(audit_radius.size),
        "collocation_max_rms_residual": float(np.max(solution.rms_residuals)),
        "scalar_eom_residual_l2": float(np.sqrt(np.trapezoid(
            measure * scalar_residual**2, audit_radius
        ))),
        "scalar_eom_residual_max": float(np.max(np.abs(scalar_residual))),
        "maxwell_eom_residual_l2": float(np.sqrt(np.trapezoid(
            measure * maxwell_residual**2, audit_radius
        ))),
        "maxwell_eom_residual_max": float(np.max(np.abs(maxwell_residual))),
        "origin_scalar_regularity_residual": float(abs(fp[0])),
        "origin_maxwell_regularity_residual": float(abs(ap[0])),
        "outer_scalar_residual": float(abs(f[-1])),
        "outer_gauge_residual": float(abs(a0[-1])),
        "profile_min": float(np.min(f)),
        "profile_max": float(np.max(f)),
        "local_frequency_min": float(np.min(local_frequency)),
        "local_frequency_max": float(np.max(local_frequency)),
        "scalar_energy": energy_scalar,
        "maxwell_electric_energy": energy_maxwell,
        "energy": energy,
        "noether_particle_charge": particle_charge,
        "electric_charge": electric_charge,
        "energy_per_particle_charge": float(energy_per_particle),
        "gauss_flux_charge": gauss_flux_charge,
        "gauss_law_residual": float(gauss_residual),
        "profile_is_stationary_solution": True,
        "maxwell_backreaction_included": True,
        "vierbein_backreaction_included": False,
        "linear_stability_established": False,
        "nonlinear_stability_established": False,
        "outer_boundary_is_finite_cutoff": True,
        "radius": audit_radius,
        "amplitude": f,
        "amplitude_derivative": fp,
        "electrostatic_potential": a0,
        "electric_radial_field": -ap,
        "local_frequency": local_frequency,
        "scalar_equation_residual": scalar_residual,
        "maxwell_equation_residual": maxwell_residual,
        "scope": (
            "Stationary electrostatic gauged-U(1) Q-ball with self-consistent "
            "Maxwell A0 on a finite flat radial domain. No magnetic field, "
            "vierbein/TEGR backreaction, or perturbative stability analysis."
        ),
    }
