"""Stationary symmetric two-field PSU–Anti-Twin Q-ball with Maxwell backreaction.

Potential:
    I = |Phi_+|^2 + |Phi_-|^2
    J = |Phi_+|^2 - |Phi_-|^2
    V = m^2 I - lam I^2 + g I^3 + alpha J^2
        - kappa (Phi_+^* Phi_- + Phi_-^* Phi_+)

Conventions: (-,+,+,+), D_mu=partial_mu-iqA_mu,
Phi_±=f_±(r) exp(i omega t), A=A0(r)dt,
Omega=omega-q A0. This solver computes the exchange-symmetric,
common-phase branch f_+=f_-=F/sqrt(2), with a shared diagonal U(1).
It is flat-space/electrostatic and does not establish stability or
vierbein/TEGR backreaction.
"""
from __future__ import annotations

import numpy as np
from scipy.integrate import solve_bvp

from src.rgs_psu_qball_bvp import qball_stationary_bvp


def psu_antitwin_qball_maxwell_bvp(
    *,
    gauge_charge: float = 0.05,
    central_amplitude: float = 0.5,
    radius_max: float = 20.0,
    points: int = 160,
    mass_scale: float = 1.0,
    quartic: float = 1.5,
    sextic: float = 1.0,
    imbalance_penalty: float = 0.2,
    exchange_coupling: float = 0.1,
    tolerance: float = 1.0e-7,
    max_nodes: int = 20000,
) -> dict:
    """Solve the exchange-symmetric two-field stationary Maxwell/Q-ball BVP.

    The field amplitudes are F/sqrt(2) on each sheet. The interaction
    parameters are explicit model inputs, not fitted or universal constants.
    The symmetric branch has m_eff^2 = m^2-kappa; alpha drops out of its
    background equation but enters the relative-sheet amplitude Hessian.

    The finite-radius boundary conditions are regularity at r=0,
    F(Rmax)=0, A0'(0)=0 and A0(Rmax)=0. The latter fixes the residual
    constant-gauge freedom. The asymptotic frequency interval is necessary,
    not sufficient, for a localized solution.
    """
    q = float(gauge_charge)
    m, lam, g = float(mass_scale), float(quartic), float(sextic)
    alpha, kappa = float(imbalance_penalty), float(exchange_coupling)
    f0, rmax, tol = float(central_amplitude), float(radius_max), float(tolerance)
    values = [q, m, lam, g, alpha, kappa, f0, rmax, tol]
    if not np.all(np.isfinite(values)):
        raise ValueError("all scalar parameters must be finite")
    if m <= 0.0 or lam <= 0.0 or g <= 0.0:
        raise ValueError("mass_scale, quartic and sextic must be positive")
    if alpha < 0.0 or kappa < 0.0:
        raise ValueError("imbalance_penalty and exchange_coupling must be nonnegative")
    if f0 <= 0.0 or rmax <= 0.0 or tol <= 0.0:
        raise ValueError("central_amplitude, radius_max and tolerance must be positive")
    if points < 40 or max_nodes < points:
        raise ValueError("require points>=40 and max_nodes>=points")

    effective_mass_sq = m * m - kappa
    lower_sq = effective_mass_sq - lam * lam / (4.0 * g)
    if lower_sq <= 0.0:
        raise ValueError(
            "symmetric branch has no positive Q-ball frequency window: "
            "require m^2-kappa > quartic^2/(4 sextic)"
        )
    effective_mass = float(np.sqrt(effective_mass_sq))
    omega_min = float(np.sqrt(lower_sq))
    omega_max = effective_mass
    split = 1.0 / np.sqrt(2.0)

    # Continue from the already-tested one-field global Q-ball at q=0,
    # using the exact symmetric-branch mass shift m_eff^2=m^2-kappa.
    seed = qball_stationary_bvp(
        central_amplitude=f0,
        radius_max=rmax,
        points=points,
        mass_scale=effective_mass,
        quartic=lam,
        sextic=g,
        tolerance=tol,
        max_nodes=max_nodes,
    )
    radius = np.linspace(0.0, rmax, int(points))
    profile = np.interp(radius, seed["radius"], seed["amplitude"])
    derivative = np.gradient(profile, radius, edge_order=2)
    shape = np.exp(-radius / max(1.0, rmax / 8.0)) * (1.0 - radius / rmax) ** 2
    a_seed = q * f0 * f0 * shape
    ap_seed = np.gradient(a_seed, radius, edge_order=2)
    initial = np.vstack((
        split * profile, split * derivative,
        split * profile, split * derivative,
        a_seed, ap_seed,
    ))

    def equations(r, y, parameter):
        omega = parameter[0]
        f_plus, fp_plus, f_minus, fp_minus, a0, ap = y
        local_frequency = omega - q * a0
        invariant_I = f_plus * f_plus + f_minus * f_minus
        invariant_J = f_plus * f_plus - f_minus * f_minus
        common = (
            (m * m - local_frequency * local_frequency)
            - 2.0 * lam * invariant_I
            + 3.0 * g * invariant_I * invariant_I
        )
        rhs_plus = (
            common * f_plus + 2.0 * alpha * invariant_J * f_plus
            - kappa * f_minus
        )
        rhs_minus = (
            common * f_minus - 2.0 * alpha * invariant_J * f_minus
            - kappa * f_plus
        )
        rhs_maxwell = -2.0 * q * local_frequency * invariant_I
        return np.vstack((
            fp_plus, rhs_plus, fp_minus, rhs_minus, ap, rhs_maxwell
        ))

    def boundary(ya, yb, parameter):
        # Seven conditions for six radial state variables plus omega.
        # The omitted f_-(R)=0 follows identically from the exchange-symmetric
        # branch f_-=f_+; imposing it again would overdetermine solve_bvp.
        return np.asarray((
            ya[1],
            ya[0] - split * f0,
            ya[3],
            ya[2] - ya[0],
            yb[0],
            ya[5],
            yb[4],
        ), dtype=float)

    singular = np.zeros((6, 6), dtype=float)
    singular[1, 1] = -2.0
    singular[3, 3] = -2.0
    singular[5, 5] = -2.0
    solution = solve_bvp(
        equations,
        boundary,
        radius,
        initial,
        p=np.asarray((seed["omega"],)),
        S=singular,
        tol=tol,
        max_nodes=int(max_nodes),
        verbose=0,
    )
    if solution.status != 0:
        raise RuntimeError(
            "two-field PSU–Anti-Twin Q-ball/Maxwell BVP did not converge: "
            f"status={solution.status}; {solution.message}"
        )

    omega = float(solution.p[0])
    if not omega_min < omega < omega_max:
        raise RuntimeError(
            f"converged frequency {omega:.12g} outside necessary "
            f"symmetric-branch window ({omega_min:.12g}, {omega_max:.12g})"
        )

    audit_radius = np.linspace(0.0, rmax, max(4 * int(points), 801))
    state = solution.sol(audit_radius)
    second = solution.sol(audit_radius, 2)
    f_plus, fp_plus, f_minus, fp_minus, a0, ap = state
    fpp_plus, _, fpp_minus, _, app, _ = second
    local_frequency = omega - q * a0
    invariant_I = f_plus**2 + f_minus**2
    invariant_J = f_plus**2 - f_minus**2

    common = (
        (m * m - local_frequency**2)
        - 2.0 * lam * invariant_I
        + 3.0 * g * invariant_I**2
    )
    rhs_plus = common * f_plus + 2.0 * alpha * invariant_J * f_plus - kappa * f_minus
    rhs_minus = common * f_minus - 2.0 * alpha * invariant_J * f_minus - kappa * f_plus
    rhs_maxwell = -2.0 * q * local_frequency * invariant_I

    lap_plus = fpp_plus.copy()
    lap_minus = fpp_minus.copy()
    lap_a0 = app.copy()
    lap_plus[1:] += 2.0 * fp_plus[1:] / audit_radius[1:]
    lap_minus[1:] += 2.0 * fp_minus[1:] / audit_radius[1:]
    lap_a0[1:] += 2.0 * ap[1:] / audit_radius[1:]
    lap_plus[0] = 3.0 * fpp_plus[0]
    lap_minus[0] = 3.0 * fpp_minus[0]
    lap_a0[0] = 3.0 * app[0]
    residual_plus = lap_plus - rhs_plus
    residual_minus = lap_minus - rhs_minus
    residual_maxwell = lap_a0 - rhs_maxwell

    potential = (
        m * m * invariant_I
        - lam * invariant_I**2
        + g * invariant_I**3
        + alpha * invariant_J**2
        - 2.0 * kappa * f_plus * f_minus
    )
    measure = 4.0 * np.pi * audit_radius**2
    scalar_energy_density = (
        local_frequency**2 * invariant_I
        + fp_plus**2 + fp_minus**2 + potential
    )
    maxwell_energy_density = 0.5 * ap**2
    scalar_energy = float(np.trapezoid(measure * scalar_energy_density, audit_radius))
    maxwell_energy = float(np.trapezoid(measure * maxwell_energy_density, audit_radius))
    energy = scalar_energy + maxwell_energy

    particle_charge = float(np.trapezoid(
        measure * 2.0 * local_frequency * invariant_I, audit_radius
    ))
    electric_charge = q * particle_charge
    gauss_flux_charge = float(-4.0 * np.pi * rmax**2 * solution.sol(rmax)[5])
    gauss_residual = abs(gauss_flux_charge - electric_charge)

    # Local amplitude Hessian in the exchange-even/odd basis. These are
    # pointwise curvature coefficients, not eigenvalues of the full
    # constrained radial fluctuation operator.
    total_profile = np.sqrt(np.maximum(invariant_I, 0.0))
    common_curvature = (
        effective_mass_sq - local_frequency**2
        - 6.0 * lam * total_profile**2
        + 15.0 * g * total_profile**4
    )
    relative_curvature = (
        m * m + kappa - local_frequency**2
        + (-2.0 * lam + 4.0 * alpha) * total_profile**2
        + 3.0 * g * total_profile**4
    )

    return {
        "schema": "RGS_PSU_ANTITWIN_TWO_FIELD_GAUGED_QBALL_BVP_V1",
        "converged": True,
        "solver_status": int(solution.status),
        "solver_message": str(solution.message),
        "gauge_charge": q,
        "mass_scale": m,
        "quartic": lam,
        "sextic": g,
        "imbalance_penalty": alpha,
        "exchange_coupling": kappa,
        "effective_mass": effective_mass,
        "omega": omega,
        "frequency_window": {
            "lower_open": omega_min,
            "upper_open": omega_max,
            "satisfied": True,
            "necessary_only": True,
        },
        "central_total_amplitude": f0,
        "central_sheet_amplitudes": [float(f_plus[0]), float(f_minus[0])],
        "exchange_symmetry_residual": float(np.max(np.abs(f_plus - f_minus))),
        "radius_max": rmax,
        "mesh_nodes": int(solution.x.size),
        "audit_points": int(audit_radius.size),
        "collocation_max_rms_residual": float(np.max(solution.rms_residuals)),
        "plus_eom_residual_max": float(np.max(np.abs(residual_plus))),
        "minus_eom_residual_max": float(np.max(np.abs(residual_minus))),
        "maxwell_eom_residual_max": float(np.max(np.abs(residual_maxwell))),
        "origin_plus_regularity_residual": float(abs(fp_plus[0])),
        "origin_minus_regularity_residual": float(abs(fp_minus[0])),
        "origin_maxwell_regularity_residual": float(abs(ap[0])),
        "outer_plus_residual": float(abs(f_plus[-1])),
        "outer_minus_residual": float(abs(f_minus[-1])),
        "outer_gauge_residual": float(abs(a0[-1])),
        "local_frequency_min": float(np.min(local_frequency)),
        "local_frequency_max": float(np.max(local_frequency)),
        "scalar_energy": scalar_energy,
        "maxwell_electric_energy": maxwell_energy,
        "energy": energy,
        "noether_particle_charge": particle_charge,
        "electric_charge": electric_charge,
        "energy_per_particle_charge": float(energy / particle_charge),
        "gauss_flux_charge": gauss_flux_charge,
        "gauss_law_residual": float(gauss_residual),
        "common_amplitude_curvature_min": float(np.min(common_curvature)),
        "relative_amplitude_curvature_min": float(np.min(relative_curvature)),
        "full_linear_stability_established": False,
        "vierbein_backreaction_included": False,
        "finite_outer_cutoff": True,
        "radius": audit_radius,
        "f_plus": f_plus,
        "f_minus": f_minus,
        "electrostatic_potential": a0,
        "electric_radial_field": -ap,
        "local_frequency": local_frequency,
        "plus_equation_residual": residual_plus,
        "minus_equation_residual": residual_minus,
        "maxwell_equation_residual": residual_maxwell,
        "scope": (
            "Exchange-symmetric two-field PSU–Anti-Twin Q-ball with common "
            "gauged U(1) and electrostatic Maxwell backreaction on a finite "
            "flat radial domain. No full fluctuation stability or TEGR "
            "backreaction is established."
        ),
    }
