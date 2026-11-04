"""Stationary three-dimensional global-U(1) Q-ball boundary-value solver.

The solver is a numerical existence/convergence tool for the ungauged
complex-scalar model. It does not claim perturbative or nonlinear stability,
and it does not solve the gauged/curved RGS extension.
"""
from __future__ import annotations

import numpy as np
from scipy.integrate import solve_bvp


def qball_stationary_bvp(
    *,
    central_amplitude: float = 1.0,
    radius_max: float = 20.0,
    points: int = 180,
    mass_scale: float = 1.0,
    quartic: float = 1.0,
    sextic: float = 0.5,
    tolerance: float = 1.0e-7,
    max_nodes: int = 12000,
) -> dict:
    """Solve the nonlinear radial Q-ball BVP, treating omega as an eigenvalue.

    Equation:
        f'' + 2 f'/r = (m²-omega²) f - 2 lambda f³ + 3 g f⁵
    Conditions:
        f(0)=central_amplitude, f'(0)=0, f(radius_max)=0.

    The singular radial term is handled by scipy.solve_bvp's regular-singular
    matrix. Fixing f(0) and solving for omega prevents convergence to the
    trivial f=0 solution. The finite outer radius approximates infinity and
    must be increased in a separate truncation-convergence study.
    """
    m = float(mass_scale)
    lam = float(quartic)
    g = float(sextic)
    f0 = float(central_amplitude)
    rmax = float(radius_max)
    tol = float(tolerance)
    if not np.all(np.isfinite([m, lam, g, f0, rmax, tol])):
        raise ValueError("all scalar parameters must be finite")
    if m <= 0.0 or lam <= 0.0 or g <= 0.0:
        raise ValueError("mass_scale, quartic and sextic must be positive")
    if f0 <= 0.0 or rmax <= 0.0 or points < 40 or tol <= 0.0 or max_nodes < points:
        raise ValueError("require f0,radius_max,tolerance>0, points>=40, max_nodes>=points")

    lower_sq = m * m - lam * lam / (4.0 * g)
    if lower_sq <= 0.0:
        raise ValueError("potential has no positive Q-ball frequency lower bound")
    omega_min = float(np.sqrt(lower_sq))
    omega_max = m
    omega_guess = float(np.sqrt(0.5 * (lower_sq + m * m)))
    decay = np.sqrt(m * m - omega_guess * omega_guess)

    radius = np.linspace(0.0, rmax, int(points))
    profile = f0 * np.exp(-decay * radius) * (1.0 + 0.15 * decay * radius)
    derivative = np.gradient(profile, radius, edge_order=2)
    initial = np.vstack((profile, derivative))

    def equation(r, y, parameter):
        omega = parameter[0]
        f = y[0]
        return np.vstack((
            y[1],
            (m * m - omega * omega) * f - 2.0 * lam * f**3 + 3.0 * g * f**5,
        ))

    def boundary(ya, yb, parameter):
        return np.asarray((ya[1], ya[0] - f0, yb[0]), dtype=float)

    singular = np.asarray(((0.0, 0.0), (0.0, -2.0)), dtype=float)
    solution = solve_bvp(
        equation,
        boundary,
        radius,
        initial,
        p=np.asarray((omega_guess,)),
        S=singular,
        tol=tol,
        max_nodes=int(max_nodes),
        verbose=0,
    )
    if solution.status != 0:
        raise RuntimeError(
            f"Q-ball BVP did not converge: status={solution.status}; {solution.message}"
        )

    omega = float(solution.p[0])
    if not omega_min < omega < omega_max:
        raise RuntimeError(
            f"converged frequency {omega:.12g} lies outside "
            f"the necessary window ({omega_min:.12g}, {omega_max:.12g})"
        )

    # Evaluate the collocation solution on a uniform audit grid.
    audit_radius = np.linspace(0.0, rmax, max(4 * int(points), 801))
    values = solution.sol(audit_radius)
    first = solution.sol(audit_radius, 1)
    second = solution.sol(audit_radius, 2)
    f = values[0]
    fp = values[1]
    fpp = second[0]
    rhs = (m * m - omega * omega) * f - 2.0 * lam * f**3 + 3.0 * g * f**5
    radial_laplacian = fpp.copy()
    radial_laplacian[1:] += 2.0 * fp[1:] / audit_radius[1:]
    radial_laplacian[0] = 3.0 * fpp[0]
    residual = radial_laplacian - rhs

    potential = m * m * f**2 - lam * f**4 + g * f**6
    measure = 4.0 * np.pi * audit_radius**2
    energy_density = omega**2 * f**2 + fp**2 + potential
    charge_density = 2.0 * omega * f**2
    energy = float(np.trapezoid(measure * energy_density, audit_radius))
    charge = float(np.trapezoid(measure * charge_density, audit_radius))
    energy_per_charge = float(energy / charge) if charge > 0.0 else float("inf")

    return {
        "schema": "RGS_PSU_GLOBAL_U1_QBALL_STATIONARY_BVP_V1",
        "converged": True,
        "solver_status": int(solution.status),
        "solver_message": str(solution.message),
        "omega": omega,
        "frequency_window": {
            "lower_open": omega_min,
            "upper_open": omega_max,
            "satisfied": True,
        },
        "central_amplitude": f0,
        "radius_max": rmax,
        "mesh_nodes": int(solution.x.size),
        "audit_points": int(audit_radius.size),
        "collocation_max_rms_residual": float(np.max(solution.rms_residuals)),
        "radial_eom_residual_l2": float(np.sqrt(np.trapezoid(
            measure * residual**2, audit_radius
        ))),
        "radial_eom_residual_max": float(np.max(np.abs(residual))),
        "origin_regular_boundary_residual": float(abs(values[1, 0])),
        "outer_dirichlet_residual": float(abs(f[-1])),
        "profile_min": float(np.min(f)),
        "profile_max": float(np.max(f)),
        "profile_nonnegative_within_tolerance": bool(np.min(f) >= -100.0 * tol),
        "energy": energy,
        "noether_charge": charge,
        "energy_per_charge": energy_per_charge,
        "below_free_particle_threshold": bool(energy_per_charge < m),
        "profile_is_stationary_solution": True,
        "linear_stability_established": False,
        "nonlinear_stability_established": False,
        "outer_boundary_is_finite_cutoff": True,
        "radius": audit_radius,
        "amplitude": f,
        "amplitude_derivative": fp,
        "equation_residual": residual,
        "scope": (
            "Converged stationary solution of the finite-radius, flat-space, "
            "ungauged global-U(1) radial BVP. Stability and gauged/curved "
            "extensions are not established."
        ),
    }
