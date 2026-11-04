"""Hilbert--Cantor / QHD / U(1) monopole / MHD coupling primitives.

This module provides finite, independently testable discretizations. It is
not a claimed continuum derivation of RGS physics. The 6-bit Q6 address is
mapped to a level-3 Hilbert-grid point and to a middle-third Cantor coordinate.
The fluid and gauge routines are periodic 1-D semidiscrete benchmarks.
Monopole charge is computed from an oriented closed-surface flux cochain.
"""
from __future__ import annotations
import numpy as np


def hilbert_d2xy(order: int, distance: int) -> tuple[int, int]:
    """Map Hilbert distance d in [0,4**order) to integer grid coordinates."""
    if order < 0 or not 0 <= distance < 4**order:
        raise ValueError("distance must lie in [0, 4**order)")
    x = y = 0
    t = int(distance)
    s = 1
    while s < (1 << order):
        rx = 1 & (t // 2)
        ry = 1 & (t ^ rx)
        if ry == 0:
            if rx == 1:
                x, y = s - 1 - x, s - 1 - y
            x, y = y, x
        x += s * rx
        y += s * ry
        t //= 4
        s *= 2
    return x, y


def q6_hilbert_cantor_table() -> dict:
    """Return all 64 Q6 addresses, Hilbert coordinates and Cantor embeddings."""
    rows = []
    for g in range(64):
        bits = tuple((g >> k) & 1 for k in range(6))
        # Standard Hilbert index is the integer encoded by the six-bit word.
        x, y = hilbert_d2xy(3, g)
        cantor = sum(2.0 * bit / (3.0 ** (k + 1)) for k, bit in enumerate(bits))
        rows.append((g, bits, x, y, cantor))
    xy = np.asarray([(r[2], r[3]) for r in rows])
    if len({tuple(v) for v in xy.tolist()}) != 64:
        raise RuntimeError("Hilbert map is not bijective")
    return {
        "vertices": 64,
        "hilbert_order": 3,
        "hilbert_grid_side": 8,
        "hilbert_bijection": True,
        "rows": rows,
        "cantor_min": float(min(r[4] for r in rows)),
        "cantor_max": float(max(r[4] for r in rows)),
        "cantor_coordinate_is_finite_prefix": True,
        "warning": "Six-bit prefixes are finite cylinders, not the infinite Cantor set."
    }


def cantor_cylinder_mass(bits, branch_probabilities=(0.5, 0.5)) -> float:
    """Bernoulli cylinder mass for a finite binary address."""
    p0, p1 = map(float, branch_probabilities)
    if min(p0, p1) <= 0 or not np.isclose(p0 + p1, 1.0):
        raise ValueError("branch probabilities must be positive and sum to one")
    value = 1.0
    for bit in bits:
        if int(bit) not in (0, 1):
            raise ValueError("bits must be binary")
        value *= p1 if int(bit) else p0
    return value


def periodic_derivative(field: np.ndarray, dx: float) -> np.ndarray:
    """Centered periodic first derivative."""
    f = np.asarray(field, dtype=float)
    if f.ndim != 1 or f.size < 3 or dx <= 0:
        raise ValueError("field must be 1-D (n>=3) and dx positive")
    return (np.roll(f, -1) - np.roll(f, 1)) / (2.0 * dx)


def periodic_laplacian(field: np.ndarray, dx: float) -> np.ndarray:
    """Centered periodic second derivative."""
    f = np.asarray(field, dtype=float)
    if f.ndim != 1 or f.size < 3 or dx <= 0:
        raise ValueError("field must be 1-D (n>=3) and dx positive")
    return (np.roll(f, -1) - 2.0 * f + np.roll(f, 1)) / dx**2


def qhd_madelung_rhs(rho, velocity, dx: float, hbar: float = 1.0,
                     mass: float = 1.0, potential=None):
    """Isothermal-free Madelung benchmark: continuity and quantum Euler RHS.

    Q = -(hbar**2/(2m)) (d2 sqrt(rho)/dx2)/sqrt(rho);
    m (dv/dt + v dv/dx) = -d(V+Q)/dx.
    """
    r = np.asarray(rho, dtype=float)
    v = np.asarray(velocity, dtype=float)
    if r.shape != v.shape or r.ndim != 1 or np.any(r <= 0) or mass <= 0:
        raise ValueError("rho,v must be matching 1-D arrays; rho,mass positive")
    root = np.sqrt(r)
    quantum_potential = -(hbar**2 / (2.0 * mass)) * periodic_laplacian(root, dx) / root
    V = np.zeros_like(r) if potential is None else np.asarray(potential, dtype=float)
    if V.shape != r.shape:
        raise ValueError("potential shape mismatch")
    drho_dt = -periodic_derivative(r * v, dx)
    dv_dt = -v * periodic_derivative(v, dx) - periodic_derivative(
        V + quantum_potential, dx
    ) / mass
    return {"drho_dt": drho_dt, "dvelocity_dt": dv_dt,
            "quantum_potential": quantum_potential}


def u1_covariant_difference(psi, links):
    """Gauge-covariant forward difference on a periodic 1-D link lattice.

    links[i] transports the fibre at i+1 to i; psi[i+1] is multiplied by it.
    Under psi_i -> g_i psi_i, U_i -> g_i U_i g_(i+1)^*, covariance is exact.
    """
    z = np.asarray(psi, dtype=complex)
    u = np.asarray(links, dtype=complex)
    if z.ndim != 1 or u.shape != z.shape or np.any(np.abs(np.abs(u) - 1.0) > 1e-10):
        raise ValueError("psi and unit-modulus links must be matching 1-D arrays")
    return u * np.roll(z, -1) - z


def u1_gauge_transform(psi, links, phases):
    z = np.asarray(psi, dtype=complex)
    u = np.asarray(links, dtype=complex)
    g = np.exp(1j * np.asarray(phases, dtype=float))
    if z.shape != u.shape or z.shape != g.shape:
        raise ValueError("psi, links and phases must have matching shapes")
    # U_i is the transport i+1 -> i.
    transformed_links = g * u * np.conjugate(np.roll(g, -1))
    return g * z, transformed_links


def monopole_chern_number(oriented_face_fluxes) -> dict:
    """Compute c1=(1/2pi) integral F on a closed oriented 2-cycle."""
    flux = np.asarray(oriented_face_fluxes, dtype=float)
    if flux.ndim != 1 or flux.size == 0 or not np.all(np.isfinite(flux)):
        raise ValueError("provide a finite nonempty oriented face-flux array")
    total = float(np.sum(flux))
    charge = total / (2.0 * np.pi)
    nearest = int(np.rint(charge))
    return {
        "total_flux": total,
        "chern_number": float(charge),
        "nearest_integer_charge": nearest,
        "quantization_residual": float(abs(charge - nearest)),
        "quantized_within_tolerance": bool(abs(charge - nearest) < 1e-10),
        "interpretation": "Requires a closed oriented 2-cycle and consistent patching; flux values alone do not construct a global monopole potential."
    }


def mhd_1d_transverse_rhs(velocity, magnetic, dx: float, density: float = 1.0,
                          magnetic_permeability: float = 1.0,
                          viscosity: float = 0.0, resistivity: float = 0.0):
    """Periodic 1-D transverse ideal/resistive MHD benchmark.

    u_t + u u_x = -(B B_x)/(mu rho) + nu u_xx
    B_t + (u B)_x = eta B_xx.
    This is a reduced transverse-field system, not full 3-D MHD.
    """
    u = np.asarray(velocity, dtype=float)
    b = np.asarray(magnetic, dtype=float)
    if u.shape != b.shape or u.ndim != 1 or density <= 0 or magnetic_permeability <= 0:
        raise ValueError("u,B must match and density, permeability be positive")
    du = (-u * periodic_derivative(u, dx)
          - b * periodic_derivative(b, dx) / (magnetic_permeability * density)
          + viscosity * periodic_laplacian(u, dx))
    db = (-periodic_derivative(u * b, dx)
          + resistivity * periodic_laplacian(b, dx))
    current = periodic_derivative(b, dx) / magnetic_permeability
    lorentz_force = -b * periodic_derivative(b, dx) / magnetic_permeability
    return {"dvelocity_dt": du, "dmagnetic_dt": db,
            "current": current, "lorentz_force": lorentz_force}


def audit_hilbert_cantor_qhd_gauge_mhd() -> dict:
    table = q6_hilbert_cantor_table()
    n = 64
    x = np.arange(n) / n
    rho = 1.0 + 0.1 * np.cos(2.0 * np.pi * x)
    velocity = 0.05 * np.sin(2.0 * np.pi * x)
    qhd = qhd_madelung_rhs(rho, velocity, 1.0 / n)
    magnetic = 0.2 + 0.03 * np.cos(2.0 * np.pi * x)
    mhd = mhd_1d_transverse_rhs(velocity, magnetic, 1.0 / n)
    psi = np.exp(1j * 2.0 * np.pi * x)
    links = np.exp(1j * 0.1 * np.cos(2.0 * np.pi * x))
    phases = 0.3 * np.sin(4.0 * np.pi * x)
    d0 = u1_covariant_difference(psi, links)
    psi_g, links_g = u1_gauge_transform(psi, links, phases)
    d1 = u1_covariant_difference(psi_g, links_g)
    expected = np.exp(1j * phases) * d0
    monopole = monopole_chern_number(np.full(4, 2.0 * np.pi / 4.0))
    fluxon = fluxon_topological_charge(2.0 * np.pi * x)
    return {
        "schema": "RGS_PSU_HILBERT_CANTOR_QHD_FLUXON_MHD_V2",
        "q6_hilbert": {k: table[k] for k in ("vertices", "hilbert_order", "hilbert_grid_side", "hilbert_bijection")},
        "cantor": {"finite_prefix": True, "uniform_6_bit_cylinder_mass": 1.0 / 64.0},
        "qhd": {"continuity_mass_residual": float(abs(np.sum(qhd["drho_dt"]))),
                "quantum_potential_finite": bool(np.all(np.isfinite(qhd["quantum_potential"])))},
        "u1_gauge": {"covariance_residual": float(np.linalg.norm(d1 - expected))},
        "monopole": monopole,
        "fluxon": fluxon,
        "mhd": {"induction_mass_residual": float(abs(np.sum(mhd["dmagnetic_dt"]))),
                "lorentz_force_finite": bool(np.all(np.isfinite(mhd["lorentz_force"])))},
        "scope": "finite periodic benchmark; no continuum RGS field-equation theorem claimed"
    }


def fluxon_gauge_invariant_differences(phases, link_phases=None) -> np.ndarray:
    """Compact gauge-invariant phase drop on each periodic link.

    With U_{i<-i+1}=exp(i a_i), the invariant phase is
    theta[i+1]-theta[i]+a_i.  Under theta_i -> theta_i+chi_i and
    a_i -> a_i+chi_i-chi_{i+1}, it is unchanged modulo 2 pi.
    """
    theta = np.asarray(phases, dtype=float)
    if theta.ndim != 1 or theta.size < 3 or not np.all(np.isfinite(theta)):
        raise ValueError("phases must be a finite 1-D array of length >= 3")
    a = np.zeros_like(theta) if link_phases is None else np.asarray(link_phases, dtype=float)
    if a.shape != theta.shape or not np.all(np.isfinite(a)):
        raise ValueError("link_phases must be finite and match phases")
    raw = np.roll(theta, -1) - theta + a
    return np.angle(np.exp(1j * raw))


def fluxon_topological_charge(phases, link_phases=None, atol: float = 1e-10) -> dict:
    """Audit integer fluxon winding on a periodic 1-D compact U(1) lattice.

    Quantization is asserted only when the total background link holonomy is
    trivial (equivalently total link phase is 0 modulo 2 pi).
    """
    delta = fluxon_gauge_invariant_differences(phases, link_phases)
    a = np.zeros_like(delta) if link_phases is None else np.asarray(link_phases, dtype=float)
    total_link_phase = float(np.sum(a))
    holonomy_residual = float(abs(np.exp(1j * total_link_phase) - 1.0))
    winding = float(np.sum(delta) / (2.0 * np.pi))
    nearest = int(np.rint(winding))
    quantized = holonomy_residual <= atol and abs(winding - nearest) <= atol
    return {
        "gauge_invariant_phase_drops": delta.tolist(),
        "total_covariant_phase": float(np.sum(delta)),
        "fluxon_number": winding,
        "nearest_integer_fluxon_number": nearest,
        "total_link_holonomy_residual": holonomy_residual,
        "quantized": bool(quantized),
        "quantization_residual": float(abs(winding - nearest)),
    }


def sine_gordon_fluxon_profile(x, center: float = 0.0,
                               josephson_length: float = 1.0) -> dict:
    """Analytic unit-kink solution of lambda_J^2 phi_xx = sin(phi).

    Its phase changes by 2 pi, and B=(Phi0/2pi) d(phi)/dx integrates to Phi0.
    """
    xx = np.asarray(x, dtype=float)
    lam = float(josephson_length)
    if xx.ndim != 1 or xx.size < 3 or not np.all(np.isfinite(xx)):
        raise ValueError("x must be a finite 1-D grid with at least 3 points")
    if not np.isfinite(center) or not np.isfinite(lam) or lam <= 0:
        raise ValueError("center must be finite and Josephson length positive")
    y = (xx - center) / lam
    phase = 4.0 * np.arctan(np.exp(np.clip(y, -350.0, 350.0)))
    derivative = 2.0 / lam / np.cosh(np.clip(y, -350.0, 350.0))
    return {"x": xx, "phase": phase, "phase_gradient": derivative,
            "normalized_magnetic_field": derivative / (2.0 * np.pi),
            "phase_jump": float(phase[-1] - phase[0])}


def sine_gordon_static_residual(phase, dx: float,
                                josephson_length: float = 1.0) -> np.ndarray:
    """Centered finite-difference residual lambda_J^2 phi_xx - sin(phi)."""
    phi = np.asarray(phase, dtype=float)
    lam = float(josephson_length)
    if phi.ndim != 1 or phi.size < 3 or dx <= 0 or lam <= 0:
        raise ValueError("phase must be 1-D (n>=3), dx and Josephson length positive")
    d2 = (phi[2:] - 2.0 * phi[1:-1] + phi[:-2]) / dx**2
    return lam**2 * d2 - np.sin(phi[1:-1])


def qball_potential_amplitude(amplitude, mass_scale: float = 1.0,
                              quartic: float = 1.0,
                              sextic: float = 0.2):
    """V(f)=m^2 f^2-lambda f^4+g f^6 for a complex scalar Phi=f exp(i omega t)."""
    f = np.asarray(amplitude, dtype=float)
    m, lam, g = float(mass_scale), float(quartic), float(sextic)
    if m <= 0 or lam < 0 or g <= 0 or not np.all(np.isfinite(f)):
        raise ValueError("mass_scale,sextic must be positive; quartic nonnegative")
    return m**2 * f**2 - lam * f**4 + g * f**6


def qball_potential_derivative(amplitude, mass_scale: float = 1.0,
                               quartic: float = 1.0,
                               sextic: float = 0.2):
    f = np.asarray(amplitude, dtype=float)
    m, lam, g = float(mass_scale), float(quartic), float(sextic)
    if m <= 0 or lam < 0 or g <= 0 or not np.all(np.isfinite(f)):
        raise ValueError("mass_scale,sextic must be positive; quartic nonnegative")
    return 2.0 * m**2 * f - 4.0 * lam * f**3 + 6.0 * g * f**5


def qball_radial_residual(radius, amplitude, frequency: float,
                          mass_scale: float = 1.0, quartic: float = 1.0,
                          sextic: float = 0.2) -> np.ndarray:
    """Residual of f''+2f'/r = (1/2)V'(f)-omega^2 f in 3 spatial dimensions."""
    r = np.asarray(radius, dtype=float)
    f = np.asarray(amplitude, dtype=float)
    if (r.ndim != 1 or f.shape != r.shape or r.size < 5
            or not np.all(np.isfinite(r)) or not np.all(np.isfinite(f))
            or r[0] < 0 or np.any(np.diff(r) <= 0) or frequency <= 0):
        raise ValueError("radius/amplitude must be finite increasing 1-D arrays; frequency positive")
    fp = np.gradient(f, r, edge_order=2)
    fpp = np.gradient(fp, r, edge_order=2)
    radial_laplacian = fpp.copy()
    radial_laplacian[1:] += 2.0 * fp[1:] / r[1:]
    if r[0] == 0.0:
        radial_laplacian[0] = 3.0 * fpp[0]
    else:
        radial_laplacian[0] += 2.0 * fp[0] / r[0]
    rhs = 0.5 * qball_potential_derivative(
        f, mass_scale=mass_scale, quartic=quartic, sextic=sextic
    ) - float(frequency)**2 * f
    return radial_laplacian - rhs


def qball_energy_charge(radius, amplitude, frequency: float,
                        mass_scale: float = 1.0, quartic: float = 1.0,
                        sextic: float = 0.2) -> dict:
    """3-D radial energy and Noether charge for an ungauged Q-ball ansatz."""
    r = np.asarray(radius, dtype=float)
    f = np.asarray(amplitude, dtype=float)
    if (r.ndim != 1 or f.shape != r.shape or r.size < 5
            or r[0] < 0 or np.any(np.diff(r) <= 0)
            or not np.all(np.isfinite(r)) or not np.all(np.isfinite(f))
            or frequency <= 0):
        raise ValueError("invalid radial profile or frequency")
    fp = np.gradient(f, r, edge_order=2)
    V = qball_potential_amplitude(
        f, mass_scale=mass_scale, quartic=quartic, sextic=sextic
    )
    density_energy = frequency**2 * f**2 + fp**2 + V
    density_charge = 2.0 * frequency * f**2
    measure = 4.0 * np.pi * r**2
    energy = float(np.trapz(measure * density_energy, r))
    charge = float(np.trapz(measure * density_charge, r))
    return {
        "energy": energy,
        "noether_charge": charge,
        "energy_per_charge": energy / charge if charge > 0 else float("inf"),
        "frequency": float(frequency),
        "mass_threshold": float(mass_scale),
        "below_free_particle_threshold": bool(charge > 0 and energy / charge < mass_scale),
        "profile_is_stationary_solution": False,
        "interpretation": "Energy/charge is a variational diagnostic; stability requires a converged stationary profile and a perturbation analysis."
    }


def qball_gaussian_trial(radius, amplitude_scale: float = 0.5,
                         width: float = 2.0, frequency: float = 0.8,
                         mass_scale: float = 1.0, quartic: float = 1.0,
                         sextic: float = 0.2) -> dict:
    """Reproducible Gaussian trial profile; not a claim of an exact Q-ball solution."""
    r = np.asarray(radius, dtype=float)
    if r.ndim != 1 or r.size < 5 or r[0] < 0 or np.any(np.diff(r) <= 0):
        raise ValueError("radius must be an increasing 1-D grid")
    if amplitude_scale <= 0 or width <= 0:
        raise ValueError("amplitude_scale and width must be positive")
    f = float(amplitude_scale) * np.exp(-(r / float(width))**2)
    result = qball_energy_charge(
        r, f, frequency, mass_scale=mass_scale,
        quartic=quartic, sextic=sextic
    )
    residual = qball_radial_residual(
        r, f, frequency, mass_scale=mass_scale,
        quartic=quartic, sextic=sextic
    )
    result.update({
        "radius": r,
        "amplitude": f,
        "radial_eom_residual_l2": float(np.sqrt(np.trapz(4.0 * np.pi * r**2 * residual**2, r))),
        "radial_eom_residual_max": float(np.max(np.abs(residual))),
        "profile_is_stationary_solution": False,
        "ansatz": "Gaussian variational trial only"
    })
    return result


def qball_frequency_window(mass_scale: float = 1.0, quartic: float = 1.0,
                           sextic: float = 0.2) -> dict:
    """Necessary thin-wall existence window from min_{f>0} V(f)/f^2 < omega^2 < m^2.

    For V=m^2 f^2-lambda f^4+g f^6, min(V/f^2)=m^2-lambda^2/(4g).
    A positive lower endpoint additionally requires lambda^2 < 4 g m^2.
    """
    m, lam, g = float(mass_scale), float(quartic), float(sextic)
    if m <= 0 or lam <= 0 or g <= 0:
        raise ValueError("mass_scale, quartic and sextic must be positive")
    lower_sq = m**2 - lam**2 / (4.0 * g)
    stable_vacuum = lower_sq > 0.0
    return {
        "omega_min_squared": float(lower_sq),
        "omega_min": float(np.sqrt(lower_sq)) if stable_vacuum else None,
        "omega_max": m,
        "stable_vacuum": bool(stable_vacuum),
        "nonempty_frequency_window": bool(stable_vacuum and lower_sq < m**2),
        "criterion": "sqrt(min(V(f)/f^2)) < omega < mass_scale; necessary Q-ball existence window, not a stability proof"
    }
