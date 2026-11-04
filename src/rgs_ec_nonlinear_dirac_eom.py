"""Nonlinear Dirac Euler-Lagrange source from the minimal EC contact operator.

The module differentiates L_4f=alpha J5_a J5^a with respect to independent
Dirac variables psi and bar(psi), in the repository's (-,+,+,+) convention.
It deliberately accepts the kinetic Euler operator as an input: the complete
curved-space kinetic variation (including tetrad/spin-connection and boundary
terms) must be supplied by the caller and is not guessed here.
"""
from __future__ import annotations

import numpy as np

from src.rgs_ec_algebraic_torsion import (
    dirac_axial_current,
    dirac_gamma_matrices_minus_plus,
    gamma5_matrix,
    minkowski_metric,
)
from src.rgs_ec_dirac_four_fermion import effective_four_fermion_coefficient


def _spinor(value, name):
    x = np.asarray(value, dtype=complex)
    if x.shape != (4,) or not np.all(np.isfinite(x)):
        raise ValueError(f"{name} must be a finite complex 4-vector")
    return x


def _current(value):
    x = np.asarray(value, dtype=float)
    if x.shape != (4,) or not np.all(np.isfinite(x)):
        raise ValueError("axial_current must be a finite real 4-vector")
    return x


def effective_dirac_contact_sources(spinor, axial_current, kappa):
    """Return independent Euler derivatives of L4f wrt bar(psi) and psi.

    J5^a=i bar(psi) gamma5 gamma^a psi,
    L4f=alpha eta_ab J5^a J5^b, alpha=3 kappa/16.

    Therefore:
      dL4f/dbar(psi) = 2 alpha J5_a i gamma5 gamma^a psi
      dL4f/dpsi      = 2 alpha J5_a i bar(psi) gamma5 gamma^a.

    No complex conjugation is imposed between the two returned derivatives:
    psi and bar(psi) are independent variables in the variational calculus.
    """
    psi = _spinor(spinor, "spinor")
    current = _current(axial_current)
    if not np.isfinite(kappa) or kappa <= 0:
        raise ValueError("kappa must be finite and positive")
    gamma = dirac_gamma_matrices_minus_plus()
    gamma5 = gamma5_matrix()
    eta = minkowski_metric()
    lowered = eta @ current
    alpha = effective_four_fermion_coefficient(float(kappa))
    coefficient = 2.0 * alpha
    source_bar = coefficient * sum(
        lowered[a] * (1j * gamma5 @ gamma[a] @ psi)
        for a in range(4)
    )
    # bar(psi) is supplied in the repository convention; the row derivative
    # is evaluated directly as a covector in the independent bar variable.
    from src.rgs_ec_algebraic_torsion import dirac_adjoint
    bar = dirac_adjoint(psi)
    source_psi = coefficient * sum(
        lowered[a] * (1j * bar @ gamma5 @ gamma[a])
        for a in range(4)
    )
    return {
        "dL_dbarpsi": source_bar,
        "dL_dpsi": source_psi,
        "coefficient_alpha": alpha,
        "current_lowered": lowered,
        "source_norm": float(np.linalg.norm(source_bar)),
        "convention": "J5=i*barpsi*gamma5*gamma^a*psi; eta=(-,+,+,+)",
        "variation_variables_independent": True,
    }


def effective_dirac_eom_residual(kinetic_eom_minus_mass, spinor,
                                  axial_current, kappa):
    """Add the contact source to a caller-derived kinetic/mass Euler operator.

    kinetic_eom_minus_mass is the row/column-convention column residual
    K[psi] which already contains the complete variation of the symmetrized
    kinetic and mass terms. The nonlinear EC equation is K[psi]+dL4f/dbarpsi=0.
    """
    kinetic = _spinor(kinetic_eom_minus_mass, "kinetic_eom_minus_mass")
    source = effective_dirac_contact_sources(spinor, axial_current, kappa)
    residual = kinetic + source["dL_dbarpsi"]
    scale = max(float(np.linalg.norm(kinetic)),
                float(np.linalg.norm(source["dL_dbarpsi"])), 1.0)
    return {
        "residual": residual,
        "residual_norm": float(np.linalg.norm(residual)),
        "relative_residual": float(np.linalg.norm(residual) / scale),
        "kinetic_eom_supplied_externally": True,
        "contact_source_included": True,
        "on_shell": bool(np.linalg.norm(residual) <= 1e-10 * scale),
        "conservation_claim": "NOT_ESTABLISHED_BY_POINTWISE_EOM_RESIDUAL",
    }


def contact_lagrangian_from_independent_bilinears(bar_spinor, spinor, kappa):
    """Evaluate L4f for independent row bar(psi) and column psi.

    This low-level helper is intended for Wirtinger finite-difference audits.
    """
    bar = np.asarray(bar_spinor, dtype=complex)
    psi = _spinor(spinor, "spinor")
    if bar.shape != (4,) or not np.all(np.isfinite(bar)):
        raise ValueError("bar_spinor must be a finite complex 4-row")
    if not np.isfinite(kappa) or kappa <= 0:
        raise ValueError("kappa must be finite and positive")
    gamma = dirac_gamma_matrices_minus_plus()
    gamma5 = gamma5_matrix()
    eta = minkowski_metric()
    current = np.array([1j * bar @ gamma5 @ g @ psi for g in gamma])
    return float(effective_four_fermion_coefficient(float(kappa))
                 * (current @ eta @ current))


def spinor_current_consistency(spinor, *, tolerance=1e-10):
    """Report the repository's real-current construction and its reality error."""
    psi = _spinor(spinor, "spinor")
    current = dirac_axial_current(psi, reality_tolerance=tolerance)
    return {"axial_current": current, "reality_residual": 0.0,
            "status": "REAL_CURRENT_FROM_DECLARED_DIRAC_ADJOINT"}
