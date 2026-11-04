"""Local torsion-free Levi-Civita spin connection from a differentiable coframe.

Index convention:
    coframe[a, mu] = e^a_mu
    coframe_derivatives[a, nu, mu] = partial_mu e^a_nu
    inverse_coframe[mu, a] = e_a^mu

This is a local continuum differential-geometric kernel. It does not infer
derivatives from a finite RGS point cloud, choose a teleparallel inertial
connection, or add an internal gauge connection.
"""
from __future__ import annotations

import numpy as np

from src.clifford_qlg import ETA, gamma_matrices


def _finite_array(value, shape, name):
    array = np.asarray(value, dtype=float)
    if array.shape != shape or not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must be a finite array with shape {shape}")
    return array


def levi_civita_spin_connection(coframe, coframe_derivatives, *, tolerance=1e-10):
    """Compute Christoffel and torsion-free spin connections at one point.

    The tetrad postulate used is
      partial_mu e^a_nu + omega^a_{b mu} e^b_nu
          - Gamma^rho_{mu nu} e^a_rho = 0.

    Returns omega_mixed[a,b,mu] = omega^a_{b mu} and
    omega_lower[a,b,mu] = eta_ac omega^c_{b mu].
    """
    e = _finite_array(coframe, (4, 4), "coframe")
    de = _finite_array(coframe_derivatives, (4, 4, 4), "coframe_derivatives")
    if not np.isfinite(tolerance) or tolerance <= 0:
        raise ValueError("tolerance must be finite and positive")
    determinant = float(np.linalg.det(e))
    if abs(determinant) <= 1e-14:
        raise ValueError("coframe must be invertible")

    metric = e.T @ ETA @ e
    inverse_metric = np.linalg.inv(metric)
    inverse_coframe = np.linalg.inv(e)

    metric_derivatives = np.empty((4, 4, 4), dtype=float)
    for mu in range(4):
        de_mu = de[:, :, mu]
        metric_derivatives[:, :, mu] = de_mu.T @ ETA @ e + e.T @ ETA @ de_mu

    christoffel = np.zeros((4, 4, 4), dtype=float)
    for rho in range(4):
        for mu in range(4):
            for nu in range(4):
                christoffel[rho, mu, nu] = 0.5 * sum(
                    inverse_metric[rho, sigma] * (
                        metric_derivatives[sigma, nu, mu]
                        + metric_derivatives[sigma, mu, nu]
                        - metric_derivatives[mu, nu, sigma]
                    )
                    for sigma in range(4)
                )

    omega_mixed = np.zeros((4, 4, 4), dtype=float)
    for a in range(4):
        for b in range(4):
            for mu in range(4):
                omega_mixed[a, b, mu] = sum(
                    (
                        sum(christoffel[rho, mu, nu] * e[a, rho]
                            for rho in range(4))
                        - de[a, nu, mu]
                    ) * inverse_coframe[nu, b]
                    for nu in range(4)
                )

    omega_lower = np.einsum("ac,cbm->abm", ETA, omega_mixed)
    antisymmetry_residual = float(np.linalg.norm(
        omega_lower + omega_lower.swapaxes(0, 1)
    ))

    tetrad_postulate = np.zeros((4, 4, 4), dtype=float)
    torsion = np.zeros((4, 4, 4), dtype=float)
    for a in range(4):
        for mu in range(4):
            for nu in range(4):
                tetrad_postulate[a, mu, nu] = (
                    de[a, nu, mu]
                    + sum(omega_mixed[a, b, mu] * e[b, nu] for b in range(4))
                    - sum(christoffel[rho, mu, nu] * e[a, rho]
                          for rho in range(4))
                )
                torsion[a, mu, nu] = (
                    de[a, nu, mu] - de[a, mu, nu]
                    + sum(
                        omega_mixed[a, b, mu] * e[b, nu]
                        - omega_mixed[a, b, nu] * e[b, mu]
                        for b in range(4)
                    )
                )

    postulate_residual = float(np.linalg.norm(tetrad_postulate))
    torsion_residual = float(np.linalg.norm(torsion))
    passed = max(antisymmetry_residual, postulate_residual, torsion_residual) <= tolerance
    return {
        "metric": metric,
        "inverse_metric": inverse_metric,
        "inverse_coframe": inverse_coframe,
        "christoffel": christoffel,
        "omega_mixed": omega_mixed,
        "omega_lower": omega_lower,
        "torsion": torsion,
        "metric_signature": "(-,+,+,+)",
        "coframe_determinant": determinant,
        "antisymmetry_residual": antisymmetry_residual,
        "tetrad_postulate_residual": postulate_residual,
        "torsion_free_residual": torsion_residual,
        "passed": bool(passed),
        "status": "LOCAL_LEVI_CIVITA_SPIN_CONNECTION_FROM_COFRAME_JET",
        "scope": (
            "torsion-free Levi-Civita connection only; no teleparallel inertial "
            "connection and no internal gauge field"
        ),
    }


def spinor_connection_from_lorentz_connection(omega_lower):
    """Return Omega_mu = (1/4) omega_ab_mu gamma^a gamma^b."""
    omega = np.asarray(omega_lower, dtype=float)
    if omega.shape != (4, 4, 4) or not np.all(np.isfinite(omega)):
        raise ValueError("omega_lower must be a finite array with shape (4,4,4)")
    gamma = gamma_matrices()
    result = np.zeros((4, 4, 4), dtype=complex)
    for mu in range(4):
        for a in range(4):
            for b in range(4):
                result[mu] += 0.25 * omega[a, b, mu] * (gamma[a] @ gamma[b])
    return result


def dirac_covariant_derivative(partial_derivative, spinor, spinor_connection):
    """Compute D_mu psi = partial_mu psi + Omega_mu psi.

    Internal gauge terms are intentionally excluded and must be added through
    a separately represented gauge connection, not folded into Omega_mu.
    """
    partial = np.asarray(partial_derivative, dtype=complex)
    psi = np.asarray(spinor, dtype=complex)
    omega_spinor = np.asarray(spinor_connection, dtype=complex)
    if partial.shape != (4, 4) or not np.all(np.isfinite(partial)):
        raise ValueError("partial_derivative[:,mu] must have shape (4,4)")
    if psi.shape != (4,) or not np.all(np.isfinite(psi)):
        raise ValueError("spinor must be a finite complex 4-vector")
    if omega_spinor.shape != (4, 4, 4) or not np.all(np.isfinite(omega_spinor)):
        raise ValueError("spinor_connection must have shape (4,4,4)")
    covariant = np.empty((4, 4), dtype=complex)
    for mu in range(4):
        covariant[:, mu] = partial[:, mu] + omega_spinor[mu] @ psi
    return covariant


def orthonormal_dirac_derivative(coordinate_derivative, inverse_coframe):
    """Convert D_mu psi to D_a psi; output column a is e_a^mu D_mu psi."""
    derivative = np.asarray(coordinate_derivative, dtype=complex)
    inverse = np.asarray(inverse_coframe, dtype=float)
    if derivative.shape != (4, 4) or not np.all(np.isfinite(derivative)):
        raise ValueError("coordinate_derivative must have shape (4,4)")
    if inverse.shape != (4, 4) or not np.all(np.isfinite(inverse)):
        raise ValueError("inverse_coframe must have shape (4,4)")
    return derivative @ inverse
