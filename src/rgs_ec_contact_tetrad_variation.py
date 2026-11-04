"""Explicit tetrad variation of the minimal EC axial contact operator.

Convention: eta=(-,+,+,+), J5^a is a local-Lorentz-frame contravariant
current, and L_4f = c4 * eta_ab J5^a J5^b.  In the tetrad variation below
the spinor components (hence J5^a) are held fixed.  This isolates the
explicit Hilbert variation of the contact operator; it is NOT the complete
on-shell stress tensor of the interacting Dirac system.
"""
from __future__ import annotations
import numpy as np

from src.rgs_ec_algebraic_torsion import minkowski_metric
from src.rgs_ec_dirac_four_fermion import effective_four_fermion_coefficient


def _finite_array(value, shape, name):
    x = np.asarray(value, dtype=float)
    if x.shape != shape or not np.all(np.isfinite(x)):
        raise ValueError(f"{name} must be a finite array with shape {shape}")
    return x


def _oriented_tetrad(tetrad):
    e = _finite_array(tetrad, (4, 4), "tetrad")
    det = float(np.linalg.det(e))
    if det <= 0.0:
        raise ValueError("tetrad must be invertible and positively oriented")
    if abs(det) < 1e-14:
        raise ValueError("tetrad is singular or numerically degenerate")
    return e, det, np.linalg.inv(e)


def contact_lagrangian_density(axial_current, kappa):
    """Return L4f=c4 J5_a J5^a in a local orthonormal frame."""
    j = _finite_array(axial_current, (4,), "axial_current")
    if not np.isfinite(kappa) or kappa <= 0:
        raise ValueError("kappa must be finite and positive")
    c4 = float(effective_four_fermion_coefficient(float(kappa)))
    return c4 * float(j @ minkowski_metric() @ j)


def contact_tetrad_variation(tetrad, axial_current, kappa):
    """Return explicit tetrad derivative and Hilbert tensor of L4f.

    With S4f=int d4x e L4f and J5^a held fixed:
      (1/e) d(e L4f)/d e^a_mu = L4f e_a^mu,
      T_mu_nu^(4f, explicit) = g_mu_nu L4f.
    The latter is only the explicit contact-operator contribution.
    """
    e, det_e, inverse_e = _oriented_tetrad(tetrad)
    lagrangian = contact_lagrangian_density(axial_current, kappa)
    eta = minkowski_metric()
    metric_cov = e.T @ eta @ e
    metric_contra = np.linalg.inv(metric_cov)
    tetrad_gradient = det_e * lagrangian * inverse_e.T
    hilbert_cov = lagrangian * metric_cov
    hilbert_contra = lagrangian * metric_contra
    return {
        "det_tetrad": det_e,
        "metric_covariant": metric_cov,
        "metric_contravariant": metric_contra,
        "contact_lagrangian": lagrangian,
        "action_density": det_e * lagrangian,
        "tetrad_action_gradient": tetrad_gradient,
        "hilbert_stress_covariant_explicit": hilbert_cov,
        "hilbert_stress_contravariant_explicit": hilbert_contra,
        "variation_scope": "EXPLICIT_CONTACT_TERM_AT_FIXED_LOCAL_SPINOR_COMPONENTS",
        "full_on_shell_dirac_stress_tensor": None,
        "fluid_equation_of_state": None,
    }


def finite_difference_tetrad_gradient(tetrad, axial_current, kappa, *,
                                      step=1e-6):
    """Central-difference audit of d(e L4f)/d e^a_mu."""
    if not np.isfinite(step) or step <= 0:
        raise ValueError("step must be finite and positive")
    e = _finite_array(tetrad, (4, 4), "tetrad").copy()
    gradient = np.empty((4, 4), dtype=float)
    for a in range(4):
        for mu in range(4):
            plus, minus = e.copy(), e.copy()
            plus[a, mu] += step
            minus[a, mu] -= step
            ep, dp, _ = _oriented_tetrad(plus)
            em, dm, _ = _oriented_tetrad(minus)
            fp = dp * contact_lagrangian_density(axial_current, kappa)
            fm = dm * contact_lagrangian_density(axial_current, kappa)
            gradient[a, mu] = (fp - fm) / (2.0 * step)
    return gradient
