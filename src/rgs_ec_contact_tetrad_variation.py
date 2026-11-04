"""Fixed-state tetrad variation of the minimal EC axial contact operator.

Convention: eta=(-,+,+,+), J5^a=i*bar(psi)*gamma5*gamma^a*psi,
L4f=(3*kappa/16) J5_a J5^a. The spinor components / local-frame
ensemble moments are held fixed during the explicit coframe variation.
This computes only that explicit contribution, not state-response terms.
"""
from __future__ import annotations
import numpy as np
from src.rgs_ec_algebraic_torsion import minkowski_metric
from src.rgs_ec_dirac_four_fermion import effective_four_fermion_coefficient


def contact_tetrad_stress_energy(second_moment, kappa):
    """Return <T_ab> from explicit tetrad variation at fixed local-frame state.

    Since J5^a is a tangent-frame bilinear, M^{ab}=<J5^a J5^b> has no
    explicit coframe dependence at fixed spinor state. Thus
      delta(e L4f)/delta(e^a_mu) = e e_a^mu <L4f>,
    and T_ab^(contact, explicit) = eta_ab <L4f>.
    """
    m = np.asarray(second_moment, dtype=float)
    if m.shape != (4, 4) or not np.all(np.isfinite(m)):
        raise ValueError("second_moment must be a finite (4,4) matrix")
    if not np.isfinite(kappa) or kappa <= 0:
        raise ValueError("kappa must be finite and positive")
    eta = minkowski_metric()
    # Only the symmetric part contributes to J^a J^b.
    m_sym = 0.5 * (m + m.T)
    alpha = effective_four_fermion_coefficient(kappa)
    mean_lagrangian = float(alpha * np.einsum("ab,ab->", eta, m_sym))
    stress = eta * mean_lagrangian
    rho = -mean_lagrangian
    pressure = mean_lagrangian
    return {
        "schema": "RGS_EC_CONTACT_TETRAD_VARIATION_V1",
        "variation_mode": "explicit_coframe_variation_at_fixed_local_spinor_state",
        "mean_contact_lagrangian": mean_lagrangian,
        "stress_energy_tangent": stress,
        "interaction_energy_density_explicit": rho,
        "interaction_pressure_explicit": pressure,
        "equation_of_state_explicit": (-1.0 if abs(rho) > 0 else None),
        "trace": float(np.einsum("ab,ab->", eta, stress)),
        "state_response_included": False,
        "cosmological_closure_status": "OPEN_STATE_RESPONSE_AND_CONSERVATION",
    }
def _finite_array(value, shape, name):
    x = np.asarray(value, dtype=float)
    if x.shape != shape or not np.all(np.isfinite(x)):
        raise ValueError(f"{name} must be a finite array with shape {shape}")
    return x


def _oriented_tetrad(tetrad):
    e = _finite_array(tetrad, (4, 4), "tetrad")
    det = float(np.linalg.det(e))
    if abs(det) < 1e-14:
        raise ValueError("tetrad is singular or numerically degenerate")
    if det <= 0.0:
        raise ValueError("tetrad must be positively oriented (det(e)>0)")
    # det(e)>0 alone does not select a time orientation in 3+1 dimensions:
    # -I has determinant +1 but reverses the chosen future-directed time leg.
    # Fix the chart's time orientation by requiring e^0{}_0 > 0.
    if e[0, 0] <= 0.0:
        raise ValueError("tetrad violates the declared future-directed time orientation")
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
