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
