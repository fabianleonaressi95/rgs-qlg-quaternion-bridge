"""On-shell symmetric Hilbert tensor for torsion-eliminated EC-Dirac matter.

Conventions are inherited from rgs_ec_algebraic_torsion:
eta=(-,+,+,+), gamma^a anti-Hermitian Clifford generators and
bar(psi)=i psi^dagger gamma^0.  The kinetic tensor is the Belinfante/Hilbert
form obtained after varying the symmetrized Dirac action and eliminating the
spin-connection variation.  The axial contact contribution is varied
explicitly at fixed local-frame state.  This module does not solve the Dirac
equation or a cosmological state equation.
"""
from __future__ import annotations

import numpy as np

from src.rgs_ec_algebraic_torsion import (
    dirac_adjoint,
    dirac_axial_current,
    dirac_gamma_matrices_minus_plus,
    minkowski_metric,
)
from src.rgs_ec_dirac_four_fermion import effective_four_fermion_coefficient


def _finite_complex(value, shape, name):
    x = np.asarray(value, dtype=complex)
    if x.shape != shape or not np.all(np.isfinite(x)):
        raise ValueError(f"{name} must be a finite array with shape {shape}")
    return x


def dirac_kinetic_hilbert_from_bilinears(right_bilinear, left_bilinear,
                                          *, reality_tolerance=1e-10):
    """Build T_ab=(i/4)[barψγ_a↔D_bψ + barψγ_b↔D_aψ].

    Arrays use R[a,b]=barψ gamma_a D_b ψ and
    L[a,b]=(D_b barψ) gamma_a ψ.  The factor i is required by the
    repository's declared gamma/bar convention.  Inputs may be complex,
    while the physical symmetric tensor must be real.
    """
    if not np.isfinite(reality_tolerance) or reality_tolerance < 0:
        raise ValueError("reality_tolerance must be finite and nonnegative")
    right = _finite_complex(right_bilinear, (4, 4), "right_bilinear")
    left = _finite_complex(left_bilinear, (4, 4), "left_bilinear")
    raw = 0.25j * (right + right.T - left - left.T)
    imag_residual = float(np.linalg.norm(raw.imag))
    scale = max(float(np.linalg.norm(raw.real)), 1.0)
    if imag_residual > reality_tolerance * scale:
        raise ArithmeticError(
            "kinetic Hilbert tensor is not real under the declared Dirac convention"
        )
    tensor = 0.5 * (raw.real + raw.real.T)
    return {
        "tensor": tensor,
        "symmetry_residual": float(np.linalg.norm(tensor - tensor.T)),
        "reality_residual": imag_residual,
        "construction": "BELINFANTE_SYMMETRIZED_DIRAC_HILBERT_ON_SHELL",
    }


def dirac_kinetic_hilbert_from_spinor(spinor, covariant_derivative_spinor,
                                       *, reality_tolerance=1e-10):
    """Compute the kinetic Hilbert tensor from D_b psi in a local frame.

    covariant_derivative_spinor[:, b] is D_b psi.  Metric-compatible
    spin connection implies D_b bar(psi)=bar(D_b psi), which is used here.
    """
    psi = np.asarray(spinor, dtype=complex)
    if psi.shape != (4,) or not np.all(np.isfinite(psi)):
        raise ValueError("spinor must be a finite complex 4-vector")
    derivative = _finite_complex(
        covariant_derivative_spinor, (4, 4), "covariant_derivative_spinor"
    )
    gamma_up = dirac_gamma_matrices_minus_plus()
    eta = minkowski_metric()
    gamma_down = tuple(sum(eta[a, c] * gamma_up[c] for c in range(4))
                       for a in range(4))
    bar = dirac_adjoint(psi)
    right = np.empty((4, 4), dtype=complex)
    left = np.empty((4, 4), dtype=complex)
    for a in range(4):
        for b in range(4):
            dpsi = derivative[:, b]
            dbar = dirac_adjoint(dpsi)
            right[a, b] = bar @ gamma_down[a] @ dpsi
            left[a, b] = dbar @ gamma_down[a] @ psi
    result = dirac_kinetic_hilbert_from_bilinears(
        right, left, reality_tolerance=reality_tolerance
    )
    result["construction"] = (
        "SPINOR_BILINEARS_WITH_METRIC_COMPATIBLE_SPIN_CONNECTION"
    )
    return result


def axial_contact_hilbert_tensor(second_moment, kappa):
    """Return T_ab^(4f)=eta_ab (3kappa/16) eta_cd <J5^c J5^d>."""
    moment = np.asarray(second_moment, dtype=float)
    if moment.shape != (4, 4) or not np.all(np.isfinite(moment)):
        raise ValueError("second_moment must be a finite (4,4) matrix")
    if not np.isfinite(kappa) or kappa <= 0:
        raise ValueError("kappa must be finite and positive")
    moment = 0.5 * (moment + moment.T)
    eta = minkowski_metric()
    alpha = effective_four_fermion_coefficient(float(kappa))
    scalar = float(alpha * np.einsum("ab,ab->", eta, moment))
    return {
        "tensor": eta * scalar,
        "contact_lagrangian": scalar,
        "contact_coefficient": alpha,
        "moment_symmetry_residual": float(np.linalg.norm(moment-moment.T)),
        "construction": "EXPLICIT_FIXED_LOCAL_STATE_CONTACT_VARIATION",
    }


def full_dirac_hilbert_tensor(kinetic_tensor, second_moment, kappa):
    """Combine on-shell kinetic and explicit contact Hilbert contributions.

    The caller must supply the kinetic tensor evaluated with a solution (or
    controlled approximate solution) of the effective Dirac equation.  The
    returned object is not a closed perfect fluid: no state-response,
    thermodynamic closure, or conservation equation is inferred.
    """
    kinetic = np.asarray(kinetic_tensor, dtype=float)
    if kinetic.shape != (4, 4) or not np.all(np.isfinite(kinetic)):
        raise ValueError("kinetic_tensor must be a finite (4,4) matrix")
    kinetic = 0.5 * (kinetic + kinetic.T)
    contact = axial_contact_hilbert_tensor(second_moment, kappa)
    total = kinetic + contact["tensor"]
    return {
        "schema": "RGS_EC_FULL_DIRAC_HILBERT_TENSOR_V1",
        "metric_signature": "(-,+,+,+)",
        "kinetic_hilbert_tensor": kinetic,
        "contact_hilbert_tensor": contact["tensor"],
        "total_hilbert_tensor": total,
        "contact_lagrangian": contact["contact_lagrangian"],
        "kinetic_symmetry_residual": float(np.linalg.norm(kinetic-kinetic.T)),
        "total_symmetry_residual": float(np.linalg.norm(total-total.T)),
        "state_response_included": False,
        "dirac_eom_enforced": False,
        "conservation_status": "REQUIRES_EFFECTIVE_DIRAC_EOM_AND_STATE_CLOSURE",
        "cosmological_fluid_status": "NOT_INFERRED",
        "torsion_propagating_dofs": 0,
    }


def pure_spinor_second_moment(spinor):
    """Return J5^a J5^b for one classical spinor sample."""
    current = dirac_axial_current(spinor)
    return np.outer(current, current)
