"""Minimal Einstein-Cartan-Dirac axial current and algebraic torsion.

Conventions: eta=(-,+,+,+), bar(psi)=i psi^dagger gamma^0,
gamma5=i gamma^0 gamma^1 gamma^2 gamma^3, J5^a=i bar(psi) gamma5 gamma^a psi,
epsilon_0123=+1 as an alternating internal symbol. The first-order Cartan
solution is T_abc=-(kappa/2) epsilon_abcd J5^d. Minimal EC torsion is algebraic.
"""
from __future__ import annotations
import itertools
import numpy as np

METRIC_SIGNATURE = "(-,+,+,+)"
DIRAC_EC_CONVENTION_ID = "Bonder-Corral-2018_eta(-+++)_bar=i-psi-dagger-gamma0_gamma5=i-gamma0123_eps0123symbol=+1"
DERIVED_CARTAN_COEFFICIENT = -0.5


def minkowski_metric():
    return np.diag([-1.0, 1.0, 1.0, 1.0])


def dirac_gamma_matrices_minus_plus():
    """Return gamma^a satisfying {gamma^a,gamma^b}=2 diag(-,+,+,+)."""
    eye = np.eye(2, dtype=complex)
    zero = np.zeros((2, 2), dtype=complex)
    sigma1 = np.array([[0, 1], [1, 0]], dtype=complex)
    sigma2 = np.array([[0, -1j], [1j, 0]], dtype=complex)
    sigma3 = np.array([[1, 0], [0, -1]], dtype=complex)
    gamma_plus = [np.block([[eye, zero], [zero, -eye]])]
    for sigma in (sigma1, sigma2, sigma3):
        gamma_plus.append(np.block([[zero, sigma], [-sigma, zero]]))
    return tuple(1j * matrix for matrix in gamma_plus)


def gamma5_matrix():
    gamma = dirac_gamma_matrices_minus_plus()
    return 1j * gamma[0] @ gamma[1] @ gamma[2] @ gamma[3]


def dirac_adjoint(spinor):
    """Use bar(psi)=i psi^dagger gamma^0."""
    psi = np.asarray(spinor, dtype=complex)
    if psi.shape != (4,) or not np.all(np.isfinite(psi)):
        raise ValueError("spinor must be a finite complex 4-vector")
    return 1j * psi.conj() @ dirac_gamma_matrices_minus_plus()[0]


def dirac_axial_current(spinor, *, reality_tolerance=1e-10):
    """Return J5^a=i bar(psi) gamma5 gamma^a psi as a real contravector."""
    psi = np.asarray(spinor, dtype=complex)
    bar = dirac_adjoint(psi)
    gamma = dirac_gamma_matrices_minus_plus()
    gamma5 = gamma5_matrix()
    values = np.array([1j * bar @ gamma5 @ ga @ psi for ga in gamma])
    if np.max(np.abs(values.imag)) > reality_tolerance:
        raise ArithmeticError("axial current is not real within tolerance")
    return values.real


def levi_civita_symbol():
    """Alternating symbol (not metric-raised tensor), epsilon_0123=+1."""
    eps = np.zeros((4, 4, 4, 4), dtype=float)
    for p in itertools.permutations(range(4)):
        inversions = sum(p[i] > p[j] for i in range(4) for j in range(i + 1, 4))
        eps[p] = -1.0 if inversions % 2 else 1.0
    return eps


def _validate_metadata(convention_id, units_id):
    if not isinstance(convention_id, str) or not convention_id.strip():
        raise ValueError("convention_id must be explicit and nonempty")
    if not isinstance(units_id, str) or not units_id.strip():
        raise ValueError("units_id must be explicit and nonempty")


def axial_torsion_from_current(
    axial_current, gravitational_coupling, *,
    cartan_coefficient=DERIVED_CARTAN_COEFFICIENT,
    convention_id=DIRAC_EC_CONVENTION_ID,
    units_id="geometrized: kappa*J5 has torsion units",
):
    """Build T_abc=c_T*kappa*epsilon_abcd*J5^d in the fixed EC convention."""
    _validate_metadata(convention_id, units_id)
    j5 = np.asarray(axial_current, dtype=float)
    if j5.shape != (4,) or not np.all(np.isfinite(j5)):
        raise ValueError("axial_current must be a finite contravariant 4-vector")
    if not np.isfinite(gravitational_coupling) or gravitational_coupling <= 0:
        raise ValueError("gravitational_coupling must be finite and positive")
    if not np.isfinite(cartan_coefficient):
        raise ValueError("cartan_coefficient must be finite")
    torsion = (float(cartan_coefficient) * float(gravitational_coupling)
               * np.einsum("abcd,d->abc", levi_civita_symbol(), j5))
    anti = max(
        float(np.linalg.norm(torsion + np.swapaxes(torsion, 0, 1))),
        float(np.linalg.norm(torsion + np.swapaxes(torsion, 1, 2))),
    )
    return {
        "torsion": torsion, "axial_current": j5.copy(),
        "gravitational_coupling": float(gravitational_coupling),
        "cartan_coefficient": float(cartan_coefficient),
        "metric_signature": METRIC_SIGNATURE,
        "convention_id": convention_id, "units_id": units_id,
        "totally_antisymmetric": bool(anti < 1e-13),
        "antisymmetry_residual": anti,
        "status": "ACTION_DERIVED_MINIMAL_EC_ALGEBRAIC_TORSION",
    }


def axial_torsion_from_dirac_spinor(
    spinor, gravitational_coupling, *,
    convention_id=DIRAC_EC_CONVENTION_ID,
    units_id="geometrized: kappa*J5 has torsion units",
):
    """Derive the axial source from a spinor and solve the Cartan relation."""
    current = dirac_axial_current(spinor)
    result = axial_torsion_from_current(
        current, gravitational_coupling,
        cartan_coefficient=DERIVED_CARTAN_COEFFICIENT,
        convention_id=convention_id, units_id=units_id,
    )
    result["spinor_axial_current"] = current.copy()
    return result


def cartan_algebraic_residual(
    torsion, axial_current, gravitational_coupling, *,
    cartan_coefficient=DERIVED_CARTAN_COEFFICIENT,
    convention_id=DIRAC_EC_CONVENTION_ID,
    units_id="geometrized: kappa*J5 has torsion units",
):
    _validate_metadata(convention_id, units_id)
    t = np.asarray(torsion, dtype=float)
    j5 = np.asarray(axial_current, dtype=float)
    if t.shape != (4, 4, 4) or j5.shape != (4,):
        raise ValueError("torsion must have shape (4,4,4), current shape (4,)")
    expected = (float(cartan_coefficient) * float(gravitational_coupling)
                * np.einsum("abcd,d->abc", levi_civita_symbol(), j5))
    residual = t - expected
    norm = float(np.linalg.norm(residual))
    return {"residual_tensor": residual, "residual_norm": norm,
            "equation_satisfied": bool(norm < 1e-12),
            "convention_id": convention_id, "units_id": units_id}


def lower_torsion_indices(torsion):
    t = np.asarray(torsion, dtype=float)
    if t.shape != (4, 4, 4):
        raise ValueError("torsion must have shape (4,4,4)")
    eta = minkowski_metric()
    return np.einsum("ai,bj,ck,ijk->abc", eta, eta, eta, t)


def torsion_quadratic_invariant(torsion):
    t = np.asarray(torsion, dtype=float)
    if t.shape != (4, 4, 4):
        raise ValueError("torsion must have shape (4,4,4)")
    return float(np.einsum("abc,abc->", lower_torsion_indices(t), t))


def audit_algebraic_axial_torsion(
    axial_current=(0.2, -0.1, 0.3, 0.4), gravitational_coupling=1.0, *,
    cartan_coefficient=DERIVED_CARTAN_COEFFICIENT,
    convention_id=DIRAC_EC_CONVENTION_ID,
    units_id="geometrized: kappa*J5 has torsion units",
):
    built = axial_torsion_from_current(
        axial_current, gravitational_coupling,
        cartan_coefficient=cartan_coefficient,
        convention_id=convention_id, units_id=units_id,
    )
    check = cartan_algebraic_residual(
        built["torsion"], axial_current, gravitational_coupling,
        cartan_coefficient=cartan_coefficient,
        convention_id=convention_id, units_id=units_id,
    )
    return {
        "schema": "RGS_EC_DIRAC_ACTION_DERIVED_TORSION_V2",
        "metric_signature": METRIC_SIGNATURE, "convention_id": convention_id,
        "units_id": units_id, "cartan_coefficient": float(cartan_coefficient),
        "torsion_antisymmetry_residual": built["antisymmetry_residual"],
        "cartan_equation_residual": check["residual_norm"],
        "torsion_quadratic_invariant": torsion_quadratic_invariant(built["torsion"]),
        "algebraic_equation_satisfied": check["equation_satisfied"],
        "propagating_torsion_dofs": 0,
    }
