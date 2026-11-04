"""Convention-explicit algebraic axial-torsion response benchmark.

Finite-dimensional constitutive benchmark, NOT a full Einstein-Cartan-Dirac
solver. The convention-dependent coefficient remains an explicit input.

Model relation:
    T_{abc} = c_T * kappa * epsilon_{abcd} * J5^d

The default c_T=-1/4 is a plumbing-test placeholder. Match it to a fully
specified action, epsilon convention and spin-current normalization before
physical use.
"""
from __future__ import annotations
import itertools
import numpy as np

METRIC_SIGNATURE = "(-,+,+,+)"
DEFAULT_CONVENTION_ID = "eta-(-+++)_eps0123plus_J5-upper_explicit-cT"


def minkowski_metric():
    return np.diag([-1.0, 1.0, 1.0, 1.0])


def levi_civita_symbol():
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
    axial_current, gravitational_coupling, *, cartan_coefficient=-0.25,
    convention_id=DEFAULT_CONVENTION_ID,
    units_id="geometrized: kappa*J5 has torsion units",
):
    """Build a totally antisymmetric response from a contravariant axial current."""
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
        "status": "PARAMETERIZED_ALGEBRAIC_RESPONSE_NOT_FULL_EC_DERIVATION",
    }


def cartan_algebraic_residual(
    torsion, axial_current, gravitational_coupling, *,
    cartan_coefficient, convention_id, units_id,
):
    """Evaluate the residual of the explicitly parameterized Cartan relation."""
    _validate_metadata(convention_id, units_id)
    t = np.asarray(torsion, dtype=float)
    j5 = np.asarray(axial_current, dtype=float)
    if t.shape != (4, 4, 4) or j5.shape != (4,):
        raise ValueError("torsion must have shape (4,4,4), current shape (4,)")
    expected = (float(cartan_coefficient) * float(gravitational_coupling)
                * np.einsum("abcd,d->abc", levi_civita_symbol(), j5))
    residual = t - expected
    norm = float(np.linalg.norm(residual))
    return {
        "residual_tensor": residual, "residual_norm": norm,
        "equation_satisfied": bool(norm < 1e-12),
        "convention_id": convention_id, "units_id": units_id,
    }


def lower_torsion_indices(torsion):
    """Lower all indices with eta=(-,+,+,+)."""
    t = np.asarray(torsion, dtype=float)
    if t.shape != (4, 4, 4):
        raise ValueError("torsion must have shape (4,4,4)")
    eta = minkowski_metric()
    return np.einsum("ai,bj,ck,ijk->abc", eta, eta, eta, t)


def torsion_quadratic_invariant(torsion):
    """Return T_abc T^abc for a contravariant rank-three tensor."""
    t = np.asarray(torsion, dtype=float)
    if t.shape != (4, 4, 4):
        raise ValueError("torsion must have shape (4,4,4)")
    return float(np.einsum("abc,abc->", lower_torsion_indices(t), t))


def audit_algebraic_axial_torsion(
    axial_current=(0.2, -0.1, 0.3, 0.4), gravitational_coupling=1.0, *,
    cartan_coefficient=-0.25, convention_id=DEFAULT_CONVENTION_ID,
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
        "schema": "RGS_EC_ALGEBRAIC_AXIAL_TORSION_BENCHMARK_V1",
        "metric_signature": METRIC_SIGNATURE, "convention_id": convention_id,
        "units_id": units_id, "cartan_coefficient": float(cartan_coefficient),
        "torsion_antisymmetry_residual": built["antisymmetry_residual"],
        "cartan_equation_residual": check["residual_norm"],
        "torsion_quadratic_invariant": torsion_quadratic_invariant(built["torsion"]),
        "algebraic_equation_satisfied": check["equation_satisfied"],
        "propagating_torsion_dofs": 0,
        "scope": "Algebraic axial response only; no kinetic term, cosmological averaging, TEGR map or observational constraint.",
    }
