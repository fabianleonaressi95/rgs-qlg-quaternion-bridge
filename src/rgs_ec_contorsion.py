"""Index-explicit contorsion reconstruction for minimal EC Dirac torsion.

Input torsion is T_{abc} with all three internal indices lowered, using
eta=(-,+,+,+) and epsilon_0123=+1. The spin-connection contorsion convention
is K^a{}_b = K^a{}_{bc} e^c and T^a=K^a{}_b wedge e^b. For totally
antisymmetric T_{abc}, this convention gives K_{abc}=-T_{abc}/2.
"""
from __future__ import annotations
import numpy as np


def minkowski_metric():
    return np.diag([-1.0, 1.0, 1.0, 1.0])


def _tensor3(value, name):
    array = np.asarray(value, dtype=float)
    if array.shape != (4, 4, 4) or not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must be a finite (4,4,4) tensor")
    return array


def total_antisymmetry_residual(tensor):
    t = _tensor3(tensor, "tensor")
    return float(max(
        np.linalg.norm(t + np.swapaxes(t, 0, 1)),
        np.linalg.norm(t + np.swapaxes(t, 1, 2)),
        np.linalg.norm(t + np.swapaxes(t, 0, 2)),
    ))


def contorsion_from_dirac_torsion(torsion_lower, *, tolerance=1e-11):
    """Return K_abc=-T_abc/2 for fully antisymmetric torsion."""
    t = _tensor3(torsion_lower, "torsion_lower")
    residual = total_antisymmetry_residual(t)
    if residual > tolerance:
        raise ValueError("the -T/2 shortcut requires totally antisymmetric torsion")
    k = -0.5 * t
    return {
        "contorsion_lower": k,
        "torsion_antisymmetry_residual": residual,
        "connection_antisymmetry_residual": float(
            np.linalg.norm(k + np.swapaxes(k, 0, 1))
        ),
        "convention": "K^a_b=K^a_bc e^c; T^a=K^a_b wedge e^b",
    }


def reconstruct_torsion_from_contorsion(contorsion_lower):
    """Reconstruct T_abc=K_acb-K_abc from the declared wedge convention."""
    k = _tensor3(contorsion_lower, "contorsion_lower")
    return np.swapaxes(k, 1, 2) - k


def raise_torsion_indices(torsion_lower):
    """Raise all indices of T_abc with eta=(-,+,+,+)."""
    t = _tensor3(torsion_lower, "torsion_lower")
    eta = minkowski_metric()
    return np.einsum("ai,bj,ck,ijk->abc", eta, eta, eta, t)


def axial_current_norm(axial_current_upper):
    j = np.asarray(axial_current_upper, dtype=float)
    if j.shape != (4,) or not np.all(np.isfinite(j)):
        raise ValueError("axial_current_upper must be a finite 4-vector")
    return float(j @ minkowski_metric() @ j)


def torsion_quadratic_invariant(torsion_lower):
    """Compute T_abc T^abc, with input indices explicitly all lowered."""
    t = _tensor3(torsion_lower, "torsion_lower")
    return float(np.einsum("abc,abc->", t, raise_torsion_indices(t)))


def audit_dirac_torsion_identities(torsion_lower, axial_current_upper, kappa):
    """Check contorsion reconstruction and T^2=-(3/2)kappa^2 J5^2."""
    if not np.isfinite(kappa) or kappa <= 0:
        raise ValueError("kappa must be finite and positive")
    t = _tensor3(torsion_lower, "torsion_lower")
    j = np.asarray(axial_current_upper, dtype=float)
    if j.shape != (4,) or not np.all(np.isfinite(j)):
        raise ValueError("axial_current_upper must be a finite 4-vector")
    built = contorsion_from_dirac_torsion(t)
    k = built["contorsion_lower"]
    reconstructed = reconstruct_torsion_from_contorsion(k)
    reconstruction_residual = float(np.linalg.norm(reconstructed - t))
    t2 = torsion_quadratic_invariant(t)
    j2 = axial_current_norm(j)
    expected_t2 = -1.5 * float(kappa)**2 * j2
    invariant_residual = abs(t2 - expected_t2)
    return {
        "contorsion_reconstruction_residual": reconstruction_residual,
        "connection_antisymmetry_residual": built["connection_antisymmetry_residual"],
        "torsion_quadratic_invariant": t2,
        "axial_current_norm": j2,
        "expected_torsion_invariant": expected_t2,
        "torsion_invariant_identity_residual": invariant_residual,
        "torsion_fully_antisymmetric": built["torsion_antisymmetry_residual"] < 1e-11,
    }
