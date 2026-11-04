"""Minimal Einstein-Cartan-Dirac effective axial four-fermion interaction.

Convention bridge:
  A^a = bar(psi) gamma^5 gamma^a psi
  J5^a = i A^a  (real current used by rgs_ec_algebraic_torsion)
The standard minimal-EC result L4f=-(3 kappa/16) A_a A^a therefore becomes
L4f=+(3 kappa/16) J5_a J5^a in this real-current convention.
This module keeps the convention conversion explicit.
"""
from __future__ import annotations
import numpy as np

from src.rgs_ec_algebraic_torsion import (
    DIRAC_EC_CONVENTION_ID, minkowski_metric,
)

RAW_AXIAL_BILINEAR_CONVENTION = "A^a=bar(psi) gamma5 gamma^a psi"
REAL_AXIAL_CURRENT_CONVENTION = "J5^a=i*A^a; real contravector"
EFFECTIVE_ACTION_CONVENTION_ID = (
    "minimal-EC-Dirac;eta(-+++);"
    "L4f=-3kappa/16*A_a*A^a=+3kappa/16*J5_a*J5^a"
)


def _validate_current(current):
    value = np.asarray(current, dtype=float)
    if value.shape != (4,) or not np.all(np.isfinite(value)):
        raise ValueError("axial current must be a finite real 4-vector")
    return value


def minkowski_inner(left, right):
    a = _validate_current(left)
    b = _validate_current(right)
    return float(a @ minkowski_metric() @ b)


def raw_to_real_axial_current(raw_current):
    """Convert A=bar(psi)gamma5 gamma psi to J5=iA; A may be complex."""
    raw = np.asarray(raw_current, dtype=complex)
    if raw.shape != (4,) or not np.all(np.isfinite(raw)):
        raise ValueError("raw axial bilinear must be a finite complex 4-vector")
    real = 1j * raw
    if np.max(np.abs(real.imag)) > 1e-10:
        raise ArithmeticError("i times raw axial bilinear is not real")
    return real.real


def real_to_raw_axial_current(real_current):
    """Inverse convention map A=-i J5."""
    return -1j * _validate_current(real_current)


def effective_four_fermion_coefficient(gravitational_coupling):
    """Coefficient multiplying J5_a J5^a for the declared real-current convention."""
    if not np.isfinite(gravitational_coupling) or gravitational_coupling <= 0:
        raise ValueError("gravitational_coupling must be finite and positive")
    return 3.0 * float(gravitational_coupling) / 16.0


def effective_four_fermion_density(axial_current, gravitational_coupling):
    """Return L_4f=(3 kappa/16) J5_a J5^a in the fixed convention."""
    current = _validate_current(axial_current)
    coefficient = effective_four_fermion_coefficient(gravitational_coupling)
    return coefficient * minkowski_inner(current, current)


def multi_flavour_axial_current(currents):
    """Sum species currents before squaring; cross-flavour terms are retained."""
    values = np.asarray(currents, dtype=float)
    if values.ndim != 2 or values.shape[1] != 4 or not np.all(np.isfinite(values)):
        raise ValueError("currents must be a finite array with shape (n_flavours,4)")
    return np.sum(values, axis=0)


def multi_flavour_effective_density(currents, gravitational_coupling):
    total = multi_flavour_axial_current(currents)
    return {
        "total_axial_current": total,
        "density": effective_four_fermion_density(total, gravitational_coupling),
        "coefficient": effective_four_fermion_coefficient(gravitational_coupling),
        "includes_cross_flavour_terms": True,
    }


def effective_eom_source(axial_current, gravitational_coupling):
    """Derivative dL4f/dJ5^a = (3kappa/8) J5_a."""
    current = _validate_current(axial_current)
    coefficient = effective_four_fermion_coefficient(gravitational_coupling)
    return 2.0 * coefficient * (minkowski_metric() @ current)


def audit_four_fermion_reduction(axial_current, gravitational_coupling):
    current = _validate_current(axial_current)
    raw = real_to_raw_axial_current(current)
    raw_norm = float(np.real(raw @ minkowski_metric() @ raw))
    real_norm = minkowski_inner(current, current)
    coefficient_real = effective_four_fermion_coefficient(gravitational_coupling)
    coefficient_raw = -3.0 * float(gravitational_coupling) / 16.0
    density_real = coefficient_real * real_norm
    density_raw = coefficient_raw * raw_norm
    return {
        "schema": "RGS_EC_DIRAC_FOUR_FERMION_AUDIT_V1",
        "metric_signature": "(-,+,+,+)",
        "source_convention_id": DIRAC_EC_CONVENTION_ID,
        "effective_action_convention_id": EFFECTIVE_ACTION_CONVENTION_ID,
        "raw_axial_bilinear_convention": RAW_AXIAL_BILINEAR_CONVENTION,
        "real_axial_current_convention": REAL_AXIAL_CURRENT_CONVENTION,
        "coefficient_for_real_J5_squared": coefficient_real,
        "coefficient_for_raw_A_squared": coefficient_raw,
        "real_current_norm": real_norm,
        "raw_bilinear_norm": raw_norm,
        "density_real_convention": density_real,
        "density_raw_convention": density_raw,
        "convention_conversion_residual": abs(density_real-density_raw),
        "eom_source": effective_eom_source(current, gravitational_coupling),
        "torsion_is_algebraic": True,
        "propagating_torsion_dofs": 0,
    }
