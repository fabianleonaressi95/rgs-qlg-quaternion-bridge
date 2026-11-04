"""Convention-safe EC Dirac spin ensemble and cosmological adapter.

The returned rank-three object is the Cartan-normalized source Sigma_abc,
defined by T_abc = kappa Sigma_abc. It is deliberately NOT labelled the
canonical Noether spin current: that identification requires a separately
declared variational normalization.
"""
from __future__ import annotations
import numpy as np

from src.rgs_ec_algebraic_torsion import (
    DIRAC_EC_CONVENTION_ID, levi_civita_symbol, minkowski_metric,
)
from src.rgs_ec_dirac_four_fermion import effective_four_fermion_coefficient


SPIN_ENSEMBLE_SCHEMA = "RGS_EC_DIRAC_SPIN_ENSEMBLE_V1"
CARTAN_SOURCE_CONVENTION = (
    "Sigma_abc=-1/2*epsilon_abcd*J5^d; T_abc=kappa*Sigma_abc"
)


def _metadata(convention_id, units_id, derivation_id):
    for name, value in (("convention_id", convention_id),
                        ("units_id", units_id),
                        ("derivation_id", derivation_id)):
        if not isinstance(value, str) or not value.strip():
            raise ValueError(f"{name} must be explicit and nonempty")


def _current_array(currents):
    x = np.asarray(currents, dtype=float)
    if x.ndim != 2 or x.shape[1] != 4 or x.shape[0] < 1:
        raise ValueError("currents must have shape (n>=1,4)")
    if not np.all(np.isfinite(x)):
        raise ValueError("currents must be finite")
    return x


def cartan_normalized_spin_source(axial_current):
    """Return Sigma_abc=-1/2 epsilon_abcd J5^d in the declared frame."""
    j = np.asarray(axial_current, dtype=float)
    if j.shape != (4,) or not np.all(np.isfinite(j)):
        raise ValueError("axial_current must be a finite contravariant 4-vector")
    return -0.5 * np.einsum("abcd,d->abc", levi_civita_symbol(), j)


def torsion_from_cartan_source(axial_current, kappa):
    if not np.isfinite(kappa) or kappa <= 0:
        raise ValueError("kappa must be finite and positive")
    sigma = cartan_normalized_spin_source(axial_current)
    return float(kappa) * sigma


def _weighted_moments(currents, weights):
    x = _current_array(currents)
    if weights is None:
        w = np.full(x.shape[0], 1.0 / x.shape[0])
    else:
        w = np.asarray(weights, dtype=float)
        if w.shape != (x.shape[0],) or not np.all(np.isfinite(w)):
            raise ValueError("weights must be a finite vector matching the samples")
        if np.any(w < 0) or np.sum(w) <= 0:
            raise ValueError("weights must be nonnegative with positive sum")
        w = w / np.sum(w)
    mean = np.einsum("n,na->a", w, x)
    centered = x - mean
    covariance = np.einsum("n,na,nb->ab", w, centered, centered)
    second = np.einsum("n,na,nb->ab", w, x, x)
    eta = minkowski_metric()
    mean_sq = float(mean @ eta @ mean)
    covariance_lorentz_trace = float(np.einsum("ab,ab->", eta, covariance))
    second_lorentz_trace = float(np.einsum("ab,ab->", eta, second))
    decomposition_residual = abs(second_lorentz_trace -
                                 mean_sq - covariance_lorentz_trace)
    return {
        "sample_count": int(x.shape[0]), "normalized_weights": w,
        "mean_axial_current": mean, "covariance": covariance,
        "second_moment": second, "mean_current_lorentz_norm": mean_sq,
        "fluctuation_lorentz_trace": covariance_lorentz_trace,
        "mean_of_lorentz_norm_squared": second_lorentz_trace,
        "moment_decomposition_residual": decomposition_residual,
    }


def spatial_isotropy_residual(second_moment, mean_current=None):
    """Test FLRW isotropy using M^{ab}=<J^a J^b>, not covariance alone."""
    m = np.asarray(second_moment, dtype=float)
    if m.shape != (4, 4) or not np.all(np.isfinite(m)):
        raise ValueError("second_moment must be a finite (4,4) matrix")
    m = 0.5 * (m + m.T)
    spatial = m[1:, 1:]
    q = float(np.trace(spatial) / 3.0)
    target = q * np.eye(3)
    scale = max(float(np.linalg.norm(spatial)), 1.0)
    mean = np.zeros(4) if mean_current is None else np.asarray(mean_current, dtype=float)
    if mean.shape != (4,) or not np.all(np.isfinite(mean)):
        raise ValueError("mean_current must be a finite 4-vector")
    # M^{ab}=<J^a J^b> contains both connected fluctuations and the
    # coherent mean: Cov^{ab}=M^{ab}-<J^a><J^b>.  Keep the two
    # time-space diagnostics distinct; FLRW isotropy below is tested against
    # the full second moment and the mean spatial polarization separately.
    covariance = m - np.outer(mean, mean)
    covariance = 0.5 * (covariance + covariance.T)
    time_space_covariance_norm = float(np.linalg.norm(covariance[0, 1:]))
    time_space_second_moment_norm = float(np.linalg.norm(m[0, 1:]))
    return {
        "spatial_second_moment_isotropic_amplitude": q,
        "spatial_isotropy_residual": float(np.linalg.norm(spatial-target)/scale),
        "time_space_covariance_norm": time_space_covariance_norm,
        "time_space_second_moment_norm": time_space_second_moment_norm,
        "mean_spatial_polarization_norm": float(np.linalg.norm(mean[1:])),
        "flrw_isotropy_residual": float(np.sqrt(
            (np.linalg.norm(spatial-target)/scale)**2
            + time_space_second_moment_norm**2 / max(np.linalg.norm(m), 1.0)**2
            + np.linalg.norm(mean[1:])**2 / max(np.linalg.norm(mean), 1.0)**2
        )),
    }


def audit_spin_ensemble(currents, kappa, *, weights=None,
                        convention_id=DIRAC_EC_CONVENTION_ID,
                        units_id="geometrized: kappa*J5 has torsion units",
                        derivation_id="minimal-EC-Dirac-Cartan"):
    _metadata(convention_id, units_id, derivation_id)
    if not np.isfinite(kappa) or kappa <= 0:
        raise ValueError("kappa must be finite and positive")
    m = _weighted_moments(currents, weights)
    mean = m["mean_axial_current"]
    mean_source = cartan_normalized_spin_source(mean)
    mean_torsion = float(kappa) * mean_source
    kappa2 = float(kappa)**2
    torsion_invariant_mean = -1.5 * kappa2 * m["mean_current_lorentz_norm"]
    torsion_invariant_average = -1.5 * kappa2 * m["mean_of_lorentz_norm_squared"]
    contact_coefficient = effective_four_fermion_coefficient(kappa)
    contact_average = contact_coefficient * m["mean_of_lorentz_norm_squared"]
    contact_of_mean = contact_coefficient * m["mean_current_lorentz_norm"]
    isotropy = spatial_isotropy_residual(m["second_moment"], mean)
    return {
        "schema": SPIN_ENSEMBLE_SCHEMA,
        "metric_signature": "(-,+,+,+)",
        "convention_id": convention_id, "units_id": units_id,
        "derivation_id": derivation_id,
        "cartan_source_convention": CARTAN_SOURCE_CONVENTION,
        **m, **isotropy,
        "mean_cartan_source": mean_source,
        "mean_torsion": mean_torsion,
        "mean_torsion_quadratic_invariant": torsion_invariant_mean,
        "ensemble_mean_torsion_quadratic_invariant": torsion_invariant_average,
        "four_fermion_coefficient": contact_coefficient,
        "ensemble_mean_four_fermion_density": contact_average,
        "four_fermion_density_of_mean_current": contact_of_mean,
        "fluctuation_contact_contribution": contact_average-contact_of_mean,
        "cartan_equation_residual": float(np.linalg.norm(
            mean_torsion - torsion_from_cartan_source(mean, kappa))),
        "rho_spin": None,
        "rho_spin_status": "NOT_INFERRED: requires metric variation and matter-state EoS",
        "canonical_spin_current_identification": "NOT_ASSUMED",
        "cosmological_beta": None,
        "cosmological_beta_status": "NOT_DERIVED",
        "torsion_propagating_dofs": 0,
    }


def build_cosmological_spin_state(currents, kappa, *, scale_factor,
                                  hubble_rate, energy_density, pressure,
                                  weights=None, convention_id=DIRAC_EC_CONVENTION_ID,
                                  units_id, derivation_id):
    """Build an adapter record; does not invent an effective rho_spin or beta."""
    for name, value in (("scale_factor", scale_factor),
                        ("hubble_rate", hubble_rate),
                        ("energy_density", energy_density),
                        ("pressure", pressure)):
        if not np.isfinite(value):
            raise ValueError(f"{name} must be finite")
    if scale_factor <= 0:
        raise ValueError("scale_factor must be positive")
    _metadata(convention_id, units_id, derivation_id)
    audit = audit_spin_ensemble(
        currents, kappa, weights=weights, convention_id=convention_id,
        units_id=units_id, derivation_id=derivation_id,
    )
    return {
        "schema": "RGS_SPIN_TORSION_COSMOLOGY_ADAPTER_V1",
        "scale_factor": float(scale_factor), "hubble_rate": float(hubble_rate),
        "matter_energy_density": float(energy_density),
        "matter_pressure": float(pressure),
        "units_id": units_id, "convention_id": convention_id,
        "derivation_id": derivation_id, "spin_ensemble": audit,
        "effective_energy_density": None,
        "effective_pressure": None,
        "torsion_lambda_beta": None,
        "status": "SPIN_SOURCE_READY_COSMOLOGICAL_REDUCTION_OPEN",
    }
