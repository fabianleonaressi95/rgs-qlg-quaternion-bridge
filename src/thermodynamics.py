"""Finite-spectrum thermodynamics for the RGS-Q6/Mobius fermion boundary."""
from __future__ import annotations
import numpy as np

def _spectrum(values,name):
    x=np.asarray(values,dtype=float)
    if x.ndim!=1 or x.size==0 or not np.all(np.isfinite(x)):
        raise ValueError(f"{name} must be a non-empty finite 1D array")
    return x

def fermion_partition(energies,temperature,degeneracies=None):
    """Grand-canonical fermion determinant and filled-negative-sea normalization."""
    e=_spectrum(energies,"energies"); t=float(temperature)
    if not np.isfinite(t) or t<=0: raise ValueError("temperature must be positive and finite")
    d=np.ones_like(e) if degeneracies is None else np.asarray(degeneracies,dtype=float)
    if d.shape!=e.shape or not np.all(np.isfinite(d)) or np.any(d<=0):
        raise ValueError("degeneracies must be positive, finite and match energies")
    beta=1/t; x=beta*e; f=np.exp(-np.logaddexp(0.0,x))
    lograw=float(np.sum(d*np.logaddexp(0.0,-x)))
    evac=float(np.sum(d*np.minimum(e,0.0)))
    logrel=float(lograw+beta*evac); uraw=float(np.sum(d*e*f)); urel=uraw-evac
    entropy=float(lograw+beta*uraw)
    cv=float(beta*beta*np.sum(d*e*e*f*(1-f)))
    return {"ensemble":"independent_fermion_modes_mu_zero",
      "vacuum_convention":"all_negative_energy_levels_filled","temperature":t,"beta":beta,
      "mode_count":int(e.size),"total_degeneracy":float(np.sum(d)),"vacuum_energy":evac,
      "logZ_grand_canonical_raw":lograw,"Z_grand_canonical_raw":float(np.exp(min(lograw,700))),
      "logZ_vacuum_normalized":logrel,"Z_vacuum_normalized":float(np.exp(min(logrel,700))),
      "mean_occupation_sum":float(np.sum(d*f)),"energy_raw":uraw,
      "energy_above_filled_vacuum":float(urel),"entropy_thermodynamic":entropy,
      "free_energy_above_vacuum":float(-t*logrel),"heat_capacity":cv}

def modular_partition(correlation_eigenvalues,endpoint_tolerance=1e-12):
    """Gaussian modular partition and binary entanglement entropy."""
    z=_spectrum(correlation_eigenvalues,"correlation_eigenvalues"); tol=float(endpoint_tolerance)
    if not np.isfinite(tol) or tol<0 or tol>=.5: raise ValueError("invalid endpoint_tolerance")
    if np.any(z < -tol) or np.any(z > 1+tol):
        raise ValueError("correlation eigenvalues must lie in [0,1] within tolerance")
    z=np.clip(z,0,1); z[z<=tol]=0; z[z>=1-tol]=1
    terms=np.zeros_like(z); m0=z>0; m1=z<1
    terms[m0]-=z[m0]*np.log(z[m0]); terms[m1]-=(1-z[m1])*np.log1p(-z[m1])
    entropy=float(np.sum(terms)); n1=int(np.count_nonzero(z==1))
    logz=None if n1 else float(-np.sum(np.log1p(-z)))
    zmod=None if n1 else float(np.exp(min(logz,700)))
    eps=[None if v==0 or v==1 else float(np.log1p(-v)-np.log(v)) for v in z]
    return {"mode_count":int(z.size),"endpoint_tolerance":tol,
      "projected_zero_modes":int(np.count_nonzero(z==0)),"projected_unit_modes":n1,
      "logZ_modular":logz,"Z_modular":zmod,"entropy_entanglement":entropy,
      "modular_energies":eps,
      "status":"MODULAR_PARTITION_DIVERGENT_AT_UNIT_EIGENVALUE" if n1 else "FINITE_MODULAR_PARTITION"}

def compare_boundary_spectra(data,endpoint_tolerance=1e-12):
    """Compare naive/Wilson thermal and modular spectra from a JSON dictionary."""
    required=("T_Hawking","E_b_naive","E_b_wilson","zeta_naive","zeta_wilson")
    missing=[k for k in required if k not in data]
    if missing: raise ValueError(f"missing required spectrum fields: {missing}")
    t=float(data["T_Hawking"])
    tn=fermion_partition(data["E_b_naive"],t); tw=fermion_partition(data["E_b_wilson"],t)
    mn=modular_partition(data["zeta_naive"],endpoint_tolerance)
    mw=modular_partition(data["zeta_wilson"],endpoint_tolerance)
    reported={"naive":data.get("S_EE_naive"),"wilson":data.get("S_EE_wilson")}
    reconstructed={"naive":mn["entropy_entanglement"],"wilson":mw["entropy_entanglement"]}
    consistency={}
    for key in ("naive","wilson"):
        v=reported[key]
        consistency[key]=None if v is None else {"reported":float(v),
          "reconstructed":reconstructed[key],"absolute_residual":float(abs(float(v)-reconstructed[key]))}
    pn=data.get("hawking_power_naive"); pw=data.get("hawking_power_wilson"); power=None
    if pn is not None and pw is not None:
        power={"naive":float(pn),"wilson":float(pw),"wilson_to_naive_ratio":float(pw)/float(pn),
          "note":"Imported diagnostic; not derived without transmission coefficients and a flux prescription."}
    return {"audit":"RGS_MOBIUS_HAWKING_THERMODYNAMICS","temperature":t,"beta":1/t,
      "boundary_length":data.get("L_boundary"),"surface_gravity":data.get("kappa"),
      "naive":{"thermodynamics":tn,"modular":mn},"wilson":{"thermodynamics":tw,"modular":mw},
      "entanglement_entropy_consistency":consistency,"reported_hawking_power":power,
      "interpretation":{"thermal_entropy_and_entanglement_entropy_are_distinct":True,
        "modular_partition_is_not_an_independent_factor_of_thermal_partition":True,
        "hawking_temperature_is_an_input_scale_not_derived_by_this_module":True,
        "continuum_hawking_flux_and_gravitational_backreaction":"OPEN"}}
