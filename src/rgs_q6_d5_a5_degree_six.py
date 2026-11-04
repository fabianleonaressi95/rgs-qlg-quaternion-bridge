"""Degree-six A5-invariant coefficient: continuous S9 pushforward vs Q6 design."""
from __future__ import annotations
import math
import numpy as np
from src.rgs_icosahedral_q6_intertwiner import icosahedral_signed_axis_representation
from src.rgs_q6_d5_sphere_intertwiner import build_q6_d5_relative_direction_intertwiner
from src.rgs_q6_relative_sphere import build_q6_a5_relative_sphere_map

def _legendre_p6(x):
    x=np.asarray(x,dtype=float)
    return (231*x**6-315*x**4+105*x**2-5)/16

def _invariant_harmonic_values(directions,rotations,axis=None):
    n=np.asarray(directions,dtype=float)
    if n.ndim!=2 or n.shape[1]!=3: raise ValueError("directions must have shape (N,3)")
    u=np.asarray([1.,2.,3.] if axis is None else axis,dtype=float)
    if u.shape!=(3,) or not np.all(np.isfinite(u)) or np.linalg.norm(u)==0:
        raise ValueError("axis must be a finite nonzero 3-vector")
    u=u/np.linalg.norm(u)
    orbit_axes=np.asarray(rotations,dtype=float)@u
    return np.mean(_legendre_p6(n@orbit_axes.T),axis=1)

def _uniform_sphere_harmonic_rms(rotations,axis=None):
    mu,weights=np.polynomial.legendre.leggauss(16)
    phi=2*math.pi*np.arange(32)/32
    directions=np.asarray([[math.sqrt(max(0.,1-m*m))*math.cos(t),
        math.sqrt(max(0.,1-m*m))*math.sin(t),m] for m in mu for t in phi])
    values=_invariant_harmonic_values(directions,rotations,axis)
    return float(math.sqrt(np.sum(np.repeat(weights/(2*32),32)*values**2)))

def _continuous_direction_batch(rng,count,intertwiner):
    x=rng.standard_normal((count,10)); q,p=x[:,:5],x[:,5:]
    pairs=[(i,j) for i in range(5) for j in range(i+1,5)]
    wedge=np.column_stack([q[:,i]*p[:,j]-q[:,j]*p[:,i] for i,j in pairs])
    vectors=wedge@intertwiner.T
    norms=np.linalg.norm(vectors,axis=1)
    if np.any(norms<=np.finfo(float).tiny): raise ArithmeticError("sample hit zero locus")
    return vectors/norms[:,None]

def build_q6_d5_degree_six_coefficient_audit(samples_per_batch=4096,batches=16,seed=20261104):
    if not isinstance(samples_per_batch,int) or samples_per_batch<256: raise ValueError("samples_per_batch >=256")
    if not isinstance(batches,int) or batches<4: raise ValueError("batches >=4")
    if not isinstance(seed,int): raise ValueError("seed must be integer")
    rotations=np.asarray(icosahedral_signed_axis_representation()["rotations"],dtype=float)
    bridge=build_q6_d5_relative_direction_intertwiner()
    if not bridge["all_checks_passed"]: raise ArithmeticError("intertwiner audit failed")
    T=np.asarray(bridge["intertwiner"],dtype=float)
    directions=np.asarray(build_q6_a5_relative_sphere_map()["unique_directions"],dtype=float)
    if directions.shape!=(32,3): raise ArithmeticError(f"expected 32 directions, got {directions.shape}")
    axis=np.asarray([1.,2.,3.]); rms=_uniform_sphere_harmonic_rms(rotations,axis)
    if not np.isfinite(rms) or rms<=1e-14: raise ArithmeticError("degree-six invariant harmonic vanished")
    finite_raw=_invariant_harmonic_values(directions,rotations,axis)
    finite_coefficient=float(np.mean(finite_raw)/rms)
    rng=np.random.default_rng(seed); block_means=[]
    for _ in range(batches):
        n=_continuous_direction_batch(rng,samples_per_batch,T)
        block_means.append(float(np.mean(_invariant_harmonic_values(n,rotations,axis))))
    block_means=np.asarray(block_means)
    raw=float(np.mean(block_means)); se=float(np.std(block_means,ddof=1)/math.sqrt(batches))
    continuous=raw/rms; coefficient_se=se/rms; difference=continuous-finite_coefficient
    return {
      "audit":"Q6_D5_A5_DEGREE_SIX_COEFFICIENT_COMPARISON_V1",
      "harmonic":"A5 Reynolds average of P6(u dot n), u=(1,2,3)/sqrt(14)",
      "normalization":"L2 norm under uniform probability measure on S2",
      "harmonic_uniform_sphere_rms":rms,
      "finite_q6_design":{"direction_count":len(directions),"raw_harmonic_mean":float(np.mean(finite_raw)),
        "normalized_degree_six_coefficient":finite_coefficient},
      "continuous_s9_pushforward":{"samples":samples_per_batch*batches,"batches":batches,
        "samples_per_batch":samples_per_batch,"seed":seed,"raw_harmonic_mean":raw,
        "normalized_degree_six_coefficient":continuous,"standard_error":coefficient_se,
        "normal_95pct_interval":[continuous-1.96*coefficient_se,continuous+1.96*coefficient_se],
        "block_means_raw":block_means.tolist()},
      "difference_continuous_minus_finite":difference,"difference_standard_error":coefficient_se,
      "difference_z_score":difference/coefficient_se if coefficient_se>0 else None,
      "finite_value_outside_continuous_95pct_interval":bool(abs(difference)>1.96*coefficient_se),
      "checks":{"source_intertwiner_passed":bool(bridge["all_checks_passed"]),
        "finite_design_has_32_directions":len(directions)==32,"harmonic_rms_positive":rms>0,
        "continuous_estimate_finite":bool(np.isfinite(continuous)),
        "finite_estimate_finite":bool(np.isfinite(finite_coefficient)),"positive_sampling_error":coefficient_se>0},
      "all_checks_passed":bool(bridge["all_checks_passed"] and len(directions)==32 and rms>0
        and np.isfinite(continuous) and np.isfinite(finite_coefficient) and coefficient_se>0),
      "interpretation":"Seeded Monte Carlo estimate; interval quantifies sampling noise only. Overlap is not proof of equality.",
      "scope_limits":["not an exact continuous integral","no model uncertainty included",
        "compact S9 phase-shell only","not physical Liouville or global D5 Bergman measure",
        "no 3+1D limit or fine-structure-constant derivation"]}

