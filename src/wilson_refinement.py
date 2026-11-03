"""Analytic refinement audit for naive and Wilson fermions on periodic tori.

The symbol is the one induced by the repository's centered spatial operator:
  D_naive(p) = (2/a) sum_mu gamma_mu sin(p_mu a)
  W(p)       = (r(a)/a) sum_mu [1-cos(p_mu a)] I.
This module audits the universal Fourier mechanism; it does not identify
momenta on the curved finite Q6 graph.
"""
from __future__ import annotations
from itertools import product
import math
import numpy as np


def wilson_symbol(momentum, spacing: float, r: float = 1.0) -> float:
    p=np.asarray(momentum,dtype=float)
    a=float(spacing); r=float(r)
    if p.ndim!=1 or not np.all(np.isfinite(p)):
        raise ValueError("momentum must be a finite one-dimensional vector")
    if not np.isfinite(a) or a<=0:
        raise ValueError("spacing must be positive and finite")
    if not np.isfinite(r) or r<0:
        raise ValueError("r must be non-negative and finite")
    return float((r/a)*np.sum(1.0-np.cos(p*a)))


def naive_symbol_energy(momentum, spacing: float) -> float:
    p=np.asarray(momentum,dtype=float); a=float(spacing)
    if p.ndim!=1 or not np.all(np.isfinite(p)):
        raise ValueError("momentum must be a finite one-dimensional vector")
    if not np.isfinite(a) or a<=0:
        raise ValueError("spacing must be positive and finite")
    return float((2.0/a)*np.sqrt(np.sum(np.sin(p*a)**2)))


def wilson_brillouin_corners(dimension: int, spacing: float, r: float = 1.0):
    """Return all 2^d corners and their Wilson masses, grouped by pi-components."""
    d=int(dimension); a=float(spacing); r=float(r)
    if d<1 or d>8:
        raise ValueError("dimension must lie in 1..8")
    if a<=0 or not np.isfinite(a) or r<0 or not np.isfinite(r):
        raise ValueError("spacing must be positive and r non-negative")
    rows=[]
    for bits in product((0,1),repeat=d):
        k=int(sum(bits))
        p=np.asarray(bits,dtype=float)*(np.pi/a)
        mass=wilson_symbol(p,a,r)
        rows.append({"corner_bits":list(bits),"pi_directions":k,
                     "naive_energy":naive_symbol_energy(p,a),
                     "wilson_mass":mass,"expected_mass":2.0*r*k/a})
    return rows


def refinement_audit(
    sizes=(16,32,64,128), dimension: int=3, physical_length: float=1.0,
    mode=(1,0,0), r0: float=1.0, r_exponent: float=0.0
):
    """Track a fixed physical Fourier mode and all Brillouin corners as a->0.

    r(a)=r0*(a/L)^alpha. For a fixed physical momentum, W=O(a^(alpha+1));
    at a corner with k pi-components, W=2*k*r(a)/a. Thus decoupling and
    vanishing low-momentum correction coexist precisely for -1<alpha<1.
    """
    d=int(dimension); L=float(physical_length); r0=float(r0); alpha=float(r_exponent)
    m=np.asarray(mode,dtype=int)
    if d<1 or d>8 or m.shape!=(d,):
        raise ValueError("mode must have one integer component per dimension")
    if L<=0 or not np.isfinite(L) or r0<=0 or not np.isfinite(r0):
        raise ValueError("physical_length and r0 must be positive")
    if not np.isfinite(alpha):
        raise ValueError("r_exponent must be finite")
    if any(int(n)<0 for n in sizes) or any(int(n)<4 for n in sizes):
        raise ValueError("each refinement size must be at least 4")
    p=2.0*np.pi*m/L
    continuum_energy=2.0*float(np.linalg.norm(p))
    rows=[]
    for n0 in sizes:
        n=int(n0); a=L/n; r=r0*(a/L)**alpha
        naive=naive_symbol_energy(p,a)
        wilson=wilson_symbol(p,a,r)
        corners=wilson_brillouin_corners(d,a,r)
        positive=[c["wilson_mass"] for c in corners if c["pi_directions"]>0]
        rows.append({
            "N_per_axis":n,"spacing":a,"r":r,
            "fixed_mode_naive_energy":naive,
            "continuum_target_energy":continuum_energy,
            "fixed_mode_naive_relative_error":(
                abs(naive-continuum_energy)/continuum_energy if continuum_energy else 0.0),
            "fixed_mode_wilson_correction":wilson,
            "fixed_mode_wilson_correction_over_a":wilson/a,
            "one_direction_doubler_mass":2.0*r/a,
            "max_corner_doubler_mass":max(positive) if positive else 0.0,
            "corner_masses_by_pi_count":{
                str(k):2.0*k*r/a for k in range(1,d+1)
            }
        })
    low_vanishes=(alpha>-1.0)
    doublers_decouple=(alpha<1.0)
    return {
        "dimension":d,"physical_length":L,"fixed_integer_mode":m.tolist(),
        "physical_momentum":p.tolist(),"continuum_target_energy":continuum_energy,
        "r_scaling":"r(a)=r0*(a/L)^alpha","r0":r0,"r_exponent":alpha,
        "low_momentum_wilson_term_vanishes":low_vanishes,
        "doubler_mass_diverges":doublers_decouple,
        "simultaneous_scaling_window":bool(low_vanishes and doublers_decouple),
        "analytic_scaling_window":"-1 < alpha < 1",
        "refinements":rows,
        "corners_at_finest_refinement":wilson_brillouin_corners(d,L/int(sizes[-1]),r0*(1.0/int(sizes[-1]))**alpha)
    }


def hermitian_scalar_sum_branches(momentum, spacing: float, r: float = 1.0):
    """Eigenvalue branches of the nonstandard Hermitian sum D0+W."""
    w=wilson_symbol(momentum,spacing,r)
    e=naive_symbol_energy(momentum,spacing)
    return float(w-e),float(w+e)


def euclidean_wilson_singular_value(momentum, spacing: float, r: float = 1.0):
    """Singular-value magnitude of i*D0+W for a Clifford plane wave."""
    w=wilson_symbol(momentum,spacing,r)
    e=naive_symbol_energy(momentum,spacing)
    return float(np.hypot(e,w))


def scalar_sum_accidental_zero_angle(r: float = 1.0) -> float:
    """Return theta in (0,pi) where W(theta)-|D0(theta)| vanishes in 1D.

    For D0=2 sin(theta)/a and W=r(1-cos(theta))/a, tan(theta/2)=2/r.
    This zero shows why lifting the Nyquist corner alone is not sufficient
    when the Hermitian scalar Wilson Laplacian is simply added to D0.
    """
    r=float(r)
    if not np.isfinite(r) or r<=0:
        raise ValueError("r must be positive and finite")
    return float(2.0*np.arctan(2.0/r))
