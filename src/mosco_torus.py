"""Fourier-form prototype for Mosco convergence on a flat periodic torus.

This is a regular-lattice reference family, not an RGS-Q6 refinement theorem.
"""
from __future__ import annotations
import numpy as np
from src.wilson_refinement import naive_symbol_energy, wilson_symbol


def euclidean_wilson_form_symbol(momentum, spacing: float, r: float = 1.0) -> float:
    """Symbol of D_E^dagger D_E = E0(p)^2 + W(p)^2."""
    e=naive_symbol_energy(momentum,spacing)
    w=wilson_symbol(momentum,spacing,r)
    return float(e*e+w*w)


def torus_mosco_sampling_audit(
    sizes=(16,32,64,128,256), dimension=3, physical_length=1.0,
    mode=(1,0,0), r0=1.0, r_exponent=0.0
):
    d=int(dimension); L=float(physical_length); m=np.asarray(mode,dtype=int)
    alpha=float(r_exponent); r0=float(r0)
    if d<1 or d>8 or m.shape!=(d,):
        raise ValueError("mode must have one integer component per dimension")
    if L<=0 or r0<=0 or not np.isfinite(L+r0+alpha):
        raise ValueError("length and r0 must be positive and all inputs finite")
    p=2*np.pi*m/L; q_limit=4*float(np.dot(p,p))
    rows=[]
    for n0 in sizes:
        n=int(n0)
        if n<4: raise ValueError("each grid size must be at least 4")
        a=L/n; r=r0*(a/L)**alpha
        q=euclidean_wilson_form_symbol(p,a,r)
        rows.append({
            "N_per_axis":n,"spacing":a,"r":r,
            "fixed_mode_form_value":q,
            "continuum_form_value":q_limit,
            "relative_form_error":abs(q-q_limit)/q_limit if q_limit else 0.0,
            "one_direction_doubler_form_value":(2*r/a)**2,
            "minimum_nonphysical_corner_mass":2*r/a
        })
    return {
        "audit":"FLAT_TORUS_EUCLIDEAN_WILSON_FORM_REFINEMENT",
        "dimension":d,"physical_length":L,"mode":m.tolist(),
        "continuum_form_value":q_limit,"r0":r0,"r_exponent":alpha,
        "mosco_low_mode_consistency_condition":"alpha > -1",
        "doubler_spectral_decoupling_condition":"alpha < 1",
        "joint_consistency_and_doubler_window":"-1 < alpha < 1",
        "fixed_mode_form_converges":alpha>-1,
        "doubler_mass_diverges":alpha<1,
        "refinements":rows,
        "theorem_scope":"Fourier-embedded flat periodic torus family; not the curved RGS-Q6 graph"
    }
