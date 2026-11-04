# Q6 x periodic Euclidean time: Wilson action and Gauss constraint

## Construction

Take the cubical Q6 complex with 64 vertices, 192 spatial links and 240
square faces, and form its periodic product with a time cycle C_Nt. At each
time slice there are 64 temporal links. For every spatial edge there is one
space-time plaquette per time step.

For a canonically oriented spatial edge e=(a,b), define spatial link angles
A_e(t) and temporal link angles A_0(v,t). The compact gauge transformations
are

    A_e(t) -> A_e(t) + lambda_b(t) - lambda_a(t)
    A_0(v,t) -> A_0(v,t) + lambda_v(t+1) - lambda_v(t)

with periodic t. The plaquette fluxes are

    F_ij(t) = D1 A(t)
    F_0e(t) = A_e(t) + A_0(b,t) - A_e(t+1) - A_0(a,t).

Both transform invariantly modulo 2 pi. The Wilson action is

    S = beta_s sum_{t,f}(1-cos F_f(t))
      + beta_t sum_{t,e}(1-cos F_0e(t)).

The implementation audits exact cell counts, gauge invariance of both
plaquette families, and gauge invariance of the action.

## Gauss law from the temporal link

Varying a temporal link A_0(v,t) gives

    dS/dA_0(v,t) = beta_t sum_e D0[e,v] sin(F_0e(t))
                 = beta_t (D0^T sin(F_0(t)))_v.

Thus the lattice Euler-Lagrange equation for the temporal link is the
nonlinear lattice Gauss constraint

    D0^T sin(F_0(t)) = 0.

This is not imposed by hand: it is the temporal-link field equation. A
central finite-difference test compares the action derivative with beta_t
times the computed Gauss residual.

## Dimensional caveat

Q6 is the six-dimensional cubical hypercube. Its product with one periodic
time direction is therefore a **6+1-dimensional Euclidean lattice**, not a
3+1-dimensional spacetime. This module establishes a gauge-invariant
temporal extension and a discrete Gauss constraint, but it does not identify
three physical spatial directions or derive a 3+1D continuum limit.

## Evidence ledger

**Algebraically established / numerically audited**
- periodic Q6 x C_Nt link and plaquette incidence;
- spatial and temporal U(1) plaquette gauge invariance;
- Wilson-action gauge invariance;
- temporal-link Euler-Lagrange equation equals the nonlinear lattice Gauss law;
- finite-difference verification of the A_0 variation.

**Still open**
- a physically justified 3D spatial subcomplex or dimensional reduction;
- deriving beta_s, beta_t and their anisotropic continuum scaling from RGS-PSU;
- Hamiltonian transfer matrix and electric momentum spectrum;
- coupling to the PSU anti-twin vierbein and deriving its Lorentzian continuation.

## Reproduction

    PYTHONPATH=. python -m pytest tests/test_rgs_q6_euclidean_time_wilson.py -q
    PYTHONPATH=. python -c "from src.rgs_q6_euclidean_time_wilson import audit_q6_euclidean_time_wilson; print(audit_q6_euclidean_time_wilson()['all_checks_passed'])"
