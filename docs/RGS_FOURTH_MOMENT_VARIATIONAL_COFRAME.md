# RGS fourth-moment variational coframe sector

## Scope

Use the audited 64-site anisotropy R_abcd and define the local algebraic
3+1 density

    L_aniso = N sqrt(h) V(Q),   V(Q) = lambda/2 R_abcd Q^ab Q^cd.

Q_ab is an independent symmetric trace-free internal tensor. The spatial
coframe e^a_i is orientation-preserving and h_ij = delta_ab e^a_i e^b_j.
R is a material anisotropy tensor: its cubic stabilizer is physical data.
Under a change of internal frame, both R and Q must transform as tensors.
Holding R fixed while rotating the frame arbitrarily breaks symmetry to
the material stabilizer O_h.

## Exact first variations

Since Q and material components of R are independent of the coframe,

    delta sqrt(h) = sqrt(h) e_a^i delta e^a_i,

so

    dL/de^a_i = N sqrt(h) V(Q) e_a^i,
    dL/dN      = sqrt(h) V(Q),
    dV/dQ^ab   = lambda R_abcd Q^cd.

This is an algebraic constitutive source in the coframe equation, not the
complete gravitational Euler equation. TEGR/Einstein terms, matter,
constraints and any derivative couplings must also be varied.

## Why this does not yet produce torsion

The action contains no derivatives of e and therefore cannot itself
generate a Weitzenbock torsion contribution. Torsion remains

    T^a = de^a + omega^a_b wedge e^b.

A spin-current or torsion response requires a derivative coupling and a
consistent variation of the flat spin connection (or a specified
covariant teleparallel formulation). No ad-hoc map is inserted from the
fully symmetric R_abcd to torsion's antisymmetric index pair.

## Reproducibility and evidence

Run:
    PYTHONPATH=. python -m pytest tests/test_rgs_fourth_moment_variation.py -q

Central finite differences test coframe, lapse and trace-free Q
Euler derivatives. A covariance test rotates R and Q simultaneously.

- NUMERICALLY VERIFIED: finite-dimensional first variations and
  simultaneous internal-frame tensor covariance.
- CANDIDATE: algebraic fourth-moment constitutive sector.
- OPEN: complete coframe field equations, derivative/spin-connection
  coupling, teleparallel torsion dynamics, Q6/Yang-Mills coupling,
  and continuum/Mosco convergence.
