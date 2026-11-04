# EC Dirac contorsion and quadratic invariant audit

This module continues the fixed first-order convention in
EC_DIRAC_ACTION_DERIVATION_AND_CONVENTIONS.md. It takes T_abc with all
internal indices lowered, signature (-+++), and defines the spin-connection
contorsion by
    K^a_b = K^a_bc e^c,
    T^a = K^a_b wedge e^b.
The wedge product gives
    T_abc = K_acb - K_abc.
For totally antisymmetric Dirac torsion, this reduces to
    K_abc = -1/2 T_abc.
The implementation reconstructs T from K rather than assuming the identity.

With the action-derived Cartan relation
    T_abc = -(kappa/2) epsilon_abcd J5^d,
and epsilon_0123=+1 as an alternating symbol, Lorentz contraction yields
    T_abc T^abc = -(3/2) kappa^2 J5_a J5^a.
The code raises indices with eta=(-+++), computes both sides independently,
and reports the residual.

Evidence:
- PROVED: the wedge reconstruction identity follows from T^a=K^a_b wedge e^b
  and total antisymmetry.
- NUMERICALLY VERIFIED: tests check the reconstruction and quadratic identity
  on nontrivial timelike/spacelike axial-current data.
- OPEN: insert K into the full EC+Dirac action, retain the boundary term
  bookkeeping, and derive the precise effective four-fermion Lagrangian and
  Hehl-Datta equation in these same conventions. This module does not assign
  an independent propagating kinetic term to torsion.
