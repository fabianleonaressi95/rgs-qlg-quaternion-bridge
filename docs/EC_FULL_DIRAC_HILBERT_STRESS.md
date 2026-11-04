# Full Hilbert tensor of torsion-eliminated Einstein–Cartan–Dirac matter

## Scope and convention

This module extends the fixed-state contact variation by adding the symmetric
Belinfante/Hilbert kinetic tensor of the symmetrized Dirac action. It uses the
repository convention eta=(-,+,+,+), gamma^a and bar(psi) from
`rgs_ec_algebraic_torsion`, and
`L_4f=(3 kappa/16) J5_a J5^a`.

For a metric-compatible spin connection, define
R_ab=bar(psi) gamma_a D_b psi and
L_ab=(D_b bar(psi)) gamma_a psi. The symmetric on-shell kinetic tensor is

    T_ab^(D) = (i/4) [R_ab+R_ba-L_ab-L_ba].

The factor i belongs to the repository's declared gamma/bar convention.
The implementation checks that the resulting tensor is real and symmetric.
The spinor-based builder forms D_b bar(psi) by applying the repository's
Dirac adjoint to D_b psi; this assumes a metric-compatible spin connection.

The eliminated algebraic torsion contributes

    L_4f = (3 kappa/16) eta_cd M^{cd},
    M^{cd}=<J5^c J5^d>,
    T_ab^(4f, explicit) = eta_ab L_4f.

The combined tensor is T_ab=T_ab^(D)+T_ab^(4f,explicit).

## What is and is not derived

- ACTION-STRUCTURE: the kinetic term is the symmetric Dirac Hilbert/Belinfante
  expression after the spin-connection variation is accounted for; the
  contact term is varied explicitly at fixed local-frame state.
- NUMERICALLY TESTABLE: symmetry, reality, contact contraction, additive
  decomposition, and pure-spinor second moments.
- NOT YET DERIVED: a discrete/continuum solution of the nonlinear effective
  Dirac equation, state response under tetrad/scale-factor variation, a
  closed homogeneous fermion state, and the complete conservation reduction.

The kinetic tensor routine returns a formula, not proof that its input
spinor derivatives solve the interacting Dirac equation. The combined
routine therefore sets `dirac_eom_enforced=False` and
`conservation_status=REQUIRES_EFFECTIVE_DIRAC_EOM_AND_STATE_CLOSURE`.
No equation of state, rho_spin, beta parameter, or Friedmann correction is
invented. Minimal Einstein–Cartan torsion remains algebraic and contributes
zero propagating torsion degrees of freedom.

## Validation

Run:

    python -m pytest -q tests/test_rgs_ec_full_dirac_stress.py
