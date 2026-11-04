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


## Nonlinear effective Dirac equation

The axial contact operator is differentiated with respect to independent
Dirac variables. With alpha=3 kappa/16 and
J5^a=i bar(psi) gamma5 gamma^a psi,

    d L4f / d bar(psi) = 2 alpha J5_a i gamma5 gamma^a psi,
    d L4f / d psi      = 2 alpha J5_a i bar(psi) gamma5 gamma^a.

The column Dirac equation therefore has the form

    K[psi] + 2 alpha J5_a i gamma5 gamma^a psi = 0,

where K[psi] is the Euler operator obtained by varying the complete
symmetrized kinetic and mass action. The new module computes the nonlinear
contact source and combines it with a caller-supplied K[psi]. It does not
pretend that a pointwise residual supplies the missing curved-space kinetic
variation, boundary conditions, or a global solution.

The independent-variable derivatives are checked against central
finite-difference variations of the contact Lagrangian. This validates the
local nonlinear source, not a solved interacting Dirac field.

## Conservation remains a separate audit

On-shell covariant conservation follows from diffeomorphism invariance and
the complete matter equations, provided tetrad and spin-connection
variations, boundary terms, and state dependence are consistently included.
A local algebraic EOM residual is necessary but not sufficient for a
numerical conservation claim. This implementation therefore makes no claim
of covariant conservation and does not infer FLRW rho, pressure, equation of
state, or TorsionLambda beta.

Run the additional tests with:

    python -m pytest -q tests/test_rgs_ec_nonlinear_dirac_eom.py


## Kinetic Euler operator from the declared symmetrized action

The repository's action convention is

    S_D = 1/2 int [bar(psi) *gamma wedge D psi
                    + D bar(psi) wedge *gamma psi].

In an oriented orthonormal frame, metric compatibility gives
`*e^a wedge e^b = -eta^{ab} vol` and
`e^b wedge *e^a = eta^{ba} vol`. Therefore

    L_D = 1/2 [D_a bar(psi) gamma^a psi
               - bar(psi) gamma^a D_a psi] vol.

Varying the independent row spinor and integrating the first term by parts,
with covariantly constant gamma matrices and discarded boundary terms,
gives

    delta S_D / delta bar(psi) = -gamma^a D_a psi.

The new `rgs_ec_dirac_kinetic_eom.py` implements this local kinetic Euler
operator. The caller supplies the complete `D_a psi`, including the
chosen tetrad and spin connection. It then adds the already audited
Einstein-Cartan axial contact source, yielding the massless local equation

    -gamma^a D_a psi + (3 kappa/8) J5_a i gamma5 gamma^a psi = 0.

The declared repository action does not specify a mass term, so none is
silently added. This is an action-derived local operator, not a global
curved-space solution: spinor boundary conditions, a solved spin connection,
tetrad equations, and covariant conservation remain separate tasks.
