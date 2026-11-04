# Fourth-moment derivative coupling and SO(3) connection current

## Constitutive operator and positivity

The RGS seed supplies a self-adjoint map on the five-dimensional
symmetric-trace-free space,

    A(Q)_ab = R_abcd Q_cd.

Its spectrum has a negative Eg sector. Therefore A itself must not be
used as a kinetic metric. We define the positive shifted constitutive
operator

    H = mu I + beta A,
    mu = max(0, -beta lambda_min(A)) + epsilon, epsilon > 0.

Thus H is strictly positive on Sym_0(3), while retaining the fourth-moment
anisotropic splitting. The shift is a stability regularization; it is not
claimed to be a measured coupling or a new fundamental constant.

## Gauge-covariant derivative and action

For an so(3)-valued spatial connection omega_i and an independent
symmetric trace-free order parameter Q,

    D_i Q = partial_i Q + [omega_i,Q],

and the positive spatial gradient density is

    L_grad = kappa/2 sum_i <D_i Q, H(D_i Q)>_F, kappa >= 0.

This is invariant under simultaneous local frame transformations of
Q, R, omega and the derivative, provided H is transformed by the
induced representation on Sym_0(3). The seed anisotropy is material data;
holding its components fixed restricts the frame symmetry to its
stabilizer. This code audits the pointwise tensor law, not a global
bundle construction.

## Connection variation and spin current

With delta(D_i Q)=[delta omega_i,Q], Frobenius variation gives

    delta L_grad = <J_i, delta omega_i>_F,
    J_i = kappa [H(D_i Q), Q] in so(3).

The commutator of two symmetric matrices is antisymmetric, so the current
has the correct Lie-algebra index symmetry. Central finite differences
of the three independent so(3) connection generators audit this formula.

## Teleparallel boundary

This is a derivative coupling to an so(3) connection. It is not yet a
complete teleparallel theory. In particular, imposing flatness
R(omega)=0 requires a constrained variation (or inertial spin-connection
parameterization); the coframe equation, lapse/shift constraints,
TEGR sector, and boundary terms have not been combined here. A current
is derived, but a solved torsion field equation is not claimed.

## Reproduction and evidence

    PYTHONPATH=. python -m pytest tests/test_rgs_fourth_moment_connection.py -q

- NUMERICALLY VERIFIED: positive shifted metric, Sym_0 preservation,
  nonnegative gradient energy, antisymmetric connection current, and
  its finite-difference variation.
- CANDIDATE: fourth-moment-induced derivative constitutive coupling.
- OPEN: flat-connection constrained variation, full TEGR/coframe equations,
  coupling to Q6/Yang-Mills and fermion spinors, continuum/Mosco limit.
