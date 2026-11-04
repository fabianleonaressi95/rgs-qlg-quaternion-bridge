# Local Levi-Civita spin connection from a coframe jet

**Status:** ACTION/GEOMETRY-DERIVED local construction with finite analytic test cases. It is not yet the RGS point-cloud spin connection or a teleparallel inertial connection.

## 1. Index contract

The implementation in src/rgs_spin_connection.py uses:

- coframe[a, mu] = e^a_mu;
- coframe_derivatives[a, nu, mu] = partial_mu e^a_nu;
- inverse_coframe[mu, a] = e_a^mu;
- omega_mixed[a,b,mu] = omega^a_(b mu);
- omega_lower[a,b,mu] = eta_ac omega^c_(b mu).

The metric is g_mu_nu = eta_ab e^a_mu e^b_nu with eta=(-,+,+,+).

## 2. Derivation

The Levi-Civita Christoffel symbols are calculated from the metric and its first derivatives:

\[
\Gamma^\rho{}_{\mu\nu}
=\frac12 g^{\rho\sigma}
(\partial_\mu g_{\sigma\nu}
+\partial_\nu g_{\sigma\mu}
-\partial_\sigma g_{\mu\nu}).
\]

The tetrad postulate is imposed:

\[
\partial_\mu e^a{}_{\nu}
+\omega^a{}_{b\mu}e^b{}_{\nu}
-\Gamma^\rho{}_{\mu\nu}e^a{}_{\rho}=0.
\]

Solving for the mixed spin connection gives:

\[
\omega^a{}_{c\mu}
=(\Gamma^\rho{}_{\mu\nu}e^a{}_{\rho}
-\partial_\mu e^a{}_{\nu})e_c{}^\nu.
\]

Metric compatibility requires omega_ab_mu = -omega_ba_mu. The returned audit independently checks this antisymmetry, the tetrad postulate and vanishing torsion.

The spinor representation is:

\[
\Omega_\mu=\frac14\omega_{ab\mu}\gamma^a\gamma^b,
\qquad
D_\mu\psi=\partial_\mu\psi+\Omega_\mu\psi.
\]

The orthonormal derivative supplied to the existing EC–Dirac Euler operator is:

\[
D_a\psi=e_a{}^\mu D_\mu\psi.
\]

## 3. Validation

The tests include:

1. constant Minkowski coframe, for which omega=0;
2. spherical-coordinate Minkowski coframe e=diag(1,1,r,r sin(theta)), for which the connection is nonzero but torsion and tetrad-postulate residuals vanish;
3. spinor-connection contraction and conversion from coordinate to orthonormal derivatives;
4. invalid/singular input rejection.

## 4. Scope and next bridge

This module accepts the first coframe jet as input. It does not estimate derivatives from the 64-site RGS coframe, and it does not decide whether the physical connection is Levi-Civita, Einstein–Cartan with independent torsion, or teleparallel with an inertial spin connection.

Next, the RGS refinement must provide consistent derivatives of the coframe (or a discrete connection with a proven continuum map). The EC branch then needs the independent connection equation and algebraic torsion solution; the TEGR branch needs its own inertial spin-connection construction and identity audit.

Internal U(1), SU(2), and SU(3) gauge connections are separate representation-valued terms and are not included in Omega_mu.
