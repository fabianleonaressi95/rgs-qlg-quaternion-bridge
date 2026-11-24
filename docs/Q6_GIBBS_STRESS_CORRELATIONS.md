# Q6 Gibbs equilibrium and stress-correlation benchmark

## Scope

This module samples only the compact SU(2) link variables of the declared
conditional face-area-weighted gauge action
[
S_g[X,U]=\beta\sum_{p=1}^{240}\frac{A_p^0}{A_p(X)}
\left(1-\frac12\operatorname{ReTr}U_p\right).
]
At fixed recovered RGS coordinates, the finite Gibbs measure is
[
d\mu_{Q6}(U)=Z^{-1}\exp[-\lambda S_g(X,U)]
\prod_{e\in Q6^{(1)}}dU_e,
\qquad U_e\in SU(2),
]
where (dU_e) is normalized Haar measure and (lambda>0) is the
dimensionless sampling action scale. Since there are finitely many compact
SU(2) variables and the action is continuous and bounded, (0<Z<\infty).

The sampler uses symmetric single-link proposals
(U_e\mapsto U_e\exp(i\sigma\cdot v)), with isotropic Gaussian (v),
and the Metropolis acceptance probability
[
a(U\to U')=\min\{1,\exp[-\lambda(S_g(U')-S_g(U))]\}.
]
The proposal is symmetric under inversion; hence detailed balance holds
with respect to (d\mu_{Q6}). The code records the gauge shear observable
(Pi_{xy}^g=-V^{-1}\partial_\gamma S_g|_{\gamma=0}).

## Correlation estimator

For the sampled sequence (Pi_n), the implementation computes
[
C(k)=\frac{\widehat{\operatorname{Cov}}(\Pi_n,\Pi_{n+k})}
{\widehat{\operatorname{Var}}(\Pi_n)},\qquad C(0)=1,
]
and the finite trapezoidal sum
[
\tau_{\rm int}^{MC}(K)=\Delta n
\left[\frac12 C(0)+\sum_{k=1}^{K-1}C(k)+\frac12C(K)\right].
]
This is an **algorithmic Monte Carlo autocorrelation**. The chain's step
index is not physical time, and the integral is not a Green--Kubo viscosity.

## What is established and what remains open

- **Exact finite construction:** the compact-link Gibbs density and
  Metropolis detailed-balance condition for the declared action.
- **Numerical benchmark:** finite-chain action, stress, acceptance rate
  and normalized autocorrelation, with deterministic seeded tests.
- **Not established:** ergodicity/mixing rate for the chosen finite proposal,
  statistical convergence from any particular run, a physical Hamiltonian
  or Langevin generator, fluctuation-dissipation normalization, physical
  time calibration, thermodynamic/refinement limit, Green--Kubo viscosity,
  or Stokes matching.

Matter fields are intentionally not integrated here. The present
massless graph Dirichlet term has covariantly constant zero-mode issues;
a normalizable matter Gibbs measure requires a specified mass/constraint
and its physical interpretation. No physical viscosity is reported.
