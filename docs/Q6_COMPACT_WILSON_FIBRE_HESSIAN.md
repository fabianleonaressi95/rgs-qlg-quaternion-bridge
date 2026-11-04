# Exact Hessian identity for the compact Wilson fibre marginal

## Setup

Let (M=D_1B\in\mathbb R^{240\times5}) be the plaquette image of the
retained five-mode chart. For an explicitly supplied finite fibre ensemble
\(\{\eta_j,a_j\}_{j=1}^J\), define
[
Z_J(q)=\sum_{j=1}^J a_j
\exp\left[-\beta_s\sum_{f=1}^{240}
(1-\cos(\eta_{jf}+(Mq)_f))\right],
\qquad a_j>0.
]
Set (S_J(q)=-\log Z_J(q)). This is a finite conditional fibre sum.

At (q=0), define the normalized conditional probabilities
[
p_j=\frac{a_j e^{-V_j(0)}}{Z_J(0)}.
]
Differentiation of the finite sum gives the exact identity
[
\boxed{
\nabla^2S_J(0)
=\mathbb E_p[\nabla^2V_j(0)]
-\operatorname{Cov}_p(\nabla V_j(0),\nabla V_j(0)).
}
]
For the compact Wilson cosine,
[
(\nabla V_j)_i
=\beta_s\sum_f\sin(\eta_{jf})(D_1B)_{fi},
]
[
(\nabla^2V_j)_{ik}
=\beta_s\sum_f\cos(\eta_{jf})
(D_1B)_{fi}(D_1B)_{fk}.
]
The covariance term is essential. Replacing the compact cosine by its
quadratic approximation sets the relevant fluctuations aside and generally
does not reproduce this finite-fibre Hessian.

## Interpretation

This identity is the correct nonlinear counterpart of the Gaussian
Schur-complement calculation: the curvature of the effective action is a
conditional mean local curvature minus a fluctuation covariance. The second
term is the precise channel through which integrating out fibre degrees of
freedom can renormalize the retained-mode stiffness.

The supplied implementation evaluates this identity stably for any
explicit finite list of fibre configurations and positive base weights.
A finite list is not the complete rank-124 integer fibre. The routine
therefore reports its scope honestly and does not claim a full Wilson
marginal, Bergman density, Hua cutoff, or prediction of alpha.

## Reproducibility

```bash
pytest -q tests/test_rgs_q6_compact_wilson_fiber_hessian.py
```
