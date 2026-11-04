# Exact five-torus representation of the Q6-D5 Wilson fiber marginal

## Main identity

Let C be the integral fundamental-cycle matrix (129 x 192), and let K_Z be
the saturated rank-124 kernel of the rational five-mode projection. Choose
an integer matrix A (5 x 129) such that
K_Z = {k in Z^129 : k A^T = 0}. Clearing denominators in the transpose of
the rational mode map supplies such an A; a nonzero common scalar does not
change its integer kernel.

Define
F(theta) = product_e exp[beta_t cos((C^T theta)_e)] / I_0(beta_t),
for theta in T^129. Its absolutely convergent Fourier series is
F(theta) = sum_{a in Z^129} w(a C) exp(i a.theta), where
w(n) = product_e I_|n_e|(beta_t)/I_0(beta_t).

Substitute theta=A^T t, multiply by exp(-i(c A^T).t), and integrate
normalized Haar measure on T^5. Character orthogonality retains exactly
a such that (a-c)A^T=0, i.e. a=c+k with k in K_Z. Therefore

W_eff([c]) =
 integral_T5 exp(beta_t sum_e cos(((A C)[:,e]).t))
              / I_0(beta_t)^192
              * exp(-i(c A^T).t) dmu_Haar(t).

This is the complete fiber sum, not a finite truncation. It is independent
of the representative c modulo K_Z.

## Why this matters

The original sum has 124 integer summation directions. The exact dual
representation is a five-dimensional compact integral. It is a Fourier
coefficient of a positive analytic function; for c=0 the integrand is
positive, which is especially useful numerically.

## Numerical implementation and epistemic status

The implementation evaluates the integral using scrambled Sobol points in
five dimensions. It removes the known maximum log scale before exponentiation
and reports independent-scramble variability.

The numerical result is explicitly labelled QMC estimate, not certified
interval. A rigorous enclosure still requires validated quadrature or a
Fourier-tail/discrepancy bound. The exact identity itself follows from
absolute convergence and Fourier character orthogonality.

This representation does not identify the Q6 Wilson marginal with the D5
Hua measure. A measure pushforward and Wilson-to-Hua cutoff remain open.
