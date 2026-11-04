# Q6 Wilson-to-D5 Maxwell matching: coefficient theorem

## Statement

Consider the anisotropic compact U(1) Wilson action on the six-dimensional
spatial cubical complex Q6, with spatial spacing a_s, temporal spacing a_t,
and continuum normalization

S_E=(1/(4g_7^2)) integral d tau d^6x F_mu_nu F_mu_nu.

The temporal character expansion and integer Gauss projector give, for
fixed integer flux and beta_t -> infinity,

E_n=-(1/a_t) log(I_n(beta_t)/I_0(beta_t))
   = n^2/(2 beta_t a_t) + O(beta_t^-2/a_t).

With beta_t=a_s^4/(g_7^2 a_t) and beta_s=a_t a_s^2/g_7^2, the quadratic
Hamiltonian is

H=(g_7^2/(2a_s^4))*||p||^2
  +(a_s^2/(2g_7^2))*<A,D1^T D1 A>.

Here D0^T A=0 is the Coulomb-gauge tangent condition and the magnetic
term is the small-plaquette expansion of the spatial Wilson potential.

## Five-mode projection

The Q6 Hodge operator D1^T D1, restricted to ker(D0^T), has a
five-dimensional coexact eigenspace with eigenvalue lambda=12. Let
B in R^(192x5) be an orthonormal basis of that eigenspace. The finite-complex
audits verify

B^T B=I_5,  D0^T B=0,  D1^T D1 B=12 B.

Set A=Bq and p=B pi. Then the restricted Hamiltonian is

H_5=1/2 * [
  (g_7^2/a_s^4) pi^T pi + (12 a_s^2/g_7^2) q^T q
].

Define

beta_eff=a_s^4/g_7^2,   lambda_eff=12/a_s^2.

It follows identically that

H_5=pi^T pi/(2 beta_eff)
    + beta_eff*lambda_eff*q^T q/2.

This is exactly the coefficient convention implemented by the D5 five-mode
Maxwell Hamiltonian. The frequency is omega=sqrt(12)/a_s, independent of
g_7^2. The temporal spacing cancels from the Hamiltonian coefficients after
the transfer logarithm is divided by a_t.

For a_s=1, this reduces to beta=1/g_7^2 and lambda=12.

## What this proves—and what it does not

**Proved algebraically under the stated Wilson-to-continuum normalization:**
the quadratic continuum Hamiltonian restricted to the real five-dimensional
coexact eigenspace has exactly the D5 Maxwell oscillator coefficients.

**Verified numerically on the finite Q6 complex:** orthonormality, Gauss-tangent
condition, Hodge eigenvalue equation, and equality of both Hamiltonian
coefficients and oscillator frequencies over multiple a_s, a_t, g_7^2
choices.

**Not proved:** the real five-dimensional eigenspace is not shown to be an
exact rank-five sublattice of the integer cycle lattice
ker(D0^T) intersect Z^192. The compact finite-beta_t transfer spectrum is a
product of Bessel ratios and is not exactly quadratic. The Maxwell
identification uses the fixed-flux large-beta_t limit and a real spectral
projection. No global Q6-to-D5 map, physical 3+1D reduction, microscopic
derivation of g_7^2, or prediction of alpha follows.

## Reproducibility

Run pytest -q tests/test_rgs_q6_u1_maxwell_matching.py.

The report function is q6_u1_maxwell_matching_report(...) in
src/rgs_q6_u1_maxwell_matching.py.
