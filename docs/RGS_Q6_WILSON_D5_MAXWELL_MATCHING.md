# Q6 Wilson-to-D5 Maxwell matching: exact spectral and coefficient theorem

## 1. The theorem and its scope

Let X_6 be the full six-dimensional cubical complex of the 6-cube: 64
vertices, 192 oriented edges and 240 square faces. Let D_0 and D_1 be its
integer incidence matrices, with conventions

- D_0: C^0(X_6) -> C^1(X_6),
- D_1: C^1(X_6) -> C^2(X_6),
- D_1 D_0 = 0.

Use the anisotropic compact U(1) Wilson action with spatial spacing a_s,
temporal spacing a_t and continuum normalization

S_E = (1/(4 g_7^2)) integral d tau d^6x F_{mu nu} F_{mu nu}.

The theorem has two logically separate parts:

**(A) Exact finite-complex spectral statement.** On the Coulomb-gauge
cochain space ker(D_0^T), the positive operator D_1^T D_1 has an eigenvalue
12 of multiplicity exactly five.

**(B) Projected Wilson–Maxwell coefficient identity.** In the fixed-integer-
flux, beta_t -> infinity Hamiltonian limit, and after real orthogonal
projection onto that five-dimensional eigenspace, the quadratic Hamiltonian
is

H_5 = (1/(2 beta_eff)) pi^T pi
      + (beta_eff lambda_eff/2) q^T q,

where

beta_eff = a_s^4/g_7^2,       lambda_eff = 12/a_s^2.

Thus the five projected coordinates have oscillator frequency
omega = sqrt(12)/a_s, independent of g_7^2 and a_t.

Part (A) is an exact statement about the cubical complex, not a fit to
numerical data. Part (B) is exact algebra within the stated continuum-
matched quadratic Hamiltonian and spectral projection. The passage from the
compact finite-beta_t theory to that quadratic Hamiltonian is a controlled
fixed-flux asymptotic limit, not an exact identity at finite beta_t.

## 2. Exact derivation of the five-dimensional coexact eigenspace

The 6-cube is the Cartesian product of six copies of the interval complex
I = (two vertices, one edge). For one interval, the degree-zero Hodge
Laplacian has eigenvalues 0 and 2, while the degree-one Hodge Laplacian has
the single eigenvalue 2.

The Hodge Laplacian on a product complex is the sum of the factor
Laplacians. A degree-one cochain on X_6 has exactly one degree-one factor,
which can be chosen in six ways; the other five factors are degree zero.
Choosing k of those five vertex factors in their eigenvalue-2 state gives

  lambda_k = 2 + 2 k,
  multiplicity in the full degree-one Hodge Laplacian = 6 C(5,k),
  k = 0,...,5.

The full degree-one Hodge Laplacian is

  Delta_1 = D_0 D_0^T + D_1^T D_1.

The exact subspace im(D_0) is invariant and isomorphic, away from constants,
to the positive-spectrum part of the degree-zero Laplacian. The vertex
Laplacian on the 6-cube has eigenvalue 2m with multiplicity C(6,m).
Consequently, at lambda_k = 2(k+1), the exact 1-form multiplicity is
C(6,k+1). Subtracting exact modes leaves the coexact multiplicity

  m_k = 6 C(5,k) - C(6,k+1).

In particular,

  lambda_5 = 12,
  m_5 = 6 C(5,5) - C(6,6) = 6 - 1 = 5.

Since D_0^T D_1^T = 0, D_1^T D_1 preserves ker(D_0^T), and on that
space Delta_1 = D_1^T D_1. This proves that the eigenvalue 12 has
multiplicity exactly five on ker(D_0^T).

For reference, the complete positive coexact spectrum is

| eigenvalue | coexact multiplicity |
|---:|---:|
| 4 | 15 |
| 6 | 40 |
| 8 | 45 |
| 10 | 24 |
| 12 | 5 |

The multiplicities sum to 129 = 192 - rank(D_0), since rank(D_0)=63.
This is consistent with H^1(X_6)=0 and rank(D_1)=129.

Let B in R^(192 x 5) have orthonormal columns spanning the lambda=12
coexact eigenspace. Then, exactly,

  B^T B = I_5,
  D_0^T B = 0,
  D_1^T D_1 B = 12 B.

The finite-matrix implementation audits these identities numerically; the
multiplicity-five result itself follows from the product-complex argument
above. The A5 representation carried by this space is a further symmetry
statement and is audited separately by the Q6/D5 representation tests.

## 3. Wilson transfer limit and continuum Hamiltonian

For a temporal plaquette coupling beta_t, the character expansion is

  exp(beta_t cos F) = sum_{n in Z} I_n(beta_t) exp(i n F).

Haar integration over temporal links imposes the exact integer Gauss law

  D_0^T n = 0,   n in Z^192.

After vacuum normalization, the temporal transfer eigenvalue for a fixed
integer flux vector n is

  Lambda_t(n) = product_e I_|n_e|(beta_t) / I_0(beta_t).

For each fixed integer n_e,

  -log(I_|n_e|(beta_t)/I_0(beta_t))
      = n_e^2/(2 beta_t) + O(beta_t^-2).

Therefore

  E_n = -log Lambda_t(n)/a_t
      = ||n||^2/(2 beta_t a_t) + O(beta_t^-2/a_t).

This is a fixed-flux asymptotic expansion. It is not uniform over arbitrary
flux vectors whose components grow with beta_t, and it does not turn the
finite-beta_t Bessel spectrum into an exactly quadratic spectrum.

Matching the anisotropic Wilson couplings to the stated continuum action in
d spatial dimensions gives

  beta_t = a_s^(d-2)/(g_(d+1)^2 a_t),
  beta_s = a_t a_s^(d-4)/g_(d+1)^2.

The electric Legendre term and the quadratic expansion of the spatial
plaquette potential are therefore

  H_E = g_(d+1)^2/(2 a_s^(d-2)) ||p||^2,

  H_B = a_s^(d-4)/(2 g_(d+1)^2)
        <A, D_1^T D_1 A>.

For Q6, d=6 and g_(d+1)^2 = g_7^2. Thus

  H = g_7^2/(2 a_s^4) ||p||^2
      + a_s^2/(2 g_7^2) <A,D_1^T D_1 A>,

on the Gauss/Coulomb tangent space D_0^T A=0. The magnetic expression is
the quadratic (small-plaquette) part of
(beta_s/a_t) sum_f [1-cos((D_1 A)_f)].

## 4. Orthogonal spectral projection and coefficient identity

Write the projected canonical variables as

  A = B q,       p = B pi,

where q, pi are in R^5. Because B^T B=I_5, this is a canonical orthonormal
mode normalization: the projected symplectic form is
sum_i d pi_i wedge d q_i. The electric term becomes

  ||p||^2 = pi^T B^T B pi = pi^T pi.

The magnetic term becomes

  A^T D_1^T D_1 A
    = q^T B^T D_1^T D_1 B q
    = 12 q^T q.

Consequently,

  H_5 = 1/2 [
          (g_7^2/a_s^4) pi^T pi
          + (12 a_s^2/g_7^2) q^T q
        ].

Set

  beta_eff = a_s^4/g_7^2,
  lambda_eff = 12/a_s^2.

Then 1/beta_eff = g_7^2/a_s^4 and
beta_eff lambda_eff = 12 a_s^2/g_7^2, hence

  H_5 = pi^T pi/(2 beta_eff)
        + beta_eff lambda_eff q^T q/2.

This is precisely the coefficient convention of the five-mode D5 Maxwell
oscillator Hamiltonian. The frequency follows from the product of the
kinetic and potential Hessians:

  omega^2 = (1/beta_eff)(beta_eff lambda_eff)
           = lambda_eff = 12/a_s^2.

The cancellation of g_7^2 is algebraic. The cancellation of a_t occurs
when the transfer logarithm is divided by a_t and the anisotropic Wilson
couplings are substituted. For a_s=1, beta_eff=1/g_7^2 and
lambda_eff=12.

## 5. Epistemic status and non-implications

**Exact finite-complex result:** the coexact lambda=12 eigenspace has
dimension five, by the tensor-product Hodge spectrum of the 6-cube.

**Exact conditional coefficient theorem:** given the continuum-matched
quadratic U(1) Hamiltonian and the orthogonal spectral projection B, both
D5 Maxwell Hessian coefficients and the oscillator frequency agree
identically, for every positive a_s and g_7^2.

**Numerically audited implementation:** orthonormality, Gauss tangency,
Hodge eigen-equation, projected Hessian, and coefficient/frequency
equalities are checked by the repository test suite over multiple parameter
choices.

**Not established by this theorem:**

1. The real five-dimensional eigenspace is not thereby an exact rank-five
   sublattice of the integer cycle lattice
   ker(D_0^T) intersect Z^192. The integer Gauss-law lattice has rank 129;
   an arbitrary real spectral subspace need not be a sublattice.
2. The finite-beta_t compact transfer spectrum is a product of Bessel
   ratios, not exactly Gaussian/quadratic.
3. The quadratic spatial magnetic approximation does not prove equality of
   the full compact Wilson theory and a free Maxwell theory at finite field.
4. The five-mode projection is not a global embedding of the compact Q6
   gauge quotient into the noncompact type-IV domain D5. The D5 statement
   concerns the ten-real-dimensional cotangent tangent sector and its
   oscillator normalization.
5. No physical 3+1D reduction, microscopic RGS–PSU determination of g_7^2,
   electric-charge normalization, renormalization prescription, or
   prediction/derivation of the fine-structure constant alpha follows.

## 6. Reproducibility

Run:

  pytest -q tests/test_rgs_q6_u1_maxwell_matching.py
  pytest -q tests/test_rgs_q6_hodge_spectral_theorem.py

The coefficient report is
q6_u1_maxwell_matching_report(...) in
src/rgs_q6_u1_maxwell_matching.py.
