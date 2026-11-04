# Q6 U(1) Wilson transfer operator in the electric-flux basis

## 1. Exact temporal kernel

For a temporal plaquette, the character expansion is

\[
e^{\beta_t\cos F_{0e}} =
\sum_{n_e\in\mathbb Z} I_{n_e}(\beta_t)e^{i n_e F_{0e}}.
\]

Integrating all compact temporal links with normalized Haar measure imposes
\(D_0^T n=0\) at every spatial vertex. Thus the temporal transfer kernel is
diagonal in the integer cycle lattice. After dividing out the zero-flux
normalization, its exact eigenvalue for a flux configuration is

\[
\lambda_{\boldsymbol n}(\beta_t)
=\prod_{e=1}^{192}\frac{I_{|n_e|}(\beta_t)}{I_0(\beta_t)},
\qquad D_0^T\boldsymbol n=0.
\]

The corresponding exact finite-\(\beta_t\) electric energy, with the vacuum
energy subtracted, is

\[
E_{\boldsymbol n}(a_t)
=-\frac{1}{a_t}\log\lambda_{\boldsymbol n}(\beta_t).
\]

## 2. Symmetric transfer operator

Let \(V_s(A)=\beta_s\sum_f[1-\cos((D_1A)_f)]\). In link-angle
representation the symmetric one-step transfer operator is

\[
T=e^{-V_s/2}\,K_t\,e^{-V_s/2},
\]

where

\[
K_t(A',A)=
\sum_{\substack{\boldsymbol n\in\mathbb Z^{192}\\D_0^T\boldsymbol n=0}}
\lambda_{\boldsymbol n}(\beta_t)
e^{i\boldsymbol n\cdot(A'-A)}.
\]

The implementation evaluates this kernel for an explicitly supplied finite
set of integer flux states. This is a Galerkin truncation of the infinite
cycle-lattice sum; it is not silently presented as the full transfer operator.
The spatial Wilson potential remains an exact multiplication factor.

## 3. Anisotropic continuum matching

Starting from

\[
S_E=\frac{1}{4g_{d+1}^2}\int d\tau\,d^dx\,
F_{\mu\nu}F_{\mu\nu}
=\frac{1}{2g_{d+1}^2}\int d\tau\,d^dx
\sum_{\mu<\nu}F_{\mu\nu}^2,
\]

and assigning spatial spacing \(a_s\) and temporal spacing \(a_t\), the
small-plaquette Wilson coefficients are

\[
\boxed{\beta_t=\frac{a_s^{d-2}}{g_{d+1}^2a_t}},
\qquad
\boxed{\beta_s=\frac{a_ta_s^{d-4}}{g_{d+1}^2}}.
\]

The electric coefficient obtained from the transfer eigenvalue is

\[
-\frac{1}{a_t}\log\frac{I_n(\beta_t)}{I_0(\beta_t)}
\sim \frac{n^2}{2\beta_ta_t}
=\frac{g_{d+1}^2}{2a_s^{d-2}}n^2
\quad (\beta_t\to\infty,\ n\text{ fixed}),
\]

while the magnetic Hamiltonian coefficient is

\[
\frac{\beta_s}{a_t}=\frac{a_s^{d-4}}{g_{d+1}^2}.
\]

For \(d=3\), these become the familiar Hamiltonian lattice-gauge
coefficients \(g^2/(2a_s)\) and \(1/(g^2a_s)\). For the actual Q6
spatial complex, \(d=6\):

\[
\beta_t=\frac{a_s^4}{g_7^2a_t},\qquad
\beta_s=\frac{a_ta_s^2}{g_7^2},
\]

and \([g_7^2]=\text{length}^3\). This dimensional distinction matters:
one cannot import the 3+1D coupling scaling into Q6 unchanged.

The Bessel-ratio asymptotic is a fixed-flux, large-\(\beta_t\) statement.
It is not by itself a proof of the full continuum limit, renormalizability,
or a 3+1D reduction.

## 4. Audits and scope

The tests cover:
- positivity and vacuum normalization of the electric transfer eigenvalues;
- exact integer Gauss restriction of the supplied flux states;
- Hermiticity of the symmetric finite-character kernel;
- gauge invariance of the spatial Wilson potential;
- the fixed-flux Bessel asymptotic;
- separate \(d=3\) and \(d=6\) anisotropic scaling laws.

The default audit uses the zero flux plus positive and negative fundamental
cycles (259 states). It is a finite Galerkin diagnostic, not a converged
approximation to the complete 129-dimensional integer cycle lattice.

Still open:
1. convergence under systematic enlargement of the flux-state set;
2. derivation of \(g_7^2\), or its RGS-PSU counterpart, from the microscopic
   RGS-PSU action;
3. a physically justified three-dimensional spatial subcomplex before any
   3+1D claim;
4. the continuum and renormalization analysis of the 6+1D model.
