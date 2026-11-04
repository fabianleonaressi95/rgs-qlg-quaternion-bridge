# Q6 Wilson spectral regulator and the D5 cutoff question

## Exact regulator supplied by the compact transfer operator

For a divergence-free integer electric flux (n\in\ker(D_0^T)\cap\mathbb Z^{192}),
the normalized temporal transfer eigenvalue is

\[
\lambda_n(\beta_t)=\prod_{e=1}^{192}
\frac{I_{|n_e|}(\beta_t)}{I_0(\beta_t)},\qquad
E_n(\beta_t,a_t)=-\frac{1}{a_t}\log\lambda_n.
\]

For any threshold (0<\tau\leq1), the spectral projector is therefore
defined exactly by either equivalent condition

\[
P_{\tau}=\mathbf 1\{\lambda_n\geq\tau\}
=\mathbf 1\{E_n\leq -a_t^{-1}\log\tau\}.
\]

This is a genuine finite-\(\beta_t\) Wilson-sector regulator. It uses no
large-\(\beta_t\) approximation and introduces no arbitrary D5 coordinate.
On the implemented finite Galerkin shell (vacuum and both orientations of
the 129 fundamental integral cycles) it gives a reproducible finite
selection. The shell is not the complete infinite lattice
\(\mathbb Z^{129}\).

## Why this does not yet determine the Hua cutoff

The integer cycle lattice has rank 129. Let (B\in\mathbb R^{192\times5})
be the orthonormal basis of the selected five-mode coexact eigenspace and
let (C\in\mathbb Z^{129\times192}) be the fundamental-cycle matrix.
The five-mode momentum projection is

\[
\Pi_5:\mathbb R^{129}\to\mathbb R^5,\qquad
c\mapsto cCB.
\]

The five-mode subspace has an exact rational spectral projector. With
\(L_1=D_1^TD_1\), its projector is
\[
P_{12}=\frac{L_1(L_1-2I)(L_1-4I)(L_1-6I)(L_1-8I)(L_1-10I)}{46080}.
\]
The exact Q6 coexact spectrum gives \(\operatorname{rank}P_{12}=5\).
Because \(L_1\) and the fundamental-cycle matrix \(C\) are integer,
\(CP_{12}\) is rational. The cycle space contains the selected eigenspace,
so the map \(c\mapsto cCP_{12}\) has rank five. Clearing denominators and
applying the integer-kernel theorem gives
\[
\operatorname{rank}_{\mathbb Z}
\ker(\,c\in\mathbb Z^{129}\mapsto cCP_{12}\,)=124.
\]
This is an exact arithmetic statement, not an inference from real
dimension counting.

In particular, there is a nonzero integer Gauss-law flux \(n=cC\) with
\(nP_{12}=0\). It has exactly the same five-mode projection as the vacuum,
but its finite-\(\beta_t\) Wilson transfer weight is strictly smaller:
\[
\lambda_n=\prod_e\frac{I_{|n_e|}(\beta_t)}{I_0(\beta_t)}
<1=\lambda_0,\qquad \beta_t>0.
\]
Thus the full Wilson weight does **not** factor through the five D5
coordinates. This proves a genuine integer-lattice obstruction.

It still does not produce a scalar D5 Hua cutoff. Rather, it proves that
the full Wilson transfer spectrum cannot be represented by the existing
five coordinates alone. An effective D5 measure remains possible only
after specifying what happens along each 124-dimensional integer fibre:
for example, a conditional measure, a fibre sum/trace of the transfer
operator, or an action-derived quotient. The resulting pushforward is
additional structure, not an identification of \(\tau\) with
\(\varepsilon\).

The full proof is recorded in
docs/Q6_INTEGER_LATTICE_D5_WILSON_OBSTRUCTION.md.
## Correct separation of the two regulators

- **Wilson spectral cutoff:** \(\lambda_n\geq\tau\), exactly defined at
  finite \(\beta_t\), with energy threshold \(-\log\tau/a_t\).
- **D5 Hua cutoff:** \(N(z)\geq\epsilon\), whose exact normalizer and
  \((5/3)\epsilon^{-4}\) divergence were derived in PR #149.
- **Not established:** a function \(\epsilon=\epsilon(\tau,\beta_t,a_t,
  a_s,g^2)\) relating these two conditions.

The Wilson threshold is therefore a derived regulator in the compact
transfer sector; the D5 epsilon remains a separate mathematical regulator.
Identifying them now would add an unsupported assumption.

## Reproducibility and status

Run:

    python scripts/run_q6_wilson_d5_regulator.py

The JSON report records the selected finite-shell states, exact Bessel
weights, energy threshold, projected rank, kernel dimension, and checks.
Tests:

    pytest -q tests/test_rgs_q6_d5_wilson_regulator.py

**Exact:** finite-beta transfer eigenvalues and their spectral threshold;
integer Gauss law; cycle rank 129; rational rank-five spectral projector;
integer kernel rank 124; existence of invisible integer fluxes with Wilson
weight different from vacuum; failure of transfer-weight factorization
through the five D5 coordinates.

**Numerically verified:** orthonormal five-mode frame and finite-shell
spectral selection.

**Open:** explicit short integer kernel generators; a canonical fibre
measure/action and induced effective D5 density; an action-derived global
quotient and physical continuum scaling.
No 3+1D reduction, physical scale, or fine-structure constant is claimed.
