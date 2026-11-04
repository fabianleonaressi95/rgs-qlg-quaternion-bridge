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

The audit verifies \(\operatorname{rank}(CB)=5\), hence the real-linear
extension has \(\dim\ker\Pi_5=129-5=124\). This proves that there is no
inverse linear map on all of \(\mathbb R^{129}\). It does **not**, by
dimension counting alone, prove that the restriction
\(\Pi_5|_{\mathbb Z^{129}}\) is non-injective: a real-linear map with a
nontrivial real kernel can still be injective on an integer lattice when
its coefficients have suitable arithmetic independence.

Accordingly, the current audit establishes a 124-dimensional real kernel,
but does not establish an integer collision, nor prove that the Wilson
weight fails to factor through the five projected coordinates on the
integer cycle lattice. Those are separate arithmetic/functional questions.
The finite-\(\beta_t\) weight is a product over all 192 edge-flux
components, and the repository currently supplies no proved factorization
through the five-mode coordinates and no global map to the D5 generic norm.
Thus no unique \(\epsilon\) is derived.

The precise obstruction currently proved is the absence of a real-linear
inverse/continuous linear reconstruction from the five modes to arbitrary
real cycle coordinates. The stronger discrete obstruction requires either
an explicit pair of distinct integer fluxes with identical projection but
different Wilson weights, or a theorem excluding a suitable factorization.
Neither is claimed here. A full action-preserving quotient or a
well-defined global map to the D5 chart is still missing.

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
integer Gauss law; cycle rank 129.

**Numerically verified:** rank-five five-mode projection and its
124-dimensional kernel; finite-shell spectral selection.

**Open:** an equivariant, action-preserving global quotient/map that retains
the discarded flux information; only after that construction can one ask
whether a Wilson spectral threshold induces a D5 generic-norm cutoff.
No 3+1D reduction, physical scale, or fine-structure constant is claimed.
