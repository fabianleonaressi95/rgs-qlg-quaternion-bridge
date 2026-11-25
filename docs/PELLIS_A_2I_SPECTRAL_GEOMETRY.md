# Pellis A — Exact spectral geometry of the Poincaré space \(S^3/2I\)

**Version:** 1.0 — 25 November 2026  
**Status:** CLOSED for the scalar Laplacian on the unit spherical quotient; spinorial, vector-valued and Pellis-resonance correspondences remain separate questions.

## 1. Geometric and group-theoretic setup

Let \(2I\subset SU(2)\) be the binary icosahedral group of order \(120\), acting freely on the unit three-sphere \(S^3\simeq SU(2)\) by left multiplication. The quotient
\[
\mathcal P=S^3/2I
\]
is the Poincaré homology sphere with its standard round quotient metric. Its Riemannian volume is
\[
\operatorname{Vol}(\mathcal P)=\frac{1}{120}\operatorname{Vol}(S^3)
=\frac{\pi^2}{60}.
\]
The scalar Hilbert space is \(L^2(\mathcal P,d\mu_{\mathcal P})\), with the quotient Riemannian measure. The non-negative scalar Laplacian is the self-adjoint operator
\[
A=-\Delta_{\mathcal P}.
\]

The action is free, so the quotient is a smooth closed Riemannian three-manifold, not an orbifold with fixed-point singularities.

## 2. Invariant polynomial ring and Molien series

For the defining two-dimensional complex representation of \(2I\), the invariant ring has Molien series
\[
M(t)=\frac{1+t^{30}}{(1-t^{12})(1-t^{20})}
=\frac{1-t^{60}}{(1-t^{12})(1-t^{20})(1-t^{30})}.
\]
Thus the basic invariant degrees are \(12,20,30\), with one relation in degree \(60\). If
\[
a_l=[t^l]M(t),
\]
then \(a_l\) is the dimension of the invariant subspace in the degree-\(l\) irreducible \(SU(2)\) representation \(\operatorname{Sym}^l(\mathbb C^2)\).

The equivalent character formula is
\[
a_l=\frac1{120}\sum_{g\in2I}\chi_l(g),
\]
where \(\chi_l\) is the character of \(\operatorname{Sym}^l(\mathbb C^2)\). This formula fixes multiplicities from the group action, not from fitted spectral weights.

## 3. Scalar Laplacian spectrum

The degree-\(l\) scalar harmonics on \(S^3\) have eigenvalue
\[
\lambda_l=l(l+2)
\]
and transform as the \((l+1)\)-dimensional irreducible representation on the left, tensored with an \((l+1)\)-dimensional multiplicity space on the right. Taking the left \(2I\)-invariants leaves
\[
\boxed{\lambda_l=l(l+2),\qquad
m_l=(l+1)a_l.}
\]
Only degrees with \(a_l>0\) occur on the quotient. The first levels are:

| \(l\) | \(a_l\) | \(\lambda_l\) | \(m_l\) |
|---:|---:|---:|---:|
| 0 | 1 | 0 | 1 |
| 12 | 1 | 168 | 13 |
| 20 | 1 | 440 | 21 |
| 24 | 1 | 624 | 25 |
| 30 | 1 | 960 | 31 |
| 32 | 1 | 1088 | 33 |
| 36 | 1 | 1368 | 37 |
| 40 | 1 | 1680 | 41 |
| 42 | 1 | 1848 | 43 |
| 44 | 1 | 2024 | 45 |
| 48 | 1 | 2400 | 49 |
| 50 | 1 | 2600 | 51 |
| 52 | 1 | 2808 | 53 |
| 54 | 1 | 3024 | 55 |
| 56 | 1 | 3248 | 57 |
| 60 | 2 | 3720 | 122 |

In particular, the first positive scalar eigenvalue is \(168\), with multiplicity \(13\).

## 4. Spectral theorem and canonical weights

Since \(\mathcal P\) is compact and has no boundary, \(A=-\Delta_{\mathcal P}\) is non-negative and self-adjoint on its standard \(H^2\) domain. Its spectral resolution is
\[
A=\int_{\sigma(A)}\lambda\,dE_A(\lambda).
\]
For a normalized state \(\psi\in L^2(\mathcal P)\), the associated spectral measure is
\[
\mu_\psi(B)=\langle\psi,E_A(B)\psi\rangle,
\qquad \mu_\psi(\mathbb R_{\ge0})=1.
\]
For a discrete eigenspace projector \(P_l\), its weight is
\[
w_l=\langle\psi,P_l\psi\rangle\ge0,\qquad \sum_lw_l=1.
\]
The operator fixes the eigenvalues and eigenspace multiplicities. The weights additionally depend on the state \(\psi\); they are not uniquely fixed until a geometrically or variationally justified state-selection rule is supplied.

The Poisson-type spectral semigroup
\[
R_s=e^{-s\sqrt A}
\]
is well-defined for \(s\ge0\). Setting \(s=2\log\varphi\) gives the mode factor
\[
e^{-2\log(\varphi)\sqrt{l(l+2)}}
=\varphi^{-2\sqrt{l(l+2)}}.
\]
This is a valid golden-weighted semigroup once \(s\) is chosen. The choice \(s=2\log\varphi\) is an additional parameter choice, not a consequence of \(2I\)-symmetry alone.

## 5. Exact result versus proposed Pellis correspondence

**Closed exact result**
- the free \(2I\)-action and quotient geometry;
- the Molien series and invariant degrees \(12,20,30\);
- the allowed scalar harmonic degrees;
- the scalar Laplace eigenvalues and quotient multiplicities;
- self-adjointness and the spectral measure of the scalar Laplacian.

**Not established by this result**
- an identification of invariant degrees \(12,20,30\) with the Pellis exponents \(2,3,5\);
- a derivation of the coefficients in Golden Spheres Resonance;
- a unique state-selection rule fixing all spectral weights;
- equivalence with a spinor, Dirac, vector or QLG operator;
- an intertwiner to the Q6/PRISM operator.

The appearance of \(\varphi\) in icosahedral geometry does not by itself derive the fine-structure constant. The spectral construction is an independent exact input for Objective I, not a completed physical derivation.

## 6. Reproducibility and source boundary

The invariant multiplicities can be reproduced from the Molien series or the finite-group character average. The scalar spectrum then follows from \(\lambda_l=l(l+2)\) and \(m_l=(l+1)a_l\). The repository's explicit quaternion implementation of \(2I\) provides a finite computational group audit; the scalar spectral formula above is the representation-theoretic statement.

The separate numerical observation that a recovered RGS connection lies close to \(2I\) is not used in this theorem: proximity of link matrices to a finite subgroup is a numerical candidate, not a proof that the connection has \(2I\) holonomy.


## 7. Explicit quaternion implementation and reproduction

The quaternion implementation is part of this repository, not merely an external dependency:
- implementation: `src/binary_icosahedral.py`;
- construction and audit specification: `docs/binary-icosahedral-2i.md`;
- deterministic construction: the 8 coordinate units, 16 half-coordinate points, and 96 even-permutation/sign points of the 600-cell vertex set;
- multiplication convention: scalar-first Hamilton quaternions (q=(w,x,y,z)), with active rotations on column vectors.

Reproduce the implementation audit with:

```bash
python -c "from src.binary_icosahedral import audit_binary_icosahedral; print(audit_binary_icosahedral())"
pytest -q
```

The finite audit checks all (120^2=14,400) products, inverses, unit norms, and the two-to-one map to 60 rotations. This is an exhaustive floating-point implementation audit of the explicitly generated set, not a proof-assistant certificate. The scalar spectral formula and Molien multiplicities remain the representation-theoretic result; the computational group audit supplies reproducibility, not a new proof of the quotient spectrum.
