# Q6–D5: A5-invariant angular pushforward measure

## Setup

The five-dimensional coexact Q6 eigenspace is an orthogonal irreducible
representation \(V_5\) of \(A_5\). The existing intertwiner
\[
T:\Lambda^2(V_5)\longrightarrow\mathbb R^3
\]
is equivariant. On the ten-dimensional cotangent phase space write
\((q,p)\in\mathbb R^5\oplus\mathbb R^5\), and define
\[
F(q,p)=T(q\wedge p),\qquad
n(q,p)=\frac{F(q,p)}{\|F(q,p)\|}
\]
where \(F\ne0\).

Use the normalized surface probability measure on the unit phase sphere
\[
S^9=\{(q,p):\|q\|^2+\|p\|^2=1\}.
\]
This is an explicit compact, rotationally invariant phase-shell measure.
It is not silently identified with an unbounded Liouville measure or with
the full nonlinear D5 Bergman volume.

## Exact symmetry theorem

For each \(g\in A_5\), the phase-space action is
\[
(q,p)\mapsto(\rho_gq,\rho_gp),\qquad \rho_g\in O(5).
\]
It preserves the normalized surface measure on \(S^9\). Equivariance of
\(T\) gives
\[
F(\rho_gq,\rho_gp)=R_gF(q,p),\qquad R_g\in SO(3).
\]
Thus, wherever the normalized map is defined, its pushforward is
\(A_5\)-invariant. If the zero locus \(F^{-1}(0)\cap S^9\) has surface
measure zero, this defines an \(A_5\)-invariant probability measure
\(\nu\) on \(S^2\).

The dimension of the fixed subspace of degree-\(\ell\) spherical
harmonics is the character average
\[
m_\ell=\frac1{60}\left[
(2\ell+1)+15\chi_\ell(\pi)+20\chi_\ell(2\pi/3)
+12\chi_\ell(2\pi/5)+12\chi_\ell(4\pi/5)
\right],
\]
where
\[
\chi_\ell(\theta)=
\frac{\sin((\ell+1/2)\theta)}{\sin(\theta/2)}
\]
and \(\chi_\ell(0)=2\ell+1\).

The exact character average gives
\[
(m_0,m_1,m_2,m_3,m_4,m_5,m_6)=(1,0,0,0,0,0,1).
\]
Consequently every \(A_5\)-invariant probability measure on \(S^2\)
has vanishing nonconstant spherical-harmonic coefficients in degrees
1 through 5. Equivalently, its polynomial moments through degree five
agree with those of the uniform sphere. This applies to the continuous
pushforward \(\nu\), conditional only on the zero locus being null.

## What this does not prove

- It does not prove that the 32-point Q6 design and \(\nu\) are the same
  measure. They share moments through degree five because both are
  \(A_5\)-invariant, but their degree-six coefficients can differ.
- Degree six is the first possible anisotropy. Its invariant harmonic
  subspace is one-dimensional; symmetry does not fix its coefficient.
- The global zero locus of \(T(q\wedge p)\) has not yet been classified.
  A generic-rank calculation is evidence, not a global measure-zero proof.
- This is a compact phase-shell construction, not a physical 3+1D limit,
  a global D5 Lie-ball embedding, a derivation of Wyler's volume, or a
  prediction of the fine-structure constant.

## Epistemic status

**Proved algebraically:** the A5 character formula and the fixed-harmonic
multiplicities through degree six.

**Numerically verified by the audit:** orthogonality of the five-mode action,
proper-rotation character of the geometric action, and the existing wedge
intertwiner residual.

**Conditional theorem:** the normalized S9 pushforward has uniform S2
moments through degree five if its zero locus has S9 measure zero.

**Open:** global zero-locus classification and comparison of the continuous
and finite degree-six coefficients.
