# Pellis D6 decomposition — Step 1: the three-mode module is 1 ⊕ 2

**Status:** PROVED (finite-dimensional representation theory). This note closes only the decomposition of the specified three-coordinate permutation representation. It makes no claim about a physical Pellis mode space.

Let
\[
C_3(x_0,x_1,x_2)=(x_2,x_0,x_1),\qquad
F_3(x_0,x_1,x_2)=(x_0,x_2,x_1),
\]
so that \(C_3^3=F_3^2=I\) and \(F_3C_3F_3=C_3^{-1}\). These generate the standard three-dimensional permutation representation of \(D_3\) (equivalently the quotient action relevant to the three-mode factor).

Set
\[
u=\frac{1}{\sqrt3}(1,1,1),\qquad
V_0=\mathbb C u,\qquad
V_2=\{x\in\mathbb C^3:x_0+x_1+x_2=0\}.
\]
Both subspaces are invariant under \(C_3\) and \(F_3\), are orthogonal, and
\[
\mathbb C^3=V_0\oplus V_2,\qquad \dim V_0=1,\quad\dim V_2=2.
\]
On \(V_0\), both generators act trivially.

To prove that \(V_2\) is irreducible over \(\mathbb C\), let \(\omega=e^{2\pi i/3}\) and define
\[
v_+=(1,\omega,\omega^2),\qquad
v_-=(1,\omega^2,\omega).
\]
They form a basis of \(V_2\), with
\[
C_3v_+=\omega v_+,\quad C_3v_-=\omega^2v_-,
\qquad F_3v_+=v_-,\quad F_3v_-=v_+.
\]
Since the two eigenvalues of \(C_3|_{V_2}\) are distinct, every one-dimensional \(C_3\)-invariant subspace is one of \(\mathbb Cv_+\) or \(\mathbb Cv_-\). The reflection \(F_3\) exchanges these two lines, so neither is invariant under both generators. A two-dimensional complex representation has no other possible proper nonzero subspace. Hence \(V_2\) is irreducible.

Therefore the representation decomposes as the direct sum of one trivial irreducible and one irreducible two-dimensional representation:
\[
\boxed{\mathbb C^3\cong\mathbf1\oplus\mathbf2.}
\]

The corresponding orthogonal projections are
\[
P_0=uu^*=\frac13
\begin{pmatrix}1&1&1\\1&1&1\\1&1&1\end{pmatrix},
\qquad P_2=I-P_0.
\]
They satisfy \(P_i^2=P_i=P_i^*\), \(P_0P_2=0\), \(P_0+P_2=I\), and commute with both representation generators.

## Scope and what this does not prove

1. This is an exact algebraic decomposition of the stated finite permutation representation; it is not a numerical spectral observation.
2. It does not establish that the physical GSR/Pellis model contains these three modes.
3. It does not yet derive the full operator commutant or any radial transfer-matrix form. Those are separate subsequent steps.
4. In particular, no physical degeneracy or equality of radial blocks is asserted here.

A focused executable check is provided by `tests/test_pellis_step1_d6_decomposition.py`. It checks only the invariant splitting and the generator action on the two complex eigenvectors; it is not a substitute for the proof above.
