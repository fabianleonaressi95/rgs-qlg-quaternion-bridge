# Refinement obstruction audit: what Q6 incidence does and does not determine

**Purpose.** Continue the trigram/Krist-All refinement programme using the canonical Q6 incidence bundle already stored in this repository. This note adds a rigorous obstruction analysis before any geometric transfer matrix is claimed.

## 1. Exact input already present

The canonical Q6 bundle, documented in data/q6/README.md, records
\[
B_1\in\mathbb Z^{64\times192},\qquad B_2\in\mathbb Z^{192\times240},
\]
with \(B_1B_2=0\), \(\operatorname{rank}B_1=63\), and \(\operatorname{rank}B_2=129\). These are oriented 0-, 1-, and 2-cell incidence matrices of the six-cube. The six-bit vertex set factors as
\[
Q_6=Q_3^{\rm lower}\square Q_3^{\rm upper}.
\]
They specify a finite chain complex, but not a parent-to-child refinement: no pair of levels, subdivision maps \(S_k\), or geometric child-cell incidence list is encoded by \(B_1,B_2\) alone.

## 2. Theorem 1 — Translation-equivariant obstruction

Let \(T=\mathbb Z_2^3\), let \(\tau_g f(t)=f(t\oplus g)\) be its regular translation action, and let \(\mathcal R:\mathbb R^T\to\mathbb R^T\) commute with every translation:
\[
\mathcal R\tau_g=\tau_g\mathcal R\qquad(g\in T).
\]
Then every Walsh character \(\chi_S\) is an eigenvector of \(\mathcal R\). In particular, on the first-order Walsh subspace \(\mathcal H_1\), the compressed matrix \(E^T\mathcal R E\) is diagonal.

**Proof.** The Walsh characters are the one-dimensional irreducible representations of the finite abelian group \(T\). A linear map commuting with the regular representation preserves each character space. Since every character space is one-dimensional, \(\mathcal R\chi_S=\rho_S\chi_S\). The three first-order characters are therefore individually invariant and the compression in their ordered basis is diagonal. \(\square\)

The tribonacci target
\[
M_\psi=\begin{pmatrix}0&0&1\\1&0&0\\0&1&1\end{pmatrix}
\]
is not diagonal. Consequently, a refinement producing this target cannot be fully equivariant under all eight translations of the trigram cube. A successful derivation must identify the physical orientation/anchoring data that break full Q3 translation symmetry.

## 3. Theorem 2 — Coordinate-permutation obstruction

Let \(S_3\) act on \(T\) by permuting the three line coordinates. If \(\mathcal R\) commutes with this action, then its compression \(M=E^T\mathcal R E\) commutes with every 3x3 permutation matrix. Hence
\[
M=\alpha I_3+\beta(J_3-I_3)
\]
for some scalars \(\alpha,\beta\), where \(J_3\) is the all-ones matrix.

**Proof.** The first-order Walsh space transforms as the three-coordinate permutation representation. Equivariance implies \(PM=MP\) for every permutation matrix \(P\). A matrix commuting with all coordinate permutations has a common diagonal entry and a common off-diagonal entry, giving the displayed form. \(\square\)

The tribonacci target is not of this form: its diagonal entries are \(0,0,1\), and its off-diagonal entries are not all equal. Thus the target requires a distinguished ordering/orientation of the three trigram lines. Full S3 symmetry cannot survive.

## 4. Exact audit of the existing Q6 graph adjacency

Write the Q6 adjacency as
\[
A_6=A_3\otimes I_8+I_8\otimes A_3.
\]
Let \(u_0=8^{-1/2}(1,\ldots,1)^T\) be the normalized constant mode of the upper trigram factor, and let \(E\) be the 8x3 first-order Walsh matrix on the lower factor. Since \(A_3E=E\) and \(A_3u_0=3u_0\),
\[
A_6(E\otimes u_0)=4(E\otimes u_0).
\]
Therefore the compression of \(A_6\) to the lower three channels with the upper factor in its constant mode is exactly
\[
\boxed{4I_3,}
\]
not \(M_\psi\). This rules out unmodified Q6 graph adjacency as the desired transfer operator, but not a different geometrically defined subdivision operator.

## 5. Theorem 3 — Incidence data alone do not identify a refinement operator

Let \(E\in\mathbb R^{8\times3}\) have orthonormal columns. For any prescribed 3x3 matrix \(M\) and scalar \(\lambda\), define
\[
\mathcal R_{M,\lambda}=E M E^T+\lambda(I_8-EE^T).
\]
Then
\[
E^T\mathcal R_{M,\lambda}E=M,\qquad
(I_8-EE^T)\mathcal R_{M,\lambda}E=0.
\]

**Proof.** Use \(E^TE=I_3\) and \(E^T(I_8-EE^T)=0\). \(\square\)

Thus the three-channel compression and closure equations alone do not select \(M_\psi\): every 3x3 matrix can be realized algebraically. The geometric incidence/refinement construction is indispensable. Defining \(\mathcal R=EM_\psi E^T+\lambda(I-EE^T)\) would be circular and is not accepted as a derivation.

## 6. Required geometric input

A valid refinement must specify at every level:
1. parent cells and their dimensions/orientations;
2. child cells and their dimensions/orientations;
3. signed subdivision coefficients \(S_k(\sigma,\tau)\);
4. geometric embedding or coordinate map of every child in its parent;
5. shared-face identifications and affine transition maps;
6. verification of \(B_k^{(n+1)}S_k=S_{k-1}B_k^{(n)}\);
7. induced action on trigram labels and its first-order Walsh compression.

Only then can one compute \(M_K=E^T\mathcal R E\) from geometry and test equality with \(M_\psi\), including zero leakage.

## 7. Status

| Claim | Status |
|---|---|
| Q6 incidence bundle defines the recorded finite cubical chain complex | Exact repository data |
| Full trigram translation equivariance forces diagonal Walsh compression | Proved |
| Full coordinate-permutation equivariance forces alpha I + beta(J-I) | Proved |
| Existing Q6 adjacency induces 4I on lower first-order modes with upper constant mode | Proved |
| Incidence matrices alone determine a parent-to-child refinement | False: subdivision maps are absent |
| A Kathara/Krist-All geometric refinement induces Mpsi | Open |
| Physical orientation/anchoring supplies the required symmetry breaking | Open |

## 8. Reproduction

~~~bash
python -m pytest tests/test_trigram_refinement_obstructions.py -q
python examples/trigram_refinement_obstruction_audit.py
~~~

The audit checks exact finite-dimensional obstruction identities. It does not fabricate or claim a geometric subdivision.
