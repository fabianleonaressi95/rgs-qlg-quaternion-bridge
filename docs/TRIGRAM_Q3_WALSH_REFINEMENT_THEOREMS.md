# Trigram Q3: Walsh decomposition, face compatibility, and refinement theorem gates

**Status:** exact finite-dimensional algebraic theorems are established below. The identification of a geometric Kathara/Krist-All refinement with the tribonacci matrix is **OPEN** until its incidence/refinement map is constructed from the actual cells.

## 1. Conventions and trigram cube

Write a trigram as \(t=(a,b,c)\in T=\mathbb Z_2^3\), with coordinates ordered from the bottom line to the top line. Addition is modulo two. The three coordinate flips are \(t\mapsto t\oplus e_j\).

The eight binary triples are the vertices of the cube graph \(Q_3\). Two vertices are adjacent exactly when they differ in one coordinate. This is a combinatorial identification; it does not by itself identify the abstract cube of trigram labels with a physical three-dimensional Kathara cell complex.

For the usual I-Ching naming convention with bit order bottom-to-top and 0 = yin, 1 = yang:

| Bits | Trigram |
|---|---|
| 000 | Kun ☷ |
| 001 | Gen ☶ |
| 010 | Kan ☵ |
| 011 | Xun ☴ |
| 100 | Zhen ☳ |
| 101 | Li ☲ |
| 110 | Dui ☱ |
| 111 | Qian ☰ |

## 2. Theorem 1 — Walsh basis and canonical three-channel subspace

Let \(\mathcal H_T=\mathbb R^T\), with normalized inner product
\[
\langle f,g\rangle_T=\frac18\sum_{t\in T}f(t)g(t).
\]
For every subset \(S\subseteq\{1,2,3\}\), define the Walsh character
\[
\chi_S(t)=(-1)^{\sum_{j\in S}t_j}.
\]
Then the eight functions \(\chi_S\) form an orthonormal basis of \(\mathcal H_T\). Consequently,
\[
\mathcal H_T=\mathcal H_0\oplus\mathcal H_1\oplus
\mathcal H_2\oplus\mathcal H_3,
\]
where
\[
\begin{aligned}
\mathcal H_0&=\operatorname{span}\{1\},\\
\mathcal H_1&=\operatorname{span}\{\chi_{\{1\}},\chi_{\{2\}},\chi_{\{3\}}\},\\
\mathcal H_2&=\operatorname{span}\{\chi_{\{1,2\}},\chi_{\{1,3\}},\chi_{\{2,3\}}\},\\
\mathcal H_3&=\operatorname{span}\{\chi_{\{1,2,3\}}\}.
\end{aligned}
\]
The dimensions are \(1,3,3,1\). In particular, the three first-order Walsh modes define an intrinsic three-dimensional channel space; they are not extra variables appended to the trigram set.

**Proof.** For subsets \(S,R\),
\[
\langle\chi_S,\chi_R\rangle_T
=\frac18\sum_{t\in T}\chi_{S\triangle R}(t).
\]
The sum is one if \(S=R\) and zero otherwise, by cancellation along any coordinate in \(S\triangle R\). Thus the eight characters are orthonormal; since \(\dim\mathcal H_T=8\), they are a basis. Grouping by subset cardinality gives the stated decomposition. \(\square\)

Let \(E\in\mathbb R^{8\times3}\) have rows indexed by \(t\in T\), in lexicographic order \(000,001,\ldots,111\), and entries
\[
E_{t,j}=\frac1{\sqrt8}\chi_{\{j\}}(t).
\]
Then \(E^TE=I_3\), and \(P_1=EE^T\) is the Euclidean orthogonal projector onto the three first-order channels. The normalized inner product above is equivalent to using the unscaled characters as an orthonormal basis.

## 3. Theorem 2 — Spectrum of the trigram-cube adjacency operator

Define
\[
(Af)(t)=\sum_{j=1}^3f(t\oplus e_j).
\]
For every Walsh character,
\[
A\chi_S=(3-2|S|)\chi_S.
\]
Therefore
\[
\operatorname{Spec}(A)=\{3,1^{(3)},-1^{(3)},-3\},
\qquad A|_{\mathcal H_1}=I_3.
\]

**Proof.** Flipping coordinate \(j\) multiplies \(\chi_S\) by \(-1\) if \(j\in S\), and leaves it unchanged otherwise. Summing the three contributions gives \(3-2|S|\). Multiplicities follow from the number \(\binom3{|S|}\) of subsets of each cardinality. \(\square\)

This theorem concerns adjacency only. In particular, adjacency does not mix the three first-order channels and does not generate a tribonacci transfer matrix.

## 4. Theorem 3 — Exact criterion for a closed three-channel refinement

Let \(\mathcal R:\mathbb R^8\to\mathbb R^8\) be any linear operator representing a proposed refinement on trigram-indexed data. Define its compression to the three channels by
\[
M_K=E^T\mathcal R E.
\]
The three-channel subspace is invariant under \(\mathcal R\) if and only if
\[
\boxed{(I_8-EE^T)\mathcal R E=0.}
\]
When this holds,
\[
\boxed{\mathcal R E=E M_K.}
\]

**Proof.** Decompose
\[
\mathcal R E=EE^T\mathcal R E+(I_8-EE^T)\mathcal R E.
\]
The second term is precisely the component orthogonal to the three-channel subspace. It vanishes exactly when the image of every channel lies in that subspace. In that case the first term equals \(E(E^T\mathcal R E)=EM_K\). The converse is immediate. \(\square\)

This is a criterion, not a construction of \(\mathcal R\). Its entries must be obtained from the oriented cell-refinement/incidence data.

## 5. Theorem 4 — Affine compatibility on a shared face

Let an affine face be
\[
F=P_0+\operatorname{span}(u_1,\ldots,u_d)
\]
and let two affine rules be
\[
\Phi_\alpha(P)=A_\alpha P+b_\alpha,\qquad
\Phi_\beta(P)=A_\beta P+b_\beta.
\]
They agree at every point of \(F\) if and only if
\[
(A_\alpha-A_\beta)P_0+(b_\alpha-b_\beta)=0
\]
and
\[
(A_\alpha-A_\beta)u_j=0,\qquad j=1,\ldots,d.
\]

**Proof.** Substitute \(P=P_0+\sum_js_ju_j\) into \(\Phi_\alpha(P)-\Phi_\beta(P)\). The result is affine in the free parameters \(s_j\); it vanishes identically if and only if its constant term and every linear coefficient vanish. \(\square\)

The theorem applies to a specified geometric shared face. A square face of the abstract label cube \(Q_3\) is a four-element set of trigram labels; it must not be identified with a physical Kathara face unless an explicit geometric incidence map makes that identification.

## 6. Corollary — Refinement must commute with boundary

Let \(C_k^{(n)}\) be chain spaces at level \(n\), \(B_k^{(n)}:C_k^{(n)}\to C_{k-1}^{(n)}\) boundary maps, and \(S_k:C_k^{(n)}\to C_k^{(n+1)}\) refinement maps. If
\[
B_k^{(n+1)}S_k=S_{k-1}B_k^{(n)}
\]
and \(B_{k-1}^{(n)}B_k^{(n)}=0\), then
\[
B_{k-1}^{(n+1)}B_k^{(n+1)}S_k=0.
\]
Thus refinement preserves the boundary-of-boundary identity on the image of the refinement map.

**Proof.**
\[
B_{k-1}^{(n+1)}B_k^{(n+1)}S_k
=B_{k-1}^{(n+1)}S_{k-1}B_k^{(n)}
=S_{k-2}B_{k-1}^{(n)}B_k^{(n)}=0.
\]
\(\square\)

## 7. Corollary — Q6 factorization into two trigram cubes

With the six bit coordinates divided into ordered lower and upper triples,
\[
Q_6=\mathbb Z_2^6\cong\mathbb Z_2^3\times\mathbb Z_2^3
=Q_3^{\mathrm{lower}}\square Q_3^{\mathrm{upper}}.
\]
Every hexagram is an ordered pair \((t_{\rm lower},t_{\rm upper})\) of trigrams, and \(8\cdot8=64\). This is an exact product decomposition of the hypercube vertex set and graph.

## 8. Tribonacci target: conditional theorem, not yet a geometric result

The proposed target transfer matrix is
\[
M_\psi=\begin{pmatrix}0&0&1\\1&0&0\\0&1&1\end{pmatrix}.
\]
Its characteristic polynomial is
\[
\det(\lambda I-M_\psi)=\lambda^3-\lambda^2-1.
\]
Its Perron eigenvalue is the positive real root
\[
\psi=1.465571231876768\ldots.
\]
If a geometrically constructed refinement \(\mathcal R\) satisfies both
\[
E^T\mathcal R E=M_\psi,\qquad
(I_8-EE^T)\mathcal R E=0,
\]
then its closed three-channel dynamics is exactly tribonacci and has asymptotic growth factor \(\psi\) for initial data with nonzero Perron component.

The spectral statement about \(M_\psi\) is exact. The existence of a Kathara/Krist-All geometric refinement satisfying these identities is **OPEN**. No coefficients of \(\mathcal R\) are prescribed here to force the desired answer.

## 9. Status ledger

| Statement | Status |
|---|---|
| Trigram labels form \(Q_3=\mathbb Z_2^3\) | Exact combinatorial identity |
| Walsh basis and \(1+3+3+1\) decomposition | Proved |
| Three canonical first-order channels | Proved |
| \(Q_3\) adjacency spectrum | Proved |
| Invariant-subspace criterion for a refinement | Proved |
| Affine shared-face compatibility criterion | Proved |
| Boundary/refinement commutation implication | Proved |
| \(Q_6=Q_3^{\rm lower}\square Q_3^{\rm upper}\) | Exact combinatorial identity |
| Tribonacci polynomial and Perron root for the target matrix | Proved |
| Actual Kathara/Krist-All refinement matrix \(\mathcal R\) | Open; requires geometric incidence data |
| Identification of abstract \(Q_3\) faces with physical Kathara faces | Open; requires an explicit map |
| Derivation of \(M_\psi\) from that geometry | Open |

## 10. Reproducibility

Run from the repository root:
\`\`\`bash
python -m pytest tests/test_trigram_q3_theorems.py -q
python examples/trigram_q3_theorem_audit.py
\`\`\`
The audit checks the exact finite matrices numerically to floating-point tolerance, including the Walsh orthogonality, adjacency spectrum, affine-face criterion on a sample face, Q6 product adjacency, and the target matrix's characteristic polynomial/eigenvalue. It does **not** claim to construct the missing geometric refinement.
