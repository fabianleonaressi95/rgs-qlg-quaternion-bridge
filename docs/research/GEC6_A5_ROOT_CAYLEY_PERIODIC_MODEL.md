# GEC-6: the A5 root-lattice periodic model

**Status:** exact algebraic and finite-graph construction. This is a new symmetry-compatible periodic model, not a graph-isomorphism claim about the original PRISM product.

## 1. The infinite periodic geometry

Let
\[
V=\{x\in\mathbb R^6:\sum_{i=0}^5x_i=0\}
\]
with the Euclidean metric inherited from \(\mathbb R^6\), and let
\[
A_5=V\cap\mathbb Z^6
=\{x\in\mathbb Z^6:\sum_i x_i=0\}.
\]
The flat torus
\[
\mathbb T_{A_5}=V/A_5
\]
is a compact, connected, five-dimensional flat torus. Its fundamental group and first integral homology are
\[
\pi_1(\mathbb T_{A_5})\cong H_1(\mathbb T_{A_5};\mathbb Z)
\cong A_5\cong\mathbb Z^5.
\]

The symmetric group \(S_6\) acts on \(V\) by permuting the six ambient coordinates. This action preserves the Euclidean metric and the lattice \(A_5\), so it descends to isometries of \(\mathbb T_{A_5}\). Its restriction gives a faithful \(A_6\) action. Thus the rank-five periodic space does realize the required six-coordinate symmetry, unlike the original Cartesian PRISM graph.

## 2. Root-step graph and finite periodic quotients

The roots of \(A_5\) are
\[
R=\{e_i-e_j:0\le i\ne j\le5\},\qquad |R|=30.
\]
They generate the lattice \(A_5\), and \(S_6\) permutes this root set.

For any integer \(m\ge3\), define the finite quotient group
\[
G_m=A_5/mA_5\cong(\mathbb Z/m\mathbb Z)^5
\]
and its undirected Cayley graph
\[
\Gamma_m=\operatorname{Cay}(G_m,\{\bar r:r\in R\}).
\]
Because the 30 roots remain distinct and nonzero modulo \(mA_5\) for \(m\ge3\), \(\Gamma_m\) is a simple 30-regular graph. It is connected because the roots generate \(A_5\). It has
\[
|V(\Gamma_m)|=m^5,\qquad |E(\Gamma_m)|=15m^5.
\]
For the smallest illustrative quotient \(m=3\), these are 243 vertices and 3,645 undirected edges.

Every coordinate permutation preserves \(R\), hence induces a graph automorphism of \(\Gamma_m\). The action is faithful for \(m\ge3\): the 30 root classes are distinct modulo \(mA_5\), and the permutation action on the root set is faithful. In particular, \(A_6\leq\operatorname{Aut}(\Gamma_m)\).

This is a finite periodic discretization of the five-dimensional \(A_5\) torus, with exact \(S_6\) (therefore \(A_6\)) symmetry.

## 3. Explicit integer coordinates

Use the simple-root basis
\[
\alpha_k=e_k-e_{k+1},\qquad k=0,\ldots,4.
\]
For \(y\in A_5\), its integer coordinates in this basis are
\[
q_k=\sum_{j=0}^{k}y_j.
\]
Conversely,
\[
y_0=q_0,\quad y_i=q_i-q_{i-1}\ (1\le i\le4),\quad y_5=-q_4.
\]
These formulas implement the lattice, the root steps, and the \(S_6\) action using only integer arithmetic. Reduction modulo \(m\) gives the finite quotient action.

## 4. Link to Q6 and the original PRISM

The already established composite map
\[
Z_1(Q_6^{(1)};\mathbb Z)\twoheadrightarrow A_5
\]
is a surjective homomorphism of \(S_6\)-modules. The target \(A_5\) is now realized as the fundamental lattice of \(\mathbb T_{A_5}\), and the finite graph \(\Gamma_m\) supplies explicit periodic quotients carrying the same symmetry.

The original PRISM product
\[
P=C_{64}\square C_6^3\square C_5
\]
also has \(H_1(P;\mathbb Z)\cong\mathbb Z^5\), so its homology lattice is abstractly isomorphic to \(A_5\). But its graph automorphism action on this homology cannot realize the faithful \(A_6\) action: its induced image has order at most \(2^5\cdot3!=192<360=|A_6|\), and simplicity of \(A_6\) forces every such homomorphism to be trivial.

Therefore the new model resolves the symmetry-realization problem by replacing the Cartesian PRISM graph with an \(A_5\)-root Cayley geometry. It does **not** prove that the old PRISM graph is isomorphic, metrically equivalent, or dynamically equivalent to \(\Gamma_m\). In particular, the illustrative \(m=3\) graph has 243 vertices and degree 30, unlike the 69,120-vertex degree-10 PRISM product.

## 5. Relation to the Q6 filled cube

The source of the quotient map is the Q6 graph cycle lattice, not the first homology of the filled cubical 6-cube. Filling all 240 square faces makes the cube contractible and kills \(H_1\). No nonzero \(H_1\) map from that filled cube is claimed.

## 6. Evidence ledger

| Claim | Status |
|---|---|
| \(A_5\) is a rank-five integral lattice | PROVED |
| \(\mathbb T_{A_5}=V/A_5\) is a compact flat 5-torus | PROVED |
| \(S_6\) acts isometrically and faithfully on the torus and on \(A_5\) | PROVED |
| \(\Gamma_m\) is connected, simple, 30-regular for \(m\ge3\) | PROVED |
| \(S_6\) acts faithfully by graph automorphisms on \(\Gamma_m\) | PROVED |
| The Q6 graph-cycle quotient maps onto the periodic lattice \(A_5\) | PROVED, using the preceding GEC-6 construction |
| The original Cartesian PRISM graph realizes this \(A_6\) action | IMPOSSIBLE under its graph-automorphism action on \(H_1\) |
| Equivalence of the new root graph and original PRISM dynamics/spectra/physics | OPEN; not claimed |

## 7. Reproduction

\`\`\`bash
PYTHONPATH=. python -m pytest tests/test_gec6_a5_root_cayley_torus.py -q
PYTHONPATH=. python examples/gec6_a5_root_cayley_audit.py
\`\`\`
