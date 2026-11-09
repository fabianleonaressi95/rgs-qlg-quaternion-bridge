# Exact dyadic subdivision benchmark for the canonical Q6 2-skeleton

**Status:** exact combinatorial/geometric reference construction. **Not** identified with Kathara or Krist-All.

## 1. Coarse complex

The source complex is the oriented 2-skeleton of the unit six-cube:
\[
(C_0,C_1,C_2)=(\mathbb Z^{64},\mathbb Z^{192},\mathbb Z^{240}),
\qquad B_1B_2=0.
\]
The source orientation is inherited from src/q6_complex.py: edges point from the smaller six-bit integer to the larger; each square is oriented by its ordered axis pair \(a<b\). This theorem uses that source constructor's ordering; byte-for-byte or array-by-array identity with the recovered historical NPZ in data/q6/ is **not asserted**. A permutation/chain isomorphism to that artifact must be audited separately before transferring the subdivision matrices to its stored array order.

## 2. Geometric midpoint subdivision

Embed the six-cube as \([0,1]^6\). Store coordinates in doubled units, so every refined coordinate belongs to \(\{0,1,2\}\), representing \(0,1/2,1\).

The refined 2-skeleton consists of:
- vertices with at most two coordinates equal to 1 and all others in \(\{0,2\}\);
- axis-parallel edges with one varying coordinate and at most one other coordinate fixed at 1;
- squares with two varying coordinates and all other coordinates fixed in \(\{0,2\}\), each subdivided into four half-size squares.

The counts are exact:
\[
|V'|=\sum_{j=0}^2\binom6j2^{6-j}=64+192+240=496,
\]
\[
|E'|=6\cdot2\left(2^5+5\cdot2^4\right)=1344,
\qquad |F'|=\binom62\,2^4\,2^2=960.
\]
The 240 face centres are new vertices; the 960 internal midline segments are new edges.

## 3. Oriented subdivision maps

Define integer maps
\[
S_0:C_0\to C'_0,\quad S_1:C_1\to C'_1,\quad S_2:C_2\to C'_2.
\]
A coarse vertex maps to the same geometric vertex. A coarse oriented edge maps to the sum of its two positively oriented half-edges. A coarse oriented square maps to the sum of its four positively oriented child squares.

With the boundary matrices \(B_k\) and \(B'_k\), the subdivision maps satisfy
\[
\boxed{B'_1S_1=S_0B_1,\qquad B'_2S_2=S_1B_2.}
\]
These are exact integer identities. They say that subdivision commutes with boundary. They are checked by exact integer matrix multiplication in the tests.

The construction also satisfies \(B'_1B'_2=0\). The map column sums are respectively 1, 2 and 4, recording the number of oriented children of each parent cell.

## 4. Direction-count transfer

The coarse Q6 has 32 edges in each of its six coordinate directions. Aggregate the columns of \(S_1\) by the parent edge direction and its rows by the child edge direction, and divide each parent-direction count by 32. The resulting direction-count transfer is
\[
\boxed{D_{ij}=2\delta_{ij},\qquad D=2I_6.}
\]
Indeed, every parent edge has exactly two children, both parallel to that same coordinate axis. Hence, if three ordered trigram lines were identified with any three fixed coordinate axes, this particular midpoint rule would induce \(2I_3\) on direction counts, not the tribonacci matrix.

This is a statement about **oriented edge-direction counts**, not the Walsh-character compression \(E^T\mathcal R E\). A geometric identification between those two representations would itself have to be defined and proved.

## 5. What this does and does not establish

**Established:** a concrete, orientation-preserving parent-to-child subdivision of the source-defined Q6 square complex, with exact chain-map identities.

**Not established:** a Kathara/Krist-All subdivision, a map from physical Kathara cells to Q6 cells, or a three-channel transfer matrix on the trigram Walsh space. The midpoint subdivision preserves the coordinate-axis directions by construction; it supplies no geometric mechanism for the asymmetric cross-channel couplings in
\[
M_\psi=\begin{pmatrix}0&0&1\\1&0&0\\0&1&1\end{pmatrix}.
\]

Therefore this benchmark is a control construction: it proves what ordinary cubical midpoint refinement does, and prevents it from being mistaken for the sought Krist-All rule.

## 6. Reproduction

\`\`\`bash
PYTHONPATH=. python -m pytest tests/test_q6_dyadic_refinement.py -q
PYTHONPATH=. python examples/q6_dyadic_refinement_audit.py
\`\`\`
