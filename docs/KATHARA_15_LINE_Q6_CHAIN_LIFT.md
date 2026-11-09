# Kathara 15-line structure and the integer lift to Q6

**Date:** 2026-11-09  
**Status:** exact combinatorial deductions plus an explicitly open geometric assignment. No numerical optimization or tests are claimed.

## 1. Why the source must retain 15 primary lines

The source transcription in `src/kathara_15_line_graph.py` contains 12 centres, 15 named primary lines and 21 unique elementary edges. Expanding the three vertical polylines gives

- L1: (1,2), (2,5), (5,8), (8,11), (11,12);
- L2: (3,6), (6,9);
- L3: (4,7), (7,10);
- L4..L15: the twelve individual cross-links already recorded in `PRIMARY_LINES`.

These 21 edges are pairwise distinct across the 15 lines. Consequently the line-edge incidence matrix
\[
A\in\{0,1\}^{15\times21},\qquad
A_{\ell e}=1\iff e\text{ belongs to }L_\ell
\]
has exactly one nonzero entry in every column:
\[
\sum_{\ell=1}^{15}A_{\ell e}=1\quad(e=1,\ldots,21).
\]
This is stronger and cleaner than a generic overlapping hyperedge model: the 15 primary lines form a partition of the 21 elementary edges, while retaining their ordered polyline structure.

For a line-constant edge weight vector \(a\in\mathbb R^{15}\), the induced elementary-edge weights are
\[
w=A^T a,\qquad w_e=a_{\ell(e)}.
\]
This formula introduces no numerical hierarchy. The source motivates distinguished categories (central vertical, two lateral verticals, cross-links), but it does not determine numerical values of \(a\). In particular, equal weights within a source-defined line are a transparent structural convention, not a physical law.

## 2. The induced reflection on all 15 lines

The planar reflection
\[
\sigma=(3\ 4)(6\ 7)(9\ 10)
\]
fixes centres 1,2,5,8,11,12 and swaps the two lateral branches. Its action on primary-line IDs is
\[
\rho=(2\ 3)(4\ 6)(5\ 7)(8\ 10)(9\ 11)(12\ 13)(14\ 15),
\]
with line 1 fixed. Thus the line-edge incidence is equivariant:
\[
A_{\rho(\ell),\,\sigma(e)}=A_{\ell e}.
\]
Here \(\sigma(e)\) denotes the unordered edge obtained by applying \(\sigma\) to both endpoints. This gives an exact test for any proposed lift: if the target action is specified, the vertex, edge-path and face-chain maps must intertwine these source actions. No target subgroup is presumed.

## 3. A strict graph embedding into Q6 is obstructed

The Kathara 1-skeleton contains 3-cycles, for example
\[
1\to2\to3\to1
\]
(using the edges (1,2), (2,3), (1,3)). The Q6 1-skeleton is bipartite, since every edge changes the parity of the six-bit Hamming weight. Therefore Q6 has no odd cycles.

It follows that there is no graph map sending every Kathara elementary edge to one Q6 edge while preserving all three edges of this triangle. In particular, an injective edge-preserving graph embedding is impossible.

This is not an obstruction to a cellular chain map. An elementary Kathara edge may map to an integral 1-chain supported on a Q6 path with the same endpoints. The existing Kathara-to-Q3 construction already uses this path-chain notion.

## 4. The correct integer chain-map problem

Let \(K\) be the chosen Kathara cellular complex (12 vertices, 21 edges, 10 filled bounded faces) and let \(Q_6^{(2)}\) be the cubical 2-skeleton (64 vertices, 192 edges, 240 square faces). Fix the canonical boundary conventions in `src/q6_complex.py`.

A candidate lift consists of integer matrices
\[
R_0\in\mathbb Z^{64\times12},\quad
R_1\in\mathbb Z^{192\times21},\quad
R_2\in\mathbb Z^{240\times10},
\]
satisfying exactly
\[
B_1^{Q6}R_1=R_0B_1^K,\qquad
B_2^{Q6}R_2=R_1B_2^K.
\]
The first equation requires every column of \(R_1\) to be an integral path-chain (or a more general 1-chain) between the images of the corresponding Kathara endpoints. The second fills each mapped Kathara face boundary by Q6 squares.

The matrix \(R_0\) must be specified separately. If it is a vertex function, each column is a standard basis vector; if distinct Kathara centres must remain distinct, that is an additional injectivity constraint. Neither condition is implied by the chain equations.

## 5. Existence and non-uniqueness of face fillings

The cubical 2-skeleton of the 6-cube has \(H_1=0\). One way to see this is that adding 3-cells to a 2-skeleton cannot change its first homology, while the full cube is contractible. Hence every integral 1-cycle in \(Q_6^{(2)}\) is an integral boundary of 2-cells.

Therefore, for any chosen vertex map \(R_0\) and any choice of edge 1-chains \(R_1\) satisfying the first chain identity, every column of \(R_1B_2^K\) is a cycle and admits an integer square filling. Thus an \(R_2\) exists.

But this does not make the filling unique. Since the target is only the 2-skeleton, it has no 3-boundary map to quotient out its 2-cycles. The ranks are
\[
\operatorname{rank}B_1^{Q6}=63,\qquad
\operatorname{rank}B_2^{Q6}=192-63=129,
\]
and hence
\[
\operatorname{rank}H_2(Q_6^{(2)})
=240-129=111.
\]
For each source 2-cell, two fillings can therefore differ by an element of the rank-111 integer cycle lattice (more precisely, by an integral 2-cycle). The ten columns of \(R_2\) must not be treated as canonical merely because they solve the boundary equation.

This is a key distinction between:
- **existence:** guaranteed once \(R_0,R_1\) satisfy the first chain identity;
- **uniqueness:** false without an additional selection principle;
- **RGS canonicity:** open until an invariant geometric/tangent/spectral criterion is defined on these lifts.

## 6. How the 15-line structure constrains, but does not yet select, a lift

The line partition supplies a natural bookkeeping constraint on \(R_1\). Write \(R_1=[r_e]_{e\in E(K)}\). For every primary line \(L_\ell=(v_0,\ldots,v_m)\), define its lifted line chain
\[
\mathcal R_\ell=\sum_{j=0}^{m-1}r_{(v_j,v_{j+1})},
\]
with signs adjusted to the fixed orientation of the Kathara edge basis. Its boundary is forced:
\[
\partial\mathcal R_\ell=R_0(v_m)-R_0(v_0).
\]
The line labels therefore distinguish the three vertical paths from the twelve cross-links even when their target chains have the same endpoints or overlap.

A legitimate next constraint is to require the three lifted vertical lines to be treated as ordered path chains, and to test the source reflection \(\sigma\) against the line permutation \(\rho\). A stronger condition—such as pairwise disjoint target supports, shortest paths, or a chosen ordering of bit flips—must be declared as an extra hypothesis, not smuggled into the definition.

## 7. What the RGS seed contributes—and what it does not

The recovered seed loader in `src/rgs_ellis_refinement.py` constructs 64 coordinates in the repository's six-bit state order and verifies the archived coordinate, edge and coframe hashes. The Q6 cell order in `src/q6_complex.py` uses a separate explicit convention (six-bit big-endian axis order). Any implementation must compare the endpoint sets and orientation signs explicitly; it must not assume that two 192-edge arrays share the same ordering.

The RGS seed provides a geometric realization of the 64 labelled Q6 states. It does not by itself define which 12 of those states are the images of the 12 Kathara centres, nor which Q6 path chains represent the 21 Kathara edges. A nearest-coordinate assignment or a threshold on the planar Kathara drawing would add an external convention. The 15-line incidence and reflection constraints can filter assignments, but no unique \(R_0,R_1\) has yet been derived from them.

## 8. Implementation gates

1. **Incidence normalization:** construct and export \(A\), the 15 line paths, the 21-edge ordering, and the exact source permutations \(\sigma,\rho\).
2. **Q6 order reconciliation:** construct an explicit bijection between the loader's 192 edge records and the canonical Q6 edge ordering, including orientation signs.
3. **Vertex assignment search:** enumerate or characterize \(R_0\) candidates using only explicitly stated RGS geometric and source-symmetry constraints. Report all minimizers/candidates; do not use arbitrary nearest-neighbour tie-breaking.
4. **Path-chain lift:** construct integer \(R_1\) and check \(B_1R_1=R_0B_1^K\) exactly. Preserve the ordered 15-line decomposition.
5. **Square filling:** solve \(B_2R_2=R_1B_2^K\) over integers and report the rank-111 filling ambiguity.
6. **Equivariance and RGS gate:** impose only a target symmetry subgroup actually verified on Q6 cells; then compare tangent/spectral functionals using projectors, never a chosen tangent frame.

## 9. Evidence ledger

| Claim | Status |
|---|---|
| 12 centres, 15 named lines, 21 elementary edges | **SOURCE-TRANSCRIBED IN REPOSITORY** |
| Each elementary edge belongs to exactly one primary line | **EXACT FROM THE TRANSCRIPTION** |
| Reflection action \(\sigma\) and induced line permutation \(\rho\) | **EXACT COMBINATORIAL DEDUCTION** |
| No edge-preserving graph embedding of the Kathara 1-skeleton into Q6 | **PROVED BY ODD-CYCLE / BIPARTITE OBSTRUCTION** |
| Integer path-chain map can exist for a chosen vertex map and compatible edge chains | **ALGEBRAICALLY CHARACTERIZED** |
| Q6 2-skeleton has \(H_1=0\), rank-111 \(H_2\) | **EXACT CUBICAL-HOMOLOGY DEDUCTION** |
| RGS-derived unique \(R_0,R_1,R_2\) | **OPEN** |
| Canonical RGS filling among the rank-111 ambiguity | **OPEN** |

No tests, numerical search, or repository execution are claimed in this note.
