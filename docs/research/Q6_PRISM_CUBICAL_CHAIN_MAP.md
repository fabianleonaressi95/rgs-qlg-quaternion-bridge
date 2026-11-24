# An explicit cubical chain map from Q6 into the PRISM corner

**Status:** exact integral chain map, with the chain identity checked on all 729 source cells. It is not a homotopy equivalence to the full PRISM product and does not imply a transport law.

## 1. Target subcomplex

The PRISM graph/product is
\[
P=C_{64}\square C_6\square C_6\square C_6\square C_5.
\]
In each cycle choose the edge with vertices \(\{0,1\}\). Their Cartesian product is a cubical subcomplex \(Q_5\subset P\), with 32 vertices and 80 edges. We use target coordinates in factor order \((C_{64},C_6,C_6,C_6,C_5)\).

## 2. Vertex map and its geometric meaning

Write a Q6 vertex as \(v=(v_0,v_1,v_2,v_3,v_4,v_5)\in\{0,1\}^6\). Define
\[
F_0(v)=(v_0\oplus v_5,v_1,v_2,v_3,v_4)\in\{0,1\}^5.
\]
This is a surjective graph homomorphism \(Q_6\to Q_5\): flipping any one source bit flips exactly one target bit. It is two-to-one on vertices. The two source directions 0 and 5 are identified as the first target direction; the remaining four directions are retained separately.

This is materially stronger than mapping all six directions to a single parity edge: five independent target directions survive.

## 3. Integral cubical chain map

Represent an oriented \(k\)-cell of \(Q_6\) by \([v;S]\), where \(S\subseteq\{0,\ldots,5\}\) is the increasing list of free coordinates and \(v_d=0\) for \(d\in S\). Its boundary is
\[
\partial[v;S]=\sum_{j=1}^{k}(-1)^{j-1}
\bigl([v+e_{s_j};S\setminus\{s_j\}]-[v;S\setminus\{s_j\}]\bigr).
\]

If \(S\) contains both 0 and 5, set \(F_k[v;S]=0\): two source directions have folded onto one target direction, so the cell loses dimension.

Otherwise, replace whichever of 0 or 5 occurs in \(S\) by target direction 0, and retain directions 1–4. The resulting target direction set is sorted increasingly. The target basepoint is the image vertex \(F_0(v)\), reset to coordinate 0 in every free target direction. The coefficient is the orientation sign:
- if direction 0 varies, multiply by \((-1)^{v_5}\);
- if direction 5 varies, multiply by \((-1)^{v_0}(-1)^{|S|-1}\), where the second factor reorders the last source direction into the first target position;
- otherwise the coefficient is \(+1\).

This defines \(F_k\) on every cell and extends linearly over \(\mathbb Z\). The accompanying implementation verifies, with exact integer arithmetic,
\[
\boxed{\partial F_k=F_{k-1}\partial}
\]
for every source cell in all dimensions. There are \(3^6=729\) cubical cells in \(Q_6\), and the check is exhaustive, not sampled.

## 4. Exact consequences

- The vertex map is 2-to-1: 64 source vertices map onto all 32 target vertices.
- All 192 source edges map to target edges; all 80 target edges are reached.
- The chain map is surjective on target cell generators in every degree.
- Every source cell containing both folded directions 0 and 5 maps to zero.
- Since both cubes are contractible, the induced homology map is an isomorphism on \(H_0\) and is trivial in positive homological degrees (both positive-degree homology groups vanish).

The kernel is a genuine chain subcomplex because \(F\) is a chain map. It records exactly which combinations of cubical cells are lost under the fold.

## 5. The limitation that must remain explicit

The map lands in the binary-corner \(Q_5\), not in a way that wraps around the full cycles \(C_{64},C_6,C_6,C_6,C_5\). Consequently it does not capture the nontrivial torus homology of the full PRISM product, its global cycle generators, its Laplacian spectrum, or its transport coefficients. It also does not produce a homotopy equivalence between Q6 and PRISM.

The next mathematically distinct problem is to construct a chain map whose image reaches the periodic directions of PRISM and then compute its induced map on \(H_1\) (and, if appropriate, on the full homology ring). That requires an explicit periodic wrapping rule; it cannot be inferred from this local fold.

## 6. Reproduction

    PYTHONPATH=. python -m pytest tests/test_q6_prism_chain_map.py -q
    PYTHONPATH=. python examples/q6_prism_chain_map_audit.py

## 7. Evidence ledger

| Claim | Status |
|---|---|
| \(F_0\) is a surjective 2-to-1 graph homomorphism \(Q_6\to Q_5\) | PROVED |
| The signed cell formula defines an integral chain map | PROVED by formula; exhaustive exact identity check |
| Chain identity on all 729 cells | NUMERICALLY VERIFIED with exact integer arithmetic |
| Surjective on all 243 target cell generators | NUMERICALLY VERIFIED |
| \(H_0\) map is an isomorphism | PROVED |
| Positive-degree homology of the cube source/target | ZERO, since both are contractible |
| Map captures full PRISM torus homology or transport | NOT ESTABLISHED |

