# Kathara-Q6 Gate D: non-degenerate integral face lift

**Status:** exact integral chain-map construction for the existing explicit Kathara vertex chart, using a declared rank-effective Q6 path-order convention. The RGS coordinate attachment is checked. The path-order convention is not RGS-canonical; RGS-equivariant canonicality and compatible quaternionic transport are not claimed.

## 1. Why the earlier map could not work

The Gate C edge map used root paths P(g(v))-P(g(u)). It is a vertex-potential difference, so every closed source edge cycle telescopes: R1 B2^K = 0. Consequently R2=0 was forced if the target filling was not otherwise specified. That map is valid but carries no source face geometry.

Gate D changes the edge map: each source edge is sent to a shortest Q6 path between its already-declared endpoint labels. The six target axes are traversed in the fixed order

`AXIS_ORDER = (1, 2, 0, 3, 4, 5)`.

This is an explicit combinatorial path convention. It is not the default numerical axis order, is not selected by the RGS metric, and must not be described as RGS-canonical. Its role is to give a deterministic non-degenerate chain map for the declared vertex chart.

## 2. Orientation of source-edge paths

The columns of R1 are indexed by the canonically oriented source edges (u,v), u<v. If a face boundary traverses such an edge backwards, its path word must be the exact reverse of the canonical column path:

[
W(v	o u)=operatorname{reverse}(W(u	o v)).
]

Recomputing a shortest path from the reversed endpoints in the fixed axis order is not equivalent. The audit uses this oriented-edge rule before concatenating each face word.

## 3. Exact cubical filling algorithm

For a source face boundary, concatenate the oriented axis words of its edge paths. Since the boundary is closed, every bit axis occurs an even number of times.

In Q6, two adjacent distinct bit flips commute across a unique square. Replace each adjacent inversion relative to AXIS_ORDER by the opposite two-edge path, recording the oriented square whose boundary is the difference between the two paths. Bubble sorting the word puts equal axes together; equal adjacent flips cancel as a backtrack. The accumulated signed squares form an integral 2-chain whose boundary is exactly the original loop.

This constructs each column of R2 over Z, without floating-point least squares or arbitrary cell weights.

## 4. Resulting chain map

With the existing explicit coordinate chart g: V_K -> V_Q6, define R0 as the one-hot vertex map, R1 as shortest axis-ordered Q6 edge paths between mapped endpoints, and R2 as the exact cubical square filling produced by the swap algorithm.

The audit checks using integer arrays:

[
B_1^{Q6}R_1=R_0B_1^K,qquad
B_2^{Q6}R_2=R_1B_2^K.
]

The face matrix has full column rank 10, certified by row reduction modulo 1,000,003. Full rank modulo this prime proves rank at least 10 over Q; with 10 columns, its rational rank is exactly 10.

The axis order is a declared combinatorial choice that makes this particular face lift non-degenerate. It is not a variational result and does not establish uniqueness or geometric naturality.

The code attaches the 12 mapped labels to the archived X_Oh coordinates through the existing hash-verifying seed loader and checks that all 21 source edges have positive chord length.

## 5. Scope and remaining geometry gate

This is a non-degenerate integral chain map, not yet a cellular embedding or an RGS-canonical map. The shortest-path and axis-order rules are declared deterministic combinatorial conventions. They have not been derived by minimizing or intertwining the RGS metric. The integer coordinate chart from Gate C is reused unchanged.

The recovered Ellis spin-connection code acts on points of the exact embedded hypersurface
F(u,n) = (sqrt(u^2+r0^2) n, r0 asinh(u/r0)) in R4.
The archived X_Oh in R3 alone does not specify (u,n), a physically selected r0, or a unique lift of each Q6 edge path to a curve on that hypersurface. Therefore this Gate D does not yet compute SU(2) links or plaquette spin holonomies. The next gate must declare and justify that geometric lift before calling ellis_edge_spin_link; it must not silently identify chordal RGS edges with Ellis geodesics.

Likewise, the 188D tangent basis Q in R192 is not itself a 3D tangent frame at each of the 64 vertices. A local-frame extraction and its equivariance must be specified before contracting SU(2) links with the three Walsh channels.

## 6. Reproduction

Run from the repository root:
```bash
PYTHONPATH=. python -m pytest tests/test_kathara_q6_nondegenerate_chain_map.py -q
```

**Evidence ledger**
- Integral chain identities: exact integer matrix equalities.
- rank_Q R2=10: certified by full rank modulo 1,000,003.
- Nonzero face fills: checked for all ten source faces.
- RGS coordinate attachment: checked using the hash-verifying loader; not a metric-canonical selection.
- Cellular embedding / RGS equivariance: open.
- Ellis/SU(2) transport on this map: open pending a justified (u,n,r0) lift and edge-path rule.
- Supergolden three-channel intertwiner: open; no target matrix is imposed.
