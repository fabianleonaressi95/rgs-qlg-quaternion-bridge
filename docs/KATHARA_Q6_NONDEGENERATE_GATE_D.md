# Kathara-Q6 Gate D: non-degenerate integral face lift

**Status:** exact integral chain-map construction from the repository's existing explicit Kathara vertex chart; RGS coordinate attachment is checked. RGS-equivariant canonicality and compatible quaternionic transport are not claimed.

## 1. Why the earlier map could not work

The Gate C edge map used root paths P(g(v))-P(g(u)). It is a vertex-potential difference, so every closed source edge cycle telescopes: R1 B2^K = 0. Consequently R2=0 was forced if the target filling was not otherwise specified. That map is valid but carries no source face geometry.

Gate D changes the edge map: each oriented source edge is sent to a shortest Q6 path between its already-declared endpoint labels, with differing bit axes traversed in increasing axis order. It is not a difference of root potentials.

## 2. Exact cubical filling algorithm

For a source face boundary, concatenate the axis words of its oriented edge paths. Since the boundary is closed, every bit axis occurs an even number of times.

In Q6, two adjacent distinct bit flips commute across a unique square. Replace each adjacent inversion a>b by b,a, recording the oriented square whose boundary is the difference between the two two-edge paths. Bubble sorting the word puts equal axes together; equal adjacent flips cancel as a backtrack. The accumulated signed squares form an integral 2-chain whose boundary is exactly the original loop.

This constructs each column of R2 over Z, without floating-point least squares or arbitrary cell weights.

## 3. Resulting chain map

With the existing explicit coordinate chart g: V_K -> V_Q6, define R0 as the one-hot vertex map, R1 as shortest axis-ordered Q6 edge paths between mapped endpoints, and R2 as the exact cubical square filling produced by the swap algorithm.

The audit checks using integer arrays:
B1^Q6 R1 = R0 B1^K, and B2^Q6 R2 = R1 B2^K.

The face matrix is certified to have full column rank 10 by row reduction modulo 1,000,003. Full rank modulo this prime proves rank at least 10 over Q; with 10 columns, its rational rank is exactly 10.

The code attaches the 12 mapped labels to the archived X_Oh coordinates through the existing hash-verifying seed loader and checks that all 21 source edges have positive chord length.

## 4. Scope and remaining geometry gate

This is a non-degenerate chain map, not yet a cellular embedding or an RGS-canonical map. The shortest-path rule is a declared deterministic combinatorial convention. It has not been derived by minimizing or intertwining the RGS metric. The integer coordinate chart from Gate C is reused unchanged.

The recovered Ellis spin-connection code acts on points of the exact embedded hypersurface
F(u,n) = (sqrt(u^2+r0^2) n, r0 asinh(u/r0)) in R4.
The archived X_Oh in R3 alone does not specify (u,n), a physically selected r0, or a unique lift of each Q6 edge path to a curve on that hypersurface. Therefore this Gate D does not yet compute SU(2) links or plaquette spin holonomies. The next gate must declare and justify that geometric lift before calling ellis_edge_spin_link; it must not silently identify chordal RGS edges with Ellis geodesics.

Likewise, the 188D tangent basis Q in R192 is not itself a 3D tangent frame at each of the 64 vertices. A local-frame extraction and its equivariance must be specified before contracting SU(2) links with the three Walsh channels.

## 5. Reproduction

Run from the repository root:
PYTHONPATH=. python -m pytest tests/test_kathara_q6_nondegenerate_chain_map.py -q

**Evidence ledger**
- Integral chain identities: proved by the exact construction and checked by integer matrix equality.
- rank_Q R2=10: certified by full rank modulo 1,000,003.
- RGS coordinate attachment: checked using the hash-verifying loader; not a metric-canonical selection.
- Cellular embedding / RGS equivariance: open.
- Ellis/SU(2) transport on this map: open pending a justified (u,n,r0) lift and edge-path rule.
- Supergolden three-channel intertwiner: open; no target matrix is imposed.
