# Kathara–Q6 gate C: line incidence, equivariance and exact obstructions

## Purpose and status

This note extends the Kathara–Q3 / RGS canonical-selection audit by retaining the source's 15 primary lines, rather than reducing Kathara immediately to its 21 elementary edges. It specifies the exact line-to-edge incidence data and the symmetry constraints that a future integer chain map must satisfy.

**Status:** the incidence and symmetry statements below follow from the repository transcription in `src/kathara_15_line_graph.py`. The Kathara-to-Q6 matrices are **not** constructed here. No canonical map is claimed and no tests are claimed to have been executed.

## 1. Fixed source conventions

The source graph has 12 centres and 15 labelled primary lines. Expanding the polylines and sorting undirected edges lexicographically gives the following fixed edge order:

| index | edge | index | edge | index | edge |
|---:|:---:|---:|:---:|---:|:---:|
| 1 | 1–2 | 8 | 4–7 | 15 | 8–10 |
| 2 | 1–3 | 9 | 5–6 | 16 | 8–11 |
| 3 | 1–4 | 10 | 5–7 | 17 | 9–11 |
| 4 | 2–3 | 11 | 5–8 | 18 | 9–12 |
| 5 | 2–4 | 12 | 6–9 | 19 | 10–11 |
| 6 | 2–5 | 13 | 7–10 | 20 | 10–12 |
| 7 | 3–6 | 14 | 8–9 | 21 | 11–12 |

The 15 primary lines expand as follows (indices refer to the table above):

- L1, central vertical: (1, 6, 11, 16, 21), the polyline 1–2–5–8–11–12.
- L2, left vertical: (7, 12), the polyline 3–6–9.
- L3, right vertical: (8, 13), the polyline 4–7–10.
- L4: (2), edge 1–3.
- L5: (4), edge 2–3.
- L6: (3), edge 1–4.
- L7: (5), edge 2–4.
- L8: (18), edge 9–12.
- L9: (17), edge 9–11.
- L10: (20), edge 10–12.
- L11: (19), edge 10–11.
- L12: (14), edge 8–9.
- L13: (15), edge 8–10.
- L14: (9), edge 5–6.
- L15: (10), edge 5–7.

Every elementary edge belongs to exactly one primary line in this transcription. Thus the line-edge incidence matrix A is a 15-by-21 binary matrix with exactly one 1 in each column. Its rows are specified completely by the index sets above; this is preferable to silently treating the 15 source lines as 15 additional edges.

## 2. Vertical hierarchy marker

Define a four-class marker V in {0,1}^{4 x 15}, with rows ordered
(central vertical, left vertical, right vertical, cross-link), and columns L1,...,L15:

- central row: 1 only in column L1;
- left row: 1 only in column L2;
- right row: 1 only in column L3;
- cross-link row: 1 in columns L4,...,L15.

Every column of V sums to one. This records the documented distinction among the three vertical lines without assigning numerical energy weights. In particular, it does **not** justify setting a_l to any particular positive values in a weighted functional.

The left/right distinction is a marking, not an asymmetry of the unmarked planar graph.

## 3. Kathara reflection and its action on line data

The reflection
`sigma = (3 4)(6 7)(9 10)`
fixes centres 1,2,5,8,11,12 and exchanges the left and right branches.

Its induced action on primary-line labels is

`L1 fixed;`
`(L2 L3)(L4 L6)(L5 L7)(L8 L10)(L9 L11)(L12 L13)(L14 L15).`

The line L1 remains distinguished. The two lateral vertical lines are exchanged. The cross-link class is preserved as a set. This gives an exact action on A and V: if P_E and P_L are the induced permutation matrices on elementary edges and line labels, then
`P_L A = A P_E`
under the convention that permutation matrices act on column-indexed objects by relabelling. Implementations must document their left/right action convention and test the identity with that convention.

## 4. Chain-map gate into Q6

Let K be the filled Kathara planar complex already defined in the repository, with 12 vertices, 21 edges and 10 bounded 2-cells. Let Q6^(2) be the cubical 2-skeleton with 64 vertices, 192 edges and 240 square faces.

An integer chain map consists of

- R0 in Z^(64 x 12),
- R1 in Z^(192 x 21),
- R2 in Z^(240 x 10),

satisfying exactly

`B1_Q6 R1 = R0 B1_K,`
`B2_Q6 R2 = R1 B2_K.`

These equations are necessary and sufficient for a degree-0,1,2 chain map between the displayed truncated chain complexes. They do not by themselves imply that the map is induced by a continuous cellular map, injective, or an embedding.

For a vertex-induced cellular map, R0 must have exactly one 1 in each column. R1 must then represent oriented target edge paths joining the images of the endpoints. R2 must fill the resulting image cycles by integer 2-chains. A general integer chain map is a broader algebraic object and should not be described as a geometric embedding without extra conditions.

### Equivariance constraint

Let P0_K, P1_K, P2_K be the source chain actions induced by sigma, and let U0_Q, U1_Q, U2_Q be a *specified* cellular involution of Q6^(2). Equivariance requires

`U0_Q R0 = R0 P0_K,`
`U1_Q R1 = R1 P1_K,`
`U2_Q R2 = R2 P2_K.`

This is the appropriate way to use the Kathara left/right hierarchy. The target involution cannot be chosen implicitly: it must be supplied by an independently established Q6 cellular action. The coordinate O_h action on the 64-point RGS seed and its action on Q6 cells are separate audits; the former does not automatically provide the latter.

## 5. Exact obstruction to a one-edge-per-edge embedding

The Kathara 1-skeleton contains triangles, for example 1–2–3–1. The Q6 1-skeleton is bipartite: each Q6 edge flips exactly one bit, hence changes the parity of the Hamming weight.

Consequently, there is no injective graph embedding of the un-subdivided Kathara 1-skeleton into the Q6 1-skeleton that sends every source edge to one target edge. Such an embedding would send a 3-cycle to a 3-cycle, which cannot exist in a bipartite graph.

This is a genuine, exact obstruction, but it is deliberately narrow. It does **not** rule out:
- a map that sends source edges to multi-edge paths;
- a map after subdividing the Kathara source;
- a non-injective cellular map;
- a general integer chain map satisfying the two boundary identities.

Any proposed gate-C construction must state which of these categories it belongs to.

## 6. Why the 15-line data do not select a Q6 bit assignment

The line-incidence matrix A and marker V are intrinsic to the transcribed Kathara combinatorics. They do not identify:
- one of the six Q6 bit coordinates with a Kathara direction;
- a Q6 vertex for any Kathara centre;
- a target cellular involution U induced by sigma;
- a unique path for each source edge;
- a unique 2-chain filling each mapped face.

The RGS seed supplies 64 embedded points indexed by Q6 labels, but the source has only 12 centres and no source-to-seed incidence correspondence. A distance minimization or nearest-seed assignment would add a metric matching rule; it is a candidate convention, not a consequence of A or V.

Likewise, the qualitative two-grid / Kristalla–Kristos language in the source does not specify a map of vertices, edges or 2-cells. It may motivate a later doubled construction, but cannot be inserted as an exact chain map without further source data.

## 7. Gate-C deliverables and acceptance criteria

A future implementation should emit, in machine-readable form:

1. the fixed vertex, edge, face and primary-line orderings;
2. A and V, together with their exact permutation identities under sigma;
3. the proposed Q6 involution U and proof/audit that it preserves vertices, edges and square faces;
4. R0, R1, R2 as integer matrices;
5. exact zero residuals for both chain identities;
6. the map category: vertex-induced cellular map, subdivision map, quotient, or abstract chain map;
7. injectivity/rank and image-support diagnostics, without conflating rank with topological embedding;
8. the provenance of every additional geometric or metric choice.

Only after these gates pass should one compare RGS tangent or D5 spectral functionals on candidate maps. The functional and its inner products must be fixed independently of the desired A/B outcome.

## Evidence ledger

| Claim | Status |
|---|---|
| 15 primary lines expand to 21 unique elementary edges | Source transcription / exact combinatorial expansion |
| Every elementary edge belongs to one primary line | Exact from the listed line expansion |
| Reflection sigma acts on line labels as listed | Exact from the source vertex permutation |
| Q6 1-skeleton is bipartite | PROVED from parity of six-bit Hamming weight |
| No injective one-edge-per-edge graph embedding K^(1) -> Q6^(1) | PROVED by the triangle/bipartite obstruction |
| A nontrivial path-valued or abstract chain map exists with desired extra properties | OPEN |
| A canonical Kathara-to-Q6 map follows from the source hierarchy | NOT ESTABLISHED |
| An RGS/D5 variational selection is defined on such a map | OPEN until the lift and functional are fixed |

## Repository entry points

- `src/kathara_15_line_graph.py`
- `src/kathara_q3_chain_map.py`
- `src/rgs_q6_oh_action.py`
- `docs/KATHARA_RGS_CANONICAL_SELECTION_AUDIT.md`

No tests or numerical computations are claimed by this note.
