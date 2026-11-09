# Kathara 12-center / 15-primary-line graph: source transcription and face audit

**Status:** finite combinatorial reconstruction of the published Level-1 Kathara diagram. The graph is a transcription of the source's numbered centers, three vertical primary lines and twelve cross-lines. It is not evidence for the source's energetic or cosmological interpretations.

## 1. Source controls

The Emerald Awakening page “Spirals of Creation” states that the Kathara grid has 12 centers, a central vertical line 1–12, AzurA at its midpoint, and a 2:1 ratio between the central column and the 6–7 width. It describes the stage-to-stage rotation/expansion and the center-9 / center-12 overlays (see sections 2 and 3 of that page).

The *Dance for Life, Love, Freedom* manual, “Structure of the Kathara Level-1 12-Tree Grid Core Scalar Template”, explicitly identifies 12 primary centers, 15 primary Kathara lines, three vertical lines (IDs 1, 2, 3), and the central vertical line (ID 1). The diagram labels the 12 remaining lines 4–15. The published diagram's line labels and endpoints are transcribed below. The source text describes the lines as conceptual Partiki-phasing connections; the graph formalizes only their displayed incidence.

## 2. Vertex coordinates

Use the planar integer embedding consistent with the displayed symmetric grid. Coordinates are in arbitrary units; the vertical central column has length 8 and the 6–7 separation has length 2, so the 2:1 proportion holds:
\[
\begin{array}{c|rrrrrrrrrrrr}
v&1&2&3&4&5&6&7&8&9&10&11&12\\\hline
x&0&0&-1&1&0&-1&1&0&-1&1&0&0\\
y&0&1&2&2&3&4&4&5&6&6&7&8
\end{array}
\]
The AzurA point is \((0,4)\), the midpoint of the 1–12 column and the 6–7 segment. These coordinates fix the straight-line embedding used for face enumeration; they are not a claim that the source supplies calibrated metric coordinates for every center.

## 3. Fifteen primary lines

The three vertical primary lines are represented as polylines:
\[
L_1=(1,2,5,8,11,12),\quad
L_2=(3,6,9),\quad
L_3=(4,7,10).
\]
The twelve remaining primary cross-lines are:
\[
\begin{array}{c|rrrrrrrrrrrr}
\ell&4&5&6&7&8&9&10&11&12&13&14&15\\\hline
\text{endpoints}&1\!-\!3&2\!-\!3&1\!-\!4&2\!-\!4&
9\!-\!12&9\!-\!11&10\!-\!12&10\!-\!11&
8\!-\!9&8\!-\!10&5\!-\!6&5\!-\!7
\end{array}
\]
This representation treats the three vertical lines as paths through intermediate centers, while lines 4–15 are individual cross-links. Expanding the three paths into elementary segments and adding the twelve cross-links gives 21 unique undirected edges.

For reproducibility, each elementary edge is oriented from the lower to the higher center label. This is a chosen orientation convention; the source does not assign arrows to these geometric lines. Reversing any edge changes incidence-matrix signs but not the underlying graph or its faces.

## 4. Directed graph and exact finite checks

Let \(V=\{1,\ldots,12\}\) and let \(E\) be the 21 elementary edges obtained above. The audit verifies:
- 12 vertices and 15 named primary lines;
- 21 unique elementary edges after expanding the three vertical paths;
- connectedness;
- no proper crossings in the straight-line embedding;
- an acyclic orientation under increasing center labels;
- a planar rotation system and its face-boundary walks.

The planar embedding has 11 faces including the unbounded exterior face, hence 10 bounded faces. Euler's identity is
\[
\boxed{|V|-|E|+|F|=12-21+11=2.}
\]

## 5. Bounded face boundaries

The ten bounded faces, written as cyclic vertex boundaries, are:
\[
\begin{aligned}
&(1,2,3),\quad (1,4,2),\\
&(2,4,7,5),\quad (2,5,6,3),\\
&(5,7,10,8),\quad (5,8,9,6),\\
&(8,10,11),\quad (8,11,9),\\
&(9,11,12),\quad (10,12,11).
\end{aligned}
\]
The remaining face is the exterior boundary
\[
(1,3,6,9,12,10,7,4).
\]
The list is a topological face decomposition of the chosen planar embedding. It is not automatically a list of physical 2-cells in the Kathara source ontology. A cellular-complex interpretation would require declaring which bounded regions are intended as 2-cells and whether any line crossings or higher-dimensional identifications are present.

## 6. Comparison with our RGS / trigram constructions

- The 12 Kathara centers are **not** the eight trigram labels. There is no canonical bijection between these sets.
- The eight trigram labels form the vertices of \(Q_3\), with 12 Hamming-adjacency edges and six square 2-faces. That is an exact algebraic cube.
- The Kathara graph here has 12 vertices, 21 elementary edges, 10 bounded planar faces and a three-path spine. It is not \(Q_3\), nor the Q6 2-skeleton.
- A possible future bridge is a map from the Kathara graph's vertices/edges/faces to trigram states/channels, but it must preserve the relevant incidence and orientations. No such map is asserted here.
- The square-similarity stage map \(T=\begin{pmatrix}1&1\\-1&1\end{pmatrix}\) concerns the inter-stage planar similarity. It is a separate operator from the adjacency/boundary operators of this 12-vertex graph.

## 7. Reproduction

\`\`\`bash
PYTHONPATH=. python -m pytest tests/test_kathara_15_line_graph.py -q
PYTHONPATH=. python examples/kathara_15_line_graph_audit.py
\`\`\`
