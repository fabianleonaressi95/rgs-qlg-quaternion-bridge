# Q3–Q6–Omega/Beneš–PRISM: exact finite bridge and obstruction

**Status:** finite combinatorial correspondence implemented; graph/dynamical equivalence to PRISM remains OPEN.

## 1. Fixed I Ching convention

A trigram is a three-bit word in **bottom-to-top line order**. The eight trigrams are the vertices of \(Q_3\cong\mathbb F_2^3\). For example, \(100\) and \(011\) are distinct trigrams and have Hamming distance three. A hexagram is an ordered pair
\[
h=(\ell,u),\qquad \ell,u\in\{0,1\}^3.
\]
The repository encoding is
\[
\operatorname{idx}(h)=\operatorname{int}_2(\ell)+8\operatorname{int}_2(u),
\]
so the lower trigram occupies Q6 coordinates 0–2 and the upper trigram coordinates 3–5. No King Wen ordering is imposed.

## 2. Q3 coordinate matchings are the switch stages

For coordinate \(d\in\{0,1,2\}\), define
\[
M_d=\{\{v,v\oplus 2^d\}:v_d=0\}.
\]
Each \(M_d\) is a perfect matching of four disjoint edges on the eight trigram vertices. The three matchings partition all twelve edges:
\[
E(Q_3)=M_0\sqcup M_1\sqcup M_2.
\]
A 2×2 switch stage acting on one trigram bit is exactly a choice, independently on each of the four edges in \(M_d\), to preserve or exchange the two wires.

Under the shuffle-normalized convention in the existing exhaustive network audit:
- Omega-8 uses dimensions \((0,1,2)\), one pass through the three Q3 coordinate matchings;
- Beneš-8 uses \((0,1,2,1,0)\), repeating the first and last coordinate matchings in reverse order.

This is the declared normalized stage model. Physical wire labels in a conventional shuffle drawing still require the interstage shuffle/relabel conversion specified in the network audit.

## 3. Q6 vertices as input/output demands

A Q6 vertex is interpreted as an ordered demand
\[
(\ell,u)=(\text{input trigram},\text{output trigram}).
\]
It is one cell of an 8×8 input/output matrix. A permutation \(\pi\in S_8\) is represented by the eight Q6 vertices
\[
D_\pi=\{(i,\pi(i)):i=0,\ldots,7\}.
\]
This is precisely a permutation matrix: every input and every output occurs once. Conversely, any eight Q6 demand vertices with distinct first and second coordinates decode uniquely to a permutation.

Consequences:
- \(8!=40,320\) possible permutation matrices correspond to 40,320 distinct eight-vertex subsets of Q6;
- the Omega enumerator's 4,096 routable permutations correspond to 4,096 such Q6 matchings;
- the Beneš enumerator's 40,320 routable permutations correspond to all permutation matrices.

These counts inherit the switch-stage and shuffle-normalized conventions in the network audit. They do **not** say that all Q6 vertices form a route, or that Q6 itself is the switching fabric.

## 4. Stage-by-stage path lift

For a permutation witness with stage masks \(s_1,\ldots,s_m\) and dimensions \(d_1,\ldots,d_m\), each input trigram \(i\) has a trace
\[
i=x_0\mapsto x_1\mapsto\cdots\mapsto x_m=\pi(i),
\qquad x_j=T_{d_j,s_j}(x_{j-1}),
\]
where \(T_{d,s}\) independently swaps or preserves the four pairs in \(M_d\). Each transition is an edge of Q3. The eight traces, taken together, realize the permutation witness. This is a path-level description of the finite switch model; it is not a path in the six-cube Q6.

## 5. What this says—and does not say—about PRISM

The PRISM product previously recorded in the project is
\[
C_{64}\square C_6^3\square C_5,
\]
with 69,120 vertices. Its \(C_{64}\) factor has 64 positions, so the 64 Q6 demand labels can be put in bijection with those positions. This is only a **label correspondence**.

It is not a graph isomorphism between Q6 and \(C_{64}\): every Q6 vertex has degree 6, whereas every vertex of cycle \(C_{64}\) has degree 2. No edge-preserving embedding, covering, quotient, chain map, or dynamics-preserving map to the PRISM product is asserted here. Such a map needs its own definition and tests against the full Cartesian-product adjacency.

The current bridge therefore ends at the permutation-matching layer. It does not manufacture an arbitrary Q6-to-PRISM edge map to make the diagrams commute.

## 6. Reproduction

    PYTHONPATH=. python -m pytest tests/test_q3_q6_network_bridge.py tests/test_interconnection_networks_8x8.py -q
    PYTHONPATH=. python examples/q3_q6_network_prism_audit.py

The audit is written to results/q3_q6_network_prism_audit.json.

## 7. Evidence ledger

| Claim | Status |
|---|---|
| Q3 has 8 vertices, 12 edges, and three coordinate perfect matchings | PROVED finite combinatorics |
| Each switch stage is a choice on one Q3 coordinate matching | PROVED under the declared stage model |
| Q6 vertices encode ordered input/output trigram demands | DEFINITION |
| Eight demands form a permutation iff they are a perfect matching of the 8×8 demand bipartition | PROVED finite combinatorics |
| Omega/Beneš routable counts and Q6 matching encodings | NUMERICALLY VERIFIED by exhaustive enumeration |
| Each witness lifts to eight stage-by-stage Q3 paths | NUMERICALLY VERIFIED; transition identity follows from stage definition |
| Q6 is graph-isomorphic to the PRISM C64 factor | DISPROVED by degree mismatch |
| A canonical Q6-to-full-PRISM graph map or dynamical equivalence exists | OPEN |
| Any physical Green–Kubo, transport, or continuum conclusion follows | NOT ESTABLISHED |

No Q6–Omega/Beneš–PRISM physical equivalence is claimed.
