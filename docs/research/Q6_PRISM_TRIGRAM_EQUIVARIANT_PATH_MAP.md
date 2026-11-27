# Trigram-block-equivariant path-valued Q6–PRISM map

**Status: PROVED constructive graph-homology result.** This replaces the arbitrary five-chord assignment by an explicit map derived from the six-bit lower/upper trigram split. It is a map of the Q6 graph, not of the filled cubical 6-cube.

## 1. I-Ching bit convention and symmetry

Write a Q6 vertex as (v=(ell,u)in{0,1}^3\times\{0,1}^3), where bits (0,1,2) are the lower trigram and bits (3,4,5) are the upper trigram. The source symmetry used here is

[
H=S_3^{\mathrm{lower}}\times S_3^{\mathrm{upper}},
]

the independent permutations of the three line positions within each trigram. The target action permutes the three equal (C_6) factors according to (S_3^{\mathrm{lower}}) and fixes the (C_{64}) and (C_5) factors; (S_3^{\mathrm{upper}}) acts trivially on the target.

This is a nontrivial I-Ching block symmetry, but it is not the full hyperoctahedral group: it does not claim invariance under bit flips, exchange of lower and upper trigrams, or all coordinate permutations (S_6).

## 2. Edge-path assignment

Let (k=|u|) be the Hamming weight of the upper trigram. All 64 source vertices map to the PRISM basepoint (p=(0,0,0,0,0)). Orient each source edge positively from bit (d=0) to bit (d=1).

For a source edge in lower-trigram direction (d\in\{0,1,2\}), assign a based target loop with winding vector, in factor order ((C_{64},C_6^{(0)},C_6^{(1)},C_6^{(2)},C_5)),

[
w_d(k)=
\bigl(\mathbf1_{k\ge2},\;k\,\delta_{0d},\;k\,\delta_{1d},\;k\,\delta_{2d},\;\mathbf1_{k=3}\bigr).
]

The loop traverses each indicated factor positively the prescribed integer number of times, in the fixed order (C_6^{(d)}), (C_{64}), (C_5). For edges in upper-trigram directions (3,4,5), use the constant path. Reverse-oriented edges use the reversed path. Each edge path starts and ends at (p), so the assignment defines a continuous map from the entire source graph.

The assignment depends only on the upper-trigram Hamming weight and the distinguished lower direction. Therefore it is equivariant under (S_3^{\mathrm{lower}}\times S_3^{\mathrm{upper}}) with the stated target action.

## 3. Five square cycles certify integral surjectivity

For a mixed square in lower direction (d) and upper direction (j), with the other two upper bits having weight (k\in\{0,1,2\}), orient the boundary as lower edge, upper edge, reversed lower edge, reversed upper edge. Its target winding is

[
q_{d,k}=
\bigl(\mathbf1_{k\ge2}-\mathbf1_{k+1\ge2},\;
-\delta_{0d},-\delta_{1d},-\delta_{2d},\;
\mathbf1_{k=3}-\mathbf1_{k+1=3}\bigr).
]

In particular, the five source squares with ((d,k)=(0,0),(1,0),(2,0),(0,1),(0,2)) map respectively to

[
-e_{C_6^{(0)}},\
-e_{C_6^{(1)}},\
-e_{C_6^{(2)}},\
-e_{C_{64}}-e_{C_6^{(0)}},\
-e_{C_5}-e_{C_6^{(0)}}.
]

These five vectors form an integral basis of (mathbb Z^5): the corresponding (5\times5) minor has determinant (\pm1). Since each is the image of an actual Q6 graph cycle, the induced homomorphism

[
H_1(Q_6^{(1)};\mathbb Z)\cong\mathbb Z^{129}
\longrightarrow H_1(P;\mathbb Z)\cong\mathbb Z^5
]

is surjective over the integers, not merely of rank five over the rationals.

## 4. What is and is not preserved

**Preserved:** the lower/upper trigram block decomposition; independent permutations of line positions within each block; the three equal (C_6) target factors transform as the lower-trigram direction triple; the induced map is integrally surjective.

**Not claimed:** preservation of Q6 graph distances, local adjacency, the individual 3-bit trigram values, bit-flip symmetries, lower/upper block exchange, the full (S_6) coordinate symmetry, or any dynamics/transport law.

The construction is still path-valued and collapses all vertices to the common target basepoint. It is more structured than the five-arbitrary-chords construction, but it is not yet a geometric embedding or a local model.

## 5. Relation to the filled Q6 cube

The source is only the Q6 1-skeleton. The graph cycle lattice has rank (192-64+1=129). If the map extended over all 240 filled square faces, its induced map would factor through the first homology of the filled cube, which is zero. Therefore this surjective graph-homology map cannot extend over the entire filled cubical 6-cube. The nontrivial images of the mixed square cycles are precisely why those faces cannot all be filled in the target.

## 6. Reproduction

```bash
PYTHONPATH=. python -m pytest tests/test_q6_prism_trigram_equivariant_map.py -q
PYTHONPATH=. python examples/q6_prism_trigram_equivariant_map_audit.py
```

The implementation also checks path-level equivariance for all 36 elements of (S_3\times S_3) and all 192 positively oriented source edges.

## 7. Evidence ledger

| Claim | Status |
|---|---|
| Edge-path assignment defines a continuous map on the Q6 graph | PROVED by construction |
| Five square cycles map to a unimodular basis of target (H_1) | PROVED |
| Induced map (\mathbb Z^{129}\to\mathbb Z^5) is surjective | PROVED |
| Equivariance under (S_3^{lower}\times S_3^{upper}) | PROVED by formula; exhaustively checked in the implementation |
| Full (S_6), bit-flip, or block-exchange equivariance | NOT CLAIMED |
| Extension over all 240 filled faces | IMPOSSIBLE |
| Metric, dynamical, or physical interpretation | OPEN |
