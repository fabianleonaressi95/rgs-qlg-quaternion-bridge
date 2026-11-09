# Kathara–Q6 integral chain-map gate (C)

**Status:** explicit, exact-integer chain-map candidate constructed; canonicality and geometric compatibility remain open. No test execution is claimed by this commit.

## 1. Source and target complexes

The source is the repository's declared Kathara cellular complex:
- 12 centres;
- 21 elementary edges, obtained by expanding the 15 primary lines;
- 10 bounded planar faces, treated as 2-cells by the existing convention.

The target is the canonical cubical 2-skeleton from `src/q6_complex.py):
- 64 vertices;
- 192 edges;
- 240 square faces;
- integer labels 0..63, six-bit big-endian;
- canonical edge and face orientations.

The 15 primary Kathara lines are retained as explicit edge-membership metadata. They are not converted into arbitrary numerical weights.

## 2. Explicit vertex code

For a Kathara centre with drawing coordinates `(x,y)`, use the exact column code
`-2 -> 0, 0 -> 1, +2 -> 2`, and retain the integer level `y in {0,...,8}`.
The Q6 label is
[
g(x,y)=16c(x)+y.
]
This uses the source's integer coordinate chart directly. It uses no threshold, nearest-neighbour search, floating-point comparison, or tunable parameter. It is injective on the 12 listed centres.

**Important limitation:** the column/level coding is a declared coordinate chart, not a consequence of the RGS metric or of an intrinsic isomorphism. It is therefore a reproducible candidate, not a canonical map selected by RGS.

## 3. Chain matrices

Let (P(g)in C_1(Q_6;mathbb Z)) be the canonical path chain from 0 to (g), flipping the six bits in axis order 0,...,5. For each oriented Kathara edge (e=(u,v)), define
[
R_0[v_K]=[g(v_K)],qquad
R_1[e]=P(g(v))-P(g(u)),qquad R_2=0.
]
Here the target canonical edge orientation is the one in `q6_complex.py`; path coefficients are signed accordingly.

The matrices have shapes
[
R_0inmathbb Z^{64	imes12},quad
R_1inmathbb Z^{192	imes21},quad
R_2inmathbb Z^{240	imes10}.
]
Since (partial_1P(g)=[g]-[0]), one obtains exactly
[
B_1^{Q6}R_1=R_0B_1^K.
]
For every source face, the edge-chain image telescopes:
[
R_1B_2^K=0=B_2^{Q6}R_2.
]
Thus this is an integral chain map. The identities are algebraic consequences of the construction; the accompanying tests are provided for repository execution.

## 4. What this map does and does not achieve

- It is injective on the 12 source vertices.
- It maps each source edge to a root-routed 1-chain in the Q6 1-skeleton.
- It collapses all ten source 2-cells: (R_2=0).
- It is not a cellular embedding and does not preserve the nontrivial face geometry.
- It uses the actual RGS seed by pulling back (X_{g(v)}) at the 12 assigned Q6 labels; the loader verifies the archived coordinate/coframe hashes.
- It does **not** prove that those 12 labels are geometrically preferred by (X_g). The RGS coordinates are attached after the exact combinatorial code; no metric minimization or symmetry selection is being smuggled in.

The construction proves that an incidence-valid integer chain map exists under an explicit coordinate convention. It does not solve the stronger problem of a nondegenerate, RGS-equivariant Kathara-to-Q6 cellular embedding.

## 5. Next mathematical obstruction

A useful Gate C successor must require at least one of:
1. a nonzero (R_2) preserving the ten face classes;
2. a cellular embedding of the source disk into the Q6 2-skeleton;
3. an RGS-equivariant assignment of the 12 source centres, with a source action specified on the 15 marked lines;
4. a metric/intertwining condition derived independently of this coordinate code.

The source complex is a disk. A map that sends its entire boundary-face structure to zero is valid as a chain map but too coarse to support a meaningful area/Hodge/Krist-All comparison. Gate D must therefore not use this degenerate map as evidence for canonical selection.

## 6. Reproduction

```bash
PYTHONPATH=. python -m pytest tests/test_kathara_q6_chain_map.py -q
```

The audit entry point is `audit_kathara_to_q6_chain_map()`. It returns the three integer matrices, exact chain-identity flags, the 15-line incidence metadata, and the RGS coordinates at the image labels.

**Evidence classification:** chain-map formula = proved algebraically from the definitions; implementation = tests added, not yet executed here; RGS canonicality = open; nondegenerate cellular lift = open.
