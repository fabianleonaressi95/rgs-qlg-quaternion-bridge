# Trigram–Krist-All–supergolden crosswalk

**Purpose:** reconcile the matrices and recurrence claims already present in this repository without treating a shared characteristic polynomial as a geometric derivation.

## 1. Two explicit three-state matrices

The trigram refinement target in `docs/TRIGRAM_Q3_WALSH_REFINEMENT_THEOREMS.md) is
[
M_T=\begin{pmatrix}0&0&1\\1&0&0\\0&1&1\end{pmatrix}.
]

The earlier Cantor–Q6 transfer model in `docs/PHASE21_CANTOR_Q6_SUPERGOLDEN_TRANSFER.md) uses
[
A_\psi=\begin{pmatrix}1&0&1\\1&0&0\\0&1&0\end{pmatrix}.
]
These matrices are permutation-similar. With
[
P=\begin{pmatrix}0&1&0\\0&0&1\\1&0&0\end{pmatrix},
qquad P^{-1}=P^T,
]
direct multiplication gives
[
oxed{M_T=P A_\psi P^{-1}.}
]
Consequently both have characteristic polynomial
[
\lambda^3-\lambda^2-1
]
and the same Perron eigenvalue (psi=1.465571231876768\ldots).

This is an exact change of ordering of the three abstract states. It is **not** an intertwiner derived from geometric cells, and it does not identify the Phase 21 symbolic states with the lower/central/upper trigram lines.

## 2. The separate Krist-All polynomial record

The master catalogue records the reported relation
[
K^3=K^2+I
]
and 11,725 Krist-All classes, but explicitly marks the claim as a reported finite algebraic identity / numerical verification pending a complete definition of (K), its coefficient field, a reproducible construction and an independent exact-arithmetic check.

If a concrete finite matrix (K) satisfies this identity, its minimal polynomial divides (x^3-x^2-1). This fact alone does not imply that (K) is similar to either 3x3 matrix above: dimension, minimal polynomial, rational canonical form and multiplicities of invariant factors must be supplied. Nor does it construct a parent-to-child map for the Q6 or Kathara cells.

## 3. What the new Q6 subdivision contributes

The midpoint benchmark in [Q6_DYADIC_SUBDIVISION_CHAIN_MAP.md](Q6_DYADIC_SUBDIVISION_CHAIN_MAP.md) constructs explicit oriented maps
[
S_0:C_0(Q_6)\to C'_0,quad
S_1:C_1(Q_6)\to C'_1,quad
S_2:C_2(Q_6)\to C'_2
]
and verifies
[
B'_1S_1=S_0B_1,qquad B'_2S_2=S_1B_2.
]
It is a real geometric subdivision of the standard cubical 2-skeleton, but it preserves coordinate-axis labels. It does not identify Kathara cells, and it does not provide a map from refined chains to the three trigram Walsh modes.

## 4. Required next input for a genuine Krist-All derivation

A geometric derivation still needs source-controlled data specifying:
1. the actual Kathara/Krist-All parent cells, with coordinates and dimensions;
2. the child-cell list for every parent, with orientations;
3. shared-face identifications and exact signed incidence coefficients;
4. the embedding maps of child cells into parent cells;
5. the map from the three ordered trigram lines to geometric directions/channels;
6. an explicit relation between the finite Krist-All operator (K) and the subdivision maps.

Only after these are specified can one assemble the subdivision chain maps, induce a three-channel operator, and test whether its matrix is (M_T) (or the permutation-conjugate (A_\psi)). Until then the geometric identification and the supergolden derivation remain **OPEN**.

## 5. Reproduction

```bash
PYTHONPATH=. python -m pytest tests/test_trigram_supergolden_crosswalk.py -q
PYTHONPATH=. python examples/trigram_supergolden_crosswalk_audit.py
```
