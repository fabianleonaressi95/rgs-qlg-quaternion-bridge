# Symmetry obstruction for Q6 graph cycles mapping onto PRISM periodic H1

**Status: PROVED under the stated equivariance hypothesis.** This is a no-go theorem, not a construction. It explains why the path-valued surjection in the preceding study cannot be made fully coordinate-symmetry-equivariant while the target remains the existing PRISM graph product.

## 1. Setup

Let (G=Q_6^{(1)}) be the 1-skeleton of the six-dimensional cube, regarded as a graph (the 240 square faces are not filled). Then
[
H_1(G;mathbb Z)=Z_1(G;mathbb Z)congmathbb Z^{129}.
]
Let
[
P=C_{64}square C_6square C_6square C_6square C_5.
]
Its cubical first homology is (H_1(P;mathbb Z)congmathbb Z^5).

The coordinate-permutation subgroup (A_6leq S_6leqoperatorname{Aut}(G)) acts on the source. We ask whether a map (f:G	o P), allowing source edges to map to finite target paths, can induce a surjection on (H_1) and be equivariant for this (A_6)-action and an action on (P) by graph automorphisms.

## 2. The target cannot carry a nontrivial A6-action on H1

The five Cartesian cycle factors have lengths (64,6,6,6,5). Graph automorphisms can reverse each factor, and can permute only the three equal (C_6) factors. Consequently the image of (operatorname{Aut}(P)) on (H_1(P;mathbb Z)) is contained in the signed-permutation group
[
(mathbb Z/2)^5times S_3,
]
whose order is (2^5cdot 3!=192).

The group (A_6) is simple of order 360. Any homomorphism from (A_6) to a group of order at most 192 is trivial: its kernel is normal, and an injective image of order 360 cannot fit in a group of order 192. Thus every action of (A_6) on the target induced by graph automorphisms acts trivially on (H_1(P;mathbb Z)).

## 3. The source cycle module has no A6-invariant linear functional

Orient every Q6 edge from its endpoint with bit (i=0) to its endpoint with bit (i=1). An (A_6)-invariant 1-chain has a coefficient (a_k) on every edge for which the other five bits contain exactly (k) ones, (k=0,ldots,5). This follows because (A_6) is transitive on the (k)-subsets of the other five coordinates.

At a vertex of Hamming weight (w), the coefficient of its boundary is
[
w,a_{w-1}-(6-w),a_w=0,
]
where terms outside (0leq kleq5) are omitted. Taking (w=0,1,ldots,6) successively gives
[
a_0=a_1=cdots=a_5=0.
]
Therefore the (A_6)-invariant subspace of the rational graph cycle space (Z_1(G;mathbb Q)) is zero.

Because (A_6) is finite and representations over (mathbb Q) are semisimple, the dual cycle space has no invariant vectors either. Equivalently,
[
operatorname{Hom}_{A_6}(H_1(G;mathbb Q),mathbb Q)=0.
]

## 4. No-go theorem

Suppose (f_*:H_1(G;mathbb Z)	o H_1(P;mathbb Z)) is induced by an (A_6)-equivariant map, with the target action induced by graph automorphisms of (P). The target action on (H_1(P;mathbb Q)) is trivial by Section 2. Hence every coordinate functional of (f_*) is an (A_6)-invariant linear functional on the source cycle space. Section 3 says each such functional is zero. Therefore
[
oxed{f_*=0.}
]
In particular, no such equivariant map can be surjective onto (mathbb Z^5).

## 5. Interpretation and escape routes

This theorem does **not** say that path-valued surjections do not exist: the preceding construction gives one, but it is basis-dependent and breaks the full coordinate symmetry.

It says that the following three requirements cannot all hold simultaneously:

1. source symmetry includes the full even coordinate-permutation group (A_6);
2. target symmetry comes from graph automorphisms of the existing PRISM product;
3. the induced map on periodic (H_1) is nonzero (hence cannot be surjective).

At least one requirement must be relaxed. Mathematically clean options are:
- choose and state a symmetry-breaking datum (a preferred coordinate, orientation, or cycle basis);
- enlarge or change the target so its (H_1)-representation can carry a nontrivial (A_6)-module;
- seek equivariance only under a proper subgroup compatible with the PRISM factor structure.

No physical, dynamical, or transport interpretation follows from this obstruction alone.

## 6. Evidence ledger

| Claim | Status |
|---|---|
| (operatorname{Aut}(P)) acts on (H_1(P)) through a group of order at most 192 | PROVED for the Cartesian product graph |
| Every homomorphism (A_6	ooperatorname{im}(operatorname{Aut}(P)	o GL_5(mathbb Z))) is trivial | PROVED |
| (Z_1(Q_6;mathbb Q)^{A_6}=0) | PROVED by the boundary recurrence |
| Every equivariant induced map (H_1(G)	o H_1(P)) is zero | PROVED |
| Existence of a canonical I-Ching-compatible surjection | OPEN |
| Equivariance under a smaller subgroup | OPEN |
