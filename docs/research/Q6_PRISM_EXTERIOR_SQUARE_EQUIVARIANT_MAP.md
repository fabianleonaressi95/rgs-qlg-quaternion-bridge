# An A6-equivariant exterior-square map from Q6 graph cycles

**Status: PROVED algebraic construction.** This construction resolves the representation obstruction by enriching the target with a 15-dimensional homology module carrying the natural coordinate-pair representation. It does not claim a surjection onto the original five-dimensional PRISM homology.

## 1. Source and enriched target

Let (G=Q_6^{(1)}) be the Q6 graph, with graph cycle lattice
[
Z_1(G;mathbb Z)congmathbb Z^{129}.
]
The original PRISM product is
[
P=C_{64}square C_6^3square C_5,qquad H_1(P;mathbb Z)=mathbb Z^5.
]

Enlarge the target by a 15-torus with one circle factor for each unordered pair of the six Q6 coordinates:
[
widehat P=P	imes (S^1)^{15},
qquad
H_1(widehat P;mathbb Z)
=mathbb Z^5oplusLambda^2(mathbb Z^6).
]
The added 15 generators are (e_iwedge e_j), (0le i<jle5). A finite graph model replaces each added circle by (C_5). The group (S_6), hence (A_6), acts on this added module by
[
sigma(e_iwedge e_j)=e_{sigma(i)}wedge e_{sigma(j)},
]
including the sign when the resulting pair is reordered.

## 2. A canonical cycle homomorphism

Orient each Q6 edge in coordinate direction (d) from bit (d=0) to bit (d=1). For a positive edge based at (x), define the (Lambda^2(mathbb Z^6))-valued 1-cochain
[
alpha_d(x)=sum_{i<d}x_i,e_iwedge e_d.
]
For a graph cycle (z=sum_e n_e e), set
[
Phi(z)=sum_e n_ealpha(e).
]
This is an integral homomorphism from the graph cycle lattice to the added target lattice.

For the oriented square in directions (i<j), based at a vertex with (x_i=x_j=0), its discrete curl is
[
alpha_j(x+e_i)-alpha_j(x)
-alpha_i(x+e_j)+alpha_i(x)
=e_iwedge e_j.
]
Thus every one of the 15 basis vectors is attained by a square boundary. There are 16 squares for each coordinate pair, and all 240 square faces give the prescribed coordinate-pair wedge. Hence (Phi) is surjective.

Since the source cycle rank is 129 and the image rank is 15,
[
operatorname{rank}kerPhi=129-15=114.
]

## 3. Equivariance

The square-boundary characterization is independent of the chosen coordinate ordering: the image of an oriented (ij)-square is the corresponding oriented wedge (e_iwedge e_j). Under any coordinate permutation (sigmain S_6), the square is sent to the square in directions (sigma(i),sigma(j)), and its image transforms by the exterior-square action. Therefore
[
Phi(sigma z)=Lambda^2(sigma)Phi(z).
]
In particular, (Phi) is (A_6)-equivariant.

The displayed edge cochain is a convenient gauge for realizing this homomorphism by paths. Different coordinate orderings change the cochain by a coboundary and therefore do not change its action on graph cycles.

## 4. Relation to the original PRISM

The induced map into the enriched target is
[
oxed{
Z_1(Q_6^{(1)};mathbb Z)
longrightarrow
H_1(widehat P;mathbb Z)
=mathbb Z^5oplusLambda^2(mathbb Z^6),
quad
zmapsto(0,Phi(z)).
}
]
It is surjective onto the added 15-dimensional summand, not onto the original PRISM (mathbb Z^5) summand. This distinction is essential: the preceding symmetry obstruction still forbids a nonzero (A_6)-equivariant projection onto the original trivial (H_1(P)) module.

## 5. Filled plaquettes and scope

The source is the Q6 graph, so its 240 square boundaries are cycles that generate the graph cycle lattice. If all 240 plaquettes are filled, these cycles become boundaries and (H_1) of the filled cube is zero. The map (Phi) therefore does not extend as a homomorphism from the first homology of the filled cubical Q6 to a nonzero target (H_1).

This construction is a canonical, symmetry-compatible graph-cycle invariant. It is not a metric embedding, not a map preserving all six-dimensional local distances, and not a physical transport law. The added 15-torus is a mathematically defined enrichment whose physical interpretation remains open.

## 6. Verification and evidence ledger

The audit checks all 240 square curvatures and all (720\times15=10800) coordinate-permutation actions on the wedge basis.

| Claim | Status |
|---|---|
| (Phi:Z_1(Q_6^{(1)};mathbb Z)	oLambda^2mathbb Z^6) is integral | PROVED |
| (Phi) is surjective | PROVED, each wedge basis vector is a square image |
| (kerPhi) has rank 114 | PROVED |
| (Phi) is (S_6)-equivariant, hence (A_6)-equivariant | PROVED |
| Projection onto original PRISM (mathbb Z^5) summand is nonzero equivariantly | IMPOSSIBLE under the preceding obstruction |
| Physical interpretation of the added 15 cycles | OPEN |
