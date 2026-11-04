# Global torus structure in Q6–D5–Wilson and its RGS–PSU interfaces

**Purpose.** Give one mathematical picture of the different tori appearing in the project. A geometric flat torus, a compact gauge-configuration torus, a harmonic/moduli torus, and the annihilator torus are related but are not interchangeable.

## 1. Four distinct meanings of “torus”

1. **Geometric torus** \(\mathbb T^d=\mathbb R^d/(2\pi\mathbb Z)^d\). It is a compact flat manifold with \(H_1(\mathbb T^d;\mathbb Z)=\mathbb Z^d\). It is the continuum reference space used by the flat finite-element and Fourier-refinement results.
2. **Link-configuration torus** \(\mathbb T^E=U(1)^E\). One angle is assigned to each oriented edge. This is a compact abelian Lie group, not a spacetime manifold.
3. **Flat-connection/moduli torus.** On a connected cell complex, compact U(1) link fields modulo vertex gauge transformations and with vanishing plaquette curvature have harmonic holonomies parametrized by \(H^1(X;U(1))\). For a torsion-free first homology group of rank \(b_1\), this is \(U(1)^{b_1}\).
4. **Annihilator torus.** Given an integer row lattice \(L\subset\mathbb Z^E\), its annihilator is
   \[
   H_L=\{x\in\mathbb T^E:\ell\cdot x=0\pmod{2\pi}\;\forall\ell\in L\}.
   \]
   It is a closed subgroup of \(\mathbb T^E\). Its dimension and disconnected components are controlled by the Smith normal form (SNF) of a basis of \(L\).

The word “toroid” in the PSU–Anti-Twin geometric picture may additionally mean an embedded toroidal surface/handlebody or a two-sheeted configuration. That is a geometric ansatz and must not be identified with any of the compact abelian groups above without an explicit map.

## 2. Q6 integer topology

For the six-dimensional cubical Q6 complex:
\[
V=64,\qquad E=192,\qquad F=240.
\]
The 1-skeleton is connected and 6-regular. Therefore
\[
\operatorname{rank}D_0=V-1=63,\qquad
\dim\ker D_0=E-63=129.
\]
The integer cycle lattice has rank 129; equivalently the first Betti number is \(b_1=129\). The 240 square faces generate the plaquette boundaries, with the chain identity \(D_1D_2=0\). The quotient of link angles by vertex gauge transformations has dimension 129; its harmonic flat-connection torus is \(\mathbb T^{129}\) when the integral homology has no torsion.

This is a topological count. It does not select five physical dimensions, and it does not turn the six-dimensional Q6 complex into 3+1-dimensional spacetime.

## 3. The five-mode projection and its integer fibre

Let \(C\in\mathbb Z^{129\times192}\) be a fundamental-cycle matrix and let \(P_{12}\) be the exact rational projector onto the selected coexact five-mode sector. The map from integer cycle coordinates to the retained sector is represented by
\[
\Pi_5(c)=cCP_{12}.
\]
The exact arithmetic audit gives rank five for this map and rank 124 for its integer kernel. After clearing denominators and taking the saturated integer kernel, let \(B\) be a rank-124 integer basis. The corresponding edge-flux row lattice is
\[
L=BC\subset\mathbb Z^{192}.
\]
The rank statement is exact; it is not inferred from floating-point dimension counting.

The annihilator is
\[
H_L=\operatorname{Ann}(L)\subset\mathbb T^{192}.
\]
If \(ULV=S\) is an SNF decomposition with nonzero diagonal entries \(d_1,\ldots,d_{124}\), then
\[
H_L\cong \mathbb T^{68}\times
\prod_{i=1}^{124}\mathbb Z/d_i\mathbb Z.
\]
For the audited Q6 lattice all 124 invariant factors are one. Thus the finite component group is trivial and
\[
\boxed{H_L\cong\mathbb T^{68}.}
\]
This is an exact connectedness and dimension result.

## 4. Why 68 = 63 + 5 is structurally meaningful

The gauge subgroup in link-angle space is
\[
G=\operatorname{im}(D_0^T:\mathbb T^{64}\to\mathbb T^{192}),
\]
with the constant vertex phase in its kernel, hence \(\dim G=63\). Every cycle pairs trivially with a gradient:
\[
C D_0^T=0.
\]
Consequently \(G\subset H_L\). Since \(H_L\) is a connected 68-torus, the quotient by this connected 63-dimensional subtorus is a connected 5-torus:
\[
H_L/G\cong\mathbb T^5.
\]
The isomorphism is noncanonical until an integral splitting or a compatible SNF basis is fixed.

This explains the arithmetic decomposition
\[
192-124=68=\underbrace{63}_{\text{vertex-gauge directions}}
+\underbrace{5}_{\text{retained quotient directions}}.
\]
It is the strongest current structural interpretation of the 68 dimensions. The last identification—showing that the quotient coordinates are exactly the chosen D5/coexact coordinates, rather than merely another five-torus—requires an explicit integral comparison of the gauge sublattice, the SNF splitting and the five-mode projector.

## 5. Haar integration, Bessel characters and Gauss law

For \(\beta>0\), the normalized one-link factor has the absolutely convergent character expansion
\[
\frac{e^{\beta\cos x}}{I_0(\beta)}
=\sum_{n\in\mathbb Z}\frac{I_{|n|}(\beta)}{I_0(\beta)}e^{inx}.
\]
On \(\mathbb T^{192}\), the product expansion is indexed by integer edge fluxes \(n\in\mathbb Z^{192}\). Haar orthogonality on \(H_L\) projects onto the annihilator of the annihilator lattice, namely the corresponding saturated integer flux lattice. This is the exact Fourier origin of the complete fibre sum.

The gauge subgroup gives a second orthogonality condition. Averaging a character over \(G\) yields zero unless
\[
D_0 n=0,
\]
which is the integer lattice Gauss constraint. Thus the 63 gauge directions are not extra physical polarizations: in the character representation their Haar average enforces charge conservation/Gauss law. They cannot simply be dropped from the original integral without carrying out this quotient or conditional averaging correctly.

**Important distinction.** The bare factor \(\prod_e e^{\beta\cos x_e}\) is not itself invariant under local vertex gauge transformations. Gauge invariance emerges in the correct character/constraint construction (or from a plaquette Wilson action), not by assertion. Therefore an integral over \(H_L\) must not be casually called an integral over the physical \(\mathbb T^5\). One must specify the observable, the gauge average, and the induced measure on the quotient.

## 6. Relation to the other torus results in the repository

- **Flat P1 torus (G0–G6):** a genuine continuum manifold \(\mathbb T_L^d\), with a proved Mosco limit for the regular flat periodic P1 family. This is a rigorous analytic prototype, not a curved RGS-Q6 limit.
- **Fourier Wilson torus:** a periodic momentum-space model with a fixed Hilbert embedding. Its symbol and scaling window belong to that specific flat translation-invariant model; they do not follow from the finite Q6 graph.
- **Q6 link and cycle tori:** finite-dimensional compact groups determined by the incidence matrices and integer homology. These are the correct setting for U(1) holonomy, integer electric flux and Wilson character sums.
- **PSU–Anti-Twin / Möbius / double-torus pictures:** possible geometric or field-configuration ansätze. To connect them to the Q6 compact groups, one needs a specified map (for example, a coframe-induced connection and its holonomy representation), not a visual resemblance or equality of dimensions.

## 7. What is proved, what is not

**Exact algebra/topology (under the stated Q6 incidence and saturated-kernel construction):**
- \(\operatorname{rank}D_0=63\), cycle rank 129;
- rank-five rational coexact projector and rank-124 integer fibre kernel;
- SNF invariants all one for the audited annihilator lattice;
- connected annihilator \(H_L\cong\mathbb T^{68}\);
- inclusion of the 63-dimensional vertex-gauge subgroup and a five-dimensional quotient torus;
- character orthogonality as the source of integer Gauss constraints and the Bessel fibre representation.

**Not yet proved:**
- that the quotient \(H_L/G\) is canonically the specific D5 coexact coordinate torus;
- that the induced quotient measure/conditional fibre sum is a D5 Hua measure;
- that the five retained modes define three physical spatial dimensions plus time;
- that the Q6 torus groups converge to a curved RGS–PSU continuum geometry;
- that the PSU–Anti-Twin embedded toroid is dynamically derived from the Q6 connection;
- any prediction of the fine-structure constant from this torus construction.

## 8. The right next theorem

The next task should not be a larger QMC run. It should be an exact integral-structure theorem:

1. Compute a unimodular basis adapted simultaneously to \(\operatorname{im}D_0^T\subset H_L\) and the saturated five-mode quotient.
2. Prove the induced map \(H_L/G\to\mathbb T^5\) is an isomorphism of compact groups, or state the precise finite kernel/cokernel if it is not.
3. Push Haar measure through this map and express the Wilson character integral as a conditional fibre average.
4. Determine whether the resulting five-variable density is invariant under the relevant D5 symmetry and whether it has any relation to the Hua density. This last step must be derived, not assumed.

Only after this theorem should numerical integration be used, with a certified error bound if numerical values are needed.
