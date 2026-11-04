# RGS fourth moment → Maxwell constitutive bridge

## 1. Input already established by the 64-site audit

For the hash-verified RGS coordinate seed (X_n\in\mathbb R^3), (n=1,\ldots,64), define
\[
M^{(4)}_{ijkl}=\sum_n X_{n,i}X_{n,j}X_{n,k}X_{n,l}.
\]
Its Frobenius-orthogonal (SO(3))-isotropic projection is
\[
M^{(4),iso}_{ijkl}=c(\delta_{ij}\delta_{kl}+\delta_{ik}\delta_{jl}
+\delta_{il}\delta_{jk}),\qquad
c=\frac{M^{(4)}_{iikk}}{3(3+2)},
\]
and (R=M^{(4)}-M^{(4),iso}) is fully symmetric, trace-free on every pair contraction, and (O_h)-invariant. The existing audit verifies these finite-seed identities and the (E_g\oplus T_{2g}) splitting.

## 2. The tensorial no-go: no direct (RFF) Maxwell term

The spatial electromagnetic field strength is antisymmetric:
\[
F_{ij}=-F_{ji}.
\]
Since (R_{ijkl}=R_{jikl}), contraction in the first pair gives
\[
R_{ijkl}F_{ij}F_{kl}=0.
\]
This is an exact index-symmetry obstruction, not a numerical cancellation. Thus the tempting direct rank-four correction (R_{ijkl}F_{ij}F_{kl}) is identically zero. The fully symmetric fourth moment cannot be identified with a Riemann-type constitutive tensor or teleparallel torsion.

## 3. Correct intermediate: an STF order parameter

Let (Q_{ij}=Q_{ji}), (Q_{ii}=0), be a separately defined internal/coframe order parameter. The fourth moment induces
\[
A_{ij}(Q)=R_{ijkl}Q_{kl}.
\]
Because (R) is symmetric and trace-free, (A(Q)\in\operatorname{Sym}_0(3)). The (O_h)-equivariance is
\[
A_{g\cdot R}(gQg^T)=gA_R(Q)g^T.
\]
The anisotropy can therefore enter the electric and magnetic response matrices:
\[
Z_E(Q)=z_E I+g_E A(Q),\qquad
Z_B(Q)=z_B I+g_B A(Q),
\]
with (z_E,z_B>0) and positivity conditions
\[
|g_E|\,\|A(Q)\|_2<z_E,\qquad
|g_B|\,\|A(Q)\|_2<z_B.
\]
The 3+1 constitutive density is
\[
\mathcal L_{EM}(Q)=\frac12 E_i(Z_E)_{ij}E_j
-\frac12 B_i(Z_B)_{ij}B_j.
\]
This is gauge invariant because it depends on (F=dA), and it is covariant under simultaneous transformations of the internal frame, (R), (Q), (E), and (B). The implementation checks positive definiteness rather than clipping negative eigenvalues.

## 4. What this does and does not imply for Wilson → Maxwell

The compact Q6 Wilson marginal and its five retained modes are already defined in the repository. To insert the RGS fourth moment, one still needs a specified map
\[
T:\mathbb R^5_{\text{Q6 retained}}\longrightarrow\operatorname{Sym}_0(3),
\qquad Q=Tq,
\]
with the correct group action and equivariance. Equal dimensions (five) do not establish such an intertwiner: the Q6 (A_5) mode representation and the RGS cubic (O_h) representation are not to be identified by dimension alone.

If such a map is derived, (Z_E(Tq)) and (Z_B(Tq)) define a field-dependent constitutive response to be combined with the Wilson effective action. The existing quadratic Gaussian Wilson marginal gives (12\beta_s I_5) and only a (q)-independent determinant; it does not generate this tensor or a Hua factor. The nonlinear compact marginal must supply the measure and dynamics for (q) (or (Q)).

There is also a useful consequence for the vacuum normalization. If (Q) is integrated with a centered, (Q\mapsto-Q) symmetric measure, then (\langle A(Q)\rangle=0): there is no first-order anisotropic correction to the vacuum Maxwell quadratic term. Integrating a Gaussian (Q) coupled linearly to (F^2) produces an (F^4) term at leading nontrivial order, not a shift of the scalar (F^2) coefficient. A nonzero background (Q_0), a non-centered measure, or Wilson-induced field dependence would have to be derived before a quadratic anisotropic response appears.

Consequently this bridge provides a well-defined **candidate constitutive sector**, but it does not yet determine the scalar (Z_F), the renormalized electric charge, or
\[
\alpha=\frac{1}{4\pi Z_F}
\]
(in rationalized conventions). No Wyler normalization or Hua cutoff is inserted.

## 5. The representation-theoretic bridge: restrict to the common tetrahedral A4

There is a nontrivial route between the two five-dimensional sectors, but it is **not** an identification of the full groups (A_5) and (O_h).

The full groups cannot be identified by an embedding: (A_5) is simple of order 60, while (O_h) has order 48. A homomorphism (A_5\to O_h) is therefore trivial (an injective image would have order 60 dividing 48). So a nonzero full-(A_5)-equivariant map into an (O_h)-representation is not the right demand.

Instead, both groups contain a tetrahedral subgroup (A_4):
- (A_4\subset A_5) as the stabilizer of one letter in the natural five-letter action;
- (A_4\subset O_h) as the proper tetrahedral rotation subgroup.

For the five-dimensional irreducible representation of (A_5), the character on the four (A_4) conjugacy classes (identity, double transpositions, and the two classes of 3-cycles) is
\[
\chi_{5\downarrow A_4}=(5,1,-1,-1).
\]
The (A_4) irreducible characters are
\[
\chi_{1}=(1,1,1,1),\quad
\chi_{1'}=(1,1,\omega,\omega^2),\quad
\chi_{1''}=(1,1,\omega^2,\omega),\quad
\chi_3=(3,-1,0,0),
\qquad \omega=e^{2\pi i/3}.
\]
Their character inner products give
\[
\boxed{5\downarrow A_4=1'\oplus1''\oplus3.}
\]
On the RGS side, the cubic decomposition (\operatorname{Sym}_0(3)=E_g\oplus T_{2g}), restricted to the same tetrahedral subgroup, has exactly the same decomposition: (E_g\downarrow A_4=1'\oplus1'') (as a real 2-plane) and (T_{2g}\downarrow A_4=3).

Therefore an (A_4)-equivariant real isomorphism
\[
T:V_5^{A_5}\longrightarrow\operatorname{Sym}_0(3)
\]
**exists after restriction to (A_4)**. It is not canonical: one must choose compatible subgroup generators and orthonormal bases, then solve the finite intertwiner equations
\[
T\rho_{Q6}(h)=\rho_{STF}(h)T,\qquad h\in\mathcal G_{A_4},
\]
and fix the independent normalization/sign on the 2- and 3-dimensional irreducible blocks. This is the concrete next algebraic construction; it replaces the invalid full-group identification with a valid common-subgroup bridge.

## 6. Representation type: the Dyson check

The common-subgroup statement is a statement about real representations, not evidence that either five-dimensional sector is quaternionic.

Over the real numbers,
`1' + 1''` is the real two-dimensional irreducible representation of the tetrahedral group (the two complex one-dimensional characters are exchanged by conjugation). The three-dimensional block is of real type. Consequently the real commutant of the restricted five-dimensional representation is
[
\operatorname{End}_{A_4}(V_5)\simeq\mathbb C\oplus\mathbb R.
]
For orthogonal intertwiners this leaves a phase on the two-plane and a sign on the three-dimensional block, once compatible invariant inner products and bases are fixed.

Dyson's Threefold Way classifies irreducible symmetry types by the real, complex, and quaternionic division algebras. It does **not** make this representation quaternionic. A quaternionic structure requires an antiunitary operator (J) commuting with the group action and satisfying (J^2=-I). No such (J) is established here. In particular, the real five-dimensional irreducible representation of (A_5) is of real type. A possible quaternionic structure on a doubled PSU–Anti-Twin space is a separate construction and remains open.

## 7. Hermann's bundle-geometric formulation

Let (P\to M) be a principal (G)-bundle and (\rho:G\to GL(V)) a representation. Its associated vector bundle is (E=P\times_\rho V). A principal connection (\omega) induces
[
\nabla=d+\rho_*(\omega),\qquad
\Omega=d\omega+\tfrac12[\omega,\omega].
]
An (A_4)-intertwiner (T:V_{Q6}\to V_{RGS}) induces a bundle map between the associated bundles when the reductions and transition functions use the same (A_4) action. It is connection-preserving only if the two induced connections satisfy
[
\nabla_{RGS}\circ\Phi=\Phi\circ\nabla_{Q6}.
]
This is an additional compatibility equation; representation equivalence alone does not imply it. The construction follows the gauge-bundle/Cartan–Ehresmann framework treated by Robert Hermann.

Thus the concrete next computation is: choose and document common (A_4) generators; extract their matrices on the five Q6 retained modes and on (E_g\oplus T_{2g}\); solve (T\rho_Q(h)=\rho_R(h)T) for generators; normalize the two blocks; then test connection compatibility and curvature intertwining. Do not label the map quaternionic unless the antiunitary (J) test above is passed.

## 8. Evidence status

- **Exact tensor identity:** for fully symmetric (R) and antisymmetric (F), (R_{ijkl}F_{ij}F_{kl}=0).
- **Exact representation-theory consequence:** the two five-dimensional real sectors restrict to the same (A_4) type (2\oplus3); an (A_4)-intertwiner exists, but is not canonical.
- **Exact Dyson classification consequence:** the real commutant is (\mathbb C\oplus\mathbb R); no quaternionic structure follows.
- **Numerically audited:** the 64-site seed decomposition, (O_h) covariance and constitutive positivity in the documented test.
- **Candidate:** the (Q)-dependent electric/magnetic constitutive matrices.
- **Open:** compute the actual intertwiner matrices; derive (Q) from the coframe/PSU–Anti-Twin geometry; derive the (Q)-measure from the complete nonlinear Wilson fibre; establish connection compatibility and the continuum limit; determine scalar Maxwell normalization and only then compare with Wyler.

## 9. Reproduction and references

\`\`\`bash
PYTHONPATH=. python -m pytest tests/test_rgs_fourth_moment_maxwell_bridge.py -q
PYTHONPATH=. python -m src.rgs_fourth_moment_maxwell_bridge
\`\`\`

The rank-four seed construction is documented in [RGS_64_FOURTH_MOMENT_BRIDGE.md](RGS_64_FOURTH_MOMENT_BRIDGE.md). The Wilson lattice-gauge and transfer-matrix framework is Kogut, *Rev. Mod. Phys.* 51 (1979), 659–713, DOI: 10.1103/RevModPhys.51.659. The general constitutive-tensor formulation is treated in Hehl & Obukhov, *Foundations of Classical Electrodynamics* (Birkhäuser, 2003), DOI: 10.1007/978-94-010-0205-8. These references support the general frameworks; they do not prove the RGS-specific coupling or its physical interpretation.
