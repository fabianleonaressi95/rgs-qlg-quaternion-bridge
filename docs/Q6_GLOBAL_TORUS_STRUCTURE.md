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

**Literature and repository audit:** see [Q6–D5 literature and repository evidence map](Q6_D5_LITERATURE_AND_REPOSITORY_EVIDENCE_MAP.md) for external references, the internal proof-file map and open obligations.

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
The **graph cycle lattice** is \(Z_1=\ker(\partial_1)\subset\mathbb Z^{192}\), where \(\partial_1\) maps oriented edges to their endpoint difference; its rank is \(E-V+1=129\). This is the cycle rank of the 1-skeleton, not the first Betti number before face boundaries are imposed. The 240 square faces define \(\partial_2:\mathbb Z^{240}\to Z_1\), with \(\partial_1\partial_2=0\), and \(H_1(Q6;\mathbb Z)=Z_1/\operatorname{im}\partial_2\). For the full six-cube, contractibility gives \(H_1=0\); equivalently the square-face boundaries generate the full integral cycle lattice and the quotient has no torsion. If one were given only the three cell counts, without the cubical boundary maps or contractibility, this conclusion would not follow.

This distinction separates graph cycles (the domain used by the integer-cycle coordinates below) from physical flat holonomies (which also quotient plaquette boundaries). It does not select five physical dimensions or turn the six-dimensional Q6 graph into 3+1-dimensional spacetime.

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
- **Q6 link and cycle tori:** finite-dimensional compact groups determined by incidence matrices and integer lattices. The graph cycle torus, the flat-connection moduli and the Wilson annihilator must be computed as distinct quotients; the latter is not obtained by merely naming the graph cycle rank as \(b_1\).
- **PSU–Anti-Twin / Möbius / double-torus pictures:** possible geometric or field-configuration ansätze. To connect them to the Q6 compact groups, one needs a specified map (for example, a coframe-induced connection and its holonomy representation), not a visual resemblance or equality of dimensions.

## 7. What is proved, what is not

**Exact algebra/topology (under the stated Q6 incidence and saturated-kernel construction):**
- \(\operatorname{rank}D_0=63\), graph cycle rank 129;
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


## 9. Harmonic-analysis interpretation and the two-point Smith example

The structural statements above are instances of Pontryagin duality. The character group of
\(\mathbb T^E\) is \(\mathbb Z^E\), with pairing
\[
\langle n,x\rangle=e^{i n\cdot x}.
\]
For a closed subgroup \(H\leq\mathbb T^E\), its annihilator
\[
H^\perp=\{n\in\mathbb Z^E:e^{in\cdot x}=1\;\forall x\in H\}
\]
is a subgroup of the integer character lattice, and \((H^\perp)^\perp=H\). Consequently normalized Haar integration is a projector on Fourier modes:
\[
\int_H e^{in\cdot x}\,d\mu_H(x)
=\begin{cases}1,&n\in H^\perp,\\0,&n\notin H^\perp.\end{cases}
\]
This is the exact mathematical operation behind the integer selection rules; it is not a numerical cancellation or a heuristic averaging argument.

For a full-row-rank integer matrix \(A\in\mathbb Z^{r\times E}\), the subgroup
\(H=\{x:Ax=0\pmod{2\pi}\}\) is classified by the Smith form
\(UAV=\operatorname{diag}(d_1,\ldots,d_r,0,\ldots,0)\).
In Smith coordinates each nonzero diagonal equation has \(d_i\) solutions on the circle, while each zero column contributes one free circle. Thus
\[
H\simeq\mathbb T^{E-r}\times\prod_i\mathbb Z/d_i\mathbb Z.
\]
The free dimension is \(E-r\); the torsion factors describe disconnected components. In particular, rank alone does not prove connectedness: one must also prove that every nonzero invariant factor is one. This is why the saturated integer lattice and its SNF, rather than a floating-point nullspace, are decisive.

The toy test in the failed QMC branch is an exact miniature of this theorem. For
\(A=(2\;0)\), the annihilator is
\[
H=\{(x_1,x_2):2x_1=0\}
=\{0,\pi\}\times\mathbb T.
\]
With normalized Haar measure and \(\beta>0\), define
\[
J_m=\int_H\frac{e^{\beta(\cos x_1+\cos x_2)}}{I_0(\beta)^2}
e^{-imx_1}\,d\mu_H.
\]
The free-circle integral gives \(I_0(\beta)\), while the two-point average gives
\[
J_0=\frac{\cosh\beta}{I_0(\beta)},\qquad
J_1=\frac{\sinh\beta}{I_0(\beta)}.
\]
Both are real because the measure and weight are invariant under complex conjugation paired with \(x\mapsto-x\). A numerical implementation must form the exponent separately at each sample; summing \(\cos x\) over the entire sample-by-coordinate array changes the mathematical integrand and causes overflow. That is an implementation error, not a failure of the character identity.

## 10. What “toroid” can mean in the PSU–Anti-Twin programme

There are at least three geometrically different objects that must remain distinguished:

- A **torus as a manifold**, such as \(S^1\times S^1\), embedded in a three-dimensional ambient space. Its genus, induced metric, second fundamental form and extrinsic curvature are geometric data.
- A **torus as a compact group**, such as \(U(1)^k\), which is a space of phases, holonomies or moduli. It has a group law and Haar measure; it is not automatically an embedded surface.
- A **toroidal field configuration**, for example a pair of coupled sheets, a Möbius-twisted bundle, or a PSU–Anti-Twin ovaloid. This is a field/geometric ansatz whose topology, metric, connection and dynamics must be specified separately.

A map between these objects can be built only after specifying the data. A useful route is a principal \(U(1)\) bundle with connection \(A\): its holonomy gives a character of closed cycles, while its curvature \(F=dA\) controls local plaquette flux. A coframe \(e^a\) and a metric-compatible connection then provide the RGS/teleparallel geometric data. To identify the Q6 phase torus with a PSU toroidal configuration, one must construct a map from the discrete connection/coframe data to the continuum bundle and prove compatibility of holonomy, curvature, gauge transformations and the limiting measure. Equality of dimensions, a shared word “torus”, or a visual resemblance is not such a proof.

For the proposed \(68=63+5\) split, the clean mathematical target is therefore a diagram of compact groups
\[
\mathbb T^{64}/\mathbb T^1\longrightarrow H_L\longrightarrow H_L/G\simeq\mathbb T^5,
\]
together with the exact character-lattice maps dual to each arrow. The first quotient removes the constant vertex phase; the second arrow is the inclusion of gauge directions into the annihilator; the last quotient is the retained five-dimensional group. The diagram must be realized integrally (unimodularly), not only at the level of real vector spaces. Its continuum counterpart would then be a separate convergence theorem for connections and measures.

## 11. Literature and scope of the import

The framework uses established results, not a new general theory of tori:

- Pontryagin duality identifies compact abelian groups with the duals of discrete character groups; annihilators exchange closed subgroups and quotient groups.
- Smith normal form classifies integer sublattices and the connected/disconnected structure of their annihilator subgroups in a torus.
- Haar orthogonality is the exact Fourier projector onto characters trivial on the subgroup.
- The Wilson character expansion converts compact-link Boltzmann weights into integer representation/flux sums; link integration imposes local flux constraints. In compact \(U(1)\), the coefficients are modified Bessel functions \(I_n(\beta)\).
- In Hamiltonian lattice gauge theory, Gauss's law selects the gauge-invariant Hilbert subspace; it is not an optional deletion of link variables.

References:

1. S. A. Morris, *Pontryagin Duality and the Structure of Locally Compact Abelian Groups*, Cambridge University Press (1977), especially the annihilator and quotient-duality theorems.
2. T. Tao, “Using the Smith normal form to manipulate lattice subgroups and closed torus subgroups” (20 July 2022), an explicit account of the integer-lattice/closed-subtorus correspondence.
3. J. B. Kogut, “An Introduction to Lattice Gauge Theory and Spin Systems,” *Rev. Mod. Phys.* **51** (1979), 659–713.
4. J. B. Kogut and L. Susskind, “Hamiltonian Formulation of Wilson's Lattice Gauge Theories,” *Phys. Rev. D* **11** (1975), 395–408, DOI: 10.1103/PhysRevD.11.395.
5. K. G. Wilson, “Confinement of Quarks,” *Phys. Rev. D* **10** (1974), 2445–2459, DOI: 10.1103/PhysRevD.10.2445.

These references justify the mathematical tools and the standard lattice-gauge mechanism. They do not establish the proposed physical identification of the Q6 quotient with a D5 Hua space, nor the PSU–Anti-Twin continuum interpretation; those remain project-specific theorems to prove.


## 12. The Liouville–Bergman law: exact statement and dynamical meaning

The five-mode Hamiltonian sector is a ten-real-dimensional symplectic manifold. In canonical coordinates
\((q,p)\in\mathbb R^5\times\mathbb R^5\),
\[
\omega_{\mathrm{can}}=\sum_{j=1}^5dq_j\wedge dp_j,\qquad
d\mu_L=\frac{\omega_{\mathrm{can}}^5}{5!}
=\prod_{j=1}^5dq_j\,dp_j.
\]
For
\[
H_5(q,p)=\frac12\bigl(\beta\lambda\,|q|^2+
\beta^{-1}|p|^2\bigr),\qquad
\Omega_H=\sqrt{\lambda},
\]
the PR #145 chart uses the positive scale \(\kappa_{145}=\sqrt{\beta\lambda}\):
\[
\Phi_{145}(q,p)=z=\frac{\sqrt{\kappa_{145}}\,q+
i\,p/\sqrt{\kappa_{145}}}{\sqrt{20}}.
\]
This is a valid real-linear volume chart, and its Jacobian is independent of the positive scale. However, for the Hamiltonian written above it is not the oscillator action-angle chart at general \(\beta\). The flow-adapted canonical scale is instead
\[
\kappa_H=\beta\sqrt{\lambda},
\qquad
\Phi_H(q,p)=\frac{\sqrt{\kappa_H}\,q+
i\,p/\sqrt{\kappa_H}}{\sqrt{20}}.
\]
Indeed Hamilton's equations give \(\dot q=\beta^{-1}p\), \(\dot p=-\beta\lambda q\), so
\(\Phi_H(t)=e^{-i\Omega_Ht}\Phi_H(0)\). The two scales agree when \(\beta=1\), or can be reconciled by an explicit canonical rescaling of \(q,p\). Both charts have the same real Jacobian \(20^{-5}\), so the exact volume pullback law below holds for either scale. Its relevant domain is
\[
\Omega_5=\Phi^{-1}(D_5),\qquad
D_5=\{z:N(z)>0,\ |z^Tz|<1\},
\]
where
\[
N(z)=1-2z^\dagger z+|z^Tz|^2.
\]
Thus the chart identifies \(\Omega_5\), not all of \(\mathbb R^{10}\), with \(D_5\).

For the type-IV Bergman metric in the normalization fixed by PR #145,
\[
\det g_B(z)=10^5N(z)^{-5},\qquad
dV_B=2^5\det(g_B)\,d^{10}z.
\]
The real Jacobian is \(|\det D\Phi|=20^{-5}\). Consequently
\[
\boxed{\Phi^*dV_B=N(\Phi(q,p))^{-5}\,d\mu_L.}
\]
This is an exact volume-form identity on \(\Omega_5\), not a fitted density. It does not by itself imply that the Bergman Kähler form equals the canonical symplectic form away from the origin; PR #145 proves that equality at the origin only.

### Liouville transport equation

Let \(\varphi_t\) be the Hamiltonian flow of any smooth \(H\) on this chart. Since
\(\mathcal L_{X_H}d\mu_L=0\), a density \(\rho\,d\mu_L\) is invariant precisely when
\[
X_H(\rho)+\rho\,\operatorname{div}_{\mu_L}X_H=0.
\]
For the Bergman pullback density \(\rho_B=N^{-5}\), this becomes
\[
\boxed{\mathcal L_{X_H}(\rho_Bd\mu_L)=0
\quad\Longleftrightarrow\quad X_HN=0.}
\]
The equivalence holds wherever \(N>0\). Thus the Bergman-weighted measure is not invariant under every Hamiltonian: it is invariant exactly when the Hamiltonian flow preserves the generic norm.

For the isotropic fivefold oscillator above, the flow-adapted chart \(\Phi_H\) gives
\[
q(t)=q(0)\cos(\Omega_Ht)+
\frac{p(0)}{\beta\Omega_H}\sin(\Omega_Ht),
\]
\[
p(t)=p(0)\cos(\Omega_Ht)-
\beta\Omega_Hq(0)\sin(\Omega_Ht),
\qquad \Omega_H=\sqrt\lambda,
\]
and therefore \(z_H(t)=e^{-i\Omega_Ht}z_H(0)\). Both
\(z_H^\dagger z_H\) and \(|z_H^Tz_H|\) are invariant under this common phase rotation. Hence \(N(z_H(t))=N(z_H(0))\), and
\[
\boxed{\varphi_t^*(N^{-5}d\mu_L)=N^{-5}d\mu_L.}
\]
This is the precise dynamical Liouville statement for the degenerate \(\lambda\)-eigenspace. For anisotropic mode frequencies, a common phase rotation is absent and invariance of \(N\) must be checked anew; it cannot be assumed.

## 13. Hua normalization and the cutoff probability law

Let \(d\nu_E\) denote normalized Euclidean volume on \(D_5\). The Hua identity in the normalization used by PR #146 is
\[
\int_{D_5}N(z)^s\,d\nu_E(z)
=\frac{60}{(s+1)(s+2)(s+3)(s+4)(s+5/2)},
\quad \Re s>-1.
\]
The induced density of \(t=N(z)\) under \(d\nu_E\) is
\[
f_N(t)=\frac{20}{3}(1-\sqrt t)^4(1+4\sqrt t+t),
\quad 0<t<1.
\]
It follows that the Bergman volume is not finite without a cutoff, because its exponent is \(s=-5\). For \(0<\varepsilon<1\),
\[
Z_5(\varepsilon)=\int_\varepsilon^1t^{-5}f_N(t)\,dt
\]
is exactly the finite expression already recorded in PR #149, with
\[
Z_5(\varepsilon)\sim\frac53\varepsilon^{-4}.
\]
The normalized cutoff measure is
\[
d\mu_{B,\varepsilon}(z)=
\frac{\mathbf1_{\{N(z)\geq\varepsilon\}}N(z)^{-5}}
{Z_5(\varepsilon)}\,d\nu_E(z).
\]
This is a genuine probability measure for each fixed \(\varepsilon>0\). Its uncut limit is not a finite probability measure: the normalizer diverges and mass escapes toward the boundary \(N=0\). No physical value of \(\varepsilon\) is fixed by Hua's theorem.

## 14. Exact Wilson-to-Liouville measure matching: the required theorem

The compact-group reduction and the Liouville law live on spaces of different dimensions and different measure types. They must be connected in two stages.

First, after the exact integral splitting of the constrained flux lattice, write
\[
\Lambda_{\rm cyc}\simeq\mathbb Z^{129},\qquad
K=\ker(P:\Lambda_{\rm cyc}\to\mathbb Z^5)
\simeq\mathbb Z^{124},\qquad
\Lambda_{\rm cyc}/K\simeq\mathbb Z^5.
\]
Use the row-vector convention of the repository: \(C\in\mathbb Z^{129\times192}\), \(B\in\mathbb Z^{5\times129}\), and \(K_0\in\mathbb Z^{124\times129}\), with the stacked matrix
\[
\begin{bmatrix}B\\K_0\end{bmatrix}\in GL(129,\mathbb Z).
\]
Every constrained cycle-coordinate flux has the unique form
\[
c=m^TB+k^TK_0,\qquad m\in\mathbb Z^5,\quad k\in\mathbb Z^{124},
\]
and its edge-flux vector is \(n=cC\in\mathbb Z^{192}\). For the factorized compact-U(1) link character weight, the exact marginal is
\[
W_\beta(m)=
\sum_{k\in\mathbb Z^{124}}
\prod_{e=1}^{192}
\frac{I_{|((m^TB+k^TK_0)C)_e|}(\beta)}{I_0(\beta)},
\]
with the understanding that the actual Q6 cycle embedding and any additional plaquette/temporal indices must be included exactly as specified by the chosen Wilson action. The absolute convergence and tail control of the fibre sum are the subject of the earlier integer-fibre and Bessel-majorant results. Define
\[
\mathcal Z_\beta=\sum_{m\in\mathbb Z^5}W_\beta(m).
\]
Then the exact quotient flux law is the discrete probability
\[
\mathbb P_\beta(m)=W_\beta(m)/\mathcal Z_\beta.
\]
The angular quotient carries normalized Haar measure \(d\theta/(2\pi)^5\) on \(\mathbb T^5\). Thus the exact reduced Wilson object is a measure on
\(\mathbb T^5\times\mathbb Z^5\), not yet a smooth measure on a ten-dimensional cotangent bundle.

Second, to compare with the five-mode Liouville law, one must specify a canonical cotangent lift and a continuum interpolation of the integer flux. At finite coupling the measure supported on integer momenta is singular with respect to
\(d^5\theta\,d^5p\); therefore it cannot equal the absolutely continuous Bergman measure. A smoothing/interpolation kernel is necessary and is additional mathematical data, not a consequence of Haar orthogonality. For a normalized kernel \(K_\sigma(p-m)\), one possible regularized measure is
\[
d\widetilde\mu_{\beta,\sigma}(\theta,p)=
\frac{d^5\theta}{(2\pi)^5}
\sum_{m\in\mathbb Z^5}\mathbb P_\beta(m)
K_\sigma(p-m)\,d^5p.
\]
A scale map from \((\theta,p)\) to canonical \((q,p_{\rm can})\), preserving
\(\sum dq_j\wedge dp_{{\rm can},j}\), and a chart into \(\Omega_5\) must then be fixed.

Let \(F_{\beta,\sigma}\) denote this map on a chosen fundamental-cell chart, and let \(r_{\beta,\sigma}(z)\,d\nu_E(z)\) be the density of its pushforward, if it is absolutely continuous. The exact finite-cutoff matching statement is the Radon–Nikodym identity
\[
\boxed{
F_{\beta,\sigma*}\widetilde\mu_{\beta,\sigma}
=d\mu_{B,\varepsilon}
\quad\Longleftrightarrow\quad
r_{\beta,\sigma}(z)=
\frac{\mathbf1_{\{N(z)\geq\varepsilon\}}N(z)^{-5}}
{Z_5(\varepsilon)}
\quad(d\nu_E\text{-a.e.}).
}
\]
If \(F_{\beta,\sigma}\) is many-to-one, the left-hand density is computed by the coarea/disintegration formula, summing the conditional Wilson weight over the full fibre with the inverse Jacobian. This is where the exact Bessel fibre sum must enter; it cannot be replaced by a five-mode projection before marginalization.

This equivalence is the precise theorem to prove for a Wilson derivation of the Liouville–Bergman law. It is not yet proved by the separate Haar, SNF, Bergman or Hua identities.

## 15. What the Bessel asymptotic does—and does not—give

For large \(\beta\) with flux \(n=O(\sqrt\beta)\), the local Gaussian regime of the modified Bessel coefficients gives, at leading order,
\[
\frac{I_n(\beta)}{I_0(\beta)}
\approx \exp\!\left(-\frac{n^2}{2\beta}\right).
\]
This is the standard large-parameter/large-order asymptotic regime of modified Bessel functions; see NIST DLMF §§10.40–10.41. It motivates a Gaussian approximation to the constrained integer-flux sum, but it is not an exact identity at finite \(\beta\).

If the quadratic approximation is uniform on the relevant flux window and the lattice local-limit theorem applies, summing over the fibre produces an effective quadratic form on the five quotient fluxes. With the Euclidean edge quadratic form and real fibre subspace \(K_{0,\mathbb R}\), its leading exponent is governed by the Schur-complement metric
\[
Q_{\rm eff}=B^T(I-\Pi_K)B,\qquad
\Pi_K=K_0(K_0^TK_0)^{-1}K_0^T.
\]
This expression is independent of the chosen section \(B\) up to the induced integral change of quotient coordinates. For a nontrivial Wilson action, the edge quadratic form must be replaced by its actual positive Hessian, and the corresponding weighted orthogonal projection must be used.

The resulting Gaussian is only a leading local/weak-coupling approximation. It does not equal the nonlinear density \(N^{-5}\) on \(D_5\) by itself. The latter contains the full type-IV geometry and its boundary singularity. To obtain it, the exact nonlinear Jacobian of the canonical chart and the non-Gaussian corrections to the Wilson marginal must combine to give the Radon–Nikodym identity in §14. This is a testable mathematical target, not a consequence of the leading Bessel asymptotic.

## 16. A dimension/topology obstruction that prevents a false identification

The quotient \(H_L/G\simeq\mathbb T^5\) is five-dimensional and compact. The D5 domain is ten-dimensional, noncompact and contractible. Therefore there is no diffeomorphism
\[
\mathbb T^5\cong D_5,
\]
and a direct pushforward of five-dimensional Haar measure cannot produce a ten-dimensional absolutely continuous Bergman volume.

The correct candidate bridge is the cotangent lift
\[
T^*\mathbb T^5\simeq\mathbb T^5\times\mathbb R^5,
\]
which is ten-dimensional and carries its canonical Liouville measure. Even this is not globally diffeomorphic to \(D_5\): \(T^*\mathbb T^5\) has fundamental group \(\mathbb Z^5\), whereas \(D_5\) is contractible. A global identification is therefore impossible. A chart obtained by cutting the torus along a fundamental-cell boundary can be locally identified with an open subset of \(\mathbb R^{10}\); only on such a chart can the audited linear map \(\Phi\) be used. The cut locus and the boundary of \(\Omega_5\) must be handled in the measure construction.

This is not a defect of the five-mode theory. It specifies the correct scope: a local Darboux/Bergman chart or a measurable pushforward may exist, while a global symplectomorphism of the full cotangent torus with D5 cannot.

## 17. Status ledger and next proof obligations

**Exact results already established in the repository**
- the five-mode chart and the exact volume pullback \(N^{-5}d\mu_L\) (PR #145);
- equality of canonical and Bergman symplectic forms at the origin, not globally (PR #145);
- the Hua Mellin law, the exact generic-norm density and finite-cutoff normalizer (PRs #146, #148, #149);
- the integer fibre, saturated lattice, Bessel majorant and annihilator-torus structure (PRs #153, #156, #157, #159, #162);
- invariance of the Bergman-weighted measure under the isotropic common-phase five-mode oscillator, by \(X_HN=0\).

**New theorem not yet established**
- construct the exact quotient character map \(P\), integral section \(B\), and fibre basis \(K_0\) in the same basis as the five-mode coexact projector;
- derive the finite-\(\beta\) marginal \(W_\beta(m)\) for the complete Wilson action, with all temporal/plaquette indices retained;
- specify and justify the canonical scaling and smoothing/interpolation of integer electric fluxes;
- prove absolute continuity of the resulting pushforward and compute its Radon–Nikodym derivative;
- prove that derivative is exactly \(N^{-5}/Z_5(\varepsilon)\), with a cutoff \(\varepsilon\) derived from the Wilson data rather than inserted by hand;
- establish convergence (including tightness and boundary control) as the coupling/refinement limit is taken.

**Conclusion.** The Liouville law is now placed at the correct point in the architecture: it is the exact target density of the five-mode Bergman chart; the Q6–Wilson torus supplies a compact quotient and an integer flux marginal; the cotangent lift and a controlled continuum limit are the missing bridge. Neither the five-torus alone nor the leading Gaussian Bessel approximation proves the Hua law.

Additional references:
- NIST Digital Library of Mathematical Functions, §§10.40–10.41, modified-Bessel large-argument and large-order asymptotics: https://dlmf.nist.gov/10.40 and https://dlmf.nist.gov/10.41.
- J. B. Kogut, “An Introduction to Lattice Gauge Theory and Spin Systems,” *Rev. Mod. Phys.* 51 (1979), 659–713, DOI: 10.1103/RevModPhys.51.659.
- Hamiltonian Liouville theorem: preservation of the canonical symplectic volume under Hamiltonian flow; see, e.g., the standard symplectic-geometry treatment in the cited graduate text.


## 18. Unified RGS–PSU–Q6–Wilson–Liouville synthesis
The architecture is now summarized as a chain of different but compatible structures:
\[
Q6\text{ gauge quotient}\to\ker D_0^T\to V_5
\to T^*V_5\to D_5\to N^{-5}d\mu_L\to d\mu_{B,\varepsilon}.
\]
In parallel, the exact character route is
\[
\mathbb Z^{129}\to\mathbb Z^{129}/K_{\mathbb Z}\simeq\mathbb Z^5,
\quad
W_\beta([c])=\sum_{k\in K_{\mathbb Z}}w_\beta(c+k),
\quad H_L\simeq\mathbb T^{68}.
\]
These routes meet only after an explicit integral splitting, cotangent construction and measure pushforward. The following subsections state the interfaces and their status.

### 18.1 Q6 and the A5 mode sector

The six-cube has 64 vertices, 192 edges and 240 square faces. The U(1) link quotient has dimension 129; this is the graph-cycle rank and Coulomb tangent dimension, not the first Betti number of the complete cubical complex. The top coexact eigenvalue 12 has multiplicity five and carries the irreducible A5 character \((5,1,-1,0,0)\). It defines \(V_5\), and its cotangent lift is the ten-real-dimensional canonical phase space. The five modes are spectral degrees of freedom, not five homology classes or five spacetime dimensions.

### 18.2 Integer Wilson quotient

For a rational top-mode basis \(V_{\mathbb Q}\), the exact map \(M=CV_{\mathbb Q}\) has rank five. Its saturated integer kernel \(K_{\mathbb Z}\) has rank 124, so the quotient lattice is free of rank five. A unimodular completion \([B\ K_0]\in GL(129,\mathbb Z)\) gives the unique decomposition \(c=Bm+K_0k\). The complete marginal is
\[
W_{\beta_t}(m)=\sum_{k\in\mathbb Z^{124}}
\prod_{e=1}^{192}\frac{I_{|((Bm+K_0k)C)_e|}(\beta_t)}{I_0(\beta_t)}.
\]
Its absolute convergence follows from the factorial Bessel majorant; the certified edge-flux tail bound controls truncations. A finite sum or Sobol estimate is not the complete marginal.

The annihilator of the embedded saturated kernel is \(H_L\simeq\mathbb T^{68}\). The 63-dimensional vertex gauge torus lies inside it, and the quotient is abstractly \(\mathbb T^5\). The canonical identification of this quotient with the selected A5 mode coordinates still requires the explicit unimodular completion and intertwiner. The 68-torus, the 129-dimensional link quotient and the five-torus are distinct objects.

### 18.3 Maxwell and Liouville normalization

For \(H_5=\frac12(\beta\lambda|q|^2+\beta^{-1}|p|^2)\), the frequency is \(\Omega=\sqrt\lambda\). The flow-adapted complex scale is \(\kappa=\beta\sqrt\lambda\):
\[
z=\frac{\sqrt\kappa q+i\,p/\sqrt\kappa}{\sqrt{20}},
\qquad H_5=10\Omega z^\dagger z,\qquad z(t)=e^{-i\Omega t}z(0).
\]
The Bergman generic norm is \(N=1-2z^\dagger z+|z^Tz|^2\). With \(\det g_B=10^5N^{-5}\) and real chart Jacobian \(20^{-5}\),
\[
\Phi^*dV_B=N^{-5}d\mu_L.
\]
The weighted measure is invariant under a Hamiltonian \(H\) exactly when \(X_HN=0\). The isotropic five-mode oscillator satisfies this condition because its flow is a common phase rotation. This is not an assertion for arbitrary or anisotropic Hamiltonians.

### 18.4 Hua cutoff and the Wilson matching theorem

The normalized Euclidean Hua Mellin law is \(60/\chi_5(s)\), where
\[
\chi_5(s)=(s+1)(s+2)(s+3)(s+4)(s+5/2),\quad \Re s>-1.
\]
The generic-norm density is
\[
f_N(t)=\frac{20}{3}(1-\sqrt t)^4(1+4\sqrt t+t).
\]
The uncut Bergman volume diverges; \(Z_5(\varepsilon)=\int_\varepsilon^1t^{-5}f_N(t)\,dt\sim(5/3)\varepsilon^{-4}\). The cutoff is exact mathematics, but its physical value is not derived.

At finite coupling, the Wilson quotient is discrete in integer flux, while the Liouville target is absolutely continuous in ten real variables. A normalized smoothing/interpolation and cotangent lift are therefore required. The theorem to prove is the Radon–Nikodym identity
\[
F_{\beta,\sigma*}\widetilde\mu_{\beta,\sigma}
=\frac{\mathbf1_{N\geq\varepsilon}N^{-5}}{Z_5(\varepsilon)}d\nu_E.
\]
It requires the full fibre sum, the exact map \(F_{\beta,\sigma}\), its Jacobian/coarea formula and a cutoff derived from the Wilson data. The Gaussian Bessel limit is only leading asymptotics and cannot replace this proof.

A global diffeomorphism \(T^*\mathbb T^5\simeq D_5\) is impossible: their fundamental groups differ. A cut-cell local chart or measurable pushforward remains possible.

### 18.5 PSU–Anti-Twin geometric interface

A toroidal surface, a compact phase torus and a two-sheeted field configuration are different mathematical objects. The connection between Q6 and a continuum U(1) field is through link holonomies \(U_e=\exp(i\int_e A)\) and plaquette curvature \(\exp(i\int_fF)\), \(F=dA\). The vierbein is independent data: with spin connection \(\omega^a{}_b\), teleparallel torsion is \(T^a=de^a+\omega^a{}_b\wedge e^b\), and in Weitzenböck gauge \(T^a=de^a\). U(1) curvature and torsion have different bundle representations and cannot be identified without a specified coupled bundle, solder form and action.

An Anti-Twin model requires a double cover, sheet-exchange involution, transformation laws for \(A,e^a,\omega^a{}_b\), and gluing conditions. A discrete Q6×Z2 construction is a candidate, not by itself a proof of a wormhole or spacetime throat. The A5 action and its binary lift \(2I\subset SU(2)\) provide a route to the quaternionic/spin sector, but the equivariant bridge to the U(1) Wilson quotient and vierbein remains open.

### 18.6 Weyl/Wyler and modal status

The finite Euclidean D5 volume, the invariant Bergman volume and the Hua-cutoff probability are distinct measures. In the Wyler normalization \(V_E(D_5;R)=R^{10}\pi^5/1920\), so the associated fourth-root expression scales as \(R^{5/2}\). The Bergman volume diverges. No alpha prediction follows until the physical scale, observable measure, coupling renormalization and continuum limit are derived.

The microhistory identity \(|S_s\cdot h|/s!=1/\prod_e|n_e|!\) is an exact orbit-stabilizer statement and supports the Bessel factorial bound. It is not Gödel's ontological proof. A formal modal theorem requires its own language, axioms, interpretation and soundness argument.

### 18.7 Final status

**Exact:** Q6 incidence ranks; top coexact five-mode eigenspace and A5 character; saturated rank-124 kernel and rank-five free quotient; SNF annihilator \(\mathbb T^{68}\); exact Bessel character/fibre identity and factorial tail bound; Maxwell Legendre transform under \(\beta=g^{-2}\); Bergman volume pullback; isotropic weighted Liouville invariance; Hua Mellin law and finite cutoff.

**Open:** deterministic integral splitting aligned to the A5 frame; full spacetime Wilson marginal; Wilson-to-D5 Radon–Nikodym theorem; physical derivation of \(\varepsilon\); tight continuum limit; coupled U(1)–Spin–vierbein Anti-Twin construction; physical 3+1 reduction; alpha derivation; formal Gödel/modal theorem.


## Literature anchors for the torus and Hodge statements

The algebraic Hodge framework for incidence matrices and cohomology is treated in Lim, *Hodge Laplacians on Graphs*, SIAM Review 62 (2020), https://doi.org/10.1137/18M1223101. The integer normal-form and invariant-factor statements are standard Smith/Hermite theory; see https://flintlib.org/doc/algorithms/linear-algebra/hermite-smith.html. These references support the general machinery; the concrete Q6 matrices, the rank-124 saturated kernel, the unit Smith factors and the inclusion of the 63-dimensional gauge torus are checked by the repository's exact arithmetic.
