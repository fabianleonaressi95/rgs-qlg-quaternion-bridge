# Integrated mathematical suite — formalization register

**Version:** 2026-11-25 · **Scope:** the 51-step sequence added by PR #227.

This file converts the README sequence into typed mathematical records: definitions, claims, evidence status, dependencies and proof obligations. It is a formalization index and working specification; it does not replace specialist proofs linked below.

## 0. Evidence and dependency rules

Statuses:
- **PROVED** — mathematical proof under explicit hypotheses.
- **FINITE-EXACT** — exact finite combinatorics/algebra, with matrices or certificates retained.
- **NUMERICAL** — floating-point computation, sampling or finite numerical audit.
- **CONDITIONAL** — true only under a declared ansatz, metric, weights or map.
- **CANDIDATE** — constructed but its decisive property is not established.
- **OPEN** — proof or physical identification remains outstanding.

A passing test establishes only the tested implementation property. It does not, by itself, prove a universal theorem, a continuum limit or a physical interpretation. Every bridge between models must specify source space, target space, map, and the property the map is required to preserve. Similar dimensions, spectra or symmetry labels are not maps.

Core objects:
- Q6 graph: vertices {0,1}^6; edges join words at Hamming distance one.
- Q6 square complex: the graph plus its 240 elementary square 2-cells.
- B1:C1→C0 and B2:C2→C1 are oriented boundary maps, with B1 B2=0.
- X:V(Q6)→R^3 is the recovered finite RGS seed; its pullback Gram matrix has rank at most three.
- PRISM denotes the five-factor Cartesian graph C64 □ C6 □ C6 □ C6 □ C5.
- A5 lattice means {x∈Z^6: Σx_i=0}; it is not the alternating group A5.
- Q6 (64 vertices), RGS-188D operator, PRISM (69,120 vertices), Q-ball boundary-value problem and Pellis/GSR resonance operator remain distinct objects until explicit maps are supplied.

## 1. Algebra and quaternionic geometry

### S01 — Hamilton quaternions
H has basis 1,i,j,k, relations i²=j²=k²=ijk=−1, conjugation q↦q̄ and norm N(q)=q q̄. Multiplicativity N(pq)=N(p)N(q) implies unit quaternions form Sp(1)≅SU(2). **Status: PROVED.** Dependency: none.

### S02 — Double cover SU(2)→SO(3)
The action ρ(q)v=qvq⁻¹ on imaginary quaternions is a surjective homomorphism with kernel {±1}; hence SU(2)/{±1}≅SO(3). **Status: PROVED.** Dependency: S01. q and −q induce the same rotation but are distinct spin-group elements.

### S03 — Binary icosahedral group 2I
Define 2I=ρ⁻¹(I), where I⊂SO(3) is the rotational icosahedral group of order 60. Then |2I|=120 and the restricted covering has kernel {±1}. **Status: group structure PROVED; explicit quaternion list/closure audit NUMERICAL or FINITE-EXACT according to arithmetic used.** Dependencies: S01–S02. Floating-point closure checks are not formal certificates.

### S04 — Intertwiners and Schur
For representations D₁,D₂ of G, T is an intertwiner iff T D₁(g)=D₂(g)T for every g. If both complex representations are irreducible, Hom_G(V₁,V₂) has dimension zero or one; a nonzero intertwiner is an isomorphism unique up to scalar. **Status: PROVED under irreducibility hypotheses.** Each application must separately verify irreducibility.

## 2. Q6 combinatorics, topology and Hodge theory

### S05 — Six-cube graph
V={0,1}⁶; adjacency is Hamming distance one. Thus |V|=64, degree 6, |E|=64·6/2=192; parity of coordinate sum gives a 32+32 bipartition. **Status: PROVED.**

### S06 — Graph Laplacian spectrum
For L=6I−A, Walsh characters χ_S(x)=(-1)^{Σ_{i∈S}x_i} satisfy Lχ_S=2|S|χ_S. Spectrum: 2j with multiplicity C(6,j), j=0,…,6. **Status: PROVED.** Dependency: S05.

### S07 — Oriented square complex
There are C(6,2)·2⁴=240 elementary square faces. With oriented incidence maps B1:C1→C0 and B2:C2→C1, every square has zero boundary-of-boundary, B1B2=0. **Status: PROVED; implementation FINITE-EXACT.** Dependency: S05.

### S08 — Graph cycle rank versus filled complex homology
For the connected graph, dim Z1=|E|−|V|+1=129. For the square complex, the 240 face boundaries have rank 129, giving Betti numbers (β0,β1,β2)=(1,0,111). **Status: PROVED when the integer rank certificate is retained; otherwise FINITE-EXACT audit pending.** Dependencies: S05–S07. The number 129 is the graph cycle-space rank, not β1 of the filled complex.

### S09 — Weighted combinatorial Hodge Laplacians
Choose positive-definite inner products Wk on cochains. Define coboundaries from Bk and adjoints using Wk. The finite-dimensional Hodge theorem identifies harmonic k-cochains with H^k; dimensions are (1,0,111). **Status: PROVED for the declared positive weights; finite matrices FINITE-EXACT.** Dependencies: S07–S08. These weights do not constitute a metric DEC Hodge star unless dual-cell geometry is supplied.

### S10 — Antipodal involution
a(x)=x+(1,…,1) mod 2. Its induced permutations Pk must satisfy B1P1=P0B1 and B2P2=P1B2, with signs fixed by cell orientations. **Status: FINITE-EXACT audit.** Dependencies: S05–S07. Preserve the orientation convention and chain-map checks.

## 3. RGS seed, metric and finite elements

### S11 — Recovered RGS seed
A finite map X:V(Q6)→R³, with fixed state ordering, source hash, units and transformation provenance. **Status: REPRODUCIBLE NUMERICAL reconstruction.** Dependency: S05. The recovered seed is not by itself a unique physical derivation.

### S12 — Octahedral symmetry Oh
Oh consists of 48 signed-permutation 3×3 matrices. Check X(g·v)=gX(v), the group law and induced equivariance on vertices, edges and faces; boundary maps must commute with the action. **Status: FINITE-EXACT for group/permutation data; NUMERICAL for floating coordinates.** Dependency: S11.

### S13 — Edge and face metrics
For e=(i,j), ℓ_e=||X_j−X_i||₂. Face geometry is encoded by side-vector Gram matrices, areas and closure residuals. **Status: NUMERICAL.** Dependencies: S11–S12. Report degeneracies, precision and tolerances.

### S14 — Rank of the induced pullback metric
For DX:R⁶→R³, G=(DX)ᵀDX is positive semidefinite and rank(G)≤3. **Status: PROVED by linear algebra.** Consequence: it is not a nondegenerate six-dimensional metric.

### S15 — Graph P1 FEM
On an edge of length ℓ>0, the one-dimensional P1 stiffness is ℓ⁻¹[[1,−1],[-1,1]]. The mass matrix depends on the declared measure/quadrature. **Status: PROVED for the stated 1D element; assembled matrices FINITE-EXACT or NUMERICAL.** Dependencies: S11–S13. It is not automatically a 3D FEM or Whitney/DEC Hodge star.

### S16 — Ellis–RGS curved mesh family
M_L=[−L,L]×S² with g=du²+(u²+r₀²)dΩ², together with the documented curved tetrahedral refinement family. **Status: FINITE CONSTRUCTION.** Dependencies: S11–S13. Finite levels do not establish uniform shape regularity.

### S17 — Flat periodic P1 Mosco theorem
For the specified regular periodic Freudenthal family, common Hilbert space, reconstruction, coercivity and compactness hypotheses, prove weak liminf and strong recovery for the discrete Dirichlet forms. **Status: PROVED under the flat benchmark hypotheses.** Reference: docs/FLAT_TORUS_G6_P1_MOSCO_PROOF.md. This does not imply the curved result.

### S18 — Curved Ellis–RGS Mosco problem
In H=L²(M_L,dμ_g), q(u)=∫|∇_g u|²dμ_g, and q_h the curved P1 forms extended through J_h by +∞ outside the discrete range, prove (i) weak liminf and (ii) strong recovery in the same H. **Status: OPEN.** Dependencies: S16 and the curved common-L² reconstruction. Obligations: uniform element-map/inverse bounds, Jacobian and metric bounds, coercivity, compactness, liminf and recovery. S17 is only a methodological benchmark.

## 4. Qutrit, SU(3) and switching networks

### S19 — Aperiodic qutrit chain
Finite Hilbert space (C³)^⊗n, explicit SU(3)-invariant interaction and Fibonacci bond sequence; specify boundary conditions, symmetry sector and diagonalization. Schmidt rank and entropy are state-dependent finite observables. **Status: FINITE MODEL + NUMERICAL diagonalization.** Not a QCD derivation.

### S20 — Gauge-covariant qutrit–Q6 embedding
Place the finite chain on a declared Gray path in Q6 and assign Uij∈SU(3). Under local gi, Uij↦giUijgj⁻¹; require H[U^g]=G H[U]G⁻¹. **Status: FINITE-EXACT/NUMERICAL according to implementation.** Dependencies: S05,S19. Not dynamical Yang–Mills.

### S21 — Q3, Omega and Beneš networks
Q3={0,1}³ encodes trigrams; Omega/Beneš are explicit compositions of 2×2 switching stages. Specify wiring, stage order and permutation convention, then enumerate reachable permutations. **Status: FINITE-EXACT if exhaustive enumeration and wiring are archived.**

### S22 — Trigram–hexagram encoding
An hexagram is an ordered pair of trigrams, hence a six-bit word. A selected set of eight input/output hexagrams is a permutation iff every input and output occurs exactly once. **Status: PROVED as a finite encoding under bijectivity; selected instance FINITE-EXACT.** Dependencies: S05,S21. No identification with PRISM follows.

## 5. Kathara, Krystal Spiral and refinement

### S23 — Kathara graph transcription
Declared source diagram gives 12 centres, 15 named primary lines and 21 elementary edges after expanding polylines. **Status: FINITE MODEL; counts/planarity are conditional on the transcription and embedding.**

### S24 — Kathara cellular chain complex
Build C2→C1→C0 with integer B2,B1; verify B1B2=0 and compute ranks/Betti numbers for the declared transcription. **Status: FINITE-EXACT.** Dependency: S23.

### S25 — Krystal Spiral square similarity
T=[[1,1],[-1,1]] satisfies TᵀT=2I. If the source geometry imposes clockwise 45° rotation, √2 scaling and the stated centre overlays, then p9=Tp6 and p12=T²p6. **Status: PROVED CONDITIONAL.** This is not a tribonacci or physical derivation.

### S26 — Kathara–Q6 chain-map gates
A cellular map is a family Fk:Ck(K)→Ck(Q6) with Bk(Q6)Fk=F(k−1)Bk(K). Gate C's vertex-injective candidate collapses source 2-cells and is not a cellular embedding. Gate D supplies an alternative path/filling construction under its stated ordering. **Status: CANDIDATE/FINITE-EXACT for tested matrix properties; metric canonicity OPEN.** Dependencies: S07,S23–S24.

### S27 — Barycentric triangle quadrature
For barycentric nodes λ(a) and weights wa, define mean m=Σwaλ(a) and covariance C=Σwa(λ(a)−m)(λ(a)−m)ᵀ, with Σwa=1. Compare exact first/second moments with the rule's moments. **Status: PROVED for the stated moment identities; implementation audit NUMERICAL/FINITE-EXACT.** Not the 6D covariance of Q6.

### S28 — Dyadic Q6 2-skeleton subdivision
For integer subdivision maps S0,S1,S2 require B'1S1=S0B1 and B'2S2=S1B2. The declared refined counts are 496 vertices, 1344 edges and 960 squares. **Status: FINITE-EXACT.** Reference benchmark only; not Kathara or a tribonacci derivation.

## 6. Q6–PRISM obstructions and A5 alternative

### S29 — Five-factor directional embedding obstruction
Under the restricted definition that the six independent Q6 coordinate directions must map to six distinct Cartesian product factors, a five-factor product cannot host such an edge-preserving embedding; six factors suffice for the standard binary construction. **Status: PROVED under that definition.** Not a no-go for arbitrary graph maps or path-valued maps.

### S30 — Cubical folding Q6→Q5
F(v0,…,v5)=(v0⊕v5,v1,…,v4), extended to cells with orientation signs. Verify ∂F=F∂ on the complete 729-cell source list. The image is a contractible Q5-type cubical region, not the periodic PRISM torus. **Status: FINITE-EXACT.**

### S31 — Periodic winding obstruction
For edge-to-edge maps compatible with square faces, a four-edge face boundary cannot create nonzero winding in a cyclic factor of length ≥5. Thus the induced map on H1 into the specified PRISM periodic target is zero. **Status: PROVED under edge-to-edge and face-compatibility hypotheses.** Does not apply to path-valued maps.

### S32 — Path-valued homology map
Allow each source edge to map to a target path. Given an integral cycle basis of rank 129, construct a surjection Z¹²⁹→Z⁵ with matrix [I₅|0]. **Status: PROVED as an existence construction depending on the chosen basis and path freedom.** It does not extend automatically to the 240 filled faces.

### S33 — A6 equivariance obstruction
Under the specified A6 source action and target action inherited from the original PRISM, the induced map on H1 is zero by the stated module decomposition. **Status: PROVED CONDITIONAL; retain the exact representations and target action.** Not a no-go for a different target.

### S34 — Exterior-square target
Construct an integral S6-equivariant surjection Z1(Q6 graph;Z)→Λ²Z⁶, with rank-114 kernel (source cycle rank 129, target rank 15). The original PRISM Z⁵ component is not the image of this map. **Status: FINITE-EXACT/PROVED when formula, integrality and equivariance are documented.**

### S35 — GEC-6 root-Cayley model on A5
A5={x∈Z⁶:Σxi=0}; its 30 roots ei−ej generate a 30-regular Cayley graph on A5/mA5≅(Z/mZ)⁵. Coordinate permutations give a faithful S6 action. **Status: PROVED algebraically; finite quotient FINITE-EXACT.** This is a new target, not the original PRISM.

### S36 — GEC-6 spectrum for m=3
For the stated Laplacian normalization on 243 vertices, Fourier characters give 0^(1),15^(12),24^(30),27^(50),33^(120),36^(30). Independent 243×243 diagonalization is a numerical cross-check. **Status: PROVED if Fourier character classification is complete; cross-check NUMERICAL.** Dependency: S35.

### S37 — General-modulus spectral gap
γm=(1/15) min_{k≠0} Σ_{i<j}(1−cos(2π⟨k,rij⟩/m)). The claimed limit is lim m²γm=2π²/3. **Status: PROVED only with a global minimum certificate and uniform remainder estimate; otherwise OPEN/NUMERICAL.** Dependencies: S35–S36. The gap closes quadratically; no uniform expander family.

## 7. Continuum A5, pentagon and De Aquino threshold

### S38 — Flat A5 torus spectrum
For the documented lattice metric and dual lattice, the minimum nonzero quadratic form is 5 with multiplicity 12, λ1=10π²/3, and lim m²γm=λ1/5. **Status: PROVED only under the declared metric/scale normalization; otherwise CONDITIONAL.** Dependencies: S35–S37.

### S39 — Algebraic De Aquino threshold
For mg/mi0=3−2√(1+ξ²), the zero occurs at ξ²=5/4. The proposed bridge to the A5 minimum 5 uses ξk²=kᵀQk/4. **Threshold: PROVED algebraically. Bridge: CANDIDATE/OPEN** until factor 1/4, physical momentum map and dynamical intertwiner are derived. Dependency: S38.

### S40 — Pentagon frame normalization
For the declared regular pentagon frame, EEᵀ=(5/2)I₂ and G5=EᵀE=(5/2)P. Under the stated balanced five-amplitude unit-vector ansatz, |Pq|²=1/2 and qᵀG5q=5/4. **Status: PROVED CONDITIONAL.** It does not show the physical dynamics selects the ansatz or derive p0=mi0c. Dependency: S39.

## 8. Gradients, shear and Green–Kubo

### S41 — Shear derivative of implemented Q6 action
For coframe E(γ), define Πxy=−V⁻¹ dS/dγ|₀. Differentiate the actual implemented action and compare with an independent central finite difference. Under the stated decomposition (Wilson term depends only on holonomies; hopping/twin/quartic terms only on coframe volume), the traceless linear shear response is confined to the fixed-frame regulator. **Status: PROVED for the declared action derivative; finite-difference audit NUMERICAL.** Not automatically the complete physical stress tensor.

### S42 — Inverse-square link-weight shear
For deformed edge vector d_e(γ)=(dx+γdy,dy,…), w_e=|d_e|⁻² gives w′e(0)=−2dxdy/|d_e|⁴. For the separately declared Hlink=Σe weEe, Πxy=−V⁻¹Σe w′eEe. **Status: PROVED CONDITIONAL** on this deformation and energy. Not an automatic derivative of the Wilson plaquette action.

### S43 — Metric-weighted plaquette response
Build face areas from RGS face Gram matrices and define a conditional face-Hodge action. Jacobi's determinant formula gives area derivatives when the Gram matrix is nonsingular; compare with finite differences. **Status: CONDITIONAL + NUMERICAL.** A complete DEC star requires compatible dual measures.

### S44 — Two-sheet covariant matter energy
Sm=(κ/2)Σe |Deψ|²/Le, with covariant edge difference and positive edge length. The coefficient 1/Le is the 1D P1 stiffness convention. **Status: PROVED for the discrete form; model choice CONDITIONAL; derivative audit NUMERICAL.** Not a Dirac action without spinor, Clifford, measure and limit structures.

### S45 — Gibbs stress benchmark
Single-link symmetric proposals with Metropolis acceptance satisfy detailed balance for the declared Gibbs measure. Finite stress autocorrelations are Monte Carlo estimates. **Detailed balance: PROVED; estimates: NUMERICAL.** Monte Carlo index is algorithmic time, not physical time; the integral is not yet a viscosity.

### S46 — Finite OU Green–Kubo integral
For stable OU dX=−AXdt+BdW with stationary covariance Σ solving AΣ+ΣAᵀ=BBᵀ, Cov(Xt,X0)=e^(−At)Σ. For Y=cᵀX, ∫₀∞Cov(Yt,Y0)dt=cᵀA⁻¹Σc when A is stable. Quadratic observables require Wick's formula and a separate integral. **Status: PROVED under stated stability assumptions; implementation FINITE-EXACT/NUMERICAL.** Physical viscosity additionally requires calibrated physical time, stress identification, FDT, equilibrium measure and appropriate limits.

## 9. Pellis: spectral geometry, two sheets and channel symmetry

### S47 — Scalar spectrum on S³/2I
For the round unit S³, scalar Laplacian eigenvalues are ℓ(ℓ+2); on the quotient retain 2I-invariant subspaces, with multiplicities computed by the Reynolds projector/Molien series for the specified action. **Status: PROVED with the dossier's action conventions; finite multiplicities FINITE-EXACT.** Golden Spheres Resonance weights do not follow automatically.

### S48 — Two-sheet self-adjoint operators
Use Hilbert space H+⊕H−, specify dense domain and block operator. The SU(2) invariant epsilon=iσ₂=[[0,1],[-1,0]] yields the fundamental intertwiner, unique up to phase by Schur. Distinguish complex-linear c:V→V̄ from antilinear coordinate conjugation. Self-adjointness additionally requires the correct adjoint block relations and invariant domain. **Status: SU(2) intertwiner PROVED; edge exchange CONDITIONAL on sheet-compatible weights; operator self-adjointness must be checked per domain/operator.**

### S49 — D6 channel intertwiner
For the declared six-state factorization C⁶≅C²_sheet⊗C³_mode, the three-mode representation decomposes as 1⊕2. Within the explicitly sheet-blind ansatz the commutant allows independent radial blocks M_parallel(E), M_perp(E); D6 alone does not force equality. **Status: factorization/commutant FINITE-EXACT; dipole ℓ=1 identification CANDIDATE; physical GSR selection OPEN.**

### S50 — RGS–PSU–PRISM–Pellis bridge discipline
Each proposed bridge must specify (Hs,Ht,T,P): source and target spaces, map T and required property P (intertwining, isometry, chain map, convergence, etc.). Keep Q6(64), RGS operator(188D), PRISM(69,120), Q-ball BVP and Pellis/GSR resonance operator separate. **Status: separation is a formal methodological constraint; physical maps remain OPEN.** Similarity of dimensions/symmetries is insufficient.

### S51 — D6 splitting and possible triplet degeneracy
D6 permits 1⊕2 and therefore two blocks in the sheet-blind commutant. If a larger group acts irreducibly on the complete ℓ=1 triplet, and the operator, domain and boundary conditions are all equivariant, Schur forces a scalar angular block. **Status: first statement PROVED; SO(3)-based degeneracy CONDITIONAL; its applicability to Pellis/GSR OPEN.** Do not infer SO(3) from a triangular arrangement alone.

## 10. Dependency graph

- S01→S02→S03; S04 supplies the general intertwiner framework.
- S05→S06; S05→S07→S08→S09; S05–S07→S10.
- S11→S12→S13→S15; S11–S13→S16→S18. S17 is an independent flat benchmark, not a proof premise sufficient for S18.
- S05+S19→S20; S05+S21→S22.
- S23→S24→S26; S25 is a conditional geometric branch; S07→S28.
- S05/S08→S29–S34; S35→S36→S37→S38→S39→S40, subject to the qualifications in each record.
- S13→S41–S44; S45 is a separate Monte Carlo model; S46 is an OU theorem, not a consequence of S45.
- S03→S47; S04→S48–S49; S49–S50→S51.

## 11. Completion criteria

1. Every record must link to a canonical mathematical specification and the actual implementation/tests.
2. FINITE-EXACT results retain integer matrices, orientation conventions, ranks and reproducible certificates.
3. NUMERICAL results record commit SHA, input, algorithm, tolerance, residual and whether the check is independent.
4. CONDITIONAL claims put every hypothesis in the statement, not only in a footnote.
5. Every inter-model bridge gives T:Hs→Ht, the commutative diagram and the exact property to prove.
6. Continuum claims require a refinement family, common space, forms/operators, uniform hypotheses and convergence proof.
7. Physical claims require observables, units, dynamics, equilibrium measure and limit procedure.

## 12. Existing canonical specifications

- Q6/Hodge: Q6_COMPLEX_FORMALIZATION.md; Q6_COMBINATORIAL_HODGE_LAPLACIANS.md; Q6_DYADIC_SUBDIVISION_CHAIN_MAP.md.
- Seed/metric/FEM: RGS_Q6_INDUCED_METRIC.md; Q6_BARYCENTRIC_P1_FEM.md; GATE_E1_BARYCENTRIC_COVARIANCE_QUADRATURE.md.
- Mosco: FLAT_TORUS_G6_P1_MOSCO_PROOF.md; RGS_MOSCO_CURVED_REFINEMENT_PROGRAM.md.
- Qutrit/network: QUTRIT_APERIODIC_ENTANGLEMENT.md; QUTRIT_Q6_GAUGE_COVARIANT_EMBEDDING.md; research/CULLER_INTERCONNECTION_NETWORKS_Q3_Q6_OMEGA_BENES.md.
- Kathara: KATHARA_15_LINE_GRAPH_FACE_AUDIT.md; KATHARA_CELLULAR_COMPLEX_AUDIT.md; KATHARA_KRYSTAL_SPIRAL_SQUARE_GEOMETRY.md; KATHARA_Q6_CHAIN_MAP_GATE_C.md; KATHARA_Q6_NONDEGENERATE_GATE_D.md.
- Q6–PRISM/GEC-6: research/Q6_PRISM_EMBEDDING_OBSTRUCTION.md; research/Q6_PRISM_CUBICAL_CHAIN_MAP.md; research/Q6_PRISM_PERIODIC_WINDING_OBSTRUCTION.md; research/Q6_PRISM_PATH_VALUED_HOMOLOGY.md; research/Q6_PRISM_SYMMETRY_OBSTRUCTION.md; research/Q6_PRISM_EXTERIOR_SQUARE_EQUIVARIANT_MAP.md; research/GEC6_A5_LAPLACIAN_SPECTRUM_M3.md; research/GEC6_A5_CONTINUUM_TORUS_SPECTRUM_AND_DEAQUINO_BRIDGE.md.
- Shear/Green–Kubo: Q6_ACTION_DERIVED_SHEAR_RESPONSE.md; Q6_METRIC_WILSON_RESPONSE.md; Q6_GIBBS_STRESS_CORRELATIONS.md; Q6_HYPERCUBE_GREEN_KUBO_SUPERGOLDEN.md.
- Pellis: PELLIS_A_2I_SPECTRAL_GEOMETRY.md; PELLIS_B_TWO_SHEET_SELF_ADJOINT_OPERATORS.md; PELLIS_C_D6_CHANNEL_INTERTWINER.md; PELLIS_D_RGS_PSU_PRISM_INTEGRATION.md.

**Maintenance note:** this register formalizes the 51-step scope and its proof obligations. It does not fabricate missing tests, replace specialist proofs, or promote open physical identifications to established results.

## Evidence matrix

The implementation/test/result traceability audit for all 51 records is maintained in [51-step evidence and implementation matrix](51_STEP_EVIDENCE_MATRIX_2026-11-26.md). File presence is not proof; test execution and theorem status must be assessed separately.
