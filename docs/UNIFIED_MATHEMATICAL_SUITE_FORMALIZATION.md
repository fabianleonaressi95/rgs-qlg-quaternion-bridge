# Unified mathematical suite: formal statements, evidence and dependencies

**Status:** canonical suite index and proof-obligation register  
**Scope:** the finite mathematical and computational constructions linked by the RGS–QLG Quaternion Bridge.  
**Companion programme:** [THEOREM](https://github.com/fabianleonaressi95/theorem).

This document replaces a narrative sequence of experiments with a typed mathematical dependency register. It does **not** promote every construction to a theorem. Each record distinguishes exact mathematics, finite computational evidence, conditional models and unresolved physical identifications. The implementation and audit files remain the source of reproducibility; this document states what they establish and what they do not.

## 1. Evidence vocabulary

| Label | Meaning | Minimum evidence |
|---|---|---|
| **PROVED** | A mathematical statement follows from stated hypotheses. | Definitions, hypotheses and proof or a linked theorem record. |
| **FINITE-VERIFIED** | A finite object or algorithm has been checked computationally. | Versioned implementation, test/audit and tolerances where numerical. |
| **CONDITIONAL** | The result holds for an explicitly declared model or ansatz. | Model definition and assumptions stated in the record. |
| **CANDIDATE** | A construction is explicit but a required property or canonical selection is unproved. | Construction plus precise missing gate. |
| **OPEN** | No proof or adequate verification is yet available. | A falsifiable question and closure criterion. |

A test passing is evidence for the tested finite statement only. It does not prove a continuum limit, a physical interpretation, or equivalence between distinct graphs, Hilbert spaces or actions.

## 2. Core definitions and conventions

### 2.1 Quaternion and spin conventions

Use scalar-first Hamilton quaternions
\[
q=a+bi+cj+dk,\qquad i^2=j^2=k^2=ijk=-1,
\]
with conjugation \(\bar q=a-bi-cj-dk\) and norm \(|q|^2=q\bar q\). Unit quaternions form \(SU(2)\). Identify \(\mathbb R^3\) with imaginary quaternions and define
\[
\rho(q)v=qvq^{-1}.
\]
Then \(\rho:SU(2)\to SO(3)\) is a surjective homomorphism with kernel \(\{1,-1\}\). Thus it is a two-fold covering; \(q\) and \(-q\) are distinct spin elements with the same spatial rotation. The binary icosahedral group \(2I\) is the inverse image of the rotational icosahedral group \(I\), so \(|2I|=120\) and \(|I|=60\).

For representations \(D_1,D_2\) of a group \(G\), an intertwiner is a linear map \(T\) satisfying \(TD_1(g)=D_2(g)T\) for every \(g\in G\). Schur's lemma gives scalar uniqueness only when the relevant complex representations are irreducible and equivalent (and the field and representation category are specified).

### 2.2 Q6 graph and square complex

Let \(Q_6\) have vertex set \(\{0,1\}^6\), with an edge between words differing in exactly one coordinate. Hence
\[
|V|=2^6=64,\quad |E|=6\,2^5=192,\quad \deg(v)=6.
\]
The graph is bipartite with 32 vertices in each class. Its graph Laplacian is \(L=6I-A\), with eigenvalues \(2j\) of multiplicity \(\binom 6j\), for \(j=0,\ldots,6\).

Attach the 240 coordinate squares to obtain the canonical cubical 2-skeleton. With consistently oriented cellular boundary matrices
\[
C_2\xrightarrow{B_2}C_1\xrightarrow{B_1}C_0,
\qquad B_1B_2=0.
\]
The graph cycle rank is \(192-64+1=129\). It is **not** the first Betti number of the filled square complex: the latter has \((\beta_0,\beta_1,\beta_2)=(1,0,111)\). The antipodal bit-complement map is \(v\mapsto v+(1,1,1,1,1,1)\pmod 2\); induced actions on oriented cells must include orientation signs.

### 2.3 RGS seed, metric and finite elements

The recovered seed is a map \(X:V(Q_6)\to\mathbb R^3\), reconstructed only through its hash-checked source routine. For an edge \(e=(i,j)\), define \(\ell_e=\|X_j-X_i\|\). A vertex pullback Gram matrix \((DX_i)^TDX_i\), when the differential maps into \(\mathbb R^3\), has rank at most three; it is not a nondegenerate six-dimensional metric. The graph P1 element stiffness on an edge of length \(\ell_e\) uses the one-dimensional coefficient \(1/\ell_e\).

The coordinate signed-permutation group \(O_h\) has order 48. Its action on the reconstructed point set and its induced action on Q6 labels are separate properties; neither should be inferred from the other.

### 2.4 Curved Ellis–RGS Mosco problem

The curved scalar problem uses \(M_L=[-L,L]\times S^2\), \(g=du^2+(u^2+r_0^2)d\Omega^2\), and the common Hilbert space \(H=L^2(M_L,d\mu_g)\). For reconstruction maps \(J_h:V_h\to H\), extend each discrete quadratic form by \(+\infty\) outside \(J_h(V_h)\). Mosco convergence requires both the weak liminf inequality and a strong recovery sequence in this same Hilbert space. A flat periodic P1 theorem does not imply the curved statement.

## 3. Canonical register: the 51 mathematical steps

The numbering preserves the suite crosswalk. The evidence label describes the claim at the stated scope, not a stronger physical interpretation.

| ID | Object / statement | Evidence | Dependency or explicit boundary |
|---:|---|---|---|
| 01 | Hamilton quaternion algebra, conjugation, norm and multiplication convention. | PROVED | Definition; fixes all subsequent spin conventions. |
| 02 | Double cover \(SU(2)\to SO(3)\), kernel \(\{\pm1\}\). | PROVED | 01; conjugation action on imaginary quaternions. |
| 03 | Binary icosahedral group \(2I\), order 120, over \(I\) of order 60. | PROVED / FINITE-VERIFIED | 02; exact group construction and exhaustive audit. |
| 04 | Representation decomposition and intertwiner equations; Schur uniqueness under irreducibility hypotheses. | PROVED | 03 or other declared group representations; hypotheses must be explicit. |
| 05 | Q6 graph counts: 64 vertices, 192 edges, degree 6, bipartition 32+32. | PROVED | Definition in §2.2. |
| 06 | Q6 graph-Laplacian spectrum \(2j\), multiplicity \(\binom6j\). | PROVED | 05; product-graph character basis. |
| 07 | Oriented square complex with 240 faces and \(B_1B_2=0\). | PROVED / FINITE-VERIFIED | 05; orientation convention and incidence constructor. |
| 08 | Graph cycle rank 129; filled-complex Betti numbers \((1,0,111)\). | PROVED / FINITE-VERIFIED | 07; do not conflate graph cycles with cellular homology. |
| 09 | Weighted combinatorial Hodge Laplacians and harmonic dimensions \((1,0,111)\). | PROVED under declared inner products / FINITE-VERIFIED | 07–08; weights are cochain inner products, not automatically metric DEC stars. |
| 10 | Antipodal bit-complement action on vertices, edges and faces. | PROVED / FINITE-VERIFIED | 05–07; orientation signs are part of the induced action. |
| 11 | Hash-verified reconstruction of the 64-point RGS seed. | FINITE-VERIFIED | Source routine, immutable hash and provenance. |
| 12 | Faithful order-48 signed-coordinate-permutation action \(O_h\). | PROVED / FINITE-VERIFIED | 11; coordinate action and Q6-label equivariance audited separately. |
| 13 | Edge lengths, face Gram matrices, areas and closure defects from the seed. | FINITE-VERIFIED | 11–12; finite geometric audit, not a continuum metric theorem. |
| 14 | Pullback Gram rank at most three for a map into \(\mathbb R^3\). | PROVED | Linear-algebraic rank bound; excludes a nondegenerate six-dimensional pullback. |
| 15 | Graph P1 mass and stiffness with edge stiffness coefficient \(1/\ell_e\). | PROVED for the declared 1D FEM model / FINITE-VERIFIED | 13–14; not a Whitney-form Hodge star. |
| 16 | Finite curved Ellis–RGS tetrahedral refinement family. | FINITE-VERIFIED | Recovered seed, exact embedding and mesh-generation conventions. |
| 17 | Mosco liminf and recovery for the regular flat periodic P1 benchmark. | PROVED under its stated flat hypotheses | Flat common-\(L^2\) reconstruction and uniform estimates; no curved transfer. |
| 18 | Mosco convergence for the curved Ellis–RGS family. | OPEN | 16 plus common-space reconstruction; needs all-level shape regularity, metric/Jacobian bounds, coercivity, compactness, liminf and recovery. |
| 19 | Fibonacci-modulated finite SU(3)-invariant qutrit chain and Schmidt diagnostics. | FINITE-VERIFIED | Finite-dimensional Hamiltonian; not QCD or continuum gauge dynamics. |
| 20 | Qutrit chain on a Gray path of Q6 with SU(3) link transporters and gauge covariance. | FINITE-VERIFIED | 05, 19; finite covariance only, not dynamical Yang–Mills. |
| 21 | Q3–Omega–Beneš routing correspondence; 4,096 Omega-8 and 40,320 Beneš-8 permutations under the declared convention. | FINITE-VERIFIED | Switch-stage convention and exhaustive enumerator; not a physical equivalence. |
| 22 | Trigram pair as a six-bit esagram/Q6 word; eight input-output demands form a permutation iff each input and output occurs once. | PROVED | Combinatorial encoding; does not identify Q6 with PRISM. |
| 23 | Kathara transcription: 12 centres, 15 named primary lines, 21 elementary edges; planar embedding with 10 bounded faces. | FINITE-VERIFIED / PROVED for declared embedding | Source transcription and chosen straight-line embedding; not a unique source-independent geometry. |
| 24 | Kathara cellular chain complex and Hodge ranks. | FINITE-VERIFIED | 23; exact incidence matrices and \(B_1B_2=0\). |
| 25 | Krystal Spiral square similarity \(T^TT=2I\), centre constraints \(p_9=Tp_6\), \(p_{12}=T^2p_6\). | PROVED conditionally | Declared Euclidean source geometry and orientation convention. |
| 26 | Kathara–Q6 Gate C obstruction and Gate D path/filling chain-map construction. | FINITE-VERIFIED / CANDIDATE | Gate C collapses 2-cells; Gate D has full-column-rank face map, but metric canonicality remains open. |
| 27 | Triangle barycentric moments, covariance and quadrature identities. | PROVED for stated polynomial degree / FINITE-VERIFIED | Local planar simplex coordinates; distinct from Q6 six-dimensional covariance. |
| 28 | Dyadic Q6 2-skeleton subdivision with 496 vertices, 1,344 edges, 960 squares and commuting subdivision maps. | PROVED / FINITE-VERIFIED | 07; exact integer chain-map identities. Not Kathara or tribonacci transfer. |
| 29 | No injective edge-preserving Q6 embedding into the existing five-factor Cartesian PRISM. | PROVED | Six independent Q6 directions require six distinct Cartesian factor directions. |
| 30 | Folding chain map \((v_0,\ldots,v_5)\mapsto(v_0\oplus v_5,v_1,\ldots,v_4)\), with \(\partial F=F\partial\). | PROVED / FINITE-VERIFIED | Cubical orientation signs; finite check on all 729 source cells. Image is a contractible Q5 corner, not the PRISM torus. |
| 31 | Edge-to-edge Q6 graph maps induce zero periodic winding in PRISM factors of length at least five. | PROVED | Q6 graph cycle lattice generated by four-edge squares; four-step closed walks cannot wind. |
| 32 | Path-valued graph map inducing \(\mathbb Z^{129}\to\mathbb Z^5\), represented by \([I_5\mid0]\). | PROVED as an existence construction / FINITE-VERIFIED | Choice of fundamental-cycle basis; graph only, not extension over filled 2-cells; not canonical. |
| 33 | Full \(A_6\)-equivariance obstructs the preceding surjection for the original PRISM target action. | PROVED | Representation action and target automorphism hypotheses as stated in the obstruction note. |
| 34 | Surjective \(S_6\)-equivariant graph-cycle map to \(\Lambda^2\mathbb Z^6\), rank 15, kernel rank 114. | PROVED / FINITE-VERIFIED | Graph cycle module; target is an enriched 15-torus, not the original PRISM \(\mathbb Z^5\). |
| 35 | GEC-6 periodic model \(A_5/mA_5\), generated by 30 roots with faithful coordinate-permutation symmetry. | PROVED / FINITE-VERIFIED | New A5-root Cayley model; not the Cartesian PRISM. |
| 36 | GEC-6 spectrum for \(m=3\): \(0^{(1)},15^{(12)},24^{(30)},27^{(50)},33^{(120)},36^{(30)}\); normalized gap \(1/2\). | PROVED by finite Fourier analysis / FINITE-VERIFIED | Independent 243-by-243 Laplacian diagonalization. |
| 37 | General-modulus gap asymptotic \(\lim_{m\to\infty}m^2\gamma_m=2\pi^2/3\). | PROVED under the stated A5-root normalization | Exact shortest-vector certificate and Taylor expansion; gap closes, so no uniform expander family. |
| 38 | Flat A5 torus first Laplacian eigenvalue \(10\pi^2/3\), multiplicity 12, and finite-gap scaling relation. | PROVED | Dual-lattice spectrum; \(\lim m^2\gamma_m=\lambda_1/5\). |
| 39 | De Aquino zero threshold \(\xi^2=5/4\); matching to A5 shortest value 5 needs \(\xi_k^2=k^TQk/4\). | PROVED algebraically / CANDIDATE bridge | The factor 1/4, physical momentum map and dynamical intertwiner are not derived. |
| 40 | Pentagon tight-frame identities \(EE^T=\frac52I_2\), \(G_5=\frac52P\), and balanced ansatz value \(q^TG_5q=5/4\). | PROVED for the stated frame and ansatz | Does not show De Aquino dynamics selects the ansatz; physical scale and S6 map remain open. |
| 41 | Analytic shear derivative of the implemented finite Q6 action. | PROVED for the declared action / FINITE-VERIFIED | Wilson holonomy term and volume-only sectors have zero first variation under traceless shear; remaining term is the fixed-frame regulator. Not full microscopic stress. |
| 42 | Inverse-square edge weight derivative and associated link-energy shear response. | PROVED conditionally | Requires separately declared \(H_{link}=\sum_e w_eE_e\); not implied by plaquette Wilson action. |
| 43 | Metric face-area-weighted plaquette shear derivative. | PROVED conditionally / FINITE-VERIFIED | Jacobi derivative for declared face-Hodge model; no circumcentric dual measure, hence not a complete DEC star. |
| 44 | Two-sheet graph matter energy with coefficient \(1/L_e\) and its shear derivative. | PROVED for declared graph model / FINITE-VERIFIED | Covariant graph Dirichlet energy; not a Dirac action. |
| 45 | SU(2) single-link Metropolis Gibbs benchmark and finite stress autocorrelation. | FINITE-VERIFIED | Detailed balance for declared proposal and action; Monte Carlo index is algorithmic time, not physical time. |
| 46 | Exact Green–Kubo integral for a finite OU process. | PROVED for specified OU generator, covariance and observable | Does not yield physical viscosity without calibrated dynamics, FDT, full equilibrium action and continuum/refinement control. |
| 47 | Scalar Laplacian and invariant spectral data on \(S^3/2I\). | PROVED | Molien invariants and quotient multiplicities; does not derive Golden Spheres Resonance weights. |
| 48 | Finite two-sheet self-adjoint operators and SU(2)-equivariant fibre map \(\varepsilon=i\sigma_2\). | PROVED for declared finite spaces / FINITE-VERIFIED | Distinguish complex-linear map to conjugate representation from antilinear coordinate conjugation; field exchange needs compatible edge weights. |
| 49 | Six-state Fibonacci orbit factorization \(\mathbb C^6\simeq\mathbb C^2\otimes\mathbb C^3\) and \(D_6\) intertwiner. | PROVED as representation statement | \(\mathbb C^3=\mathbf1\oplus\mathbf2\); dipole \(\ell=1\) realization is a candidate, not a selected physical sector. |
| 50 | RGS–PSU–PRISM–Pellis integration map and operator correspondence. | OPEN / CANDIDATE | Q6(64), RGS(188D), PRISM(69,120), Q-ball and GSR spaces require explicit maps and compatible operators. |
| 51 | Conditional \(M(E)\otimes I_6\) lift and radial degeneracy. | PROVED conditionally | \(D_6\) alone permits independent parallel and transverse blocks; equality requires larger symmetry or a proved dynamical selection. |

## 4. Dependency structure

The suite is a directed graph, not a chain of automatic implications.

1. **Algebraic base:** 01 → 02 → 03 → 04.
2. **Cubical base:** 05 → 06; 05 → 07 → 08 → 09; 05–07 → 10.
3. **Geometric/FEM base:** 11 → 12 → 13 → 14 → 15; 11–15 → 16. The flat theorem 17 is a separate benchmark; it does not close 18.
4. **Finite quantum/network models:** 05 + 19 → 20; 21–22 are combinatorial/network branches and do not imply PRISM dynamics.
5. **Kathara/refinement branch:** 23 → 24; 25 is conditional on its source geometry; 23 + 05 → 26; 27 is a local simplex result; 07 → 28.
6. **Q6–PRISM branch:** 05 + PRISM definition → 29–31; 32 is a graph-level path-valued existence construction; 33 is its symmetry obstruction; 34 changes the target representation; 35–38 concern a distinct A5-root periodic model.
7. **De Aquino branch:** 38 → 39; 40 supplies a frame/ansatz identity but not the missing physical map.
8. **Response branch:** 13–15 + declared action → 41–44; 45 is an algorithmic equilibrium benchmark; 46 is an OU theorem. No arrow from 45 or 46 alone to physical viscosity is valid.
9. **Pellis branch:** 03 → 47; 02 + declared fibre representation → 48; finite Fibonacci representation → 49 → 51 conditionally; 50 remains the cross-model integration problem.

An edge in this list means “uses as input or hypothesis”, not “is equivalent to”. In particular:
- Q6 graph homology, Q6 filled-complex homology and PRISM periodic homology are different objects.
- The original Cartesian PRISM, the enriched \(\Lambda^2\) target and the GEC-6 A5-root Cayley torus are different graphs/groups.
- The flat Mosco theorem does not imply curved Mosco convergence.
- A finite Gibbs autocorrelation or OU integral is not by itself a physical Green–Kubo viscosity.
- Representation equivalence does not select a physical operator or identify spectra across models.

## 5. Priority proof obligations

| Priority | Open gate | Required closure |
|---|---|---|
| P1 | Curved Ellis–RGS Mosco (18) | Uniform mesh-map and inverse bounds, metric/Jacobian equivalence, coercivity/compactness, then liminf and recovery in the common \(L^2\) space. |
| P1 | Canonical geometric Kathara–Q6 lift (26) | A nondegenerate, incidence-compatible and symmetry-compatible map, with a stated metric-selection principle. |
| P1 | Q6–PRISM physical/dynamical bridge (29–38) | A specified map of spaces and operators preserving the claimed structures; no inference from graph counts or symmetry alone. |
| P1 | Physical De Aquino normalization (39–40) | Derive the momentum map, factor 1/4, length calibration and intertwiner from a source-defined dynamical model. |
| P1 | Physical Green–Kubo viscosity (41–46) | Full metric-dependent action, physical-time dynamics, equilibrium measure, FDT, stress observable and controlled refinement/continuum limit. |
| P1 | Pellis/GSR operator identification (47–51) | Source-defined GSR operator, canonical map to the finite models, sector selection and proof of any channel degeneracy. |
| P2 | Curved spinorial/Wilson convergence | Uniform connection and mesh estimates, then strong-resolvent/Mosco-type control and doubler analysis. |
| P2 | Metric DEC/Hodge completion | Derive primal/dual measures and a consistent discrete Hodge star from the recovered geometry. |

## 6. Canonical implementation and evidence entry points

The detailed derivations, code and tests remain in their focused records. The main entry points are:

- Q6 square complex: [Q6_COMPLEX_FORMALIZATION.md](Q6_COMPLEX_FORMALIZATION.md), `src/q6_complex.py`, `tests/`.
- RGS seed and symmetry: [RGS_Q6_OH_EQUIVARIANCE.md](RGS_Q6_OH_EQUIVARIANCE.md), `src/rgs_q6_oh_action.py`.
- Q6 Hodge: [Q6_COMBINATORIAL_HODGE_LAPLACIANS.md](Q6_COMBINATORIAL_HODGE_LAPLACIANS.md).
- Q6 shear/OU: [Q6_HYPERCUBE_GREEN_KUBO_SUPERGOLDEN.md](Q6_HYPERCUBE_GREEN_KUBO_SUPERGOLDEN.md), [Q6_ACTION_DERIVED_SHEAR_RESPONSE.md](Q6_ACTION_DERIVED_SHEAR_RESPONSE.md).
- Curved Ellis refinement: [RGS_ELLIS_CURVED_REFINEMENT_G1.md](RGS_ELLIS_CURVED_REFINEMENT_G1.md), [RGS_ELLIS_CURVED_G4_COMMON_L2.md](RGS_ELLIS_CURVED_G4_COMMON_L2.md), [RGS_MOSCO_CURVED_REFINEMENT_PROGRAM.md](RGS_MOSCO_CURVED_REFINEMENT_PROGRAM.md).
- Kathara and quadrature: [KATHARA_15_LINE_GRAPH_FACE_AUDIT.md](KATHARA_15_LINE_GRAPH_FACE_AUDIT.md), [KATHARA_KRYSTAL_SPIRAL_SQUARE_GEOMETRY.md](KATHARA_KRYSTAL_SPIRAL_SQUARE_GEOMETRY.md), [GATE_E1_BARYCENTRIC_COVARIANCE_QUADRATURE.md](GATE_E1_BARYCENTRIC_COVARIANCE_QUADRATURE.md).
- Q6–PRISM and GEC-6: `docs/research/Q6_PRISM_*.md`, `docs/research/GEC6_*.md`.
- De Aquino: [DEAQUINO_PENTAGON_A5_MOMENTUM_NORMALIZATION.md](research/DEAQUINO_PENTAGON_A5_MOMENTUM_NORMALIZATION.md).
- Pellis: [PELLIS_A_2I_SPECTRAL_GEOMETRY.md](PELLIS_A_2I_SPECTRAL_GEOMETRY.md), [PELLIS_B_TWO_SHEET_SELF_ADJOINT_OPERATORS.md](PELLIS_B_TWO_SHEET_SELF_ADJOINT_OPERATORS.md), [PELLIS_C_D6_CHANNEL_INTERTWINER.md](PELLIS_C_D6_CHANNEL_INTERTWINER.md), [PELLIS_D_RGS_PSU_PRISM_INTEGRATION.md](PELLIS_D_RGS_PSU_PRISM_INTEGRATION.md).
- Programme-wide theorem records and proof statuses: [THEOREM](https://github.com/fabianleonaressi95/theorem).

## 7. Maintenance rule

When a new result is added:
1. define its mathematical object and conventions;
2. state the exact proposition and all hypotheses;
3. assign one evidence label from §1;
4. link the implementation, tests, audit output and immutable revision;
5. state what the result does **not** establish;
6. add only justified dependency edges;
7. update the short repository README, not by pasting test logs or duplicating the full proof.

No result may be upgraded from CONDITIONAL, CANDIDATE or OPEN to PROVED solely because a finite test passes. No physical identification may be inferred solely from equal dimensions, matching spectra, shared symmetry labels or a convenient normalization.
