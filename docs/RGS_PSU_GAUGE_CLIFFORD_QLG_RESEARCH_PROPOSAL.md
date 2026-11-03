# Research Proposal: RGS–PSU Gauge–Clifford QLG

**A coordinate-independent, gauge-equivariant and Clifford-steerable framework for discrete geometric structures**

## Abstract

This proposal extends the repository's quaternionic bridge from validated SU(2)/SO(3) conventions and lattice holonomies toward a reproducible framework for geometric fields on discrete complexes. It combines (1) local-frame and gauge-equivariant transport on geometric graphs and meshes, (2) representation-theoretic construction of steerable kernels and intertwiners, (3) Clifford-algebra-valued feature fields, including a Lorentzian target signature, and (4) discrete-to-continuum analysis for operators intended to approximate geometric differential operators.

The test structures are the Q6 graph, its square faces, a specified Möbius subcomplex, and the RGS geometric and quaternionic data when their exact provenance is available.

This is a mathematical and computational research proposal. It does not assume that broader physical interpretations of RGS–PSU—including claims about the origin of mass, light, gravity or spacetime—have been established. Results will be labelled PROVED, NUMERICALLY VERIFIED, CANDIDATE or OPEN.

## 1. Motivation and relation to this repository

The repository establishes a focused foundation: quaternion conventions, the SU(2) double cover of SO(3), rotation composition, and the distinction between quaternion equality and equality of represented rotations. It also identifies lattice parallel transport and plaquette holonomy as target applications.

The next question is:

> Can local geometric data, discrete connection transport and representation-valued fields be combined into operators whose outputs are independent of arbitrary choices of local frame, with symmetry constraints that can be tested independently?

This question connects to five research lines:

- **Gauge Equivariant Mesh CNNs** (de Haan, Weiler, Cohen and Welling, ICLR 2020): anisotropic convolution and parallel transport on geometric graphs.
- **A Wigner–Eckart Theorem for Group Equivariant Convolution Kernels** (Lang and Weiler, ICLR 2021): representation-theoretic characterization of equivariant kernels.
- **Coordinate Independent Convolutional Networks** (Weiler, Forré, Verlinde and Welling, 2021): isometry- and gauge-equivariant operators on Riemannian manifolds.
- **Steerable Partial Differential Operators for Equivariant Neural Networks** (Jenner and Weiler, ICLR 2022): symmetry constraints for differential operators.
- **Clifford-Steerable Convolutional Neural Networks** (Zhdanov et al., ICML 2024): Clifford-valued features and steerable operations for geometric and pseudo-Euclidean settings.

These works provide mathematical methods, not evidence for the physical validity of RGS–PSU. Our contribution must be a precise adaptation, independently reproducible tests, and—if successful—new mathematical results for the selected discrete complexes.

## 2. Research questions

**Q1 — Local-frame independence.** Given a discrete geometric complex with a local frame at each vertex, can feature transport and message-passing operators be defined so that changing local frames changes only feature coordinates, not the represented geometric output?

**Q2 — Quaternionic transport and holonomy.** Can edge transports be represented consistently by SO(3) matrices and, where appropriate, by SU(2) or unit-quaternion lifts? Can gauge covariance, inverse-edge consistency and plaquette holonomy be verified independently?

**Q3 — Möbius transitions and orientation reversal.** Can the Möbius subcomplex be equipped with transition functions that correctly encode non-orientability? Which feature representations are compatible with the reflection transition? If a spinorial lift is required, which Pin-type structure and conventions are appropriate?

**Q4 — Clifford-valued fields.** Can scalar, vector, bivector and higher-grade fields be represented in a Clifford algebra and coupled by kernels satisfying the appropriate group-intertwining constraints? Can a Lorentzian Cl(1,3) feature model be formulated without confusing the Lorentz frame group with the internal SU(2) connection?

**Q5 — Continuum limit.** Under what geometric, analytic and scaling assumptions do discrete operators converge to continuum differential operators? Equivariance alone is not a convergence theorem. Consistency, stability and compactness require separate analysis.

## 3. Discrete geometric test beds

### 3.1 Q6 graph

The Q6 graph is the six-dimensional hypercube: 64 vertices, 192 undirected edges, degree 6, and 240 elementary square faces. Vertex labels, bit-order convention and face orientations must be fixed in the data specification and tested. The incidence matrices must satisfy the chain-complex identity

\[
B_1 B_2 = 0.
\]

All incidence matrices, orderings and orientations must be exported with reproducible metadata.

### 3.2 Möbius subcomplex

The selected Möbius test complex has 12 vertices, 18 edges, 6 square faces, one boundary component, Euler characteristic chi = 0 and first Betti number beta_1 = 1.

These six faces are a subcomplex used to test orientation-reversing chart transitions. They must not be conflated with the 240 square faces of the full Q6 graph. The incidence convention must be explicit, with B_1 B_2 = 0 verified. Non-orientability must be tested by consistent face-orientation propagation, not inferred from beta_1 alone.

### 3.3 RGS geometric data

The repository may consume RGS vertex coordinates, tangent bases, vierbein and connection data produced by the broader project. Every imported dataset must include source filename and hash, array names/shapes/dtypes, coordinate and index conventions, units and normalization constants, construction provenance, and validation residuals. Missing arrays must never be silently replaced by synthetic data.

## 4. Mathematical framework

### 4.1 Local frames and gauge transformations

Let K = (V,E,F) be a discrete complex. Associate a feature space F_i with each vertex i. A local-frame change g_i in a structure group G acts through a representation rho:

\[
f_i' = \rho(g_i) f_i.
\]

For an oriented edge i to j, let P_ij transport features from the fiber at j to the fiber at i. Gauge covariance requires

\[
P_{ij}' = \rho(g_i) P_{ij} \rho(g_j)^{-1}.
\]

It follows that

\[
P_{ij}' f_j' = \rho(g_i)P_{ij}f_j.
\]

This identity is the basic local gauge-covariance test. The implementation must also enforce P_ji = P_ij^{-1} and test it numerically.

The spatial frame group and internal gauge group are distinct unless an explicit bundle construction identifies them. In particular, Lorentz/spatial frame transformations must not be silently identified with internal SU(2) gauge transformations.

### 4.2 Gauge-equivariant message passing

A candidate operator is

\[
(\mathcal D_\theta f)_i =
\sum_{j\sim i} K_\theta(\xi_{ij})P_{ij}f_j,
\]

where xi_ij denotes geometric edge data in a specified local frame. The kernel must satisfy the appropriate steerability constraint:

\[
K_\theta(g\cdot\xi) =
\rho_{\rm out}(g)K_\theta(\xi)\rho_{\rm in}(g)^{-1}.
\]

The action on xi, representations and kernel domain must be specified for each structure group. This is a design constraint, not a claim that an arbitrary learned kernel is automatically equivariant.

### 4.3 Quaternionic edge transports

For rotational transports, use a unit quaternion q_ij in SU(2) and its represented rotation R_ij in SO(3). Preserve the double-cover distinction:

\[
R(q_{ij})=R(-q_{ij}), \qquad q_{ij}^{-1}=q_{ij}^{*}.
\]

The sign of a quaternion lift is not interchangeable with its inverse. Edge reversal, path composition and loop ordering must be documented and tested.

For a closed face f, define the ordered holonomy

\[
U_f = \prod_{(i,j)\in\partial f} U_{ij}^{\epsilon_{ij}},
\]

with product order and orientation fixed in code. Tests must distinguish exact group membership, numerical distance from identity, conjugacy under gauge transformation, and projection or re-normalization error. A nontrivial plaquette holonomy is a discrete curvature indicator for the chosen connection; it is not by itself evidence of physical curvature or torsion.

### 4.4 Möbius reflection transitions

The Möbius strip has no globally consistent orientation. Its local frame transitions therefore cannot be represented globally by orientation-preserving rotations alone. For a two-dimensional tangent-frame description, the relevant orthogonal group is generally O(2). Any reduction must be demonstrated from the geometry.

For a reflection r, candidate feature representations include

\[
\rho_0(r)=1,\qquad
\rho_-(r)=-1,\qquad
\rho_{\rm reg}(r)=
\begin{pmatrix}0&1\\1&0\end{pmatrix}.
\]

Their compatibility with all transition functions must be verified. The binary icosahedral group 2I, a subgroup of SU(2), covers the rotational icosahedral group; it does not automatically include reflections. If reflection transitions are lifted to spinorial degrees of freedom, an appropriate Pin-type construction must be specified and tested. Non-orientability alone does not imply a particular SU(2) holonomy such as -I.

### 4.5 Clifford-valued feature fields

For a Lorentzian target, use

\[
\eta=\operatorname{diag}(-1,+1,+1,+1)
\]

and the corresponding Clifford algebra Cl(1,3), with the Clifford-relation convention explicitly recorded. A multivector feature is decomposed by grade:

\[
\Psi=\psi^{(0)}+\psi^{(1)}+\psi^{(2)}+\psi^{(3)}+\psi^{(4)}.
\]

The implementation must specify basis ordering and signature, Clifford product convention, grade projections, group action on each grade, input/output representation spaces, intertwiner basis and normalization, and numerical tolerances.

A Lorentzian tetrad may define the metric through

\[
g_{\mu\nu}=\eta_{AB}E^A{}_\mu E^B{}_\nu,
\]

but the tetrad is input geometric data here. This proposal does not derive the tetrad dynamically or establish a gravitational field equation.

## 5. Representation theory and intertwiners

For each symmetry group, decompose feature spaces into irreducible representations. Equivariant linear maps are intertwiners:

\[
\operatorname{Hom}_{G}(V_{\rm in},V_{\rm out})
=
\{T:T\rho_{\rm in}(g)=\rho_{\rm out}(g)T,\ \forall g\in G\}.
\]

For the binary icosahedral group 2I, the project should construct or import the 120 unit-quaternion elements with a documented multiplication convention; verify identity, inverses and closure; construct relevant representations; compute intertwiner spaces by solving the linear constraints; verify each intertwiner against every group element; and compare numerical dimensions with representation-theoretic predictions.

The Wigner–Eckart framework motivates this decomposition and can reduce independent kernel parameters. It does not determine physical coupling constants or mass scales.

## 6. Proposed software architecture

The following is a target architecture, to be reconciled with the actual repository tree before implementation:

    rgs_qlg/
      algebra/
        quaternions.py
        so3.py
        clifford.py
        representations.py
        intertwiners.py
      geometry/
        complexes.py
        q6.py
        mobius.py
        frames.py
        vierbein.py
      connections/
        edge_transport.py
        gauge_transform.py
        holonomy.py
        pin_lifts.py
      operators/
        steerable_kernels.py
        gauge_message_passing.py
        discrete_pdo.py
      continuum/
        consistency.py
        convergence.py
      audits/
        audit_quaternion_group.py
        audit_gauge_covariance.py
        audit_mobius_transitions.py
        audit_clifford_equivariance.py
        audit_intertwiners.py
      tests/
        test_quaternions.py
        test_edge_transport.py
        test_holonomy.py
        test_q6_incidence.py
        test_mobius_topology.py
        test_gauge_covariance.py
        test_clifford_algebra.py
        test_intertwiners.py
      results/
        json/
        matrices/
        figures/
      docs/

Existing files and APIs must be preserved. New modules should be added incrementally rather than replacing the validated quaternion bridge.

## 7. Verification plan and acceptance criteria

### Milestone 1 — Algebraic foundations

Test quaternion norm, inverse and multiplication; SU(2) closure and SO(3) orthogonality/determinant; R(q)=R(-q); composition order against matrix multiplication; and edge reversal.

**Acceptance:** deterministic and seeded randomized tests pass, with residuals and tolerances recorded.

### Milestone 2 — Complexes and topology

Reproduce Q6 vertex, edge and square-face counts; verify B_1 B_2 = 0; reproduce Möbius counts, boundary components and Betti numbers; verify non-orientability through orientation propagation.

**Acceptance:** exact integer incidence identities and independently checked topology.

### Milestone 3 — Gauge-covariant transport

Generate arbitrary vertex gauge transformations; transform edge transports using the stated law; verify covariance of transported features; verify gauge conjugation of plaquette holonomies; test inverse-edge and path-composition consistency.

**Acceptance:** normalized residuals below declared floating-point tolerances across randomized gauge transformations.

### Milestone 4 — Möbius feature transport

Implement reflection transitions; test scalar, sign-changing and regular representations; transport features around loops crossing the orientation-reversing seam; compare results in different local-frame choices; investigate Pin lifts only after the orthogonal transition system is correct.

**Acceptance:** local gauge equivariance and seam compatibility pass independently. No claim of torsion follows from this milestone.

### Milestone 5 — Clifford-steerable operators

Implement the selected Clifford algebra and grade decomposition; construct group actions on feature grades; compute kernel intertwiners; test the intertwining identity for every group generator; compare direct and representation-decomposed implementations.

**Acceptance:** algebra identities and equivariance tests pass, with basis conventions and numerical conditioning documented.

### Milestone 6 — Discrete differential operators and continuum analysis

Define the discrete operator and its inner product; verify symmetry or formal-adjoint properties where required; test consistency on sampled smooth fields; test stability and spectral behavior; formulate and prove—or clearly leave open—the compactness and convergence assumptions needed for a Mosco-type limit.

**Acceptance:** no continuum-limit claim without a stated theorem, hypotheses and reproducible numerical evidence.

## 8. Evidence and claim policy

Every result must be labelled:

- **PROVED:** follows from a mathematical derivation with stated assumptions.
- **NUMERICALLY VERIFIED:** reproduced by code with disclosed inputs, tolerances and tests.
- **CANDIDATE:** specified construction whose properties are not yet fully established.
- **OPEN:** unresolved mathematical or physical question.

In particular, gauge equivariance does not imply physical gauge dynamics; a nontrivial holonomy does not establish a physical particle or field; Möbius non-orientability does not imply torsion; a Clifford representation does not derive a spacetime metric; a discrete spectral pattern does not establish a mass spectrum; and a numerical convergence plot does not replace a continuum convergence proof.

## 9. Expected contributions

If completed, the project will provide:

1. a reproducible implementation of quaternionic and representation-valued transport on discrete complexes;
2. a verified gauge-equivariant operator for Q6 and a Möbius test complex;
3. a representation-theoretic construction of steerable kernels, including a Clifford-valued extension;
4. independent tests separating rotational, reflection and internal gauge structures;
5. a mathematically explicit starting point for discrete-to-continuum analysis.

The principal contribution is not a claim of physical unification. It is a carefully specified and falsifiable mathematical-computational framework in which broader interpretations must meet explicit geometric, algebraic and analytic tests.

## 10. Selected references

1. P. de Haan, M. Weiler, T. Cohen, M. Welling, **Gauge Equivariant Mesh CNNs: Anisotropic convolutions on geometric graphs**, ICLR 2020. https://arxiv.org/abs/2003.05425
2. L. Lang, M. Weiler, **A Wigner-Eckart Theorem for Group Equivariant Convolution Kernels**, ICLR 2021. https://arxiv.org/abs/2010.10990
3. M. Weiler, P. Forré, E. Verlinde, M. Welling, **Coordinate Independent Convolutional Networks: Isometry and Gauge Equivariant Convolutions on Riemannian Manifolds**, 2021. https://arxiv.org/abs/2106.06020
4. E. Jenner, M. Weiler, **Steerable Partial Differential Operators for Equivariant Neural Networks**, ICLR 2022. https://arxiv.org/abs/2106.10163
5. M. Zhdanov, D. Ruhe, M. Weiler, A. Lucic, J. Brandstetter, P. Forré, **Clifford-Steerable Convolutional Neural Networks**, ICML 2024. https://arxiv.org/abs/2402.14730

Bibliographic metadata and versions should be checked against publisher or arXiv records before formal publication.

## 11. Immediate next action

First audit the repository tree and map existing modules, tests and data files to the milestones above. Preserve all validated quaternion conventions and tests. Add new modules incrementally, beginning with Q6 incidence and gauge-covariance tests; do not replace the existing quaternion bridge with an unverified end-to-end model.
