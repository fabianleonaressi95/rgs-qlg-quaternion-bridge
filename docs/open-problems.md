# Open Problems and Research Roadmap

## 1. Purpose and epistemic status

This document records the principal open mathematical and computational questions of the RGS–QLG Quaternion Bridge.

The repository aims to establish a reproducible framework for quaternionic rotations, the double covering

\[
SU(2)\longrightarrow SO(3),
\]

discrete gauge connections, and lattice holonomies.

It does not assume that these mathematical structures already establish a physical theory.

We use the following status labels:

- **PROVED** — established by an analytical derivation under explicitly stated assumptions.
- **NUMERICALLY VERIFIED** — tested computationally with reproducible inputs and tolerances.
- **CANDIDATE** — a mathematically motivated construction that still requires proof or validation.
- **OPEN** — a question for which no sufficient result is currently available.
- **PHYSICAL HYPOTHESIS** — a proposed physical interpretation not established by the mathematical construction alone.

Historical numerical results mentioned in this repository are not automatically reproduced by the present codebase.

---

## 2. Global consistency of SU(2) lifts

### 2.1 The local lifting problem

For every rotation \(R_{gh}\in SO(3)\), there are two quaternionic or SU(2) representatives:

\[
q_{gh},\qquad -q_{gh}.
\]

They represent the same SO(3) rotation:

\[
\rho(q_{gh})=\rho(-q_{gh})=R_{gh}.
\]

For a graph with oriented edges, choosing one lift on every edge is a local assignment. It does not by itself establish global consistency around every closed cycle.

### 2.2 Cycle consistency

For a closed oriented cycle

\[
\gamma=(v_0,v_1,\ldots,v_n=v_0),
\]

define the SU(2) holonomy

\[
H_\gamma=
U_{v_0v_1}U_{v_1v_2}\cdots U_{v_{n-1}v_0}.
\]

The corresponding SO(3) holonomy is

\[
\rho(H_\gamma)=
R_{v_0v_1}R_{v_1v_2}\cdots R_{v_{n-1}v_0}.
\]

The projection to SO(3) loses the central sign:

\[
\rho(H_\gamma)=\rho(-H_\gamma).
\]

**OPEN:** determine which global sign assignments are compatible with the chosen graph, its cycles, its plaquettes, and the intended bundle structure.

A correct analysis must distinguish:

1. changing the lift of an individual edge;
2. performing a vertex-based gauge transformation;
3. changing the global lift class;
4. changing the orientation of a path.

These operations are not interchangeable.

### 2.3 Cohomological formulation

A possible formulation uses the central subgroup

\[
\mathbb Z_2=\{+1,-1\}\subset SU(2).
\]

The edgewise sign choices may be represented by a \(\mathbb Z_2\)-valued cochain on the graph.

However, the interpretation of these signs as a cohomology class requires explicit definitions of:

- the graph or cell complex;
- the allowed edge transformations;
- the face constraints;
- the equivalence relation between assignments;
- the relevant bundle and transition functions.

**CANDIDATE:** a cohomological classification of admissible lift assignments.

**OPEN:** establish the precise hypotheses under which such a classification is valid and determine whether it captures all global obstructions.

---

## 3. Quaternion inverse and central sign

The group inverse of a unit quaternion is

\[
q^{-1}=q^*.
\]

In SU(2),

\[
U(q)^{-1}=U(q)^\dagger.
\]

The other lift of the same inverse SO(3) rotation is

\[
-q^*.
\]

But

\[
qq^*=+1,
\qquad
q(-q^*)=-1.
\]

Therefore, \(-q^*\) is not the group inverse.

For an oriented edge, the natural reverse-link convention is

\[
U_{hg}=U_{gh}^{-1}=U_{gh}^{\dagger},
\]

or equivalently,

\[
q_{hg}=q_{gh}^{*}.
\]

**PROVED:** these identities follow from the unit-quaternion and SU(2) group structures.

**OPEN:** determine whether any additional central signs are required by a particular global lattice construction. Such signs must be introduced explicitly and audited rather than silently inserted into the inverse convention.

---

## 4. The two-sheet construction

A two-sheet construction has been considered in which

\[
U^{(-)}_{gh}=\overline{U^{(+)}_{gh}}.
\]

Complex conjugation is a well-defined operation on the matrix entries of the chosen representation.

It must not automatically be identified with:

- Hermitian conjugation;
- group inversion;
- the antipodal lift \(U\mapsto -U\);
- spatial parity;
- time reversal;
- charge conjugation.

These operations have different mathematical definitions and potentially different physical meanings.

**OPEN:** establish the transformation properties of the two-sheet construction under gauge transformations, path reversal, and changes of representation.

**PHYSICAL HYPOTHESIS:** interpreting the two sheets as physical sectors requires an independently specified dynamical and observational framework.

---

## 5. Gauge covariance of lattice holonomies

Let \(G_g\in SU(2)\) be a gauge transformation at vertex \(g\).

An oriented link transforms as

\[
U_{gh}\mapsto U'_{gh}
=
G_gU_{gh}G_h^{-1}.
\]

For a path from \(g_0\) to \(g_n\), its parallel transport transforms as

\[
P_\gamma\mapsto
G_{g_0}P_\gamma G_{g_n}^{-1}.
\]

For a closed loop based at \(g_0\),

\[
H_\gamma\mapsto
G_{g_0}H_\gamma G_{g_0}^{-1}.
\]

Consequently,

\[
\operatorname{Tr}(H_\gamma)
\]

is gauge invariant.

**PROVED:** these transformation laws follow by cancellation of the intermediate vertex transformations.

**NUMERICALLY VERIFIED:** the repository should test these identities for random SU(2) links and random vertex gauge transformations.

Required tests include:

1. unitarity of every link;
2. determinant equal to one;
3. reverse-link consistency;
4. path covariance;
5. closed-loop conjugation;
6. invariance of the Wilson trace;
7. invariance under equivalent SU(2) lifts when projected to SO(3).

---

## 6. Continuum limit and convergence

A finite graph does not by itself define a continuum gauge theory.

A continuum limit requires a controlled sequence of discrete structures.

Let

\[
\mathcal G_h=(V_h,E_h,F_h)
\]

denote a family of graphs or cell complexes indexed by a scale parameter \(h\).

A continuum-limit construction must specify:

- the geometric realization of each graph;
- the scaling of edge lengths and face areas;
- the discrete connection variables;
- the interpolation between discrete and continuum fields;
- the discrete measure and inner product;
- the domain of the limiting operator;
- the boundary conditions;
- the normalization of the discrete action or quadratic form.

### 6.1 Discrete connection

For a smooth continuum connection \(A_\mu(x)\), a link variable may be approximated by parallel transport:

\[
U_{gh}
=
\mathcal P\exp\left(
-\int_{\gamma_{gh}} A
\right),
\]

where the sign and representation conventions must be fixed consistently.

For a short link of length \(h\), the leading expansion depends on the chosen parametrization and representation.

A numerical comparison must therefore use the same conventions on both sides.

### 6.2 Quadratic forms and Mosco convergence

Let \(\mathcal E_h\) be a sequence of discrete quadratic forms and \(\mathcal E\) a proposed continuum form.

The relevant research question is whether

\[
\mathcal E_h
\xrightarrow{\mathrm{Mosco}}
\mathcal E.
\]

Mosco convergence requires two complementary conditions:

1. a lower-bound condition for weakly convergent sequences;
2. a recovery-sequence condition for every element in the limiting form domain.

The precise spaces, embeddings, measures, and convergence topology must be specified.

**OPEN:** construct the discrete-to-continuum framework and prove Mosco convergence for the relevant geometric and gauge-covariant operators.

**IMPORTANT:** numerical convergence of a finite list of eigenvalues is not a proof of Mosco convergence.

### 6.3 Spectral convergence

A successful continuum construction should investigate:

- convergence of low-lying eigenvalues;
- convergence of spectral projectors;
- convergence of resolvents;
- stability of multiplicities;
- control of spurious modes;
- dependence on graph refinement and boundary conditions.

These are separate requirements and should be reported separately.

---

## 7. Relation between the rank-2 QLG connection and the 188-dimensional RGS tangent structure

The QLG construction uses SU(2) link variables, which define a rank-2 complex representation.

The RGS–PSU tangent analysis contains a 188-dimensional tangent space and associated operators.

These are different mathematical objects.

The following identification has not been established:

\[
\mathbb C^2
\longleftrightarrow
\mathbb R^{188}.
\]

Nor has a verified reduction of the 188-dimensional tangent connection to an SU(2) connection been established.

**OPEN:** construct, if possible, an explicit equivariant map or reduction linking the two structures.

Any proposed map must specify:

- its domain and codomain;
- its rank and kernel;
- its transformation law;
- its compatibility with the relevant group actions;
- its effect on the connection and curvature;
- its numerical residuals.

A dimension count or spectral resemblance is not sufficient to establish such a map.

---

## 8. Candidate SU(2) to U(1) reduction

A possible reduction can be formulated using an adjoint field

\[
n_g=n_g^aT_a,
\qquad
|n_g|=1,
\]

at each vertex.

For a link \(g\to h\), define the covariant mismatch

\[
\mathcal R_{gh}
=
n_g-U_{gh}n_hU_{gh}^{-1}.
\]

A covariantly transported adjoint field satisfies

\[
\mathcal R_{gh}=0
\]

on every relevant link.

If such a field exists globally and satisfies the required regularity and compatibility conditions, it may define a preferred U(1) subgroup stabilizing the transported direction.

**CANDIDATE:** use this condition to investigate an SU(2) to U(1) reduction.

**OPEN:** establish existence, uniqueness or classification, global obstructions, and stability under refinement.

A small numerical residual alone does not prove a global reduction.

---

## 9. Yang–Mills interpretation

A lattice SU(2) connection permits the construction of plaquette holonomies and gauge-invariant quantities such as Wilson loops.

This is a mathematical starting point for lattice gauge theory.

It does not by itself establish:

- equivalence with continuum Yang–Mills theory;
- the physical identification of the gauge field with gluons;
- the emergence of QCD;
- confinement;
- a nonzero mass gap;
- a physical interpretation of the two-sheet construction.

**OPEN:** derive a controlled action, establish the continuum limit, and compare its observables with the appropriate Yang–Mills theory.

**PHYSICAL HYPOTHESIS:** any identification with QCD or with a specific microscopic interpretation of RGS–PSU requires additional derivations and independent tests.

No mass-gap result is claimed by this repository.

---

## 10. Historical QLG numerical baseline

Previous project work reported the following numerical quantities:

| Quantity | Reported value |
|---|---:|
| SO(3) edge rotations | 192 |
| SU(2) edge lifts | 192 |
| Plaquette holonomies | 240 |
| SU(2) unitarity residual | \(3.21\times10^{-16}\) |
| SU(2) to SO(3) residual | \(9.49\times10^{-16}\) |
| SO(3) orthogonality residual | \(1.17\times10^{-15}\) |
| Wilson angle RMS | 1.82391 |
| Wilson angle maximum | 5.90535 |

These are historical reported values, not results independently reproduced by this repository.

They should be promoted to **NUMERICALLY VERIFIED** here only after:

1. the exact input data or a deterministic generator is available;
2. the construction is implemented in the repository;
3. the tests reproduce the quantities;
4. the numerical tolerances are documented;
5. the output is saved with enough information to reproduce the run.

---

## 11. Research roadmap

### Phase A — Algebraic foundations

- [ ] Implement quaternion multiplication and inversion.
- [ ] Test associativity and norm multiplicativity.
- [ ] Implement the SU(2) representation.
- [ ] Verify determinant and unitarity.
- [ ] Verify the SU(2) to SO(3) map.
- [ ] Test the two-to-one property.

### Phase B — Lattice connection

- [ ] Implement oriented links and reverse links.
- [ ] Implement path transport.
- [ ] Implement plaquette holonomies.
- [ ] Test gauge covariance.
- [ ] Test Wilson-loop invariance.
- [ ] Audit global lift-sign assignments.

### Phase C — Global and topological questions

- [ ] Define the cell complex and its cycle structure.
- [ ] Classify admissible central sign assignments.
- [ ] Determine the relevant cohomological formulation.
- [ ] Separate gauge equivalence from lift ambiguity.
- [ ] Analyze the two-sheet construction.

### Phase D — Continuum analysis

- [ ] Define a sequence of graph refinements.
- [ ] Specify measures, embeddings, and interpolation maps.
- [ ] Define the discrete quadratic forms.
- [ ] Establish lower-bound and recovery-sequence conditions.
- [ ] Investigate Mosco convergence.
- [ ] Test spectral and resolvent convergence.

### Phase E — RGS–PSU and physical extensions

- [ ] Define a candidate map between the QLG and RGS tangent structures.
- [ ] Test its equivariance and rank.
- [ ] Investigate the adjoint-field reduction.
- [ ] Derive a gauge-invariant action.
- [ ] Investigate the Yang–Mills continuum identification.
- [ ] Keep physical interpretations separate from established mathematical results.

---

## 12. Scope discipline

The broader RGS–PSU programme includes quaternionic geometry, the binary icosahedral group \(2I\), the geometry of \(S^3/2I\), golden icosahedral coordinates, and possible cosmological extensions.

These subjects may become relevant to future extensions of the quaternion bridge.

They are not assumed as prerequisites for the algebraic and computational results in this repository.

Any future extension should state precisely:

1. which mathematical structure is added;
2. which theorem or identity is being used;
3. which new assumptions are introduced;
4. which numerical tests are performed;
5. which physical interpretation remains hypothetical.

The objective is a framework in which every claim can be traced to a definition, derivation, test, or explicitly identified open problem.
