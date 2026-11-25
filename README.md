# RGS–QLG Quaternion Bridge

**Research codebase for quaternionic geometry, discrete lattice structures, spectral operators and reproducible computational validation.**

This repository is part of the RGS–PSU–QLG research programme. It develops explicit mathematical constructions and computational tests involving unit quaternions, SU(2)/SO(3), discrete incidence geometry, lattice holonomies, symmetry representations and Dirac-type operators.

> **Evidence policy:** mathematical identities, numerical verification, candidate constructions and open physical hypotheses are reported separately. Passing a finite numerical test does not establish a continuum limit or a physical theory.


## Unified programme index

The canonical 27-chapter crosswalk and preservation register is maintained in the primary repository: [UNIFIED_PROGRAM_27_CHAPTER_CROSSWALK_2026-11-24.md](https://github.com/fabianleonaressi95/theorem/blob/main/research/UNIFIED_PROGRAM_27_CHAPTER_CROSSWALK_2026-11-24.md). The shared two-repository theorem graph is synchronized with the 239-record inventory; thematic-index edges remain provisional and are not proof dependencies.

## Current open problem for independent mathematical review

**One focused problem is being put forward for external scrutiny: establish (or identify an obstruction to) Mosco convergence of the curved Ellis–RGS P1 energies on the refinement family already specified in this repository.** This is a concrete scalar geometric problem. It is not a request to validate the whole RGS–PSU–QLG programme.

### Problem statement

Fix the compact manifold with boundary
`M_L=[-L,L]\times S^2`, with `L=3r_0`, `r_0=0.8292`, and Riemannian metric
`g=du^2+(u^2+r_0^2)d\Omega^2`.
For each refinement level `h\to0`, let `K_h` be the nested angular-radial tetrahedral mesh generated from the recovered RGS 64B seed (56 distinct spherical directions, 108 initial angular faces; angular faces are subdivided fourfold and radial intervals are bisected at each level).

The common Hilbert space and reconstruction are
`H=L^2(M_L,d\mu_g)`, `d\mu_g=(u^2+r_0^2)\,du\,d\Omega`, and the explicitly defined continuous P1 map `J_h:V_h\to H^1(M_L)\subset H`. The discrete mass and stiffness forms are assembled from the curved element metric induced by the exact Ellis embedding, not from fitted cube-sphere weights.

The target quadratic form is
`q(u)=\int_{M_L}|\nabla_g u|_g^2\,d\mu_g`, with form domain `H^1(M_L)` (the natural Neumann form; any alternative boundary condition must be stated explicitly). Let `q_h` be the corresponding curved P1 stiffness form on `V_h`, extended to `H` through `J_h` by assigning `+\infty` outside `J_h(V_h)`.

**Question.** Do `q_h` Mosco-converge to `q` in the common space `H` as `h\to0`? If so, prove both:

- **Liminf:** whenever `J_hu_h\rightharpoonup u` weakly in `H`, `q(u)\le\liminf_{h\to0}q_h(u_h)`.
- **Recovery:** for every `u\in H^1(M_L)`, construct `u_h\in V_h` with `J_hu_h\to u` strongly in `H` and `\limsup_{h\to0}q_h(u_h)\le q(u)`.

A valid answer may instead identify a counterexample or a missing hypothesis that prevents this statement.

### What is already available

- **G1:** a reproducible nested curved geometric refinement family, with the recovered seed, exact Ellis embedding and finite-level geometric audits.
- **G4:** an explicit common-`L^2` reconstruction `J_h` and curved-metric P1 mass/stiffness integrals on the same finite domain.
- **Flat benchmark:** a separate Mosco proof for a regular periodic flat P1 family. It is methodological precedent only; it does not prove the curved statement.

### What remains open

The central missing step is a uniform all-level analysis: shape regularity of the curved element maps and their inverses, uniform metric/Jacobian bounds, coercivity and compactness in the common `L^2` space, followed by the liminf and recovery arguments. Finite-level numerical audits and convergence plots are not substitutes for these estimates.

This question concerns a **Riemannian scalar Dirichlet form**. It does not assert Mosco convergence for a Lorentzian metric, nor convergence of tetrads, teleparallel field equations, spinors, Dirac/overlap operators, Yang–Mills theory or QCD. Those are separate problems and are not part of this review request.

### Minimal entry points and reproducibility

- [Curved refinement G1](docs/RGS_ELLIS_CURVED_REFINEMENT_G1.md)
- [Common-`L^2` embedding G4](docs/RGS_ELLIS_CURVED_G4_COMMON_L2.md)
- [Open proof programme G5–G6](docs/RGS_MOSCO_CURVED_REFINEMENT_PROGRAM.md)
- [Flat periodic P1 Mosco theorem](docs/FLAT_TORUS_G6_P1_MOSCO_PROOF.md) — comparison benchmark, not the curved result.

Run the curved finite-level audit and tests:

```bash
PYTHONPATH=. python examples/rgs_ellis_curved_refinement_audit.py
PYTHONPATH=. python -m pytest tests/test_rgs_ellis_refinement.py -v
```

**Requested contribution:** a mathematically precise proof strategy, a counterexample, or a reference theorem whose hypotheses can be checked against this exact mesh family and these exact forms. Please identify the assumptions used; do not infer a continuum theorem from the finite audits.

---

## Research scope

### Quaternionic and group geometry
- Hamilton quaternion conventions and composition order.
- SU(2) and its double covering of SO(3).
- Binary icosahedral group 2I, order 120, and selected representations.
- Representation intertwiners and symmetry-sector checks.

### Discrete geometry
- Q6: 64 vertices, 192 edges and 240 square faces.
- Integer incidence operators and chain-complex identities.
- Möbius-complex test constructions and local-frame conventions.
- Graph Laplacians, cochain operators and topology.

### Spectral and spinorial operators
- Clifford algebra conventions.
- Discrete Dirac-type and Dirac–Kähler constructions.
- Wilson-kernel and overlap-operator experiments.
- Spectral residuals, multiplicities and symmetry checks.

### Wilson spectral audit: periodic reference and curved finite graphs

The Wilson audit is split into two explicitly different settings:

- [Periodic versus curved-Q6 spectral audit](docs/WILSON_PERIODIC_VS_CURVED_Q6_SPECTRAL_AUDIT.md) derives the standard Euclidean Wilson Fourier symbol, enumerates Brillouin-corner masses and records the scaling window for the stated periodic model.
- [Finite curved-Q6 Wilson spectral audit](docs/CURVED_Q6_WILSON_FINITE_SPECTRAL_AUDIT.md) reports finite-graph Wilson-kernel eigenvalues, low-spectrum multiplicities, Euclidean Wilson–Dirac singular values and gamma5-Hermiticity residuals.

Run the focused tests:
```bash
PYTHONPATH=. python -m pytest tests/test_wilson_periodic_symbol.py tests/test_curved_wilson_spectral_audit.py -q
python scripts/run_wilson_periodic_symbol_audit.py --dimension 3 --spacing 1 --wilson-r 1 --alpha 0
```

The periodic Brillouin-corner count is not transferred to the curved graph. The finite curved audit does not establish a refinement-uniform gap, curved Mosco convergence, or strong-resolvent convergence; these remain OPEN.

The [recovered RGS/Q6 embedded-graph baseline](docs/RGS_Q6_WILSON_EMBEDDED_GRAPH_BASELINE.md) applies the finite spectral audit to the hash-checked 64-coordinate/192-edge seed with explicitly declared identity local frames and identity spin links. It isolates nonuniform embedded edge lengths only; it is not a curved tetrad or continuum spin-connection result. Run `PYTHONPATH=. python examples/rgs_q6_wilson_spectral_baseline.py`.

The [Ellis-RGS Levi-Civita spin connection](docs/ELLIS_RGS_LEVI_CIVITA_SPIN_CONNECTION.md) now constructs edgewise SU(2) links by integrating the exact ambient parallel-transport equation on the isometrically embedded Ellis hypersurface. This is a geometric finite-mesh connection; uniform refinement estimates and spinorial Mosco/resolvent convergence remain OPEN. Run `PYTHONPATH=. python -m pytest tests/test_ellis_spin_connection.py -q`.

The [intrinsic Ellis-RGS spinorial Wilson-Dirac operator](docs/ELLIS_RGS_INTRINSIC_SPIN_DIRAC.md) assembles the naive spatial Dirac block and positive Wilson graph term from local tangent Clifford directions, intrinsic Ellis edge lengths, lumped P1 vertex volumes and the Levi-Civita SU(2) links. Its finite-mesh audits include gamma5 hermiticity and local spin-frame gauge covariance; curved refinement convergence and doubler control remain OPEN. Run `PYTHONPATH=. python -m pytest tests/test_ellis_rgs_spin_dirac.py -q`; the full level-0 mesh report is generated by `PYTHONPATH=. python examples/run_ellis_rgs_spin_wilson_dirac_audit.py --level 0 --radial-base 2 --steps 12`.


A separate [multi-level mesh shape-regularity audit](examples/run_ellis_mesh_shape_regularity_audit.py) reports Gram conditioning, minimum Gram eigenvalues, chord diameters and Ellis volume convergence for levels 0–2. It is numerical evidence only, not a proof of uniform shape regularity. Run `PYTHONPATH=. python examples/run_ellis_mesh_shape_regularity_audit.py --levels 0 1 2 --radial-base 8`.
The conditional geometric criterion and remaining proof obligations are documented in [ELLIS_RGS_MESH_SHAPE_REGULARITY.md](docs/ELLIS_RGS_MESH_SHAPE_REGULARITY.md). The same audit document describes the sparse comparison between the curved common-L2 P1 mass and the Wilson operator's lumped vertex measure; run `PYTHONPATH=. python examples/run_ellis_lumped_mass_equivalence_audit.py --levels 0 1 --radial-base 8 --quadrature-order 4`. The finite-level Neumann Poincare gap check is available as `PYTHONPATH=. python examples/run_ellis_neumann_poincare_audit.py --levels 0 1 --radial-base 8 --quadrature-order 4`.
### Geometric field structures
- Lattice holonomies and gauge covariance.
- Discrete vierbein and connection constructions.
- Torsion and teleparallel research directions.

## Mathematical starting point: SU(2) to SO(3)

Unit quaternions form a group isomorphic to SU(2). The standard map

    rho: SU(2) -> SO(3)

is a double covering:

    rho(q) = rho(-q).

Thus q and -q represent the same spatial rotation, while remaining distinct elements of SU(2). This distinction is essential in spinorial calculations and in the treatment of the central element of 2I.

## Q6 reference structure

The six-dimensional hypercube graph has

    |V| = 64, |E| = 192, degree(v) = 6.

With the six-bit vertex convention and the oriented incidence matrices B1 and B2, the chain-complex condition is

    B1 B2 = 0.

Its graph Laplacian has exact eigenvalues

    lambda_j = 2j, multiplicity = binomial(6,j), j = 0,...,6.

These identities provide exact finite-dimensional reference tests. They do not, by themselves, establish a continuum geometric limit.

- **[Unified theorem dependency graph](docs/UNIFIED_THEOREM_DEPENDENCY_GRAPH.md)** — 320 typed edges across all 236 theorem records and 14 companion catalogue entries to theorem records in `fabianleonaressi95/theorem`; machine-readable [edge registry](docs/UNIFIED_THEOREM_EDGES.json).

## Repository map

- [Q6–D5 literature and repository evidence map](docs/Q6_D5_LITERATURE_AND_REPOSITORY_EVIDENCE_MAP.md) — external sources, internal proof files, repository inventory and open proof obligations.
- [Research priorities and collaboration map](docs/RESEARCH_PRIORITIES_AND_COLLABORATION_MAP.md) — ordered open problems, targeted researchers, public contacts and external software benchmarks.

- src/ — mathematical and computational implementations.
- tests/ — unit and structural tests.
- [experiments/README.md](experiments/README.md) — navigation, 311-file experiment inventory, exact-duplicate checks and safe variant-management rules.
- [experiments/EXPERIMENTS_CATALOG.json](experiments/EXPERIMENTS_CATALOG.json) — machine-readable catalog with path, size and Git blob SHA.
- [docs/EXPERIMENT_DEPENDENCY_GRAPH_AUDIT.md](docs/EXPERIMENT_DEPENDENCY_GRAPH_AUDIT.md) — static dependency graph findings, parse blockers and safe consolidation plan.
- docs/RGS_PSU_GAUGE_CLIFFORD_QLG_RESEARCH_PROPOSAL.md — research direction.
- docs/IMPLEMENTATION_ROADMAP.md — implementation roadmap.
- docs/STEM_INTEROPERABILITY_PROFILE.md — proposed interoperability profile.
- [Issue #3](https://github.com/fabianleonaressi95/rgs-qlg-quaternion-bridge/issues/3) — tracked implementation discussion.

The precise modules and callable functions used for a benchmark must always be identified by repository commit SHA and recorded in its provenance manifest.


### RGS–PSU wormhole parameter recovery

- [Wormhole parameter provenance and Ellis benchmark audit](docs/RGS_PSU_WORMHOLE_PARAMETER_RECOVERY_AUDIT.md) — records the historical values \(r_*,R_0,r_1,r_2\), the exact Ellis benchmark identities, the distinct repository backgrounds, and the unresolved uppercase \(R_1\). It does not claim the original RGS–PSU BVP derivation has been recovered.

### Kathara–Q6 chain-map gate

- [Kathara–Q6 gate C: line incidence, equivariance and exact obstructions](docs/KATHARA_Q6_GATE_C_LINE_INCIDENCE.md) — exact 15-line/21-edge incidence, the left/right reflection action, chain-map constraints and the bipartite obstruction to a one-edge-per-edge embedding. The nontrivial path-valued map remains open.
- [Integer Kathara–Q6 chain lift and filling ambiguity](docs/KATHARA_15_LINE_Q6_CHAIN_LIFT.md) — proves the rank-111 (H_2) ambiguity of the Q6 2-skeleton and separates existence of path-chain fillings from canonical selection.

## RGS–STEM interoperability

The [RGS–STEM Interoperability Profile](docs/STEM_INTEROPERABILITY_PROFILE.md) proposes a conservative sidecar-based exchange format for finite algebraic, geometric and spectral data.

It is a proposal for technical discussion, **not an approved Ontos extension or an endorsement by Stem Framework**. The current Ontos stream specification remains unchanged unless its maintainers agree to a formal extension.

## Evidence and reproducibility standards

Every reported result should state:

1. the mathematical object and conventions;
2. the implementation and exact repository revision;
3. the test method and expected invariant;
4. numerical residuals and tolerances, where relevant;
5. whether verification uses the same code, a separate module or an independent implementation;
6. limitations and unresolved assumptions.

Claims are classified as:

- **PROVED** — established by a mathematical proof.
- **NUMERICALLY VERIFIED** — supported by reproducible computation with stated tolerances.
- **CANDIDATE** — explicitly constructed but not fully established.
- **OPEN** — requiring further mathematical, numerical or empirical work.

In particular, finite Q6 or 2I computations do not automatically prove continuum quantum field theory, general relativity, Einstein–Cartan dynamics, antimatter phenomenology or cosmological claims.

## TOE research program: gradients, barycentres and unified mechanics

The [RGS–PSU–PRISM–QLG TOE research program](docs/TOE_GRADIENT_BARYCENTRIC_UNIFICATION_PROGRAM.md) consolidates the gradient/barycentre principle, Q6 and PRISM cochain structures, vierbein and teleparallel geometry, spinor and gauge currents, SU(3) flavour octets/decuplet, QCD, MHD, De Aquino source audit, and controlled limits of classical, relativistic, quantum and statistical mechanics. It is a research roadmap, not a claim that a complete TOE has been established. Evidence labels and proof gates are explicit in the document.

## TOE computational core (barycentric geometry + SU(3) flavour)

The TOE programme now has two reusable numerical modules:
- `src/barycentric_geometry.py`: affine simplex barycentric coordinates, weighted barycentres/covariance, P1 mass and Dirichlet stiffness matrices, and invariant audit.
- `src/su3_flavour.py`: Gell-Mann generators and a three-fundamental SU(3) flavour audit using permutation projectors and total quadratic Casimir. It verifies the 10 + 8 + 8 + 1 dimension content without claiming QCD dynamics.
- `tests/test_toe_geometric_su3_core.py`: tests affine reconstruction, degenerate-simplex rejection, FEM matrix properties, barycentre identities and SU(3) representation dimensions.

Run: `python -m pytest tests/test_toe_geometric_su3_core.py -q`.

## Collaboration

Independent mathematical review, alternative implementations, code review and reproducibility audits are welcome. The most useful contributions are precise, testable and accompanied by explicit conventions.

For the broader research overview, see the [research profile](docs/GITHUB_PROFILE_README.md).

## License

See the repository license file. If no license file is present, reuse permissions should not be inferred from this README.

## CI failure audit

- [2026-11-09 CI failure and theorem-coverage audit](docs/CI_FAILURE_AUDIT_2026-11-09.md) — records the compileall root cause, corrected escaped-newline defects, workflow trigger separation, T055 coverage and the new T063 theorem. It distinguishes failed steps from steps skipped because compilation stopped the job.

## Master catalogue of laws and theorems

The consolidated [RGS–PSU–PRISM–QLG Master Catalogue](docs/MASTER_CATALOGUE_LAWS_THEOREMS.md) records mathematical theorems, exact finite identities, numerical results, model definitions, conjectures and open physical hypotheses with explicit evidence labels. It is versioned and deliberately does not present unproved physical interpretations as established laws.

The [Phase 6–10 stable-ID register](docs/MASTER_CATALOGUE_PHASE6_10_ADDENDUM.md) expands the catalogue with theorem statements, hypotheses, numerical residuals, source files and explicit limits for Q6 incidence, curved Clifford fields, discrete Dirac, Wilson, overlap/Ginsparg–Wilson, Dirac–Kähler, the flat-torus Mosco prototype and Möbius finite-spectrum thermodynamics.

Machine-readable index: [theorem registry JSON](results/catalogue/laws_theorems_registry.json). The next analytical gate is documented in the [curved RGS Mosco refinement program](docs/RGS_MOSCO_CURVED_REFINEMENT_PROGRAM.md); its status remains OPEN until a curved refinement family and liminf/recovery proof are supplied.

The [flat-torus G0–G3 refinement benchmark](docs/FLAT_TORUS_G0_G3_AUDIT.md) now provides a periodic Freudenthal mesh generator, incidence matrices, P1 stiffness, positive lumped 0-form Hodge star, refinement diagnostics and dedicated CI tests. It is a flat reference benchmark, not a curved RGS Mosco result.


The [flat-torus G4 common-L2 embedding](docs/FLAT_TORUS_G4_COMMON_L2_EMBEDDING.md) now defines the P1 reconstruction map and exact consistent mass inner product. The [flat-torus G5 stability audit](docs/FLAT_TORUS_G5_STABILITY.md) proves coercivity and Rellich compactness for that regular flat family only; curved RGS G6 liminf/recovery remains OPEN.

The [flat periodic P1 G6 proof](docs/FLAT_TORUS_G6_P1_MOSCO_PROOF.md) establishes Mosco liminf and recovery for the regular Freudenthal finite-element family under the G4-G5 hypotheses. This theorem is separate from, and does not settle, the curved RGS-Q6 Mosco problem.

The [curved Ellis–RGS refinement (G1)](docs/RGS_ELLIS_CURVED_REFINEMENT_G1.md) now uses the recovered 64B/O_h coordinates and the RGS spherical-branch metric. Its finite-domain geometric family and P1 audits are distinct from the still-open curved Mosco G4–G6 proof and from a full RGS-EMS field-equation solution.

The [curved Ellis–RGS G4 common-L2 embedding](docs/RGS_ELLIS_CURVED_G4_COMMON_L2.md) now defines J_h into the single Hilbert space L2(M_L,dmu_g) and assembles mass/stiffness using the exact curved Ellis metric. Curved G5 uniform estimates and G6 Mosco convergence remain separate proof gates.


## Q6 square-complex formalization

The canonical six-cube 2-skeleton is formalized in
[Q6_COMPLEX_FORMALIZATION.md](docs/Q6_COMPLEX_FORMALIZATION.md).
The existing incidence constructor in src/q6_complex.py is complemented by
src/q6_spectral_audit.py, which checks the hypercube Laplacian spectrum,
the exact chain identity, the bit-complement automorphism and the cellular
homology rank. The reproducible exporter is
examples/q6_formalization_audit.py.

A necessary distinction is recorded: the graph cycle rank is 129, whereas
the cellular first Betti number after attaching all 240 square faces is 0;
the second cellular Betti number is 111. The abstract Q6-to-RGS X_Oh
geometric identification and SU(2) dynamics on this complex remain OPEN.

 
## Recovered RGS seed and induced O_h action (Phase 18)

`src/rgs_q6_oh_action.py` reconstructs the original 64-point RGS seed
through the hash-checked source routine and derives all 48 signed
coordinate permutations. It verifies the faithful group action,
composition law, coordinate antipode, and induced permutations on Q6
edges and square faces. Coordinate O_h validity and Q6-label
equivariance are reported as separate gates; the latter is not assumed.
See [RGS_Q6_OH_EQUIVARIANCE.md](docs/RGS_Q6_OH_EQUIVARIANCE.md).

Run:
`PYTHONPATH=. python -m pytest tests/test_rgs_q6_oh_action.py -q`
and
`PYTHONPATH=. python examples/rgs_q6_oh_equivariance_audit.py`.

The recovered 64 points form four full O_h orbits of sizes 8, 24, 24,
and 8 (represented by eight coordinate-sign blocks). This corrects the
older shorthand "eight O_h orbits" in the Phase 15-era description.

## Qutrit entanglement on an aperiodic chain

`src/qutrit_aperiodic.py` adds a three-level local Hilbert space (qutrit), the maximally entangled qutrit Bell state, Schmidt/von-Neumann diagnostics, and a finite SU(3)-invariant permutation chain whose bonds follow the Fibonacci substitution sequence. The finite-chain ground-state entanglement is computed by exact diagonalisation. See [QUTRIT_APERIODIC_ENTANGLEMENT.md](docs/QUTRIT_APERIODIC_ENTANGLEMENT.md).

Run: `PYTHONPATH=. python -m pytest tests/test_qutrit_aperiodic.py -q`. This is a finite quantum-spin model; it is not a QCD colour calculation or a continuum-limit result.


## Qutrit–Q6 gauge-covariant embedding

The finite Fibonacci-modulated qutrit sector is embedded on a genuine Gray path of Q6, with local SU(3) edge transporters and the covariance check H[U^g]=G H[U]G†. The implementation and its limits are documented in [QUTRIT_Q6_GAUGE_COVARIANT_EMBEDDING.md](docs/QUTRIT_Q6_GAUGE_COVARIANT_EMBEDDING.md).

Run:
PYTHONPATH=. python -m pytest tests/test_qutrit_q6_coupling.py -q
PYTHONPATH=. python examples/qutrit_q6_gauge_audit.py

This is a finite lattice quantum-spin construction, not yet dynamical Yang–Mills or QCD; coupling to PSU sheets and the 188D vierbein remains OPEN.


## Q6 weighted combinatorial Hodge Laplacians

`src/q6_hodge_laplacian.py` constructs positive weighted Hodge Laplacians on 0-, 1- and 2-cochains of the canonical 64/192/240 Q6 square complex. The finite audit checks positive semidefiniteness and harmonic dimensions `(1, 0, 111)`, matching the cellular Betti numbers. See [the mathematical specification](docs/Q6_COMBINATORIAL_HODGE_LAPLACIANS.md), `tests/test_q6_hodge_laplacian.py` and `examples/q6_hodge_laplacian_audit.py`. The weights are explicit cochain inner products, not a metric-derived dual-cell construction or a physical continuum model.


## Q6 hypercube law, Green–Kubo and supergolden scaling

The integrated finite bridge is implemented in
src/q6_hypercube_gk_supergolden.py, with tests, a reproducible example and
the mathematical evidence ledger in
[Q6_HYPERCUBE_GREEN_KUBO_SUPERGOLDEN.md](docs/Q6_HYPERCUBE_GREEN_KUBO_SUPERGOLDEN.md).
It combines the 240-plaquette cosine identity, the full-cube closed-2-cochain
projector, finite OU Green–Kubo integrals, the shear-stress derivative and
the exact supergolden geometric-series theorem. Physical PSU viscosity,
continuum Yang–Mills convergence and universal supergolden spectral scaling
remain open and are not asserted by the finite tests.

Run:
PYTHONPATH=. python -m pytest tests/test_q6_hypercube_gk_supergolden.py -q

## RGS–Q6 induced metric and graph P1 audit

The recovered 64-point RGS seed is now audited as a map
X: {0,1}^6 -> R^3. The new
[finite induced-metric audit](docs/RGS_Q6_INDUCED_METRIC.md) checks all
192 edge lengths, 240 face Gram matrices and closure defects, the
vertexwise rank-at-most-three pullback Gram, and metric-weighted 1D P1
mass/stiffness assembled on the Q6 graph. It deliberately does not call
the singular 6x6 pullback a nondegenerate six-dimensional metric or
identify graph P1 matrices with a Whitney 1-form Hodge star.

Run:
```bash
PYTHONPATH=. python -m pytest tests/test_rgs_q6_induced_metric.py -q
PYTHONPATH=. python -m src.rgs_q6_induced_metric --output results/rgs_q6_induced_metric
```

## Golden 6D→3D cut-and-project (RGS–PSU–Q6)

The [cut-and-project construction](docs/RGS_PSU_GOLDEN_CUT_AND_PROJECT.md)
implements the recovered golden icosahedral 3+3 split of Z^6 using
phi and its algebraic conjugate phi*=-1/phi. The canonical acceptance
window is the projection of the centred 6-cube, audited as a 32-vertex,
30-facet rhombic triacontahedron. It generates finite-radius patches of
the corresponding model set and projects the 64 Q6 states into physical
and internal coordinates.

The same module includes the previously specified torsion-coupled internal
time-axis estimator (cross-covariance SVD, 5,000-row permutation null,
2,000-row bootstrap by default). It requires a separately derived 64xm
torsion trace; the archived edge coframe is not silently treated as torsion.

Run:
```bash
PYTHONPATH=. python -m pytest tests/test_rgs_psu_cut_project.py -q
PYTHONPATH=. python -m src.rgs_psu_cut_project --output results/rgs_psu_cut_project --physical-radius 2.5
```

## Repository navigation

- [Repository map and file-placement rules](docs/REPOSITORY_MAP.md)
- [Root artifact inventory and safe relocation plan](docs/ROOT_ARTIFACT_TRIAGE.md)
- [Implementation roadmap](docs/IMPLEMENTATION_ROADMAP.md)
- [Physics architecture and compatibility boundaries](docs/PHYSICS_ARCHITECTURE_COMPATIBILITY.md)

The root artifact triage records duplicate blob identities and provenance cautions. Legacy root reports are not moved until code, documentation and external-path references have been checked.

## RGS-derived canonical-selection gate for Kathara–Q3

The local A/B path-length tie is not a global no-go theorem. The new
[repository-level RGS canonical-selection audit](docs/KATHARA_RGS_CANONICAL_SELECTION_AUDIT.md)
connects the recovered hash-checked 64-point seed, its induced O_h action,
the A14 188D tangent projector, the Q6 Galerkin compression, and the
five-dimensional A5 coexact spectral sector with its Lambda^2-to-R3
intertwiner. Gate C now has an explicit integral chain-map candidate
[documented here](docs/KATHARA_Q6_CHAIN_MAP_GATE_C.md) and implemented in
`src/kathara_q6_chain_map.py`. The candidate is injective on the 12 vertices
but collapses all ten source 2-cells; it is therefore not a cellular embedding
and is not RGS-canonical. A nondegenerate, equivariant lift is still missing,
so an RGS variational energy selecting A/B remains undefined.

## Kathara 12-center / 15-line graph and source comparison

The published Level-1 Kathara diagram is transcribed as 15 named primary lines on 12 centers. Expanding the three vertical polylines gives 21 elementary edges. The chosen straight-line embedding is planar and has 10 bounded faces plus one exterior face; Euler's formula is checked exactly.

- [Graph and face audit](docs/KATHARA_15_LINE_GRAPH_FACE_AUDIT.md)
- src/kathara_15_line_graph.py
- tests/test_kathara_15_line_graph.py
- examples/kathara_15_line_graph_audit.py
- [Source comparison of all three supplied pages](docs/THREE_SOURCE_COMPARATIVE_AUDIT_KATHARA_BAGUA.md)
- [Central theorem T058](https://github.com/fabianleonaressi95/theorem/blob/main/theorems/T058-kathara-15-line-planar-face-audit.md)
- [Cellular chain complex and Hodge audit](src/kathara_cellular_complex.py)
- [Cellular-complex tests](tests/test_kathara_cellular_complex.py)
- [Kathara to Q3 chain-map implementation](src/kathara_q3_chain_map.py)
- [Kathara to Q3 chain-map tests](tests/test_kathara_q3_chain_map.py)
- [Central theorem T060](https://github.com/fabianleonaressi95/theorem/blob/main/theorems/T060-kathara-q3-chain-map.md)
- [Cellular-complex audit documentation](docs/KATHARA_CELLULAR_COMPLEX_AUDIT.md)
- examples/kathara_cellular_complex_audit.py — JSON exporter for B1, B2 and the three Hodge matrices
- [Central theorem T059](https://github.com/fabianleonaressi95/theorem/blob/main/theorems/T059-kathara-cellular-chain-hodge-audit.md)

Run:
```bash
PYTHONPATH=. python -m pytest tests/test_kathara_15_line_graph.py -q
PYTHONPATH=. python examples/kathara_15_line_graph_audit.py
```

## Kathara / Krystal Spiral square geometry

The source-described 45° clockwise rotation and sqrt(2) expansion are now formalized as the exact integer similarity T = [[1,1],[-1,1]]. The center-9 / center-12 / center-6 overlay constraints determine p9 = T p6 and T^2 p6 = p12. This is a conditional Euclidean theorem about the stated source geometry, not a physical claim and not a tribonacci derivation.

- [Square-geometry specification](docs/KATHARA_KRYSTAL_SPIRAL_SQUARE_GEOMETRY.md)
- src/kathara_square_similarity.py
- tests/test_kathara_square_similarity.py
- examples/kathara_square_similarity_audit.py
- [Central theorem T057](https://github.com/fabianleonaressi95/theorem/blob/main/theorems/T057-kathara-square-similarity-center-overlay.md)

Run:
```bash
PYTHONPATH=. python -m pytest tests/test_kathara_square_similarity.py -q
PYTHONPATH=. python examples/kathara_square_similarity_audit.py
```

## Exact dyadic subdivision of the Q6 2-skeleton

The canonical oriented Q6 square complex now has an explicit midpoint-subdivision benchmark:
- [Mathematical specification](docs/Q6_DYADIC_SUBDIVISION_CHAIN_MAP.md)
- [Trigram–Krist-All–supergolden crosswalk](docs/TRIGRAM_KRISTALL_SUPERGOLDEN_CROSSWALK.md) — exact permutation conjugacy with the existing Phase 21 transfer matrix; no geometric intertwiner is claimed.
- `src/q6_dyadic_refinement.py`: constructs the 496-vertex / 1344-edge / 960-square refined 2-skeleton and integer subdivision maps S0, S1, S2.
- `tests/test_q6_dyadic_refinement.py`: exact chain-map and cell-count tests.
- `examples/q6_dyadic_refinement_audit.py`: reproducible JSON audit.

The defining identities are B1' S1 = S0 B1 and B2' S2 = S1 B2. This is a genuine orientation-preserving cubical refinement of the existing Q6 2-skeleton. It is a **reference benchmark only**: it is not identified with Kathara/Krist-All and does not derive the tribonacci transfer matrix. The central mathematical record is [T055 in THEOREM](https://github.com/fabianleonaressi95/theorem/blob/main/theorems/T055-q6-dyadic-subdivision-chain-map.md).

Run:
```bash
PYTHONPATH=. python -m pytest tests/test_q6_dyadic_refinement.py -q
PYTHONPATH=. python examples/q6_dyadic_refinement_audit.py
```

## Relationship to the central THEOREM research repository

This repository is the implementation- and experiment-oriented QLG/quaternionic component of the broader **RGS–PSU–PRISM–Q6** research programme. Its constructions, scripts, numerical audits and open mathematical questions are connected to the programme's central theorem and research index:

- **[THEOREM — RGS–PSU–PRISM–Q6 Mathematical Research Programme](https://github.com/fabianleonaressi95/theorem)** — central repository for theorem records, mathematical assumptions and proofs, reproducible audits, model definitions, bibliography, and the programme-wide register of open problems.
- [Central THEOREM README](https://github.com/fabianleonaressi95/theorem/blob/main/README.md) — overview of the research strands, theorem register, evidence policy and current priorities.
- [QLG–Dirac Lorentz-chirality result (T049D)](https://github.com/fabianleonaressi95/theorem/blob/main/theorems/T049D-lorentz-chiral-qlg-embedding.md) — Lorentz-covariant Weyl embeddings and the fact that Lorentz covariance alone does not select a chirality sign, when this record is available on the central repository's default branch.
- [Quaternionic SU(2)–SO(3) double-cover theorem (T063)](https://github.com/fabianleonaressi95/theorem/blob/main/theorems/T063-quaternion-su2-so3-double-cover.md) — fixes the scalar-first Hamilton convention, the defining SU(2) matrix representation, and the two-to-one conjugation projection onto SO(3).
- [Finite SU(2) lattice holonomy theorem (T064)](https://github.com/fabianleonaressi95/theorem/blob/main/theorems/T064-finite-su2-lattice-holonomy-gauge-covariance.md) — proves path-transport covariance, closed-loop conjugation and Wilson-trace gauge invariance.
- [Exact Q6 dyadic subdivision theorem (T055)](https://github.com/fabianleonaressi95/theorem/blob/main/theorems/T055-q6-dyadic-subdivision-chain-map.md) — exact integer chain maps for the 496/1344/960 midpoint-refined Q6 2-skeleton.
- [Exact cubical subdivision theorem (T065)](https://github.com/fabianleonaressi95/theorem/blob/main/theorems/T065-cubical-subdivision-chain-map-and-naturality.md) — general integral cubical subdivision and signed-permutation naturality.
- [Herasymov parity-face conformity theorem (T066)](https://github.com/fabianleonaressi95/theorem/blob/main/theorems/T066-herasymov-parity-affine-shared-face-conformity.md) — necessary and sufficient affine shared-face condition under its explicit geometric interpretation.
- [Cubical Stokes commutation theorem (T067)](https://github.com/fabianleonaressi95/theorem/blob/main/theorems/T067-stokes-commuting-cubical-face-integral-interpolation.md) — exact face-integral interpolation commutes with the cellular exterior derivative.
- [de Rham benchmark CI evidence, run #92](https://github.com/fabianleonaressi95/theorem/blob/main/audits/DE_RHAM_BENCHMARK_CI_2026-11-09.md) — finite checks and their proof-status boundaries.


The two repositories have complementary roles: **this repository develops and audits concrete quaternionic, Q6 and QLG constructions; the theorem repository organizes the corresponding mathematical statements, proof obligations, evidence levels and cross-project results.** A computational result in this repository should therefore be cited with its implementation path and commit, while any associated theorem-level claim should point to the relevant record in the central repository.

The central repository does not automatically certify every experiment or physical interpretation in this codebase. Each connection remains subject to the assumptions, proof status and limitations stated in the corresponding record. The cubical de Rham/FEEC benchmark and its theorem-level records are maintained in THEOREM; this bridge does not rerun that separate benchmark as a prerequisite. The bridge CI remains responsible for its own quaternion, QLG, Q6, Wilson and curved-Ellis implementation tests, which are not replaced by the de Rham smoke checks.

The central landing page has now been updated on `main` (README revision: [45841e0](https://github.com/fabianleonaressi95/theorem/commit/45841e0)). The update provides the programme-wide English overview and synchronizes the theorem/result index through T049D, including the QLG–Dirac chirality distinction, alongside the Einstein–Cartan, Mosco-convergence and cellular-dynamics research strands. Start with the [updated central README](https://github.com/fabianleonaressi95/theorem/blob/main/README.md) for the current programme map; use this repository for the implementation-level material and its reproducible audits.



## Kathara–Q6 Gate D: non-degenerate face lift

The new [Gate D construction](docs/KATHARA_Q6_NONDEGENERATE_GATE_D.md) replaces root-routed vertex-potential edge chains with shortest axis-ordered Q6 paths and fills every Kathara face by exact cubical square swaps. It constructs an integral chain map with a full-column-rank 240×10 face map. The RGS seed coordinates are attached through the existing hash-verifying loader. RGS metric canonicality and a justified Ellis/SU(2) lift remain open.

Run: `PYTHONPATH=. python -m pytest tests/test_kathara_q6_nondegenerate_chain_map.py -q`.
## Gate E.1 — barycentric–covariance quadrature

The local planar-triangle quadrature theorem is implemented in `src/barycentric_triangle_quadrature.py`, tested in `tests/test_barycentric_triangle_quadrature.py`, and specified in [the Gate E.1 record](docs/GATE_E1_BARYCENTRIC_COVARIANCE_QUADRATURE.md). It reuses the exact Q6 barycentric FEM results only as a separate combinatorial/gradient reference: the 6D Q6 vertex covariance is not identified with a triangle-cell covariance. The companion theorem-level statement is [T068](https://github.com/fabianleonaressi95/theorem/blob/main/theorems/T068-barycentric-covariance-quadrature-law.md). Gate E.2 (geometric maps, closed surface incidence and refinement compatibility) remains open.


## Interconnection-network topology research

The [Culler–Q3/Q6–Omega–Beneš research dossier](docs/research/CULLER_INTERCONNECTION_NETWORKS_Q3_Q6_OMEGA_BENES.md) records the Berkeley CS258 topology lectures, foundational Clos/Beneš literature, and the Q3/Q6 hypercube connection. Its bibliography is maintained in [interconnection_networks.bib](docs/research/interconnection_networks.bib). The historical Omega/Beneš 8×8 counts have now been independently reproduced under an explicit shuffle-normalized topology convention; see the exhaustive audit below.


## Exhaustive Omega-8 / Beneš-8 routing audit

The recovered historical 8×8 routing counts have now been independently reimplemented under an explicit shuffle-normalized switch-stage convention. See the [exhaustive audit](docs/research/OMEGA_BENES_8X8_EXHAUSTIVE_AUDIT.md), [enumerator](src/interconnection_networks_8x8.py), [tests](tests/test_interconnection_networks_8x8.py), and [audit runner](examples/interconnection_networks_8x8_audit.py). The model reaches 4,096 distinct permutations for Omega-8 and all 40,320 permutations for Beneš-8; bit reversal is blocked by Omega-8 and routed by Beneš-8. These are finite combinatorial results for the declared topology, not a Q6/Green–Kubo physical equivalence.


## Exact Q3–Q6–Omega/Beneš bridge and PRISM boundary

The [Q3–Q6–Omega/Beneš–PRISM bridge](docs/research/Q3_Q6_OMEGA_BENES_PRISM_BRIDGE.md) fixes the I Ching trigram convention (three bits, bottom-to-top), identifies each Q3 coordinate matching with one 2×2 switch stage, and encodes each Q6 vertex as an ordered input/output trigram demand. Eight such demands form a permutation matrix exactly when they use every input and output once. The implementation, tests and audit are in `src/q3_q6_network_bridge.py`, `tests/test_q3_q6_network_bridge.py`, and `results/q3_q6_network_prism_audit.json`.

The PRISM (C_{64}) factor admits a 64-label correspondence only; it is not graph-isomorphic to Q6 (degree 2 versus degree 6). No edge-preserving Q6-to-PRISM map or physical/dynamical equivalence is claimed.

Run: `PYTHONPATH=. python -m pytest tests/test_q3_q6_network_bridge.py tests/test_interconnection_networks_8x8.py -q`.


## Q6–PRISM embedding obstruction (proved)

The [Q6–PRISM embedding obstruction](docs/research/Q6_PRISM_EMBEDDING_OBSTRUCTION.md) proves that Q6 cannot embed injectively and edge-preservingly into the existing five-factor product `C64 □ C6^3 □ C5`: every pair of incident Q6 directions spans a square, and its image must use two distinct Cartesian factor directions, forcing six distinct factors. The theorem is sharp: a six-factor cycle product contains an induced Q6 as its binary corner subgraph. A parity-collapse homomorphism into the existing product is recorded only as a non-injective baseline, not as a geometric bridge.

Implementation and tests: `src/q6_prism_embedding_obstruction.py`, `tests/test_q6_prism_embedding_obstruction.py`; reproducible audit: `examples/q6_prism_embedding_audit.py` and `results/q6_prism_embedding_audit.json`. This is a graph-theoretic result; no physical or transport equivalence is claimed.

## Explicit Q6–PRISM cubical chain map

The [Q6–PRISM cubical chain-map note](docs/research/Q6_PRISM_CUBICAL_CHAIN_MAP.md) gives an integral chain map from cubical Q6 onto the binary-corner Q5 subcomplex of the existing five-factor PRISM product. On vertices it folds `(v0,...,v5)` to `(v0 XOR v5,v1,...,v4)`, identifying two Q6 directions while retaining five independent directions. Cells containing both folded directions map to zero, with orientation signs specified explicitly. The identity `∂F = F∂` is checked exactly on all 729 source cells.

Implementation/tests: `src/q6_prism_chain_map.py`, `tests/test_q6_prism_chain_map.py`; audit runner/result: `examples/q6_prism_chain_map_audit.py`, `results/q6_prism_chain_map_audit.json`. The image is a contractible corner subcomplex, not the full PRISM torus; this does not establish a homotopy equivalence, spectral correspondence, or transport theorem.


## Q6–PRISM periodic-winding obstruction

The [periodic-winding obstruction note](docs/research/Q6_PRISM_PERIODIC_WINDING_OBSTRUCTION.md) proves that any edge-to-edge graph homomorphism from the Q6 1-skeleton into the PRISM graph induces the zero map from the Q6 graph cycle lattice to the periodic homology of the full cubical PRISM, H1(PRISM; Z) = Z^5. The reason is exact: Q6 graph cycles are generated by square plaquettes, and a four-edge closed walk cannot wind around any PRISM factor of length 5 or greater. The shortest possible nonzero periodic winding requires five edges.

Code/tests/audit: `src/q6_prism_periodic_obstruction.py`, `tests/test_q6_prism_periodic_obstruction.py`, `examples/q6_prism_periodic_obstruction_audit.py`, `results/q6_prism_periodic_obstruction_audit.json`. A path-valued map with globally consistent winding remains open; subdividing edges is necessary to escape the four-step bound, but is not sufficient by itself.


## Path-valued Q6–PRISM homology realization

The [path-valued homology note](docs/research/Q6_PRISM_PATH_VALUED_HOMOLOGY.md) constructs a continuous map from the Q6 **graph** (not its filled cubical complex) to the PRISM 1-skeleton by sending five fundamental-cycle chords to the five periodic factor loops and all other source edges to constant paths. Its induced map on periodic first homology is the surjection `Z^129 -> Z^5` represented by `[I_5 | 0]`. This is a basis-dependent existence construction, not a canonical geometric or dynamical correspondence. It cannot extend over all 240 filled Q6 plaquettes.

Implementation/tests/audit: `src/q6_prism_path_valued_homology.py`, `tests/test_q6_prism_path_valued_homology.py`, `examples/q6_prism_path_valued_homology_audit.py`, `results/q6_prism_path_valued_homology_audit.json`.


## Q6–PRISM symmetry obstruction

The [symmetry-obstruction note](docs/research/Q6_PRISM_SYMMETRY_OBSTRUCTION.md) proves a no-go result for the path-valued graph setting: if the source action contains the even coordinate-permutation group A6 and the target action is induced by automorphisms of the existing PRISM product, then the induced map on periodic H1 must be zero. The target H1 symmetry image has order at most 192, forcing every A6 action on it to be trivial; the Q6 graph cycle module has no A6-invariant linear functional. Thus the path-valued surjection from Z^129 to Z^5 necessarily breaks this full coordinate symmetry.

Implementation/tests/audit: `src/q6_prism_symmetry_obstruction.py`, `tests/test_q6_prism_symmetry_obstruction.py`, `examples/q6_prism_symmetry_obstruction_audit.py`, `results/q6_prism_symmetry_obstruction_audit.json`. This is a graph-theoretic representation obstruction, not a physical or dynamical claim.


## A6-equivariant Q6 exterior-square cycle map

The [exterior-square construction](docs/research/Q6_PRISM_EXTERIOR_SQUARE_EQUIVARIANT_MAP.md) resolves the representation obstruction by enriching the target with a 15-torus indexed by unordered pairs of Q6 coordinates. The induced map is the natural surjection from the Q6 graph cycle lattice (rank 129) onto \\(\\Lambda^2\\mathbb Z^6\\) (rank 15), and is equivariant under the full coordinate permutation group S6. The original PRISM periodic summand \\(\\mathbb Z^5\\) remains present but is not hit by this equivariant map.

Implementation/tests/audit: `src/q6_prism_exterior_square_map.py`, `tests/test_q6_prism_exterior_square_map.py`, `examples/q6_prism_exterior_square_audit.py`, `results/q6_prism_exterior_square_audit.json`. The source is the Q6 graph, not the filled cubical complex; physical interpretation of the added 15 cycle directions remains open.

## GEC-6: A5-root periodic model with faithful S6 symmetry

The [GEC-6 A5-root Cayley torus](docs/research/GEC6_A5_ROOT_CAYLEY_PERIODIC_MODEL.md) constructs a compact flat five-torus from the A5 root lattice and finite periodic quotients
\`A5/mA5 ≅ (Z/mZ)^5\`. The 30 roots \`e_i-e_j\` define a connected 30-regular Cayley graph; coordinate permutations act faithfully by graph automorphisms, so the action includes A6. For \`m=3\`, the exact finite model has 243 vertices and 3,645 undirected edges.

Implementation/tests/audit: \`src/gec6_a5_root_cayley_torus.py\`, \`tests/test_gec6_a5_root_cayley_torus.py\`, \`examples/gec6_a5_root_cayley_audit.py\`, \`results/gec6_a5_root_cayley_audit.json\`.

This is a new A5-root periodic model, not the original Cartesian PRISM graph. Its symmetry realization does not establish graph isomorphism, spectral equivalence, dynamical equivalence or a physical identification with the 69,120-vertex degree-10 PRISM product. The Q6 source remains its graph cycle lattice, not the filled cubical complex.

### GEC-6 finite spectral audit

The [exact m=3 Laplacian spectrum](docs/research/GEC6_A5_LAPLACIAN_SPECTRUM_M3.md) is derived by finite-abelian Fourier analysis and independently checked by diagonalizing the full 243×243 integer Laplacian. The spectrum is \(0^{(1)},15^{(12)},24^{(30)},27^{(50)},33^{(120)},36^{(30)}\), with normalized gap \(1/2\). The comparison with Q6 (normalized gap \(1/3\)) and the Cartesian PRISM (approximately \(0.000963054666\)) is finite-graph only; no dynamical or physical equivalence is asserted.

Implementation/tests/audit: \`src/gec6_a5_spectral.py\`, \`tests/test_gec6_a5_spectral.py\`, \`examples/gec6_a5_spectral_audit.py\`, \`results/gec6_a5_laplacian_spectrum_m3.json\`.

Run:
\`\`\`bash
PYTHONPATH=. python -m pytest tests/test_gec6_a5_spectral.py -q
PYTHONPATH=. python examples/gec6_a5_spectral_audit.py
\`\`\`



## GEC-6: general-modulus spectral gap

The m=3 exact Fourier spectrum is supplemented by exhaustive character-sum gap audits for m=3,...,10. For each quotient, all \(m^5\) characters are evaluated using the 15 positive A5 roots. The normalized gap is
\[
\gamma_m=\frac1{15}\min_{k\ne0}\sum_{0\le i<j\le5}
\left(1-\cos\frac{2\pi\langle k,r_{ij}\rangle}{m}\right).
\]
The finite-modulus values are numerical evaluations of the exact Fourier formula, not exact algebraic-number claims.

In simple-root coordinates, \(Q=\sum(r_{ij}r_{ij}^T)\) has the matrix recorded in the spectral audit. Its shortest nonzero integer value is 5; an exact finite certificate uses \(Q\succ\frac32I\), reducing any hypothetical vector with value below 5 to the finite set \(\|k\|^2\le3\). Taylor expansion then gives
\[
\boxed{\lim_{m\to\infty}m^2\gamma_m=\frac{2\pi^2}{3}.}
\]
Thus the normalized gap vanishes quadratically; this family is not an expander family. No spectral, dynamical or physical equivalence with Q6 or the original Cartesian PRISM is implied.

Implementation: src/gec6_a5_spectral.py; tests: tests/test_gec6_a5_spectral.py; exhaustive finite audit: results/gec6_a5_gap_scaling_m3_m10.json.


## GEC-6 continuum A5 spectrum and De Aquino threshold

The [continuum A5 torus spectral result and De Aquino bridge](docs/research/GEC6_A5_CONTINUUM_TORUS_SPECTRUM_AND_DEAQUINO_BRIDGE.md) computes the exact dual-lattice spectrum of the flat torus V/A5. The shortest nonzero dual quadratic value is 5, giving first Laplacian eigenvalue 10*pi^2/3 with multiplicity 12. It also identifies the precise scaling relation between the finite Cayley gap and the continuum first eigenvalue: lim(m^2 gamma_m) = lambda_1(T_A5)/5 = 2*pi^2/3.

The De Aquino formula m_g/m_i0 = 3 - 2*sqrt(1+xi^2) has zero threshold xi^2=5/4. Matching this with the A5 shortest value 5 requires the proposed normalization xi_k^2=(k^T Q k)/4. This is explicitly a CANDIDATE NORMALIZATION only: the physical momentum map, the factor 1/4, and an operator intertwiner remain OPEN. The separate artificial A6/Dirac grid experiment is not evidence of such an identification.

Implementation/tests/result: src/gec6_a5_continuum_torus.py, tests/test_gec6_a5_continuum_torus.py, results/gec6_a5_continuum_spectrum.json.
Run: PYTHONPATH=. python -m pytest tests/test_gec6_a5_continuum_torus.py -q


## De Aquino–pentagon–A5 physical momentum normalization (candidate bridge)

The [momentum-normalization dossier](docs/research/DEAQUINO_PENTAGON_A5_MOMENTUM_NORMALIZATION.md) carries forward the earlier De Aquino source convention xi=Delta p/(m_i0 c)=n_r U/(m_i0 c^2), the regular-pentagon tight frame EE^T=(5/2)I_2, and its rank-two projector G5=E^T E=(5/2)P. The previously proposed balanced variational ansatz on unit five-channel amplitudes gives ||Pq||^2=1/2, hence q^T G5 q=5/4; this is proved for that ansatz, not derived from De Aquino dynamics.

The physical conversion is dimensionally explicit: Delta p_pent=p0 ||E q||, so xi_pent=(p0/(m_i0 c))*sqrt(q^T G5 q). At the balanced point the De Aquino threshold requires p0=m_i0 c. Independently, an A5 shortest dual mode on a torus of physical length L reaches the same threshold only at the conditional calibration L=4*pi*lambdaBar_C/sqrt(6). The length and the map from A5 modes to pentagonal channel amplitudes are not fixed by the abstract lattice. A canonical full-S6-equivariant pentagon-plane map is obstructed because a rank-two projector cannot commute with the irreducible standard S6 action on the five-dimensional A5 space. The physical bridge and common action remain OPEN.

## Q6 action-derived shear response and Green–Kubo gate

The finite Euclidean Q6 action now has an explicit coframe-shear response
observable, computed from the existing analytic action gradient:

- `src/q6_action_stress.py` — defines the symmetric shear generator,
  strained coframe, analytic \\( -V^{-1} dS/d\\gamma \\) response and an
  independent central-difference evaluator.
- `tests/test_q6_action_stress.py` — compares analytic and finite-difference
  derivatives and checks the sectors that are metric-independent under the
  declared strain.
- [Mathematical scope and evidence ledger](docs/Q6_ACTION_DERIVED_SHEAR_RESPONSE.md).

Important result: in the current action, the Wilson term depends only on
link holonomies, while hopping/twin/quartic terms depend on coframe volume
only. For a symmetric traceless shear their first variations vanish. The
remaining response comes from the fixed-frame coframe-difference regularizer.
Thus this is an exact derivative of the implemented finite action, but not
(yet) the full microscopic shear stress. Metric-dependent edge kinetics and
a metric-dependent gauge/Hodge term, followed by an equilibrium dynamics,
are still required before a Green–Kubo viscosity can be claimed.

Run: `PYTHONPATH=. python -m pytest tests/test_q6_action_stress.py -q`.

A separate geometric inverse-square edge convention is also audited:
`w_e=|d_e|^{-2}` gives the exact shear derivative
`w'_e(0)=-2 d_{e,x} d_{e,y}/|d_e|^4`; for a declared
`H_link=sum_e w_e E_e` this yields
`Pi_xy^Q6=-(1/V) sum_e w'_e E_e`. This conditional edge-energy response
is kept distinct from the plaquette-only Wilson action unless the edge
energies are explicitly derived from it.
## Metric-dependent Q6 plaquette response (conditional face-Hodge model)

`src/q6_metric_wilson_response.py` derives the shear response of a separately declared face-area-weighted plaquette action from the hash-recovered RGS embedding. It uses the actual 64-point coordinates and canonical Q6 faces; no arbitrary Q6-to-R3 assignment or edge weights are introduced. The face-area derivative follows Jacobi's formula and is tested against central differences.

The action uses fixed reference face areas and inverse current face areas. Because the recovered seed does not supply circumcentric dual measures, this is a conditional finite face-Hodge model, not a complete DEC Hodge star or continuum Yang–Mills action. It remains distinct from `unified_discrete_action.py` until the dual measure and scaling are derived.

Run: `PYTHONPATH=. python -m pytest tests/test_q6_metric_wilson_response.py -q`.

The same module now also includes a metric-derived two-sheet graph matter sector, `S_m=(kappa/2) sum_e ||D_e psi||^2/L_e`, using the exact inverse-length stiffness coefficient of the existing 1D P1 edge FEM. Its shear derivative is analytically computed and finite-difference tested. This is a covariant graph Dirichlet energy, not a Dirac action; continuum stress normalization remains open.

## Q6 Gibbs stress-correlation benchmark

`src/q6_gibbs_dynamics.py` samples the compact SU(2) links of the conditional metric face-action with a symmetric single-link Metropolis chain, records the gauge shear observable and computes finite normalized stress autocorrelations. This defines a finite Gibbs benchmark and satisfies detailed balance for the declared proposal. Its Monte Carlo index is algorithmic, not physical time: the integrated correlation is not a Green–Kubo viscosity. Matter fields are excluded pending a normalizable mass/constraint prescription.

Run: `PYTHONPATH=. python -m pytest tests/test_q6_gibbs_dynamics.py -q`.


## Pellis correspondence dossier (25 November 2026)

The four technical notes are now separated by evidential status:

- [A — Exact \(2I\) scalar spectral geometry](docs/PELLIS_A_2I_SPECTRAL_GEOMETRY.md): CLOSED for the scalar Laplacian on \(S^3/2I\), including Molien invariants, allowed degrees and quotient multiplicities. The derivation of Golden Spheres Resonance weights remains open.
- [B — Two-sheet self-adjoint operators](docs/PELLIS_B_TWO_SHEET_SELF_ADJOINT_OPERATORS.md): CLOSED for the defined finite Q6 scalar and fibre graph operators. The canonical cross-model map and continuum limit remain open.
- [C — Exact \(D_6\) channel intertwiner](docs/PELLIS_C_D6_CHANNEL_INTERTWINER.md): the six-state Fibonacci orbit is explicitly and unitarily factored as two sheets times three modes. This is a proved representation intertwiner; a candidate realization by two copies of the spherical l=1 dipole sector is now explicit; the physical Pellis operator identification remains conditional on that sector being selected by GSR.
- [D — RGS–PSU–PRISM–Pellis integration](docs/PELLIS_D_RGS_PSU_PRISM_INTEGRATION.md): integrated research map and remaining theorem obligations.

The exact intertwiner is implemented in src/fibonacci_d6_channel_intertwiner.py and tested in tests/test_fibonacci_d6_channel_intertwiner.py.

Run:
    PYTHONPATH=. python -m pytest tests/test_fibonacci_d6_channel_intertwiner.py -q

The SU(2)-equivariant fibre map c = epsilon (epsilon = i sigma_2), its Schur-lemma uniqueness, and the covariant edge-field exchange J_E are formalized in Pellis B. The implementation and focused tests are `src/pellis_su2_fibre.py` and `tests/test_pellis_su2_fibre.py`.

Run: `PYTHONPATH=. python -m pytest tests/test_pellis_su2_fibre.py -q`.

Important source distinction: the Stergios Pellis correspondence concerns Pellis Functions / Golden Spheres Resonance and its prospective self-adjoint resonance operator. The separate radial \(2\times2\) transfer/neck matching model is not identified with that resonance operator in the correspondence emails. The dossier keeps the two tracks separate until an explicit source-level identification exists.


The conditional radial D6 lift is additionally audited in \`src/pellis_d6_operator_gate.py\`: it verifies covariance of \`M(E) ⊗ I_6\`, the central sheet-parity boundary operator, and the predicted three-dimensional parity kernels. This is a finite tensor-product theorem, not a physical identification with GSR.

Run: \`PYTHONPATH=. python -m pytest tests/test_pellis_d6_operator_gate.py tests/test_fibonacci_d6_channel_intertwiner.py -q\`. tests/test_barycentric_triangle_quadrature.py::test_triangle_barycentre_is_vertex_average PASSED [  0%]
tests/test_barycentric_triangle_quadrature.py::test_triangle_covariance_matches_barycentric_second_moments PASSED [  0%]
tests/test_barycentric_triangle_quadrature.py::test_constant_and_affine_triangle_integrals_are_barycentrically_exact PASSED [  0%]
tests/test_barycentric_triangle_quadrature.py::test_quadratic_triangle_integral_is_exact PASSED [  0%]
tests/test_barycentric_triangle_quadrature.py::test_cubic_remainder_respects_one_sixth_bound_and_scales_as_h5 PASSED [  0%]
tests/test_barycentric_triangle_quadrature.py::test_degenerate_triangle_and_invalid_bound_are_rejected PASSED [  0%]
tests/test_binary_icosahedral.py::test_binary_icosahedral_exhaustive_audit PASSED [  0%]
tests/test_binary_icosahedral.py::test_quaternion_inverse_and_identity PASSED [  0%]
tests/test_binary_icosahedral.py::test_double_cover_q_and_minus_q_same_rotation PASSED [  0%]
tests/test_burgers_reconnection.py::test_equal_burgers_head_exchange_preserves_vertex_boundary_exactly PASSED [  1%]
tests/test_burgers_reconnection.py::test_unequal_burgers_head_exchange_is_rejected PASSED [  1%]
tests/test_burgers_reconnection.py::test_junction_conservation_uses_vector_sum_not_magnitudes PASSED [  1%]
tests/test_burgers_reconnection.py::test_junction_dimension_mismatch_is_rejected PASSED [  1%]
tests/test_burgers_reconnection.py::test_boundary_flux_rejects_self_loops PASSED [  1%]
tests/test_burgers_reconnection.py::test_boundary_flux_rejects_dimension_mismatch_on_disjoint_edges PASSED [  1%]
tests/test_cantor_supergolden_transfer.py::test_supergolden_identities PASSED [  1%]
tests/test_cantor_supergolden_transfer.py::test_parry_markov_data PASSED [  1%]
tests/test_cantor_supergolden_transfer.py::test_golden_probabilities PASSED [  1%]
tests/test_cantor_supergolden_transfer.py::test_pressure_special_values PASSED [  1%]
tests/test_cantor_supergolden_transfer.py::test_generalized_dimensions PASSED [  2%]
tests/test_cantor_supergolden_transfer.py::test_cylinder_probability PASSED [  2%]
tests/test_cantor_supergolden_transfer.py::test_q6_spectrum PASSED       [  2%]
tests/test_cantor_supergolden_transfer.py::test_audit_keeps_physical_claim_open PASSED [  2%]
tests/test_cantor_supergolden_transfer.py::test_invalid_phi[1.0] PASSED  [  2%]
tests/test_cantor_supergolden_transfer.py::test_invalid_phi[0.0] PASSED  [  2%]
tests/test_cantor_supergolden_transfer.py::test_invalid_phi[nan] PASSED  [  2%]
tests/test_clifford_qlg.py::test_clifford_anticommutator_signature_minus_plus_plus_plus PASSED [  2%]
tests/test_clifford_qlg.py::test_grade_projections_reconstruct_scalar_vector_bivector_and_pseudoscalar PASSED [  2%]
tests/test_clifford_qlg.py::test_2i_spinor_action_is_unitary_and_covariant_on_spatial_gamma_matrices PASSED [  3%]
tests/test_clifford_qlg.py::test_spinor_steerable_kernel_basis_is_four_dimensional PASSED [  3%]
tests/test_clifford_qlg.py::test_spinor_lift_is_block_diagonal_copy_of_fundamental PASSED [  3%]
tests/test_clifford_qlg.py::test_invalid_grade_projection_input_is_rejected PASSED [  3%]
tests/test_clifford_qlg.py::test_clifford_blade_basis_has_six_bivectors PASSED [  3%]
tests/test_curved_wilson_spectral_audit.py::test_curved_graph_wilson_spectral_audit_reports_finite_checks PASSED [  3%]
tests/test_curved_wilson_spectral_audit.py::test_curved_graph_audit_rejects_invalid_cutoff PASSED [  3%]
tests/test_dirac_kahler.py::test_q6_dirac_kahler_is_hermitian_and_degree_odd PASSED [  3%]
tests/test_dirac_kahler.py::test_q6_dirac_kahler_square_is_discrete_hodge_laplacian PASSED [  3%]
tests/test_discrete_cartan_geometry.py::test_metric_compatible_transport_and_reverse_are_exact_to_roundoff PASSED [  3%]
tests/test_discrete_cartan_geometry.py::test_q6_edge_connection_is_metric_compatible_on_all_192_edges PASSED [  4%]
tests/test_discrete_cartan_geometry.py::test_cartan_torsion_and_curvature_vanish_for_exact_flat_embedding PASSED [  4%]
tests/test_discrete_cartan_geometry.py::test_discrete_cartan_data_transform_covariantly_under_local_frame_changes PASSED [  4%]
tests/test_discrete_cartan_geometry.py::test_reflection_frame_changes_do_not_create_a_pin_lift PASSED [  4%]
tests/test_discrete_dirac.py::test_static_tetrad_and_curved_clifford_relations PASSED [  4%]
tests/test_discrete_dirac.py::test_curved_gamma_uses_inverse_coframe PASSED [  4%]
tests/test_discrete_dirac.py::test_polar_spin_links_are_unitary PASSED   [  4%]
tests/test_discrete_dirac.py::test_discrete_spatial_dirac_is_hermitian_and_has_real_spectrum PASSED [  4%]
tests/test_discrete_dirac.py::test_dirac_operator_scales_under_coordinate_dilation PASSED [  4%]
tests/test_discrete_dirac.py::test_nonpositive_vierbein_rejected PASSED  [  5%]
tests/test_discrete_dirac.py::test_curved_clifford_for_nonsymmetric_coframe PASSED [  5%]
tests/test_discrete_dirac.py::test_compact_kernel_field_transform_is_gauge_covariant PASSED [  5%]
tests/test_discrete_dirac.py::test_wilson_laplacian_is_hermitian_positive_and_scale_covariant PASSED [  5%]
tests/test_discrete_dirac.py::test_wilson_laplacian_is_locally_gauge_covariant PASSED [  5%]
tests/test_discrete_dirac.py::test_wilson_term_lifts_pi_mode_on_periodic_1d_ring PASSED [  5%]
tests/test_discrete_dirac.py::test_hermitian_scalar_perturbation_breaks_exact_chiral_anticommutation PASSED [  5%]
tests/test_discrete_dirac.py::test_assemble_spatial_dirac_rejects_negative_wilson_parameter PASSED [  5%]
tests/test_discrete_dirac.py::test_wilson_laplacian_rejects_nonunitary_spin_links PASSED [  5%]
tests/test_discrete_dirac.py::test_euclidean_wilson_dirac_has_gamma5_hermiticity PASSED [  5%]
tests/test_discrete_spin_geometry.py::test_euclidean_gamma_matrices_satisfy_clifford_algebra PASSED [  6%]
tests/test_discrete_spin_geometry.py::test_spin_lift_intertwines_a_nontrivial_so4_rotation PASSED [  6%]
tests/test_discrete_spin_geometry.py::test_reverse_transport_is_inverse PASSED [  6%]
tests/test_discrete_spin_geometry.py::test_q6_builder_returns_all_192_edge_transports PASSED [  6%]
tests/test_discrete_spin_geometry.py::test_orientation_reversing_or_degenerate_coframes_are_rejected PASSED [  6%]
tests/test_discrete_spin_geometry.py::test_metric_mismatch_is_reported_not_hidden PASSED [  6%]
tests/test_ellis_rgs_spin_dirac.py::test_intrinsic_ellis_operator_gamma5_hermiticity_and_wilson_positivity PASSED [  6%]
tests/test_ellis_rgs_spin_dirac.py::test_intrinsic_ellis_operator_is_local_spin_frame_gauge_covariant PASSED [  6%]
tests/test_ellis_rgs_spin_dirac.py::test_intrinsic_ellis_zero_wilson_and_positive_wilson_are_distinct PASSED [  6%]
tests/test_ellis_rgs_spin_dirac.py::test_intrinsic_ellis_assembler_rejects_incomplete_links PASSED [  7%]
tests/test_ellis_spin_connection.py::test_ellis_frame_is_orthonormal_and_tangent PASSED [  7%]
tests/test_ellis_spin_connection.py::test_radial_parallel_transport_is_identity_in_adapted_frame PASSED [  7%]
tests/test_ellis_spin_connection.py::test_edge_transport_is_reversible PASSED [  7%]
tests/test_ellis_spin_connection.py::test_spin_link_is_unitary_and_lift_matches_so3_rotation PASSED [  7%]
tests/test_ellis_spin_connection.py::test_ellis_mesh_links_are_built_on_refinement_edges PASSED [  7%]
tests/test_ellis_spin_connection.py::test_antipodal_edge_is_rejected_as_ambiguous PASSED [  7%]
tests/test_ellis_spin_connection.py::test_ellis_spin_holonomy_detects_nonzero_curvature_at_throat PASSED [  7%]
tests/test_ellis_spin_connection.py::test_rk4_edge_transport_converges_under_step_refinement PASSED [  7%]
tests/test_fibonacci_d6_channel_intertwiner.py::test_fibonacci_hexagon_has_exact_two_sheet_three_mode_factorization PASSED [  7%]
tests/test_fibonacci_d6_channel_intertwiner.py::test_half_turn_is_sheet_exchange_and_second_iterate_is_c3_rotation PASSED [  8%]
tests/test_fibonacci_d6_channel_intertwiner.py::test_sheet_parity_projectors_have_rank_three PASSED [  8%]
tests/test_fibonacci_d6_channel_intertwiner.py::test_triangle_modes_are_equivalent_to_spherical_l1_dipole_modes PASSED [  8%]
tests/test_fibonacci_d6_channel_intertwiner.py::test_selected_q6_hexagram_orbit_is_not_invariant_under_full_hypercube_laplacian PASSED [  8%]
tests/test_flat_torus_embedding.py::test_g4_embedding_has_exact_positive_l2_mass[2-4] PASSED [  8%]
tests/test_flat_torus_embedding.py::test_g4_embedding_has_exact_positive_l2_mass[2-6] PASSED [  8%]
tests/test_flat_torus_embedding.py::test_g4_embedding_has_exact_positive_l2_mass[3-3] PASSED [  8%]
tests/test_flat_torus_embedding.py::test_jh_reconstructs_nodal_values_and_is_periodic PASSED [  8%]
tests/test_flat_torus_embedding.py::test_jh_is_piecewise_affine_on_each_simplex[2-4] PASSED [  8%]
tests/test_flat_torus_embedding.py::test_jh_is_piecewise_affine_on_each_simplex[3-3] PASSED [  9%]
tests/test_flat_torus_embedding.py::test_discrete_inner_product_is_exact_p1_l2_product PASSED [  9%]
tests/test_flat_torus_embedding.py::test_embedding_rejects_wrong_nodal_dimension_and_points PASSED [  9%]
tests/test_flat_torus_mosco_p1.py::test_p1_mosco_recovery_sequence_for_smooth_periodic_mode PASSED [  9%]
tests/test_flat_torus_refinement.py::test_periodic_freudenthal_mesh_has_valid_chain_complex[2-4] PASSED [  9%]
tests/test_flat_torus_refinement.py::test_periodic_freudenthal_mesh_has_valid_chain_complex[2-7] PASSED [  9%]
tests/test_flat_torus_refinement.py::test_periodic_freudenthal_mesh_has_valid_chain_complex[3-3] PASSED [  9%]
tests/test_flat_torus_refinement.py::test_periodic_freudenthal_mesh_has_valid_chain_complex[3-4] PASSED [  9%]
tests/test_flat_torus_refinement.py::test_2d_torus_expected_periodic_counts PASSED [  9%]
tests/test_flat_torus_refinement.py::test_refinement_audit_runs_three_levels_and_energy_error_decreases PASSED [  9%]
tests/test_flat_torus_refinement.py::test_bad_mesh_parameters_are_rejected PASSED [ 10%]
tests/test_flat_torus_stability.py::test_g5_uniform_shape_regular_and_scaled_diameter_on_three_levels PASSED [ 10%]
tests/test_flat_torus_stability.py::test_g5_generalized_eigenvalue_respects_continuum_poincare[2-6] PASSED [ 10%]
tests/test_flat_torus_stability.py::test_g5_generalized_eigenvalue_respects_continuum_poincare[3-3] PASSED [ 10%]
tests/test_flat_torus_stability.py::test_g5_rejects_invalid_refinement_sequence PASSED [ 10%]
tests/test_geb_formal_layer.py::test_bach_supergolden_matrix_identity PASSED [ 10%]
tests/test_geb_formal_layer.py::test_escher_triangulated_torus_chain_complex PASSED [ 10%]
tests/test_geb_formal_layer.py::test_escher_involution PASSED            [ 10%]
tests/test_geb_formal_layer.py::test_connected_cover_requires_nonzero_class PASSED [ 10%]
tests/test_geb_formal_layer.py::test_godel_hash_is_canonical_but_not_incompleteness PASSED [ 11%]
tests/test_geb_formal_layer.py::test_strange_loop_residual PASSED        [ 11%]
tests/test_geb_formal_layer.py::test_audit_status_and_open_proofs PASSED [ 11%]
tests/test_geb_formal_layer.py::test_invalid_cover_class PASSED          [ 11%]
tests/test_gec6_a5_continuum_torus.py::test_q_is_symmetric_and_shortest_value_is_five PASSED [ 11%]
tests/test_gec6_a5_continuum_torus.py::test_first_eigenvalue_has_twelve_shortest_vectors PASSED [ 11%]
tests/test_gec6_a5_continuum_torus.py::test_de_aquino_bridge_is_not_overstated PASSED [ 11%]
tests/test_gec6_a5_root_cayley_torus.py::test_root_coordinate_round_trip PASSED [ 11%]
tests/test_gec6_a5_root_cayley_torus.py::test_all_30_roots_are_distinct_and_generate_a5 PASSED [ 11%]
tests/test_gec6_a5_root_cayley_torus.py::test_s6_preserves_root_system PASSED [ 11%]
tests/test_gec6_a5_root_cayley_torus.py::test_s6_action_is_faithful_on_root_set PASSED [ 12%]
tests/test_gec6_a5_root_cayley_torus.py::test_s6_action_remains_faithful_modulo_three PASSED [ 12%]
tests/test_gec6_a5_root_cayley_torus.py::test_finite_quotient_is_simple_30_regular PASSED [ 12%]
tests/test_gec6_a5_root_cayley_torus.py::test_finite_quotient_audit_counts PASSED [ 12%]
tests/test_gec6_a5_root_cayley_torus.py::test_invalid_modulus_rejected PASSED [ 12%]
tests/test_gec6_a5_spectral.py::test_exact_fourier_laplacian_spectrum_m3 PASSED [ 12%]
tests/test_gec6_a5_spectral.py::test_fourier_multiplicities_sum_to_number_of_vertices PASSED [ 12%]
tests/test_gec6_a5_spectral.py::test_trivial_character_gives_zero_laplacian_eigenvalue PASSED [ 12%]
tests/test_gec6_a5_spectral.py::test_numerical_matrix_diagonalization_matches_fourier PASSED [ 12%]
tests/test_gec6_a5_spectral.py::test_normalized_gaps_are_recorded_correctly PASSED [ 13%]
tests/test_gec6_a5_spectral.py::test_general_modulus_gap_values_m3_to_m8 PASSED [ 13%]
tests/test_gec6_a5_spectral.py::test_gap_scaling_audit_approaches_quadratic_limit PASSED [ 13%]
tests/test_gec6_a5_spectral.py::test_a5_quadratic_form_shortest_vector_certificate PASSED [ 13%]
tests/test_holonomies.py::test_link_matrices_are_in_su2 PASSED           [ 13%]
tests/test_holonomies.py::test_reverse_link_is_hermitian_conjugate PASSED [ 13%]
tests/test_holonomies.py::test_path_transport_uses_column_vector_order PASSED [ 13%]
tests/test_holonomies.py::test_reverse_path_transport_is_inverse PASSED  [ 13%]
tests/test_holonomies.py::test_closed_path_holonomy_is_in_su2 PASSED     [ 13%]
tests/test_holonomies.py::test_wilson_trace_is_real_and_normalized PASSED [ 13%]
tests/test_holonomies.py::test_gauge_transformation_preserves_su2_links PASSED [ 14%]
tests/test_holonomies.py::test_closed_holonomy_transforms_by_conjugation PASSED [ 14%]
tests/test_holonomies.py::test_wilson_trace_is_gauge_invariant PASSED    [ 14%]
tests/test_holonomies.py::test_gauge_transform_requires_all_vertices PASSED [ 14%]
tests/test_holonomies.py::test_invalid_link_matrix_is_rejected PASSED    [ 14%]
tests/test_holonomies.py::test_non_closed_path_is_rejected_as_holonomy PASSED [ 14%]
tests/test_holonomies.py::test_path_with_unknown_vertex_is_rejected PASSED [ 14%]
tests/test_integrated_audit_json_serialization.py::test_json_default_preserves_complex_real_and_imaginary_parts PASSED [ 14%]
tests/test_integrated_audit_json_serialization.py::test_json_default_serializes_numpy_arrays_and_scalars PASSED [ 14%]
tests/test_integrated_audit_json_serialization.py::test_json_default_rejects_unknown_objects PASSED [ 15%]
tests/test_integrated_mobius_psu_audit.py::test_mobius_is_six_face_subcomplex_not_the_full_q6_240_face_complex PASSED [ 15%]
tests/test_integrated_mobius_psu_audit.py::test_seeded_randomized_open_path_and_closed_loop_gauge_covariance[7] PASSED [ 15%]
tests/test_integrated_mobius_psu_audit.py::test_seeded_randomized_open_path_and_closed_loop_gauge_covariance[17] PASSED [ 15%]
tests/test_integrated_mobius_psu_audit.py::test_seeded_randomized_open_path_and_closed_loop_gauge_covariance[31] PASSED [ 15%]
tests/test_integrated_mobius_psu_audit.py::test_seeded_randomized_open_path_and_closed_loop_gauge_covariance[20261103] PASSED [ 15%]
tests/test_integrated_mobius_psu_audit.py::test_q6_to_prism_gray_map_is_injective_and_preserves_a_hamiltonian_cycle PASSED [ 15%]
tests/test_integrated_mobius_psu_audit.py::test_discrete_de_rham_hodge_operators_on_mobius_and_q6 PASSED [ 15%]
tests/test_integrated_mobius_psu_audit.py::test_sphere_haar_moment_reference_values_in_four_dimensions PASSED [ 15%]
tests/test_integrated_mobius_psu_audit.py::test_finite_2i_moment_audit_is_labelled_as_finite_not_equdistributive_proof PASSED [ 15%]
tests/test_integrated_mobius_psu_audit.py::test_multifractal_uniform_binary_partitions_have_dimension_one PASSED [ 16%]
tests/test_integrated_mobius_psu_audit.py::test_finite_metric_correspondence_and_mosco_diagnostic_do_not_claim_theorems PASSED [ 16%]
tests/test_integrated_mobius_psu_audit.py::test_clifford_grade_action_preserves_grades_under_2i_spin_rotations PASSED [ 16%]
tests/test_integrated_mobius_psu_audit.py::test_integrated_audit_json_safe_preserves_complex_components PASSED [ 16%]
tests/test_interconnection_networks_8x8.py::test_exact_stage_architectures PASSED [ 16%]
tests/test_interconnection_networks_8x8.py::test_omega_exhaustive_reachable_permutation_count PASSED [ 16%]
tests/test_interconnection_networks_8x8.py::test_benes_exhaustive_reachable_permutation_count PASSED [ 16%]
tests/test_interconnection_networks_8x8.py::test_benes_realizes_every_permutation PASSED [ 16%]
tests/test_interconnection_networks_8x8.py::test_bit_reversal_blocking_and_rearrangeability PASSED [ 16%]
tests/test_interconnection_networks_8x8.py::test_returned_settings_reproduce_requested_permutation PASSED [ 17%]
tests/test_interconnection_networks_8x8.py::test_omega_witnesses_reproduce_permutations PASSED [ 17%]
tests/test_interconnection_networks_8x8.py::test_invalid_permutation_rejected PASSED [ 17%]
tests/test_kathara_15_line_graph.py::test_source_line_inventory_and_center_count PASSED [ 17%]
tests/test_kathara_15_line_graph.py::test_oriented_graph_is_connected_planar_and_acyclic PASSED [ 17%]
tests/test_kathara_15_line_graph.py::test_face_count_and_euler_formula PASSED [ 17%]
tests/test_kathara_15_line_graph.py::test_expected_bounded_face_boundaries PASSED [ 17%]
tests/test_kathara_cellular_complex.py::test_chain_dimensions_and_boundary_squared_zero PASSED [ 17%]
tests/test_kathara_cellular_complex.py::test_exact_ranks_betti_and_euler_characteristic PASSED [ 17%]
tests/test_kathara_cellular_complex.py::test_hodge_laplacian_dimensions PASSED [ 17%]
tests/test_kathara_q3_chain_map.py::test_chain_map_commutes_with_both_boundary_operators PASSED [ 18%]
tests/test_kathara_q3_chain_map.py::test_map_is_integral_and_has_expected_shapes PASSED [ 18%]
tests/test_kathara_q3_chain_map.py::test_projection_is_explicitly_coarse_and_uses_one_square PASSED [ 18%]
tests/test_kathara_q6_chain_map.py::test_kathara_q6_chain_identities_are_exact PASSED [ 18%]
tests/test_kathara_q6_chain_map.py::test_map_dimensions_and_vertex_injectivity PASSED [ 18%]
tests/test_kathara_q6_chain_map.py::test_map_does_not_claim_to_embed_kathara_2_cells PASSED [ 18%]
tests/test_kathara_q6_chain_map.py::test_fifteen_primary_lines_are_retained_as_metadata PASSED [ 18%]
tests/test_kathara_q6_nondegenerate_chain_map.py::test_reverse_source_edge_reverses_the_whole_q6_path_word PASSED [ 18%]
tests/test_kathara_q6_nondegenerate_chain_map.py::test_nondegenerate_kathara_q6_chain_identities PASSED [ 18%]
tests/test_kathara_q6_nondegenerate_chain_map.py::test_face_map_is_nonzero_and_injective_on_face_chain_space PASSED [ 19%]
tests/test_kathara_q6_nondegenerate_chain_map.py::test_dimensions_and_vertex_chart_are_preserved PASSED [ 19%]
tests/test_kathara_q6_nondegenerate_chain_map.py::test_each_face_has_a_nonzero_cubical_fill PASSED [ 19%]
tests/test_kathara_square_similarity.py::test_stage_map_is_square_diagonal_similarity PASSED [ 19%]
tests/test_kathara_square_similarity.py::test_reported_center_overlays_determine_center9 PASSED [ 19%]
tests/test_kathara_square_similarity.py::test_two_stages_are_quarter_turn_and_double_scale PASSED [ 19%]
tests/test_kathara_square_similarity.py::test_lattice_image_is_checkerboard_index_two PASSED [ 19%]
tests/test_kathara_square_similarity.py::test_invalid_width_rejected PASSED [ 19%]
tests/test_koopman_rgs_psu.py::test_cantor_parry_q6_chain_is_stochastic_and_stationary PASSED [ 19%]
tests/test_koopman_rgs_psu.py::test_markov_koopman_and_weighted_adjoint_audit PASSED [ 19%]
tests/test_koopman_rgs_psu.py::test_markov_audit_rejects_nonstochastic_matrix PASSED [ 20%]
tests/test_koopman_rgs_psu.py::test_edmd_recovers_exact_affine_koopman_on_linear_system PASSED [ 20%]
tests/test_koopman_rgs_psu.py::test_edmd_reports_nonclosure_for_insufficient_dictionary PASSED [ 20%]
tests/test_koopman_rgs_psu.py::test_barycentric_feature_map_has_exact_centered_moire_observables PASSED [ 20%]
tests/test_koopman_rgs_psu.py::test_constant_frequency_koopman_generator_is_skew_and_semigroup_consistent PASSED [ 20%]
tests/test_koopman_rgs_psu.py::test_invalid_edmd_inputs_are_rejected PASSED [ 20%]
tests/test_koopman_rgs_psu.py::test_q6_rgs_lookup_reindexes_archived_big_endian_rows_exactly PASSED [ 20%]
tests/test_koopman_rgs_psu.py::test_q6_rgs_lookup_rejects_duplicate_or_missing_states PASSED [ 20%]
tests/test_microtubule_pipeline_contracts.py::DensityCoordinateContractTests::test_definition_and_calls_have_matching_arity 
tests/test_microtubule_pipeline_contracts.py::DensityCoordinateContractTests::test_definition_and_calls_have_matching_arity PASSED [ 20%]
tests/test_microtubule_pipeline_io.py::DownloadToPathTests::test_checksum_mismatch_preserves_existing_file PASSED [ 21%]
tests/test_microtubule_pipeline_io.py::DownloadToPathTests::test_gzip_stream_is_validated PASSED [ 21%]
tests/test_microtubule_pipeline_io.py::DownloadToPathTests::test_invalid_timeout_is_rejected PASSED [ 21%]
tests/test_microtubule_pipeline_io.py::DownloadToPathTests::test_streams_and_publishes_atomically PASSED [ 21%]
tests/test_mobius_complex.py::test_mobius_cardinalities_and_cell_structure PASSED [ 21%]
tests/test_mobius_complex.py::test_mobius_chain_complex_identity_is_exact PASSED [ 21%]
tests/test_mobius_complex.py::test_mobius_euler_characteristic_and_betti_numbers PASSED [ 21%]
tests/test_mobius_complex.py::test_orientation_propagation_detects_nonorientability PASSED [ 21%]
tests/test_mobius_complex.py::test_faces_have_cyclic_order_and_four_distinct_vertices PASSED [ 21%]
tests/test_mobius_complex.py::test_mobius_vertices_retain_big_endian_q6_labels PASSED [ 21%]
tests/test_mobius_complex.py::test_incidence_matrices_are_integer_valued PASSED [ 22%]
tests/test_mobius_gauge_audit.py::test_randomized_mobius_gauge_audit_passes_with_fixed_seed PASSED [ 22%]
tests/test_mobius_gauge_audit.py::test_randomized_audit_is_reproducible PASSED [ 22%]
tests/test_mobius_gauge_audit.py::test_gauge_transformation_preserves_su2_links_and_reverse_edges PASSED [ 22%]
tests/test_mobius_gauge_audit.py::test_gauge_audit_json_export_contains_reproducibility_metadata PASSED [ 22%]
tests/test_mobius_representations.py::test_six_shared_edge_transitions_and_nontrivial_mobius_holonomy PASSED [ 22%]
tests/test_mobius_representations.py::test_each_representation_obeys_z2_homomorphism[scalar] PASSED [ 22%]
tests/test_mobius_representations.py::test_each_representation_obeys_z2_homomorphism[sign] PASSED [ 22%]
tests/test_mobius_representations.py::test_each_representation_obeys_z2_homomorphism[regular] PASSED [ 22%]
tests/test_mobius_representations.py::test_scalar_sign_and_regular_cycle_actions PASSED [ 23%]
tests/test_mobius_representations.py::test_regular_feature_transport_swaps_channels_on_orientation_reversing_transition PASSED [ 23%]
tests/test_mobius_representations.py::test_local_frame_change_transforms_transition_and_preserves_covariance PASSED [ 23%]
tests/test_mobius_representations.py::test_invalid_transition_and_feature_shapes_are_rejected PASSED [ 23%]
tests/test_mobius_representations.py::test_audit_reports_orientation_obstruction_without_claiming_torsion PASSED [ 23%]
tests/test_mobius_representations.py::test_reflection_transitions_are_explicit_od_matrices_not_pin_lifts PASSED [ 23%]
tests/test_mosco_torus.py::test_flat_torus_wilson_form_converges_for_fixed_physical_mode PASSED [ 23%]
tests/test_mosco_torus.py::test_flat_torus_doubler_mass_diverges_for_alpha_inside_window PASSED [ 23%]
tests/test_mosco_torus.py::test_flat_torus_scaling_endpoints_do_not_decouple_doublers_or_low_modes PASSED [ 23%]
tests/test_mosco_torus.py::test_euclidean_wilson_form_has_only_the_physical_corner_zero PASSED [ 23%]
tests/test_naressi_golden_bridge.py::test_capsule_golden_parameters PASSED [ 24%]
tests/test_naressi_golden_bridge.py::test_geometric_coupling_ratio_is_constant PASSED [ 24%]
tests/test_naressi_golden_bridge.py::test_qat_corrected_hamiltonian_is_hermitian_and_unitary PASSED [ 24%]
tests/test_naressi_golden_bridge.py::test_golden_decay_dictionary_matches_q6_spectral_ratio PASSED [ 24%]
tests/test_naressi_golden_bridge.py::test_sparse_tfim_matches_dense_small_chain PASSED [ 24%]
tests/test_naressi_golden_bridge.py::test_tfim_gap_and_entropy_are_finite PASSED [ 24%]
tests/test_naressi_golden_bridge.py::test_bridge_does_not_claim_geometric_derivation PASSED [ 24%]
tests/test_naressi_golden_bridge.py::test_invalid_phi_is_rejected[1.0] PASSED [ 24%]
tests/test_naressi_golden_bridge.py::test_invalid_phi_is_rejected[0.0] PASSED [ 24%]
tests/test_naressi_golden_bridge.py::test_invalid_phi_is_rejected[nan] PASSED [ 25%]
tests/test_npy_manifest.py::test_manifest_records_shape_dtype_hash_and_duplicates PASSED [ 25%]
tests/test_npy_manifest.py::test_manifest_marks_semantics_as_unassigned PASSED [ 25%]
tests/test_oriented_three_channel_incidence.py::test_target_matrix_and_integer_characteristic_polynomial PASSED [ 25%]
tests/test_oriented_three_channel_incidence.py::test_target_type_graph_has_only_identity_automorphism PASSED [ 25%]
tests/test_oriented_three_channel_incidence.py::test_target_is_recognized_up_to_relabelling PASSED [ 25%]
tests/test_oriented_three_channel_incidence.py::test_relabelled_copy_is_recognized PASSED [ 25%]
tests/test_oriented_three_channel_incidence.py::test_antitwin_class_preserving_action_passes PASSED [ 25%]
tests/test_oriented_three_channel_incidence.py::test_antitwin_class_exchange_is_rejected PASSED [ 25%]
tests/test_overlap_dirac.py::test_finite_overlap_operator_satisfies_ginsparg_wilson_relation PASSED [ 25%]
tests/test_overlap_dirac.py::test_overlap_rejects_kernel_without_gamma5_hermiticity PASSED [ 26%]
tests/test_overlap_dirac.py::test_overlap_index_trace_resolves_a_single_negative_chirality_zero_mode PASSED [ 26%]
tests/test_overlap_dirac.py::test_overlap_spectral_flow_detects_sampled_sign_count_change PASSED [ 26%]
tests/test_overlap_dirac.py::test_overlap_locality_profile_uses_graph_shortest_path_shells PASSED [ 26%]
tests/test_peierls45_dependencies.py::test_peierls45_resolves_helpers_from_peierls46 PASSED [ 26%]
tests/test_peierls45_dependencies.py::test_peierls45_small_lattice_hamiltonian_is_hermitian PASSED [ 26%]
tests/test_peierls45_dependencies.py::test_peierls45_gap_selector_returns_valid_occupied_count PASSED [ 26%]
tests/test_pellis_d6_commutant.py::test_mode_projectors_are_orthogonal_rank_one_and_rank_two PASSED [ 26%]
tests/test_pellis_d6_commutant.py::test_d6_covariance_allows_distinct_axis_and_plane_transfers PASSED [ 26%]
tests/test_pellis_d6_commutant.py::test_threefold_degeneracy_is_the_special_equal_transfer_case PASSED [ 27%]
tests/test_pellis_d6_commutant.py::test_nonfinite_and_wrong_shape_matrices_are_rejected PASSED [ 27%]
tests/test_pellis_d6_operator_gate.py::test_lifted_radial_transfer_commutes_with_full_d6_channel_action PASSED [ 27%]
tests/test_pellis_d6_operator_gate.py::test_parity_boundary_operator_has_three_dimensional_branch_kernels PASSED [ 27%]
tests/test_pellis_d6_operator_gate.py::test_sheet_parity_projectors_are_orthogonal_rank_three_projectors PASSED [ 27%]
tests/test_pellis_d6_operator_gate.py::test_invalid_transfer_shapes_and_nonfinite_data_are_rejected PASSED [ 27%]
tests/test_pellis_su2_fibre.py::test_epsilon_is_unitary_skew_and_squares_to_minus_identity PASSED [ 27%]
tests/test_pellis_su2_fibre.py::test_su2_equivariance[0.0-0.0] PASSED    [ 27%]
tests/test_pellis_su2_fibre.py::test_su2_equivariance[0.2-0.7] PASSED    [ 27%]
tests/test_pellis_su2_fibre.py::test_su2_equivariance[1.1-2.4] PASSED    [ 27%]
tests/test_pellis_su2_fibre.py::test_su2_equivariance[1.5707963267948966--1.2] PASSED [ 28%]
tests/test_pellis_su2_fibre.py::test_fibre_exchange_is_self_adjoint_unitary_involution PASSED [ 28%]
tests/test_pellis_su2_fibre.py::test_edge_covariance_on_random_fibre_data PASSED [ 28%]
tests/test_pellis_su2_fibre.py::test_doubled_edge_exchange_covariance PASSED [ 28%]
tests/test_pellis_su2_fibre.py::test_su2_parameter_validation PASSED     [ 28%]
tests/test_pin_mobius_lift.py::test_pin_plus_generator_squares_to_positive_identity PASSED [ 28%]
tests/test_pin_mobius_lift.py::test_pin_twisted_adjoint_projects_to_orientation_reflection PASSED [ 28%]
tests/test_pin_mobius_lift.py::test_chirality_implements_grade_automorphism PASSED [ 28%]
tests/test_pin_mobius_lift.py::test_clifford_transition_is_an_algebra_automorphism_and_z2_action PASSED [ 28%]
tests/test_pin_mobius_lift.py::test_mobius_dual_cycle_has_consistent_pin_plus_holonomy_and_cocycle PASSED [ 29%]
tests/test_pin_mobius_lift.py::test_invalid_pin_transition_is_rejected PASSED [ 29%]
tests/test_projected_koopman_killing.py::test_projected_q6_graph_diffusion_is_symmetric_negative_semidefinite PASSED [ 29%]
tests/test_projected_koopman_killing.py::test_projected_q6_graph_diffusion_rejects_nonorthonormal_frame PASSED [ 29%]
tests/test_projected_koopman_killing.py::test_discrete_euclidean_rotations_are_edgewise_killing PASSED [ 29%]
tests/test_projected_koopman_killing.py::test_discrete_killing_rejects_invalid_edge_endpoint PASSED [ 29%]
tests/test_psu_two_sheet_dirac.py::test_identical_sheets_split_exactly PASSED [ 29%]
tests/test_psu_two_sheet_dirac.py::test_independent_sheets_and_unitary_transporter PASSED [ 29%]
tests/test_psu_two_sheet_dirac.py::test_nonunitary_transporter_rejected PASSED [ 29%]
tests/test_psu_two_sheet_dirac.py::test_scalar_sheet_coupling_is_mass_like_for_chirality PASSED [ 29%]
tests/test_psu_two_sheet_dirac.py::test_nonhermitian_input_rejected PASSED [ 30%]
tests/test_q3_q6_network_bridge.py::test_all_trigrams_and_q6_demands_round_trip PASSED [ 30%]
tests/test_q3_q6_network_bridge.py::test_q3_edges_are_three_disjoint_coordinate_matchings PASSED [ 30%]
tests/test_q3_q6_network_bridge.py::test_permutation_matrix_is_q6_perfect_matching PASSED [ 30%]
tests/test_q3_q6_network_bridge.py::test_reject_nonmatching_q6_demands PASSED [ 30%]
tests/test_q3_q6_network_bridge.py::test_all_network_permutations_have_distinct_q6_matching_encodings PASSED [ 30%]
tests/test_q3_q6_network_bridge.py::test_route_trace_replays_witness_to_each_q6_demand PASSED [ 30%]
tests/test_q3_q6_network_bridge.py::test_audit_records_prism_degree_obstruction PASSED [ 30%]
tests/test_q3_q6_network_bridge.py::test_iching_word_order_is_preserved_in_binary_labels PASSED [ 30%]
tests/test_q6_2i_intertwiners.py::test_exterior_signed_permutation_is_orthogonal_and_preserves_vacuum PASSED [ 31%]
tests/test_q6_2i_intertwiners.py::test_central_character_forbids_untwisted_doublet_intertwiner PASSED [ 31%]
tests/test_q6_2i_intertwiners.py::test_doublet_twist_has_explicit_equivariant_vacuum_embedding PASSED [ 31%]
tests/test_q6_2i_intertwiners.py::test_canonical_embedding_has_correct_shape_and_rank PASSED [ 31%]
tests/test_q6_2i_intertwiners.py::test_exterior_action_rejects_non_signed_permutations PASSED [ 31%]
tests/test_q6_2i_intertwiners.py::test_exterior_action_respects_signed_permutation_composition PASSED [ 31%]
tests/test_q6_action_stress.py::test_symmetric_shear_is_traceless_and_default_axes_are_spatial PASSED [ 31%]
tests/test_q6_action_stress.py::test_analytic_action_shear_response_matches_central_difference PASSED [ 31%]
tests/test_q6_action_stress.py::test_metric_independent_wilson_and_volume_only_terms_have_zero_shear_part PASSED [ 31%]
tests/test_q6_action_stress.py::test_bad_reference_volume_is_rejected PASSED [ 31%]
tests/test_q6_action_stress.py::test_inverse_square_edge_weight_derivative_matches_finite_difference PASSED [ 32%]
tests/test_q6_barycentric_fem.py::test_q6_hypercube_embedding_and_fem_matrix_invariants PASSED [ 32%]
tests/test_q6_barycentric_fem.py::test_q6_coordinate_fields_have_expected_edge_gradients_and_energy PASSED [ 32%]
tests/test_q6_barycentric_fem.py::test_translation_invariance_and_anisotropic_embedding PASSED [ 32%]
tests/test_q6_barycentric_fem.py::test_degenerate_q6_embedding_is_rejected PASSED [ 32%]
tests/test_q6_bisheet_entanglement.py::test_q6_bisheet_maximally_entangled_at_zero_tau PASSED [ 32%]
tests/test_q6_bisheet_entanglement.py::test_q6_bisheet_entropy_matches_exact_spectrum PASSED [ 32%]
tests/test_q6_bisheet_entanglement.py::test_entanglement_decreases_with_heat_kernel_scale PASSED [ 32%]
tests/test_q6_bisheet_entanglement.py::test_schmidt_probabilities_are_normalized PASSED [ 32%]
tests/test_q6_bisheet_entanglement.py::test_invalid_tau_rejected[-1.0] PASSED [ 33%]
tests/test_q6_bisheet_entanglement.py::test_invalid_tau_rejected[nan] PASSED [ 33%]
tests/test_q6_bisheet_entanglement.py::test_invalid_tau_rejected[inf] PASSED [ 33%]
tests/test_q6_complex.py::test_q6_cardinalities_and_regular_degree PASSED [ 33%]
tests/test_q6_complex.py::test_q6_chain_complex_identity_is_exact PASSED [ 33%]
tests/test_q6_complex.py::test_q6_boundary_matrices_have_expected_shapes PASSED [ 33%]
tests/test_q6_complex.py::test_every_square_face_has_four_distinct_vertices PASSED [ 33%]
tests/test_q6_complex.py::test_big_endian_axis_convention_is_fixed PASSED [ 33%]
tests/test_q6_complex.py::test_incidence_matrices_are_integer_valued PASSED [ 33%]
tests/test_q6_complex.py::test_canonical_q6_bundle_export_preserves_source_defined_ordering PASSED [ 33%]
tests/test_q6_dyadic_refinement.py::test_q6_dyadic_refinement_cell_counts_and_map_shapes PASSED [ 34%]
tests/test_q6_dyadic_refinement.py::test_q6_dyadic_refinement_is_an_exact_chain_map PASSED [ 34%]
tests/test_q6_dyadic_refinement.py::test_each_parent_cell_has_the_expected_oriented_children PASSED [ 34%]
tests/test_q6_dyadic_refinement.py::test_audit_reports_zero_chain_map_residuals PASSED [ 34%]
tests/test_q6_dyadic_refinement.py::test_midpoint_subdivision_preserves_each_axis_direction PASSED [ 34%]
tests/test_q6_gec6_a5_quotient.py::test_all_240_q6_square_curls_are_a5_roots PASSED [ 34%]
tests/test_q6_gec6_a5_quotient.py::test_root_vectors_lie_in_a5 PASSED    [ 34%]
tests/test_q6_gec6_a5_quotient.py::test_s6_equivariance_of_root_assignment PASSED [ 34%]
tests/test_q6_gec6_a5_quotient.py::test_edge_cocycle_values_have_zero_coordinate_sum PASSED [ 34%]
tests/test_q6_gec6_a5_quotient.py::test_audit_counts PASSED              [ 35%]
tests/test_q6_geometric_hodge.py::test_full_six_cube_cell_counts_and_betti_numbers PASSED [ 35%]
tests/test_q6_geometric_hodge.py::test_full_cube_one_skeleton_matches_canonical_q6_edge_orientation PASSED [ 35%]
tests/test_q6_geometric_hodge.py::test_boundary_squared_vanishes_exactly_in_all_degrees PASSED [ 35%]
tests/test_q6_geometric_hodge.py::test_whitney_mass_matrices_are_symmetric_positive_definite PASSED [ 35%]
tests/test_q6_geometric_hodge.py::test_full_cube_fills_the_111_harmonic_two_cochains_of_the_two_skeleton PASSED [ 35%]
tests/test_q6_geometric_hodge.py::test_dimension_generic_small_cubes[1] PASSED [ 35%]
tests/test_q6_geometric_hodge.py::test_dimension_generic_small_cubes[2] PASSED [ 35%]
tests/test_q6_geometric_hodge.py::test_dimension_generic_small_cubes[3] PASSED [ 35%]
tests/test_q6_geometric_hodge.py::test_dimension_generic_small_cubes[4] PASSED [ 35%]
tests/test_q6_geometric_hodge.py::test_invalid_dimension_rejected PASSED [ 36%]
tests/test_q6_gibbs_dynamics.py::test_q6_gibbs_sampler_returns_finite_trace_and_valid_acceptance PASSED [ 36%]
tests/test_q6_gibbs_dynamics.py::test_q6_gibbs_sampler_is_reproducible_for_fixed_seed PASSED [ 36%]
tests/test_q6_gibbs_dynamics.py::test_sampler_rejects_wrong_link_graph PASSED [ 36%]
tests/test_q6_gibbs_dynamics.py::test_autocorrelation_is_normalized_and_finite PASSED [ 36%]
tests/test_q6_gibbs_dynamics.py::test_integrated_correlation_is_algorithmic_not_physical PASSED [ 36%]
tests/test_q6_gibbs_dynamics.py::test_constant_trace_has_no_defined_normalized_autocorrelation PASSED [ 36%]
tests/test_q6_gibbs_dynamics.py::test_face_count_and_link_count_match_q6_contract PASSED [ 36%]
tests/test_q6_hodge_laplacian.py::test_unit_weights_match_incidence_hodge_formulas PASSED [ 36%]
tests/test_q6_hodge_laplacian.py::test_harmonic_dimensions_match_q6_cellular_betti_numbers PASSED [ 37%]
tests/test_q6_hodge_laplacian.py::test_positive_nonuniform_weights_preserve_harmonic_dimensions PASSED [ 37%]
tests/test_q6_hodge_laplacian.py::test_invalid_weights_are_rejected[kwargs0] PASSED [ 37%]
tests/test_q6_hodge_laplacian.py::test_invalid_weights_are_rejected[kwargs1] PASSED [ 37%]
tests/test_q6_hodge_laplacian.py::test_invalid_weights_are_rejected[kwargs2] PASSED [ 37%]
tests/test_q6_hodge_laplacian.py::test_invalid_weights_are_rejected[kwargs3] PASSED [ 37%]
tests/test_q6_hypercube_gk_supergolden.py::test_hypercube_cosine_identity_and_taylor_expansion PASSED [ 37%]
tests/test_q6_hypercube_gk_supergolden.py::test_closed_cochain_projector_is_orthogonal_rank_129 PASSED [ 37%]
tests/test_q6_hypercube_gk_supergolden.py::test_green_kubo_modal_sum PASSED [ 37%]
tests/test_q6_hypercube_gk_supergolden.py::test_green_kubo_detects_nondecaying_zero_mode PASSED [ 37%]
tests/test_q6_hypercube_gk_supergolden.py::test_green_kubo_zero_mode_harmless_without_stress_projection PASSED [ 38%]
tests/test_q6_hypercube_gk_supergolden.py::test_supergolden_polynomial_and_geometric_series PASSED [ 38%]
tests/test_q6_hypercube_gk_supergolden.py::test_spectral_ratio_audit_is_diagnostic PASSED [ 38%]
tests/test_q6_hypercube_gk_supergolden.py::test_shear_stress_is_negative_energy_derivative PASSED [ 38%]
tests/test_q6_hypercube_gk_supergolden.py::test_unitary_equivalence_transports_all_green_kubo_data PASSED [ 38%]
tests/test_q6_metric_wilson_response.py::test_face_area_log_derivative_matches_finite_difference PASSED [ 38%]
tests/test_q6_metric_wilson_response.py::test_metric_wilson_stress_matches_action_finite_difference PASSED [ 38%]
tests/test_q6_metric_wilson_response.py::test_reference_configuration_has_unit_face_weights PASSED [ 38%]
tests/test_q6_metric_wilson_response.py::test_shear_generator_is_symmetric_traceless PASSED [ 38%]
tests/test_q6_metric_wilson_response.py::test_metric_p1_matter_stress_matches_action_finite_difference PASSED [ 39%]
tests/test_q6_metric_wilson_response.py::test_combined_candidate_stress_is_sum_of_gauge_and_matter_sectors PASSED [ 39%]
tests/test_q6_prism_chain_map.py::test_vertex_map_is_surjective_two_to_one PASSED [ 39%]
tests/test_q6_prism_chain_map.py::test_each_q6_edge_maps_to_a_prism_corner_edge PASSED [ 39%]
tests/test_q6_prism_chain_map.py::test_every_cell_satisfies_chain_map_identity_over_integers PASSED [ 39%]
tests/test_q6_prism_chain_map.py::test_map_is_surjective_on_all_target_cells PASSED [ 39%]
tests/test_q6_prism_chain_map.py::test_cells_with_both_folded_directions_are_collapsed PASSED [ 39%]
tests/test_q6_prism_chain_map.py::test_audit_and_prism_factor_order PASSED [ 39%]
tests/test_q6_prism_embedding_obstruction.py::test_existing_prism_has_five_factors_and_cannot_host_q6_embedding PASSED [ 39%]
tests/test_q6_prism_embedding_obstruction.py::test_six_factor_extension_contains_induced_q6 PASSED [ 39%]
tests/test_q6_prism_embedding_obstruction.py::test_parity_collapse_is_a_homomorphism_but_not_an_embedding PASSED [ 40%]
tests/test_q6_prism_embedding_obstruction.py::test_audit_states_exact_scope PASSED [ 40%]
tests/test_q6_prism_exterior_square_map.py::test_all_240_square_curvatures_are_coordinate_wedges PASSED [ 40%]
tests/test_q6_prism_exterior_square_map.py::test_all_s6_coordinate_permutations_act_naturally_on_wedge_basis PASSED [ 40%]
tests/test_q6_prism_exterior_square_map.py::test_edge_cochain_has_expected_support PASSED [ 40%]
tests/test_q6_prism_exterior_square_map.py::test_audit PASSED            [ 40%]
tests/test_q6_prism_path_valued_homology.py::test_q6_graph_counts_and_fundamental_basis PASSED [ 40%]
tests/test_q6_prism_path_valued_homology.py::test_five_selected_fundamental_cycles_map_to_standard_basis PASSED [ 40%]
tests/test_q6_prism_path_valued_homology.py::test_remaining_fundamental_cycles_map_to_zero PASSED [ 40%]
tests/test_q6_prism_path_valued_homology.py::test_each_source_edge_has_a_based_target_path PASSED [ 41%]
tests/test_q6_prism_path_valued_homology.py::test_selected_paths_are_exact_factor_cycles PASSED [ 41%]
tests/test_q6_prism_path_valued_homology.py::test_audit_rank_and_scope PASSED [ 41%]
tests/test_q6_prism_periodic_obstruction.py::test_all_closed_four_step_walks_have_zero_winding PASSED [ 41%]
tests/test_q6_prism_periodic_obstruction.py::test_no_four_step_walk_can_wind_any_factor PASSED [ 41%]
tests/test_q6_prism_periodic_obstruction.py::test_five_steps_can_wind_the_shortest_factor PASSED [ 41%]
tests/test_q6_prism_periodic_obstruction.py::test_closed_walk_must_close_mod_each_factor PASSED [ 41%]
tests/test_q6_prism_periodic_obstruction.py::test_q6_graph_cycle_rank_and_audit PASSED [ 41%]
tests/test_q6_prism_symmetry_obstruction.py::test_a6_order_exceeds_target_h1_symmetry_image_bound PASSED [ 41%]
tests/test_q6_prism_symmetry_obstruction.py::test_invariant_cycle_recurrence_has_only_zero_solution PASSED [ 41%]
tests/test_q6_prism_symmetry_obstruction.py::test_audit_excludes_equivariant_surjection PASSED [ 42%]
tests/test_q6_spectral_audit.py::test_q6_laplacian_exact_hypercube_spectrum PASSED [ 42%]
tests/test_q6_spectral_audit.py::test_graph_cycle_rank_differs_from_cellular_first_betti PASSED [ 42%]
tests/test_q6_spectral_audit.py::test_antipodal_complement_is_graph_automorphism PASSED [ 42%]
tests/test_q6_spectral_audit.py::test_combined_audit_passes_but_keeps_geometric_identification_open PASSED [ 42%]
tests/test_q6_spectral_audit.py::test_graph_laplacian_is_symmetric_positive_semidefinite PASSED [ 42%]
tests/test_q6_whitney_hodge.py::test_whitney_gram_matrices_are_spd_and_have_expected_blocks PASSED [ 42%]
tests/test_q6_whitney_hodge.py::test_whitney_gram_entries_match_tensor_product_hat_integrals PASSED [ 42%]
tests/test_q6_whitney_hodge.py::test_whitney_hodge_is_symmetric_positive_and_matches_betti_numbers PASSED [ 42%]
tests/test_q6_whitney_hodge.py::test_whitney_hodge_builder_accepts_explicit_metric_object PASSED [ 43%]
tests/test_quantum_dataset_io.py::DatasetPathTests::test_finds_dataset_next_to_module PASSED [ 43%]
tests/test_quantum_dataset_io.py::DatasetPathTests::test_invalid_explicit_path_does_not_fall_back PASSED [ 43%]
tests/test_quantum_dataset_io.py::DatasetPathTests::test_missing_dataset_has_actionable_error PASSED [ 43%]
tests/test_quantum_dataset_io.py::DatasetPathTests::test_prefers_explicit_environment_path PASSED [ 43%]
tests/test_quantum_model_orchestrator.py::OrchestratorTests::test_cleanup_runs_even_if_benchmark_fails PASSED [ 43%]
tests/test_quantum_model_orchestrator.py::OrchestratorTests::test_compose_file_is_valid_json_and_uses_requested_size PASSED [ 43%]
tests/test_quantum_model_orchestrator.py::OrchestratorTests::test_compose_has_64_nodes_and_unprivileged_networking PASSED [ 43%]
tests/test_quantum_model_orchestrator.py::OrchestratorTests::test_total_nodes_must_be_at_least_two PASSED [ 43%]
tests/test_quaternion_algebra.py::test_identity_element PASSED           [ 43%]
tests/test_quaternion_algebra.py::test_zero_element PASSED               [ 44%]
tests/test_quaternion_algebra.py::test_conjugate PASSED                  [ 44%]
tests/test_quaternion_algebra.py::test_conjugate_product_reverses_order PASSED [ 44%]
tests/test_quaternion_algebra.py::test_inverse_left_and_right PASSED     [ 44%]
tests/test_quaternion_algebra.py::test_inverse_of_unit_quaternion_is_conjugate PASSED [ 44%]
tests/test_quaternion_algebra.py::test_zero_quaternion_has_no_inverse PASSED [ 44%]
tests/test_quaternion_algebra.py::test_zero_quaternion_cannot_be_normalized PASSED [ 44%]
tests/test_quaternion_algebra.py::test_norm_squared PASSED               [ 44%]
tests/test_quaternion_algebra.py::test_norm_multiplicativity PASSED      [ 44%]
tests/test_quaternion_algebra.py::test_hamilton_product_is_associative PASSED [ 45%]
tests/test_quaternion_algebra.py::test_hamilton_product_is_not_commutative PASSED [ 45%]
tests/test_quaternion_algebra.py::test_dot_product PASSED                [ 45%]
tests/test_quaternion_algebra.py::test_normalization PASSED              [ 45%]
tests/test_quaternion_algebra.py::test_scalar_multiplication PASSED      [ 45%]
tests/test_quaternion_algebra.py::test_from_iterable_requires_four_components PASSED [ 45%]
tests/test_quaternion_algebra.py::test_non_finite_components_are_rejected PASSED [ 45%]
tests/test_quaternion_algebra.py::test_opposite_lifts_are_distinct_quaternions PASSED [ 45%]
tests/test_quaternion_algebra.py::test_inverse_is_not_the_opposite_lift_of_conjugate PASSED [ 45%]
tests/test_quaternion_algebra.py::test_functional_hamilton_product[q0] PASSED [ 45%]
tests/test_quaternion_algebra.py::test_functional_hamilton_product[q1] PASSED [ 46%]
tests/test_quaternion_algebra.py::test_functional_hamilton_product[q2] PASSED [ 46%]
tests/test_quaternion_algebra.py::test_functional_hamilton_product[q3] PASSED [ 46%]
tests/test_qutrit_aperiodic.py::test_fibonacci_substitution_prefix_and_finite_period_check PASSED [ 46%]
tests/test_qutrit_aperiodic.py::test_qutrit_maximally_entangled_bell_state PASSED [ 46%]
tests/test_qutrit_aperiodic.py::test_fibonacci_qutrit_chain_is_hermitian_and_entangled PASSED [ 46%]
tests/test_qutrit_aperiodic.py::test_combined_qutrit_aperiodic_audit PASSED [ 46%]
tests/test_qutrit_aperiodic.py::test_invalid_chain_size_and_couplings_rejected PASSED [ 46%]
tests/test_qutrit_q6_coupling.py::test_gray_path_is_embedded_in_q6_and_fibonacci_modulated PASSED [ 46%]
tests/test_qutrit_q6_coupling.py::test_qutrit_q6_local_su3_gauge_covariance PASSED [ 47%]
tests/test_qutrit_q6_coupling.py::test_q6_embedded_qutrit_spectral_entanglement_audit PASSED [ 47%]
tests/test_qutrit_q6_coupling.py::test_invalid_q6_embedding_and_link_support_rejected PASSED [ 47%]
tests/test_reconnection_energy.py::test_downhill_event_is_accepted_without_event_cost PASSED [ 47%]
tests/test_reconnection_energy.py::test_event_cost_must_be_paid_by_energy_drop PASSED [ 47%]
tests/test_reconnection_energy.py::test_uphill_event_is_rejected_without_declared_budget PASSED [ 47%]
tests/test_reconnection_energy.py::test_tolerance_is_only_a_comparison_allowance PASSED [ 47%]
tests/test_reconnection_energy.py::test_invalid_energy_inputs_are_rejected[kwargs0] PASSED [ 47%]
tests/test_reconnection_energy.py::test_invalid_energy_inputs_are_rejected[kwargs1] PASSED [ 47%]
tests/test_reconnection_energy.py::test_invalid_energy_inputs_are_rejected[kwargs2] PASSED [ 47%]
tests/test_reconnection_energy.py::test_invalid_energy_inputs_are_rejected[kwargs3] PASSED [ 48%]
tests/test_representations_2i.py::test_fundamental_representation_matches_su2_lift_convention PASSED [ 48%]
tests/test_representations_2i.py::test_three_explicit_representations_are_homomorphisms PASSED [ 48%]
tests/test_representations_2i.py::test_central_minus_one_distinguishes_spinorial_and_vector_representations PASSED [ 48%]
tests/test_representations_2i.py::test_intertwiner_dimensions_match_selected_irreducible_representations PASSED [ 48%]
tests/test_representations_2i.py::test_fundamental_commutant_is_complex_scalar_identity PASSED [ 48%]
tests/test_representations_2i.py::test_invalid_group_element_shape_is_rejected PASSED [ 48%]
tests/test_rgs_d5_scale_identifiability.py::test_fixed_kernel_length_breaks_naive_scale_covariance PASSED [ 48%]
tests/test_rgs_d5_scale_identifiability.py::test_coframe_is_invariant_when_all_length_scales_are_scaled_together PASSED [ 48%]
tests/test_rgs_d5_scale_identifiability.py::test_current_q6_embedding_cannot_induce_a_positive_definite_d5_metric PASSED [ 49%]
tests/test_rgs_d5_scale_identifiability.py::test_scale_probe_must_be_finite_positive[0.0] PASSED [ 49%]
tests/test_rgs_d5_scale_identifiability.py::test_scale_probe_must_be_finite_positive[-1.0] PASSED [ 49%]
tests/test_rgs_d5_scale_identifiability.py::test_scale_probe_must_be_finite_positive[inf] PASSED [ 49%]
tests/test_rgs_d5_scale_identifiability.py::test_scale_probe_must_be_finite_positive[nan] PASSED [ 49%]
tests/test_rgs_ec_algebraic_torsion.py::test_minus_plus_gamma_matrices_obey_clifford_algebra PASSED [ 49%]
tests/test_rgs_ec_algebraic_torsion.py::test_dirac_axial_current_is_real_and_spinor_derived PASSED [ 49%]
tests/test_rgs_ec_algebraic_torsion.py::test_action_derived_cartan_coefficient_and_spinor_to_torsion_pipeline PASSED [ 49%]
tests/test_rgs_ec_algebraic_torsion.py::test_torsion_invariant_finite_and_zero_for_zero_current PASSED [ 49%]
tests/test_rgs_ec_algebraic_torsion.py::test_missing_convention_or_units_are_rejected[kwargs0] PASSED [ 49%]
tests/test_rgs_ec_algebraic_torsion.py::test_missing_convention_or_units_are_rejected[kwargs1] PASSED [ 50%]
tests/test_rgs_ec_contact_tetrad_variation.py::test_fixed_state_tetrad_variation_gives_metric_times_contact_lagrangian PASSED [ 50%]
tests/test_rgs_ec_contact_tetrad_variation.py::test_spacelike_second_moment_has_explicit_vacuum_form_stress PASSED [ 50%]
tests/test_rgs_ec_contact_tetrad_variation.py::test_contact_tetrad_variation_rejects_invalid_inputs PASSED [ 50%]
tests/test_rgs_ec_contact_tetrad_variation.py::test_tetrad_variation_matches_central_finite_difference PASSED [ 50%]
tests/test_rgs_ec_contact_tetrad_variation.py::test_explicit_hilbert_tensor_is_metric_times_contact_lagrangian PASSED [ 50%]
tests/test_rgs_ec_contact_tetrad_variation.py::test_lorentz_scalar_is_frame_invariant_under_spatial_rotation PASSED [ 50%]
tests/test_rgs_ec_contact_tetrad_variation.py::test_rejects_singular_or_wrong_orientation_tetrads PASSED [ 50%]
tests/test_rgs_ec_contact_tetrad_variation.py::test_spin_ensemble_isotropy_contract_exposes_time_space_covariance PASSED [ 50%]
tests/test_rgs_ec_contorsion.py::test_contorsion_reconstructs_totally_antisymmetric_torsion PASSED [ 50%]
tests/test_rgs_ec_contorsion.py::test_quadratic_torsion_invariant_matches_axial_current_identity PASSED [ 51%]
tests/test_rgs_ec_contorsion.py::test_reject_non_antisymmetric_torsion_for_dirac_shortcut PASSED [ 51%]
tests/test_rgs_ec_dirac_four_fermion.py::test_raw_and_real_current_conventions_have_same_density PASSED [ 51%]
tests/test_rgs_ec_dirac_four_fermion.py::test_spinor_current_to_torsion_and_four_fermion_pipeline PASSED [ 51%]
tests/test_rgs_ec_dirac_four_fermion.py::test_effective_source_is_gradient_of_density PASSED [ 51%]
tests/test_rgs_ec_dirac_four_fermion.py::test_multiflavour_current_is_summed_before_squaring PASSED [ 51%]
tests/test_rgs_ec_dirac_four_fermion.py::test_raw_real_round_trip_and_bad_kappa PASSED [ 51%]
tests/test_rgs_ec_dirac_kinetic_eom.py::test_kinetic_euler_operator_is_gamma_contraction_with_action_sign PASSED [ 51%]
tests/test_rgs_ec_dirac_kinetic_eom.py::test_kinetic_operator_rejects_invalid_derivative_shape PASSED [ 51%]
tests/test_rgs_ec_dirac_kinetic_eom.py::test_full_local_ec_dirac_residual_adds_contact_source PASSED [ 52%]
tests/test_rgs_ec_dirac_kinetic_eom.py::test_zero_derivative_and_zero_spinor_give_zero_local_residual PASSED [ 52%]
tests/test_rgs_ec_dirac_kinetic_eom.py::test_kinetic_euler_operator_is_linear_in_covariant_derivative PASSED [ 52%]
tests/test_rgs_ec_full_dirac_stress.py::test_bilinear_kinetic_hilbert_tensor_is_symmetric_and_real PASSED [ 52%]
tests/test_rgs_ec_full_dirac_stress.py::test_bilinear_kinetic_rejects_nonreal_physical_tensor PASSED [ 52%]
tests/test_rgs_ec_full_dirac_stress.py::test_spinor_builder_zero_derivative_has_zero_kinetic_tensor PASSED [ 52%]
tests/test_rgs_ec_full_dirac_stress.py::test_contact_tensor_matches_fixed_state_variation PASSED [ 52%]
tests/test_rgs_ec_full_dirac_stress.py::test_full_tensor_adds_kinetic_and_contact_without_claiming_closure PASSED [ 52%]
tests/test_rgs_ec_full_dirac_stress.py::test_pure_spinor_second_moment_is_outer_product_of_axial_current PASSED [ 52%]
tests/test_rgs_ec_nonlinear_dirac_eom.py::test_contact_sources_match_independent_barpsi_finite_difference PASSED [ 52%]
tests/test_rgs_ec_nonlinear_dirac_eom.py::test_contact_sources_match_independent_psi_finite_difference PASSED [ 53%]
tests/test_rgs_ec_nonlinear_dirac_eom.py::test_nonlinear_eom_residual_adds_contact_source_and_does_not_claim_conservation PASSED [ 53%]
tests/test_rgs_ec_nonlinear_dirac_eom.py::test_zero_current_has_zero_contact_euler_source PASSED [ 53%]
tests/test_rgs_ec_spin_ensemble.py::test_cartan_normalized_source_reproduces_action_derived_torsion PASSED [ 53%]
tests/test_rgs_ec_spin_ensemble.py::test_unpolarized_pair_has_zero_mean_but_nonzero_spin_fluctuation_contact PASSED [ 53%]
tests/test_rgs_ec_spin_ensemble.py::test_isotropic_spatial_covariance_diagnostic PASSED [ 53%]
tests/test_rgs_ec_spin_ensemble.py::test_time_space_covariance_is_connected_second_moment PASSED [ 53%]
tests/test_rgs_ec_spin_ensemble.py::test_weighted_moments_and_flavour_averaging PASSED [ 53%]
tests/test_rgs_ec_spin_ensemble.py::test_adapter_refuses_implicit_units_and_does_not_invent_beta PASSED [ 53%]
tests/test_rgs_ec_spin_ensemble.py::test_rejects_bad_weights_and_nonpositive_scale_factor PASSED [ 54%]
tests/test_rgs_ec_spin_ensemble.py::test_flrw_isotropy_uses_full_second_moment_and_mean_polarization PASSED [ 54%]
tests/test_rgs_ellis_refinement.py::test_recovered_rgs_seed_and_q6_coframe_are_reconstructed PASSED [ 54%]
tests/test_rgs_ellis_refinement.py::test_ellis_embedding_induces_the_target_warped_product_metric PASSED [ 54%]
tests/test_rgs_ellis_refinement.py::test_curved_mesh_chain_complex_and_p1_operators PASSED [ 54%]
tests/test_rgs_ellis_refinement.py::test_curved_rgs_refinement_converges_in_volume_and_diameter PASSED [ 54%]
tests/test_rgs_ellis_refinement.py::test_curved_refinement_rejects_invalid_parameters_and_levels PASSED [ 54%]
tests/test_rgs_ellis_refinement.py::test_curved_g4_common_l2_embedding_and_exact_metric_mass PASSED [ 54%]
tests/test_rgs_ellis_refinement.py::test_curved_g4_embedding_rejects_invalid_nodal_or_barycentric_data PASSED [ 54%]
tests/test_rgs_ellis_refinement.py::test_curved_lumped_mass_is_spectrally_comparable_to_common_l2_mass PASSED [ 54%]
tests/test_rgs_ellis_spin_connection.py::test_ellis_analytic_coframe_reconstructs_metric_and_jet PASSED [ 55%]
tests/test_rgs_ellis_spin_connection.py::test_ellis_spin_connection_is_torsion_free_and_metric_compatible PASSED [ 55%]
tests/test_rgs_ellis_spin_connection.py::test_ellis_chart_rejects_coordinate_axis_singularity PASSED [ 55%]
tests/test_rgs_fourth_moment_bridge.py::test_recovered_seed_has_second_moment_80_identity PASSED [ 55%]
tests/test_rgs_fourth_moment_bridge.py::test_fourth_moment_projection_is_trace_free_and_orthogonal PASSED [ 55%]
tests/test_rgs_fourth_moment_bridge.py::test_fourth_moment_has_nonzero_cubic_harmonic PASSED [ 55%]
tests/test_rgs_fourth_moment_bridge.py::test_anisotropic_response_splits_Eg_and_T2g_without_mixing PASSED [ 55%]
tests/test_rgs_fourth_moment_bridge.py::test_full_fourth_moment_bridge_audit_passes PASSED [ 55%]
tests/test_rgs_fourth_moment_connection.py::test_shifted_fourth_moment_metric_is_positive PASSED [ 55%]
tests/test_rgs_fourth_moment_connection.py::test_covariant_derivative_remains_symmetric_tracefree PASSED [ 56%]
tests/test_rgs_fourth_moment_connection.py::test_gradient_energy_is_nonnegative PASSED [ 56%]
tests/test_rgs_fourth_moment_connection.py::test_spin_current_is_antisymmetric_and_matches_connection_variation PASSED [ 56%]
tests/test_rgs_fourth_moment_connection.py::test_zero_gradient_has_zero_spin_current PASSED [ 56%]
tests/test_rgs_fourth_moment_maxwell_bridge.py::test_direct_contraction_with_antisymmetric_field_is_zero PASSED [ 56%]
tests/test_rgs_fourth_moment_maxwell_bridge.py::test_induced_response_is_stf_and_positive PASSED [ 56%]
tests/test_rgs_fourth_moment_maxwell_bridge.py::test_nonpositive_constitutive_matrix_is_rejected PASSED [ 56%]
tests/test_rgs_fourth_moment_maxwell_bridge.py::test_tensor_covariance_under_orthogonal_frame_change PASSED [ 56%]
tests/test_rgs_fourth_moment_variation.py::test_coframe_variation_matches_central_difference PASSED [ 56%]
tests/test_rgs_fourth_moment_variation.py::test_lapse_and_q_variations PASSED [ 56%]
tests/test_rgs_fourth_moment_variation.py::test_simultaneous_internal_tensor_rotation_invariance PASSED [ 57%]
tests/test_rgs_fourth_moment_variation.py::test_complete_variational_audit PASSED [ 57%]
tests/test_rgs_icosahedral_q6_intertwiner.py::test_proper_icosahedral_group_has_60_signed_axis_actions PASSED [ 57%]
tests/test_rgs_icosahedral_q6_intertwiner.py::test_signed_axis_action_preserves_c6_and_is_group_representation PASSED [ 57%]
tests/test_rgs_icosahedral_q6_intertwiner.py::test_signed_icosahedral_actions_are_q6_graph_automorphisms PASSED [ 57%]
tests/test_rgs_icosahedral_q6_intertwiner.py::test_q6_vertex_action_rejects_invalid_inputs PASSED [ 57%]
tests/test_rgs_multiscale_relational_audit.py::test_q6_exact_graph_and_spectrum PASSED [ 57%]
tests/test_rgs_multiscale_relational_audit.py::test_finite_relational_clock_constraint PASSED [ 57%]
tests/test_rgs_multiscale_relational_audit.py::test_markov_equilibrium_and_driven_cycle PASSED [ 57%]
tests/test_rgs_multiscale_relational_audit.py::test_acoustic_metric_and_pure_vortex PASSED [ 58%]
tests/test_rgs_multiscale_relational_audit.py::test_discrete_scale_toy_is_not_misrepresented_as_paper_replication PASSED [ 58%]
tests/test_rgs_multiscale_relational_audit.py::test_dimensionless_gravity_scale_diagnostic PASSED [ 58%]
tests/test_rgs_physics_architecture.py::test_core_metric_clifford_and_q6_conventions_are_compatible PASSED [ 58%]
tests/test_rgs_physics_architecture.py::test_architecture_registry_keeps_distinct_physical_sectors_separate PASSED [ 58%]
tests/test_rgs_physics_architecture.py::test_architecture_audit_does_not_claim_physical_unification PASSED [ 58%]
tests/test_rgs_physics_architecture.py::test_audit_rejects_invalid_tolerance PASSED [ 58%]
tests/test_rgs_psu_antitwin_qball_bvp.py::test_two_field_zero_gauge_limit_recovers_mass_shifted_global_qball PASSED [ 58%]
tests/test_rgs_psu_antitwin_qball_bvp.py::test_two_field_gauged_branch_solves_both_scalar_equations_and_gauss_law PASSED [ 58%]
tests/test_rgs_psu_antitwin_qball_bvp.py::test_imbalance_penalty_changes_relative_curvature_not_symmetric_background PASSED [ 58%]
tests/test_rgs_psu_antitwin_qball_bvp.py::test_two_field_qball_rejects_invalid_frequency_window_and_parameters PASSED [ 59%]
tests/test_rgs_psu_core_parse_and_moments.py::test_rgs_psu_core_sources_parse PASSED [ 59%]
tests/test_rgs_psu_core_parse_and_moments.py::test_600_cell_coordinate_count_and_unit_norms PASSED [ 59%]
tests/test_rgs_psu_core_parse_and_moments.py::test_600_cell_finite_moments_match_s3_haar PASSED [ 59%]
tests/test_rgs_psu_core_parse_and_moments.py::test_600_cell_quotient_claim_is_not_misidentified PASSED [ 59%]
tests/test_rgs_psu_cut_project.py::test_golden_conjugates_and_orthogonal_three_plus_three_split PASSED [ 59%]
tests/test_rgs_psu_cut_project.py::test_projection_matrix_is_the_recovered_icosahedral_convention PASSED [ 59%]
tests/test_rgs_psu_cut_project.py::test_projected_cube_is_the_rhombic_triacontahedral_window PASSED [ 59%]
tests/test_rgs_psu_cut_project.py::test_finite_model_set_patch_is_inside_both_cuts PASSED [ 59%]
tests/test_rgs_psu_cut_project.py::test_torsion_coupled_internal_time_axis_is_recoverable_and_stable PASSED [ 60%]
tests/test_rgs_psu_cut_project.py::test_time_axis_rejects_missing_or_constant_torsion_data PASSED [ 60%]
tests/test_rgs_psu_cut_project.py::test_model_set_window_scaling_changes_the_finite_selection PASSED [ 60%]
tests/test_rgs_psu_cut_project_window.py::test_rhombic_triacontahedral_window_has_30_rhombic_facets PASSED [ 60%]
tests/test_rgs_psu_cut_project_window.py::test_q6_edges_project_to_equal_length_icosahedral_vectors PASSED [ 60%]
tests/test_rgs_psu_fibred_intertwiner.py::test_q6_base_sign_fibre_reconstructs_archived_rgs_seed_exactly PASSED [ 60%]
tests/test_rgs_psu_fibred_intertwiner.py::test_pullback_edge_metric_and_graph_fem_are_positive PASSED [ 60%]
tests/test_rgs_psu_hilbert_qhd_mhd.py::test_q6_hilbert_curve_is_bijective PASSED [ 60%]
tests/test_rgs_psu_hilbert_qhd_mhd.py::test_cantor_cylinder_measure_is_normalized_at_each_level PASSED [ 60%]
tests/test_rgs_psu_hilbert_qhd_mhd.py::test_u1_covariant_difference_is_gauge_covariant PASSED [ 60%]
tests/test_rgs_psu_hilbert_qhd_mhd.py::test_monopole_chern_flux_quantization PASSED [ 61%]
tests/test_rgs_psu_hilbert_qhd_mhd.py::test_integrated_periodic_benchmark_residuals PASSED [ 61%]
tests/test_rgs_psu_hilbert_qhd_mhd.py::test_fluxon_winding_is_quantized_and_gauge_invariant PASSED [ 61%]
tests/test_rgs_psu_hilbert_qhd_mhd.py::test_sine_gordon_fluxon_has_one_flux_quantum_and_small_equation_residual PASSED [ 61%]
tests/test_rgs_psu_hilbert_qhd_mhd.py::test_fluxon_charge_is_not_declared_quantized_with_nontrivial_background_holonomy PASSED [ 61%]
tests/test_rgs_psu_hilbert_qhd_mhd.py::test_qball_frequency_window_matches_scalar_potential_minimum PASSED [ 61%]
tests/test_rgs_psu_hilbert_qhd_mhd.py::test_qball_gaussian_trial_has_finite_noether_charge_energy_and_residual PASSED [ 61%]
tests/test_rgs_psu_intertwiner.py::test_archived_rgs_and_golden_q6_are_not_affinely_equivalent PASSED [ 61%]
tests/test_rgs_psu_intertwiner.py::test_rgs_golden_no_go_is_translation_and_rotation_invariant PASSED [ 61%]
tests/test_rgs_psu_multifractal_bundle.py::test_seed_induced_rank_three_connection_is_metric_compatible_and_flat PASSED [ 62%]
tests/test_rgs_psu_multifractal_bundle.py::test_local_gl3_gauge_transform_preserves_flatness_and_metric_compatibility PASSED [ 62%]
tests/test_rgs_psu_multifractal_bundle.py::test_uniform_nested_measure_has_one_dimensional_generalized_dimensions PASSED [ 62%]
tests/test_rgs_psu_multifractal_bundle.py::test_binomial_cascade_has_nonconstant_generalized_dimensions PASSED [ 62%]
tests/test_rgs_psu_qball_bvp.py::test_stationary_qball_bvp_converges_to_nontrivial_regular_profile PASSED [ 62%]
tests/test_rgs_psu_qball_bvp.py::test_stationary_qball_frequency_and_energy_converge_under_mesh_refinement PASSED [ 62%]
tests/test_rgs_psu_qball_bvp.py::test_stationary_qball_rejects_potential_without_positive_frequency_window PASSED [ 62%]
tests/test_rgs_psu_qball_maxwell_bvp.py::test_zero_gauge_charge_recovers_global_qball_limit PASSED [ 62%]
tests/test_rgs_psu_qball_maxwell_bvp.py::test_gauged_qball_solves_scalar_and_maxwell_equations_and_gauss_law PASSED [ 62%]
tests/test_rgs_psu_qball_maxwell_bvp.py::test_gauged_qball_rejects_invalid_inputs PASSED [ 62%]
tests/test_rgs_q6_a4_stf_intertwiner.py::test_q6_a4_stf_linear_intertwiner PASSED [ 63%]
tests/test_rgs_q6_a5_character_theorem.py::test_q6_top_coexact_five_modes_are_irreducible_a5_five PASSED [ 63%]
tests/test_rgs_q6_clifford_module.py::test_full_q6_cube_has_all_cell_dimensions_and_boundary_squared_zero PASSED [ 63%]
tests/test_rgs_q6_clifford_module.py::test_six_dimensional_exterior_clifford_relations_and_chirality PASSED [ 63%]
tests/test_rgs_q6_clifford_module.py::test_complete_clifford_audit PASSED [ 63%]
tests/test_rgs_q6_compact_wilson_fiber_hessian.py::test_compact_wilson_fibre_hessian_matches_finite_difference PASSED [ 63%]
tests/test_rgs_q6_compact_wilson_fiber_hessian.py::test_single_fibre_point_reduces_to_cosine_hessian PASSED [ 63%]
tests/test_rgs_q6_cubical_refinement.py::test_m1_recovers_the_q6_cell_counts_and_lambda_12 PASSED [ 63%]
tests/test_rgs_q6_cubical_refinement.py::test_fundamental_joint_block_has_exact_coexact_rank_five[1] PASSED [ 63%]
tests/test_rgs_q6_cubical_refinement.py::test_fundamental_joint_block_has_exact_coexact_rank_five[2] PASSED [ 64%]
tests/test_rgs_q6_cubical_refinement.py::test_fundamental_joint_block_has_exact_coexact_rank_five[3] PASSED [ 64%]
tests/test_rgs_q6_cubical_refinement.py::test_fundamental_joint_block_has_exact_coexact_rank_five[4] PASSED [ 64%]
tests/test_rgs_q6_cubical_refinement.py::test_fundamental_joint_block_has_exact_coexact_rank_five[8] PASSED [ 64%]
tests/test_rgs_q6_cubical_refinement.py::test_fundamental_joint_block_has_exact_coexact_rank_five[16] PASSED [ 64%]
tests/test_rgs_q6_cubical_refinement.py::test_fundamental_joint_block_has_exact_coexact_rank_five[32] PASSED [ 64%]
tests/test_rgs_q6_cubical_refinement.py::test_fundamental_joint_block_has_exact_coexact_rank_five[64] PASSED [ 64%]
tests/test_rgs_q6_cubical_refinement.py::test_rescaled_fundamental_eigenvalue_converges_to_six_pi_squared PASSED [ 64%]
tests/test_rgs_q6_cubical_refinement.py::test_helmert_basis_is_the_five_dimensional_sum_zero_space PASSED [ 64%]
tests/test_rgs_q6_cubical_refinement.py::test_interval_spectrum_and_cell_counts_validate_inputs PASSED [ 64%]
tests/test_rgs_q6_d5_a5_angular_measure.py::test_a5_character_average_has_no_low_degree_angular_harmonics PASSED [ 65%]
tests/test_rgs_q6_d5_a5_angular_measure.py::test_so3_spherical_harmonic_character_identity_and_validation PASSED [ 65%]
tests/test_rgs_q6_d5_a5_angular_measure.py::test_q6_d5_phase_shell_pushforward_has_a5_symmetry_theorem PASSED [ 65%]
tests/test_rgs_q6_d5_a5_degree_six.py::test_reynolds_degree_six_harmonic_is_a5_invariant PASSED [ 65%]
tests/test_rgs_q6_d5_a5_degree_six.py::test_uniform_sphere_harmonic_rms_is_positive PASSED [ 65%]
tests/test_rgs_q6_d5_a5_degree_six.py::test_degree_six_coefficient_audit_is_reproducible_and_finite PASSED [ 65%]
tests/test_rgs_q6_d5_degree_six_sobol.py::test_sobol_direction_batch_is_unit_and_reproducible PASSED [ 65%]
tests/test_rgs_q6_d5_degree_six_sobol.py::test_sobol_crosscheck_is_reproducible_and_finite PASSED [ 65%]
tests/test_rgs_q6_d5_degree_six_sobol.py::test_sobol_crosscheck_rejects_invalid_controls[kwargs0] PASSED [ 65%]
tests/test_rgs_q6_d5_degree_six_sobol.py::test_sobol_crosscheck_rejects_invalid_controls[kwargs1] PASSED [ 66%]
tests/test_rgs_q6_d5_degree_six_sobol.py::test_sobol_crosscheck_rejects_invalid_controls[kwargs2] PASSED [ 66%]
tests/test_rgs_q6_d5_fiber_bessel_majorant.py::test_factorial_bessel_majorant_is_finite_and_tail_is_geometric PASSED [ 66%]
tests/test_rgs_q6_d5_fiber_bessel_majorant.py::test_representative_changes_only_the_counting_bound PASSED [ 66%]
tests/test_rgs_q6_d5_hua_cutoff.py::test_exact_factorized_hua_density_and_endpoints PASSED [ 66%]
tests/test_rgs_q6_d5_hua_cutoff.py::test_density_mellin_transform_is_hua_ratio[-0.9] PASSED [ 66%]
tests/test_rgs_q6_d5_hua_cutoff.py::test_density_mellin_transform_is_hua_ratio[-0.5] PASSED [ 66%]
tests/test_rgs_q6_d5_hua_cutoff.py::test_density_mellin_transform_is_hua_ratio[0.0] PASSED [ 66%]
tests/test_rgs_q6_d5_hua_cutoff.py::test_density_mellin_transform_is_hua_ratio[1.0] PASSED [ 66%]
tests/test_rgs_q6_d5_hua_cutoff.py::test_density_mellin_transform_is_hua_ratio[2.0] PASSED [ 66%]
tests/test_rgs_q6_d5_hua_cutoff.py::test_density_mellin_transform_is_hua_ratio[5.0] PASSED [ 67%]
tests/test_rgs_q6_d5_hua_cutoff.py::test_exact_cutoff_normalizer_matches_quadrature[0.5] PASSED [ 67%]
tests/test_rgs_q6_d5_hua_cutoff.py::test_exact_cutoff_normalizer_matches_quadrature[0.2] PASSED [ 67%]
tests/test_rgs_q6_d5_hua_cutoff.py::test_exact_cutoff_normalizer_matches_quadrature[0.1] PASSED [ 67%]
tests/test_rgs_q6_d5_hua_cutoff.py::test_exact_cutoff_normalizer_matches_quadrature[0.05] PASSED [ 67%]
tests/test_rgs_q6_d5_hua_cutoff.py::test_exact_cutoff_normalizer_matches_quadrature[0.01] PASSED [ 67%]
tests/test_rgs_q6_d5_hua_cutoff.py::test_cutoff_measure_is_finite_and_normalized PASSED [ 67%]
tests/test_rgs_q6_d5_hua_cutoff.py::test_cutoff_asymptotic_and_report PASSED [ 67%]
tests/test_rgs_q6_d5_hua_cutoff.py::test_invalid_cutoff_rejected[0.0] PASSED [ 67%]
tests/test_rgs_q6_d5_hua_cutoff.py::test_invalid_cutoff_rejected[1.0] PASSED [ 68%]
tests/test_rgs_q6_d5_hua_cutoff.py::test_invalid_cutoff_rejected[-0.1] PASSED [ 68%]
tests/test_rgs_q6_d5_hua_cutoff.py::test_invalid_cutoff_rejected[nan] PASSED [ 68%]
tests/test_rgs_q6_d5_hua_cutoff.py::test_invalid_cutoff_rejected[inf] PASSED [ 68%]
tests/test_rgs_q6_d5_hua_normalization.py::test_type_iv5_hua_polynomial_and_normalized_integral PASSED [ 68%]
tests/test_rgs_q6_d5_hua_normalization.py::test_hua_integral_positive_and_finite_above_threshold[-0.999] PASSED [ 68%]
tests/test_rgs_q6_d5_hua_normalization.py::test_hua_integral_positive_and_finite_above_threshold[-0.9] PASSED [ 68%]
tests/test_rgs_q6_d5_hua_normalization.py::test_hua_integral_positive_and_finite_above_threshold[-0.5] PASSED [ 68%]
tests/test_rgs_q6_d5_hua_normalization.py::test_hua_integral_positive_and_finite_above_threshold[0.0] PASSED [ 68%]
tests/test_rgs_q6_d5_hua_normalization.py::test_hua_integral_positive_and_finite_above_threshold[1.0] PASSED [ 68%]
tests/test_rgs_q6_d5_hua_normalization.py::test_hua_integral_positive_and_finite_above_threshold[4.0] PASSED [ 69%]
tests/test_rgs_q6_d5_hua_normalization.py::test_hua_integral_rejects_nonconvergent_or_nonfinite_exponent[-1.0] PASSED [ 69%]
tests/test_rgs_q6_d5_hua_normalization.py::test_hua_integral_rejects_nonconvergent_or_nonfinite_exponent[-2.0] PASSED [ 69%]
tests/test_rgs_q6_d5_hua_normalization.py::test_hua_integral_rejects_nonconvergent_or_nonfinite_exponent[-5.0] PASSED [ 69%]
tests/test_rgs_q6_d5_hua_normalization.py::test_hua_integral_rejects_nonconvergent_or_nonfinite_exponent[nan] PASSED [ 69%]
tests/test_rgs_q6_d5_hua_normalization.py::test_hua_integral_rejects_nonconvergent_or_nonfinite_exponent[inf] PASSED [ 69%]
tests/test_rgs_q6_d5_hua_normalization.py::test_hua_boundary_residue_and_bergman_divergence PASSED [ 69%]
tests/test_rgs_q6_d5_hua_normalization.py::test_hua_mellin_transform_is_product_of_beta_factors PASSED [ 69%]
tests/test_rgs_q6_d5_hua_normalization.py::test_cutoff_bergman_divergence_has_exact_leading_coefficient PASSED [ 69%]
tests/test_rgs_q6_d5_integer_fiber_quotient.py::test_exact_rational_mode_map_and_integer_kernel PASSED [ 70%]
tests/test_rgs_q6_d5_integer_fiber_quotient.py::test_kernel_translation_changes_pointwise_wilson_weight PASSED [ 70%]
tests/test_rgs_q6_d5_integer_fiber_quotient.py::test_partial_fiber_weight_is_positive_and_translation_consistent PASSED [ 70%]
tests/test_rgs_q6_d5_integer_fiber_quotient.py::test_invalid_fiber_inputs PASSED [ 70%]
tests/test_rgs_q6_d5_liouville_measure.py::test_origin_metric_and_volume_match_liouville PASSED [ 70%]
tests/test_rgs_q6_d5_liouville_measure.py::test_exact_density_and_maxwell_hamiltonian_law[0.2-12.0] PASSED [ 70%]
tests/test_rgs_q6_d5_liouville_measure.py::test_exact_density_and_maxwell_hamiltonian_law[1.0-12.0] PASSED [ 70%]
tests/test_rgs_q6_d5_liouville_measure.py::test_exact_density_and_maxwell_hamiltonian_law[3.7-8.0] PASSED [ 70%]
tests/test_rgs_q6_d5_liouville_measure.py::test_exact_density_and_maxwell_hamiltonian_law[11.0-4.0] PASSED [ 70%]
tests/test_rgs_q6_d5_liouville_measure.py::test_type_iv_determinant_formula PASSED [ 70%]
tests/test_rgs_q6_d5_liouville_measure.py::test_exact_hua_cutoff_volume_in_liouville_units PASSED [ 71%]
tests/test_rgs_q6_d5_liouville_measure.py::test_invalid_inputs_rejected PASSED [ 71%]
tests/test_rgs_q6_d5_liouville_measure.py::test_maxwell_flow_is_common_phase_and_preserves_bergman_liouville_weight[0.2-12.0] PASSED [ 71%]
tests/test_rgs_q6_d5_liouville_measure.py::test_maxwell_flow_is_common_phase_and_preserves_bergman_liouville_weight[1.0-12.0] PASSED [ 71%]
tests/test_rgs_q6_d5_liouville_measure.py::test_maxwell_flow_is_common_phase_and_preserves_bergman_liouville_weight[3.7-8.0] PASSED [ 71%]
tests/test_rgs_q6_d5_liouville_measure.py::test_maxwell_flow_is_common_phase_and_preserves_bergman_liouville_weight[11.0-4.0] PASSED [ 71%]
tests/test_rgs_q6_d5_liouville_measure.py::test_chart_uses_maxwell_canonical_reciprocal_scale PASSED [ 71%]
tests/test_rgs_q6_d5_maxwell_hamiltonian.py::test_maxwell_legendre_coefficients_and_beta_independent_frequency PASSED [ 71%]
tests/test_rgs_q6_d5_maxwell_hamiltonian.py::test_maxwell_normalized_d5_hamiltonian_audit_passes PASSED [ 71%]
tests/test_rgs_q6_d5_maxwell_hamiltonian.py::test_d5_chart_energy_is_hamiltonian_for_multiple_beta PASSED [ 72%]
tests/test_rgs_q6_d5_maxwell_hamiltonian.py::test_invalid_maxwell_parameters_rejected[0.0-12.0] PASSED [ 72%]
tests/test_rgs_q6_d5_maxwell_hamiltonian.py::test_invalid_maxwell_parameters_rejected[-1.0-12.0] PASSED [ 72%]
tests/test_rgs_q6_d5_maxwell_hamiltonian.py::test_invalid_maxwell_parameters_rejected[1.0-0.0] PASSED [ 72%]
tests/test_rgs_q6_d5_phase_space.py::test_q6_top_coexact_sector_is_five_dimensional_A5_irrep PASSED [ 72%]
tests/test_rgs_q6_d5_phase_space.py::test_cotangent_lift_is_ten_dimensional_kahler_and_symplectic PASSED [ 72%]
tests/test_rgs_q6_d5_phase_space.py::test_type_iv_domain_contains_origin_and_small_tangent_vectors PASSED [ 72%]
tests/test_rgs_q6_d5_phase_space.py::test_phase_space_to_d5_tangent_has_expected_complex_dimension PASSED [ 72%]
tests/test_rgs_q6_d5_phase_space.py::test_invalid_beta_is_rejected[0.0] PASSED [ 72%]
tests/test_rgs_q6_d5_phase_space.py::test_invalid_beta_is_rejected[-1.0] PASSED [ 72%]
tests/test_rgs_q6_d5_phase_space.py::test_invalid_beta_is_rejected[inf] PASSED [ 73%]
tests/test_rgs_q6_d5_phase_space.py::test_invalid_beta_is_rejected[nan] PASSED [ 73%]
tests/test_rgs_q6_d5_saturated_kernel_bounds.py::test_smith_hnf_returns_saturated_integer_kernel PASSED [ 73%]
tests/test_rgs_q6_d5_saturated_kernel_bounds.py::test_certified_factorial_tail_bound_is_finite_and_decreases PASSED [ 73%]
tests/test_rgs_q6_d5_saturated_kernel_bounds.py::test_tail_bound_handles_invalid_inputs PASSED [ 73%]
tests/test_rgs_q6_d5_saturated_kernel_bounds.py::test_combined_audit_reports_proved_scope PASSED [ 73%]
tests/test_rgs_q6_d5_sphere_intertwiner.py::test_exterior_square_contains_geometric_three_dimensional_intertwiner PASSED [ 73%]
tests/test_rgs_q6_d5_sphere_intertwiner.py::test_cotangent_wedge_direction_is_A5_equivariant PASSED [ 73%]
tests/test_rgs_q6_d5_sphere_intertwiner.py::test_cotangent_direction_is_invariant_under_oscillator_SO2_phase PASSED [ 73%]
tests/test_rgs_q6_d5_sphere_intertwiner.py::test_collinear_cotangent_pair_has_zero_projected_direction_or_is_rejected_cleanly PASSED [ 74%]
tests/test_rgs_q6_d5_sphere_intertwiner.py::test_direction_map_is_defined_for_generic_q_and_has_codimension_three_zero_fibre PASSED [ 74%]
tests/test_rgs_q6_d5_wilson_fiber_torus.py::test_smith_parameterization_tracks_disconnected_annihilator PASSED [ 74%]
tests/test_rgs_q6_d5_wilson_fiber_torus.py::test_annihilator_integral_reproduces_exact_bessel_fiber_sums PASSED [ 74%]
tests/test_rgs_q6_d5_wilson_fiber_torus.py::test_smith_parameterization_for_primitive_lattice PASSED [ 74%]
tests/test_rgs_q6_d5_wilson_fiber_torus.py::test_invalid_lattice_rank_rejected PASSED [ 74%]
tests/test_rgs_q6_d5_wilson_fiber_torus.py::test_scaled_wilson_integrand_factorization_is_exact_and_bounded[0.1] PASSED [ 74%]
tests/test_rgs_q6_d5_wilson_fiber_torus.py::test_scaled_wilson_integrand_factorization_is_exact_and_bounded[1.25] PASSED [ 74%]
tests/test_rgs_q6_d5_wilson_fiber_torus.py::test_scaled_wilson_integrand_factorization_is_exact_and_bounded[2.0] PASSED [ 74%]
tests/test_rgs_q6_d5_wilson_fiber_torus.py::test_scaled_wilson_integrand_factorization_is_exact_and_bounded[20.0] PASSED [ 74%]
tests/test_rgs_q6_d5_wilson_regulator.py::test_exact_transfer_spectral_cutoff_equivalence PASSED [ 75%]
tests/test_rgs_q6_d5_wilson_regulator.py::test_q6_d5_projection_obstruction_audit PASSED [ 75%]
tests/test_rgs_q6_d5_wilson_regulator.py::test_invalid_spectral_cutoff_parameters[0.0-0.05-0.2] PASSED [ 75%]
tests/test_rgs_q6_d5_wilson_regulator.py::test_invalid_spectral_cutoff_parameters[80.0-0.0-0.2] PASSED [ 75%]
tests/test_rgs_q6_d5_wilson_regulator.py::test_invalid_spectral_cutoff_parameters[80.0-0.05-0.0] PASSED [ 75%]
tests/test_rgs_q6_d5_wilson_regulator.py::test_invalid_spectral_cutoff_parameters[80.0-0.05-1.1] PASSED [ 75%]
tests/test_rgs_q6_d5_wilson_regulator.py::test_invalid_spectral_cutoff_parameters[nan-0.05-0.2] PASSED [ 75%]
tests/test_rgs_q6_equivariant_cochains.py::test_q6_cubical_cochain_dimensions_and_boundary_squared_zero PASSED [ 75%]
tests/test_rgs_q6_equivariant_cochains.py::test_icosahedral_action_lifts_to_orthogonal_cochain_representations PASSED [ 75%]
tests/test_rgs_q6_equivariant_cochains.py::test_full_equivariant_cochain_audit PASSED [ 76%]
tests/test_rgs_q6_euclidean_time_wilson.py::test_periodic_time_cell_counts_and_boundary_identity PASSED [ 76%]
tests/test_rgs_q6_euclidean_time_wilson.py::test_plaquettes_and_action_are_gauge_invariant PASSED [ 76%]
tests/test_rgs_q6_euclidean_time_wilson.py::test_temporal_link_variation_is_discrete_gauss_law PASSED [ 76%]
tests/test_rgs_q6_euclidean_time_wilson.py::test_full_euclidean_time_gauss_audit PASSED [ 76%]
tests/test_rgs_q6_euclidean_time_wilson.py::test_invalid_time_extent_rejected[0] PASSED [ 76%]
tests/test_rgs_q6_euclidean_time_wilson.py::test_invalid_time_extent_rejected[1] PASSED [ 76%]
tests/test_rgs_q6_euclidean_time_wilson.py::test_invalid_time_extent_rejected[1.5] PASSED [ 76%]
tests/test_rgs_q6_gauged_full_clifford.py::test_full_cubical_gauged_differential_has_all_degrees_and_is_hermitian PASSED [ 76%]
tests/test_rgs_q6_gauged_full_clifford.py::test_flatness_and_curvature_obstruction_to_covariant_d_squared PASSED [ 76%]
tests/test_rgs_q6_gauged_full_clifford.py::test_gauged_full_clifford_transport_audit PASSED [ 77%]
tests/test_rgs_q6_gauged_full_clifford.py::test_covariant_nilpotency_is_equivalent_to_trivial_q6_plaquette_holonomy PASSED [ 77%]
tests/test_rgs_q6_hodge_spectral_theorem.py::test_q6_coexact_hodge_spectrum_matches_product_complex_theorem PASSED [ 77%]
tests/test_rgs_q6_induced_metric.py::test_rgs_q6_metric_uses_hash_verified_seed_and_canonical_cells PASSED [ 77%]
tests/test_rgs_q6_induced_metric.py::test_local_six_dimensional_pullback_is_explicitly_degenerate PASSED [ 77%]
tests/test_rgs_q6_induced_metric.py::test_face_gram_and_quadrilateral_area_are_consistent PASSED [ 77%]
tests/test_rgs_q6_induced_metric.py::test_metric_weighted_graph_p1_mass_and_stiffness PASSED [ 77%]
tests/test_rgs_q6_induced_metric.py::test_audit_keeps_coordinate_oh_separate_from_q6_equivariance PASSED [ 77%]
tests/test_rgs_q6_integer_a5_d5_measure.py::test_exact_rational_q6_mode_space_has_five_integer_flux_generators PASSED [ 77%]
tests/test_rgs_q6_integer_a5_d5_measure.py::test_full_a5_action_preserves_integer_q6_mode_intersection PASSED [ 78%]
tests/test_rgs_q6_integer_a5_d5_measure.py::test_a5_preserves_d5_type_iv_norm_and_bergman_density PASSED [ 78%]
tests/test_rgs_q6_integer_a5_d5_measure.py::test_d5_bergman_density_rejects_points_outside_domain PASSED [ 78%]
tests/test_rgs_q6_integral_gauge_quotient.py::test_q6_integral_gauge_quotient_has_unimodular_five_mode_complement PASSED [ 78%]
tests/test_rgs_q6_oh_action.py::test_signed_permutation_group_has_order_48 PASSED [ 78%]
tests/test_rgs_q6_oh_action.py::test_all_oh_matrices_act_as_unique_coordinate_permutations PASSED [ 78%]
tests/test_rgs_q6_oh_action.py::test_induced_oh_action_is_faithful_and_respects_group_product PASSED [ 78%]
tests/test_rgs_q6_oh_action.py::test_q6_equivariance_is_reported_separately_from_coordinate_oh PASSED [ 78%]
tests/test_rgs_q6_relative_sphere.py::test_q6_vertices_map_equivariantly_to_finite_sphere_design PASSED [ 78%]
tests/test_rgs_q6_relative_sphere.py::test_relative_direction_is_unit_and_vertex_indexed PASSED [ 78%]
tests/test_rgs_q6_relative_sphere.py::test_32_point_direction_set_is_spherical_five_design_but_not_six PASSED [ 79%]
tests/test_rgs_q6_relative_sphere.py::test_icosahedral_signed_axis_frame_is_tight_and_four_radial_shells PASSED [ 79%]
tests/test_rgs_q6_spin_glass.py::test_q6_quenched_bonds_are_reproducible_and_have_expected_shape PASSED [ 79%]
tests/test_rgs_q6_spin_glass.py::test_frustration_and_energy_are_z2_gauge_invariant[gaussian] PASSED [ 79%]
tests/test_rgs_q6_spin_glass.py::test_frustration_and_energy_are_z2_gauge_invariant[bimodal] PASSED [ 79%]
tests/test_rgs_q6_spin_glass.py::test_single_spin_flip_energy_increment_matches_hamiltonian PASSED [ 79%]
tests/test_rgs_q6_spin_glass.py::test_exact_ground_state_works_on_small_induced_graph PASSED [ 79%]
tests/test_rgs_q6_spin_glass.py::test_exact_ground_state_refuses_full_q6_enumeration PASSED [ 79%]
tests/test_rgs_q6_spin_glass.py::test_parallel_tempering_returns_finite_overlap_observables PASSED [ 79%]
tests/test_rgs_q6_su2_connection.py::test_random_links_are_su2_and_invalid_links_rejected PASSED [ 80%]
tests/test_rgs_q6_su2_connection.py::test_covariant_incidence_and_dirac_kaehler_hermiticity PASSED [ 80%]
tests/test_rgs_q6_su2_connection.py::test_gauge_covariance_and_wilson_loop_invariance PASSED [ 80%]
tests/test_rgs_q6_su2_connection.py::test_full_su2_connection_audit PASSED [ 80%]
tests/test_rgs_q6_su2_wilson.py::test_recovered_symmetry_subgroups_have_expected_orders PASSED [ 80%]
tests/test_rgs_q6_su2_wilson.py::test_connection_and_plaquettes_have_canonical_dimensions PASSED [ 80%]
tests/test_rgs_q6_su2_wilson.py::test_wilson_bridge_audit_passes PASSED  [ 80%]
tests/test_rgs_q6_u1_character_electric_flux.py::test_character_expansion_reconstructs_wilson_boltzmann_factor PASSED [ 80%]
tests/test_rgs_q6_u1_character_electric_flux.py::test_bessel_coefficients_are_positive_and_symmetric_in_integer_order PASSED [ 80%]
tests/test_rgs_q6_u1_character_electric_flux.py::test_fundamental_integer_cycle_basis_has_exact_gauss_law PASSED [ 80%]
tests/test_rgs_q6_u1_character_electric_flux.py::test_integer_gauss_law_rejects_noninteger_flux PASSED [ 81%]
tests/test_rgs_q6_u1_character_electric_flux.py::test_integer_gauss_law_accepts_cycle_flux PASSED [ 81%]
tests/test_rgs_q6_u1_character_electric_flux.py::test_character_audit_passes PASSED [ 81%]
tests/test_rgs_q6_u1_character_electric_flux.py::test_invalid_beta_rejected[0.0] PASSED [ 81%]
tests/test_rgs_q6_u1_character_electric_flux.py::test_invalid_beta_rejected[-1.0] PASSED [ 81%]
tests/test_rgs_q6_u1_character_electric_flux.py::test_invalid_beta_rejected[nan] PASSED [ 81%]
tests/test_rgs_q6_u1_character_electric_flux.py::test_invalid_fourier_cutoff_rejected[-1] PASSED [ 81%]
tests/test_rgs_q6_u1_character_electric_flux.py::test_invalid_fourier_cutoff_rejected[1.5] PASSED [ 81%]
tests/test_rgs_q6_u1_maxwell_matching.py::test_wilson_projected_hamiltonian_matches_d5_maxwell_coefficients[1.0-0.05-1.0-6] PASSED [ 81%]
tests/test_rgs_q6_u1_maxwell_matching.py::test_wilson_projected_hamiltonian_matches_d5_maxwell_coefficients[0.2-0.01-0.003-6] PASSED [ 82%]
tests/test_rgs_q6_u1_maxwell_matching.py::test_wilson_projected_hamiltonian_matches_d5_maxwell_coefficients[0.7-0.03-0.4-3] PASSED [ 82%]
tests/test_rgs_q6_u1_maxwell_matching.py::test_wilson_projected_hamiltonian_matches_d5_maxwell_coefficients[1.3-0.2-2.1-4] PASSED [ 82%]
tests/test_rgs_q6_u1_maxwell_matching.py::test_matching_explicitly_does_not_claim_integer_lattice_embedding PASSED [ 82%]
tests/test_rgs_q6_u1_maxwell_matching.py::test_matching_rejects_invalid_parameters[kwargs0] PASSED [ 82%]
tests/test_rgs_q6_u1_maxwell_matching.py::test_matching_rejects_invalid_parameters[kwargs1] PASSED [ 82%]
tests/test_rgs_q6_u1_maxwell_matching.py::test_matching_rejects_invalid_parameters[kwargs2] PASSED [ 82%]
tests/test_rgs_q6_u1_maxwell_matching.py::test_matching_rejects_invalid_parameters[kwargs3] PASSED [ 82%]
tests/test_rgs_q6_u1_quotient.py::test_q6_quotient_dimensions_and_rank PASSED [ 82%]
tests/test_rgs_q6_u1_quotient.py::test_tree_gauge_fix_is_orbit_invariant PASSED [ 82%]
tests/test_rgs_q6_u1_quotient.py::test_gauge_and_icosahedral_audits PASSED [ 83%]
tests/test_rgs_q6_u1_quotient.py::test_wilson_action_has_quadratic_small_field_limit PASSED [ 83%]
tests/test_rgs_q6_u1_quotient.py::test_normalized_haar_coordinates_and_lift PASSED [ 83%]
tests/test_rgs_q6_u1_quotient.py::test_invalid_beta_rejected[0.0] PASSED [ 83%]
tests/test_rgs_q6_u1_quotient.py::test_invalid_beta_rejected[-1.0] PASSED [ 83%]
tests/test_rgs_q6_u1_quotient.py::test_invalid_beta_rejected[inf] PASSED [ 83%]
tests/test_rgs_q6_u1_quotient.py::test_invalid_beta_rejected[nan] PASSED [ 83%]
tests/test_rgs_q6_u1_transfer_operator.py::test_q6_integer_flux_transfer_is_diagonal_positive_and_gauss_restricted PASSED [ 83%]
tests/test_rgs_q6_u1_transfer_operator.py::test_electric_hamiltonian_approaches_quadratic_flux_energy PASSED [ 83%]
tests/test_rgs_q6_u1_transfer_operator.py::test_q6_spatial_wilson_potential_is_gauge_invariant PASSED [ 84%]
tests/test_rgs_q6_u1_transfer_operator.py::test_symmetric_transfer_kernel_is_hermitian PASSED [ 84%]
tests/test_rgs_q6_u1_transfer_operator.py::test_anisotropic_couplings_have_correct_d3_and_d6_scaling PASSED [ 84%]
tests/test_rgs_q6_u1_transfer_operator.py::test_transfer_operator_audit_passes PASSED [ 84%]
tests/test_rgs_q6_u1_transfer_operator.py::test_transfer_inputs_reject_nonpositive_couplings PASSED [ 84%]
tests/test_rgs_q6_wilson_fiber_hessian.py::test_q6_wilson_quadratic_fibre_marginalization PASSED [ 84%]
tests/test_rgs_q6_wilson_microhistory_modal_bridge.py::test_exact_orbit_stabilizer_and_factorial_weight_identity PASSED [ 84%]
tests/test_rgs_q6_wilson_microhistory_modal_bridge.py::test_signed_tokens_keep_orientation_in_the_quotient PASSED [ 84%]
tests/test_rgs_q6_wilson_microhistory_modal_bridge.py::test_orbit_frame_is_s5_and_box_means_invariance_on_orbit PASSED [ 84%]
tests/test_rgs_q6_wilson_microhistory_modal_bridge.py::test_audit_and_validation PASSED [ 84%]
tests/test_rgs_smilga_two_particle_audit.py::test_transform_inverse_and_jacobian PASSED [ 85%]
tests/test_rgs_smilga_two_particle_audit.py::test_mass_shell_and_casimir_are_lorentz_covariant PASSED [ 85%]
tests/test_rgs_smilga_two_particle_audit.py::test_constrained_mass_shell_measure_jacobians PASSED [ 85%]
tests/test_rgs_smilga_two_particle_audit.py::test_two_body_phase_space_normalization PASSED [ 85%]
tests/test_rgs_smilga_two_particle_audit.py::test_wyler_numeric_is_separate_from_phase_space_claim PASSED [ 85%]
tests/test_rgs_smilga_two_particle_audit.py::test_fixed_mass_shell_chart_has_five_independent_coordinates PASSED [ 85%]
tests/test_rgs_spin_connection.py::test_constant_minkowski_coframe_has_zero_spin_connection PASSED [ 85%]
tests/test_rgs_spin_connection.py::test_spherical_minkowski_coframe_has_torsion_free_connection PASSED [ 85%]
tests/test_rgs_spin_connection.py::test_spinor_derivative_uses_spin_connection_and_frame_conversion PASSED [ 85%]
tests/test_rgs_spin_connection.py::test_spin_connection_rejects_singular_coframe PASSED [ 86%]
tests/test_rgs_spin_connection.py::test_spinor_derivative_rejects_wrong_shapes PASSED [ 86%]
tests/test_rgs_spin_connection.py::test_frame_derivative_feeds_existing_ec_dirac_euler_operator PASSED [ 86%]
tests/test_rgs_unified_multiscale_operator.py::test_unified_product_operator_dimensions_and_edges PASSED [ 86%]
tests/test_rgs_unified_multiscale_operator.py::test_generator_is_markov_and_uniform_stationary PASSED [ 86%]
tests/test_rgs_unified_multiscale_operator.py::test_relational_clock_commutes_with_unified_generator PASSED [ 86%]
tests/test_rgs_unified_multiscale_operator.py::test_reversible_and_driven_entropy_production_controls PASSED [ 86%]
tests/test_rgs_unified_multiscale_operator.py::test_zero_drive_recovers_reversible_limit PASSED [ 86%]
tests/test_rgs_unified_multiscale_operator.py::test_generator_rejects_negative_rates PASSED [ 86%]
tests/test_rgs_vierbein_reconstruction.py::test_reconstructs_archived_64_site_seed_and_positive_spatial_frames PASSED [ 86%]
tests/test_rgs_vierbein_reconstruction.py::test_static_lorentzian_lift_reproduces_site_metrics PASSED [ 87%]
tests/test_rgs_vierbein_reconstruction.py::test_ellis_tetrad_reproduces_wormhole_metric PASSED [ 87%]
tests/test_rgs_vierbein_reconstruction.py::test_kerr_schild_tetrad_reproduces_exact_kerr_schild_metric PASSED [ 87%]
tests/test_rgs_weyl_c6_q6.py::test_directional_gradient_has_six_hypercube_edges_per_vertex PASSED [ 87%]
tests/test_rgs_weyl_c6_q6.py::test_positive_c6_constitutive_matrix_uses_offdiagonal_channels PASSED [ 87%]
tests/test_rgs_weyl_c6_q6.py::test_q6_operator_is_symmetric_positive_and_has_constant_kernel PASSED [ 87%]
tests/test_rgs_weyl_c6_q6.py::test_energy_is_nonnegative_and_zero_on_constants PASSED [ 87%]
tests/test_rgs_weyl_c6_q6.py::test_axis_relabelling_is_covariant_when_mapping_is_relabelled PASSED [ 87%]
tests/test_rgs_weyl_c6_q6.py::test_invalid_axis_map_and_nonpositive_parameters_rejected PASSED [ 87%]
tests/test_rgs_weyl_phase_closure.py::test_massless_weyl_helicity_projectors[momentum0] PASSED [ 88%]
tests/test_rgs_weyl_phase_closure.py::test_massless_weyl_helicity_projectors[momentum1] PASSED [ 88%]
tests/test_rgs_weyl_phase_closure.py::test_massless_weyl_helicity_projectors[momentum2] PASSED [ 88%]
tests/test_rgs_weyl_phase_closure.py::test_helicity_is_undefined_at_zero_momentum PASSED [ 88%]
tests/test_rgs_weyl_phase_closure.py::test_icosahedral_c6_quadratic_identity_and_zero_diagonal PASSED [ 88%]
tests/test_rgs_weyl_phase_closure.py::test_original_q6_weight_rule_is_uniform_not_c6_anisotropic PASSED [ 88%]
tests/test_rgs_wyler_domain_audit.py::test_wyler_unit_radius_reproduces_published_expression PASSED [ 88%]
tests/test_rgs_wyler_domain_audit.py::test_euclidean_radius_scaling_is_dimensionally_consistent PASSED [ 88%]
tests/test_rgs_wyler_domain_audit.py::test_codata_matching_radius_is_not_misreported_as_prediction PASSED [ 88%]
tests/test_rgs_wyler_domain_audit.py::test_radius_must_be_finite_and_positive[0.0] PASSED [ 88%]
tests/test_rgs_wyler_domain_audit.py::test_radius_must_be_finite_and_positive[-1.0] PASSED [ 89%]
tests/test_rgs_wyler_domain_audit.py::test_radius_must_be_finite_and_positive[inf] PASSED [ 89%]
tests/test_rgs_wyler_domain_audit.py::test_radius_must_be_finite_and_positive[nan] PASSED [ 89%]
tests/test_shared_io.py::SharedIOTests::test_atomic_write_replaces_complete_file PASSED [ 89%]
tests/test_shared_io.py::SharedIOTests::test_download_checks_size_and_sha_before_publish PASSED [ 89%]
tests/test_shared_io.py::SharedIOTests::test_download_mismatch_preserves_existing_file PASSED [ 89%]
tests/test_shared_io.py::SharedIOTests::test_gzip_download_is_validated_before_publish PASSED [ 89%]
tests/test_shared_io.py::SharedIOTests::test_zip_extracts_selected_file PASSED [ 89%]
tests/test_shared_io.py::SharedIOTests::test_zip_size_limit_is_enforced PASSED [ 89%]
tests/test_shared_io.py::SharedIOTests::test_zip_traversal_is_rejected PASSED [ 90%]
tests/test_unified_discrete_gradients.py::test_full_gradients_match_full_action_directional_difference PASSED [ 96%]
tests/test_unified_strange_loop_operator.py::test_operator_is_symmetric_psd_and_constant_mode_kernel PASSED [ 96%]
tests/test_unified_strange_loop_operator.py::test_torus_gray_embedding_and_lift PASSED [ 96%]
tests/test_unified_strange_loop_operator.py::test_unified_audit_and_prism_constant_intertwiner PASSED [ 97%]
tests/test_unified_strange_loop_operator.py::test_invalid_weights PASSED [ 97%]
tests/test_unified_strange_loop_operator.py::test_cantor_emissions_are_coupled_to_parry_q6_and_stochastic PASSED [ 97%]
tests/test_unified_strange_loop_operator.py::test_doppler_moire_phase_separation_and_relativistic_factor PASSED [ 97%]
tests/test_unified_strange_loop_operator.py::test_doppler_moire_extended_operator_is_hermitian PASSED [ 97%]
tests/test_unified_strange_loop_operator.py::test_barycenter_first_moment_covariance_and_exact_moire_identity PASSED [ 97%]
tests/test_unified_strange_loop_operator.py::test_barycentric_quadratic_expansion_converges_for_small_wavevector_gap PASSED [ 97%]
tests/test_unified_strange_loop_operator.py::test_moving_barycenter_includes_weight_transport PASSED [ 97%]
tests/test_unified_strange_loop_operator.py::test_barycenter_rejects_indefinite_metric_and_negative_weights PASSED [ 97%]
tests/test_unified_strange_loop_operator.py::test_q6_projection_recovers_affine_rgs_seed_exactly PASSED [ 98%]
tests/test_unified_strange_loop_operator.py::test_q6_projection_reports_non_affine_geometry_instead_of_hiding_it PASSED [ 98%]
tests/test_unified_strange_loop_operator.py::test_q6_projection_rejects_wrong_seed_shape_and_bad_weights PASSED [ 98%]
tests/test_wilson_periodic_symbol.py::test_three_dimensional_wilson_corners_have_expected_masses PASSED [ 98%]
tests/test_wilson_periodic_symbol.py::test_corner_shell_multiplicities_are_binomial PASSED [ 98%]
tests/test_wilson_periodic_symbol.py::test_fixed_momentum_symbol_has_continuum_kinetic_limit PASSED [ 98%]
tests/test_wilson_periodic_symbol.py::test_joint_scaling_window[-0.5] PASSED [ 98%]
tests/test_wilson_periodic_symbol.py::test_joint_scaling_window[0.0] PASSED [ 98%]
tests/test_wilson_periodic_symbol.py::test_joint_scaling_window[0.7] PASSED [ 98%]
tests/test_wilson_periodic_symbol.py::test_scaling_window_boundaries_are_strict[-1.0] PASSED [ 98%]
tests/test_wilson_periodic_symbol.py::test_scaling_window_boundaries_are_strict[1.0] PASSED [ 99%]
tests/test_wilson_periodic_symbol.py::test_scaling_window_boundaries_are_strict[1.5] PASSED [ 99%]
tests/test_wilson_periodic_symbol.py::test_invalid_symbol_parameters_rejected PASSED [ 99%]
tests/test_wilson_refinement.py::test_fixed_physical_momentum_has_second_order_naive_dispersion_error PASSED [ 99%]
tests/test_wilson_refinement.py::test_three_dimensional_wilson_corners_lift_all_seven_doublers PASSED [ 99%]
tests/test_wilson_refinement.py::test_wilson_parameter_scaling_window_is_minus_one_to_one PASSED [ 99%]
tests/test_wilson_refinement.py::test_wilson_symbol_has_expected_small_momentum_and_corner_limits PASSED [ 99%]
tests/test_wilson_refinement.py::test_invalid_refinement_inputs_are_rejected PASSED [ 99%]
tests/test_wilson_refinement.py::test_hermitian_scalar_sum_has_an_accidental_non_nyquist_zero PASSED [ 99%]
tests/test_wilson_refinement.py::test_euclidean_wilson_symbol_has_no_accidental_zero_at_that_angle PASSED [100%]
============== 951 passed, 6 subtests passed in 199.16s (0:03:19) ==============
0s
Run python examples/lattice_rotation_demo.py
====================================================================
SU(2) LATTICE HOLONOMY DEMONSTRATION
====================================================================
Plaquette path: [0, 1, 2, 0]
Holonomy:
[[ 0.88025765-0.32746999j -0.34085644-0.04155423j]
 [ 0.34085644-0.04155423j  0.88025765+0.32746999j]]
Holonomy checks:
  SU(2): True
  Unitarity residual: 3.140203083200158e-16
  Determinant: (1.0000000000000002+5.551115123125784e-17j)
Normalized Wilson trace: 0.8802576528784465
Gauge covariance:
  Holonomy conjugation residual: 4.25032032655802e-16
  Original Wilson trace: 0.8802576528784465
  Transformed Wilson trace: 0.8802576528784467
  Wilson invariance residual: 2.220446049250313e-16
====================================================================
ALL DEMONSTRATION CHECKS PASSED
====================================================================
