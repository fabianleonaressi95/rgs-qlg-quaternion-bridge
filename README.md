# RGS–QLG Quaternion Bridge

**Research codebase for quaternionic geometry, discrete lattice structures, spectral operators and reproducible computational validation.**

This repository is part of the RGS–PSU–QLG research programme. It develops explicit mathematical constructions and computational tests involving unit quaternions, SU(2)/SO(3), discrete incidence geometry, lattice holonomies, symmetry representations and Dirac-type operators.

> **Evidence policy:** mathematical identities, numerical verification, candidate constructions and open physical hypotheses are reported separately. Passing a finite numerical test does not establish a continuum limit or a physical theory.


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

The [intrinsic Ellis-RGS spinorial Wilson-Dirac operator](docs/ELLIS_RGS_INTRINSIC_SPIN_DIRAC.md) assembles the naive spatial Dirac block and positive Wilson graph term from local tangent Clifford directions, intrinsic Ellis edge lengths, lumped P1 vertex volumes and the Levi-Civita SU(2) links. Its finite-mesh audits include gamma5 hermiticity and local spin-frame gauge covariance; curved refinement convergence and doubler control remain OPEN. Run `PYTHONPATH=. python -m pytest tests/test_ellis_rgs_spin_dirac.py -q`.

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

## Relationship to the central THEOREM research repository

This repository is the implementation- and experiment-oriented QLG/quaternionic component of the broader **RGS–PSU–PRISM–Q6** research programme. Its constructions, scripts, numerical audits and open mathematical questions are connected to the programme's central theorem and research index:

- **[THEOREM — RGS–PSU–PRISM–Q6 Mathematical Research Programme](https://github.com/fabianleonaressi95/theorem)** — central repository for theorem records, mathematical assumptions and proofs, reproducible audits, model definitions, bibliography, and the programme-wide register of open problems.
- [Central THEOREM README](https://github.com/fabianleonaressi95/theorem/blob/main/README.md) — overview of the research strands, theorem register, evidence policy and current priorities.
- [QLG–Dirac Lorentz-chirality result (T049D)](https://github.com/fabianleonaressi95/theorem/blob/main/theorems/T049D-lorentz-chiral-qlg-embedding.md) — Lorentz-covariant Weyl embeddings and the fact that Lorentz covariance alone does not select a chirality sign, when this record is available on the central repository's default branch.

The two repositories have complementary roles: **this repository develops and audits concrete quaternionic, Q6 and QLG constructions; the theorem repository organizes the corresponding mathematical statements, proof obligations, evidence levels and cross-project results.** A computational result in this repository should therefore be cited with its implementation path and commit, while any associated theorem-level claim should point to the relevant record in the central repository.

The central repository does not automatically certify every experiment or physical interpretation in this codebase. Each connection remains subject to the assumptions, proof status and limitations stated in the corresponding record.

The central landing page has now been updated on `main` (README revision: [45841e0](https://github.com/fabianleonaressi95/theorem/commit/45841e0)). The update provides the programme-wide English overview and synchronizes the theorem/result index through T049D, including the QLG–Dirac chirality distinction, alongside the Einstein–Cartan, Mosco-convergence and cellular-dynamics research strands. Start with the [updated central README](https://github.com/fabianleonaressi95/theorem/blob/main/README.md) for the current programme map; use this repository for the implementation-level material and its reproducible audits.

