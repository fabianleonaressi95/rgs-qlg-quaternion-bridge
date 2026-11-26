# RGS–QLG Quaternion Bridge

**Implementation and reproducibility repository** for the quaternionic, Q6, RGS, lattice-gauge and finite-operator components of the RGS–PSU–QLG programme.

This repository develops finite mathematical constructions and computational audits within the RGS–PSU–QLG programme. It covers quaternionic groups, Q6 cubical geometry, RGS-induced metrics, finite-element operators, lattice gauge covariance, spectral representation theory and selected continuum-limit questions.

> **Evidence policy:** PROVED, FINITE-EXACT, NUMERICAL, CONDITIONAL, CANDIDATE and OPEN are distinct statuses. A finite test is not a continuum theorem or a physical derivation.

## Canonical suite formalization

The 51-step sequence formerly expanded in this README is now maintained as a dedicated mathematical register:

- **[Integrated mathematical suite — 51 formal records](docs/INTEGRATED_MATHEMATICAL_SUITE_FORMALIZATION_2026-11-25.md)** — definitions, precise claims, evidence status, dependencies and proof obligations for every step from quaternionic foundations through Q6, RGS/FEM, qutrit networks, Kathara, Q6–PRISM, A5/De Aquino, shear/Green–Kubo and Pellis.
- [Unified theorem dependency graph](docs/UNIFIED_THEOREM_DEPENDENCY_GRAPH.md) and [machine-readable edge registry](docs/UNIFIED_THEOREM_EDGES.json).
- [Master catalogue of laws and theorems](docs/MASTER_CATALOGUE_LAWS_THEOREMS.md), [Phase 6–10 addendum](docs/MASTER_CATALOGUE_PHASE6_10_ADDENDUM.md) and [machine-readable catalogue](results/catalogue/laws_theorems_registry.json).
- [27-chapter programme crosswalk](https://github.com/fabianleonaressi95/theorem/blob/main/research/UNIFIED_PROGRAM_27_CHAPTER_CROSSWALK_2026-11-24.md).

The 51 records are not one theorem. Their dependencies are typed, and open physical bridges remain open until a concrete map/operator and proof are supplied.

## Main mathematical entry points

### Quaternionic and representation theory
- [SU(2), 2I and Q6 intertwiner compatibility](docs/Q6_2I_INTERTWINER_COMPATIBILITY.md)
- [RGS icosahedral–Q6 equivariant intertwiner](docs/RGS_ICOSAHEDRAL_Q6_EQUIVARIANT_INTERTWINER.md)
- [Pellis A — spectral geometry of S³/2I](docs/PELLIS_A_2I_SPECTRAL_GEOMETRY.md)
- [Pellis B — two-sheet self-adjoint operators](docs/PELLIS_B_TWO_SHEET_SELF_ADJOINT_OPERATORS.md)
- [Pellis C — D6 channel intertwiner](docs/PELLIS_C_D6_CHANNEL_INTERTWINER.md)
- [Pellis D — RGS–PSU–PRISM integration](docs/PELLIS_D_RGS_PSU_PRISM_INTEGRATION.md)

### Q6, cochains and geometric operators
- Q6 has 64 vertices, 192 edges and 240 square faces. The graph cycle rank is 129; the filled square complex has Betti numbers (1,0,111).
- [Q6 complex formalization](docs/Q6_COMPLEX_FORMALIZATION.md)
- [Weighted combinatorial Hodge Laplacians](docs/Q6_COMBINATORIAL_HODGE_LAPLACIANS.md)
- [RGS–Q6 induced metric and graph P1 audit](docs/RGS_Q6_INDUCED_METRIC.md)
- [Barycentric P1 FEM](docs/Q6_BARYCENTRIC_P1_FEM.md)
- [Barycentric covariance/quadrature gate](docs/GATE_E1_BARYCENTRIC_COVARIANCE_QUADRATURE.md)

### Ellis–RGS refinement and Mosco convergence
The flat periodic P1 Mosco theorem is proved under its stated hypotheses. The curved Ellis–RGS Mosco problem is a separate open problem.
- [Curved refinement G1](docs/RGS_ELLIS_CURVED_REFINEMENT_G1.md)
- [Curved common-L² embedding G4](docs/RGS_ELLIS_CURVED_G4_COMMON_L2.md)
- [Curved Mosco proof programme](docs/RGS_MOSCO_CURVED_REFINEMENT_PROGRAM.md)
- [Flat periodic P1 Mosco proof](docs/FLAT_TORUS_G6_P1_MOSCO_PROOF.md)

### Q6–PRISM, Kathara and A5/GEC-6
- [Kathara graph and cellular audits](docs/KATHARA_15_LINE_GRAPH_FACE_AUDIT.md)
- [Kathara–Q6 chain-map gate C](docs/KATHARA_Q6_CHAIN_MAP_GATE_C.md)
- [Kathara nondegenerate gate D](docs/KATHARA_Q6_NONDEGENERATE_GATE_D.md)
- [Q6–PRISM embedding obstruction](docs/research/Q6_PRISM_EMBEDDING_OBSTRUCTION.md)
- [Q6–PRISM path-valued homology](docs/research/Q6_PRISM_PATH_VALUED_HOMOLOGY.md)
- [GEC-6 A5 spectrum for m=3](docs/research/GEC6_A5_LAPLACIAN_SPECTRUM_M3.md)
- [A5 continuum torus and De Aquino bridge](docs/research/GEC6_A5_CONTINUUM_TORUS_SPECTRUM_AND_DEAQUINO_BRIDGE.md)

### Gradients, shear and Green–Kubo
- [Action-derived Q6 shear response](docs/Q6_ACTION_DERIVED_SHEAR_RESPONSE.md)
- [Metric Wilson response](docs/Q6_METRIC_WILSON_RESPONSE.md)
- [Gibbs stress correlations](docs/Q6_GIBBS_STRESS_CORRELATIONS.md)
- [Q6 hypercube, finite Green–Kubo and supergolden scaling](docs/Q6_HYPERCUBE_GREEN_KUBO_SUPERGOLDEN.md)

The finite OU Green–Kubo integral is an exact result for the declared stochastic process. A physical viscosity still requires a physical-time dynamics, equilibrium/FDT identification and suitable limits.

## Focused open problem for external review: curved Ellis–RGS Mosco convergence

We invite independent proof, counterexample or missing-hypothesis analysis for the scalar Riemannian P1 forms on
M_L=[−L,L]×S², with g=du²+(u²+r₀²)dΩ² and the nested curved tetrahedral refinement family specified in the linked documents.

The target is q(u)=∫|∇_g u|²dμ_g on the common space H=L²(M_L,dμ_g). The required result is both weak liminf and strong recovery for the discrete forms q_h via the stated reconstruction J_h. Uniform shape regularity, Jacobian/metric bounds, coercivity and compactness remain proof obligations. Finite-level audits and the flat theorem do not settle this curved problem.

Reproduce the finite-level checks:
\`\`\`bash
PYTHONPATH=. python examples/rgs_ellis_curved_refinement_audit.py
PYTHONPATH=. python -m pytest tests/test_rgs_ellis_refinement.py -v
\`\`\`

This is a scalar Riemannian convergence question, not a claim about Lorentzian dynamics, spinorial convergence, teleparallel field equations, Yang–Mills or QCD.

## Reproducibility and evidence requirements

Every benchmark or claim must identify:
1. mathematical object, conventions and hypotheses;
2. exact repository revision and input provenance;
3. implementation, tests and expected invariant;
4. numerical tolerance/residual where applicable;
5. whether validation is exact, floating-point, sampled or independently implemented;
6. limitations and unresolved assumptions.

Status labels:
- **PROVED** — mathematical proof under stated hypotheses.
- **FINITE-EXACT** — exact finite algebra/combinatorics with reproducible data.
- **NUMERICAL** — computational evidence with tolerances.
- **CONDITIONAL** — valid only under explicit assumptions/ansatz.
- **CANDIDATE** — proposed construction not fully established.
- **OPEN** — outstanding proof or identification.

## Repository navigation

- src/ — mathematical and computational implementations.
- tests/ — unit, exact-identity and structural tests.
- [experiments/README.md](experiments/README.md) — experiment inventory and provenance rules.
- [Experiment dependency audit](docs/EXPERIMENT_DEPENDENCY_GRAPH_AUDIT.md).
- [Research priorities and collaboration map](docs/RESEARCH_PRIORITIES_AND_COLLABORATION_MAP.md).
- [Q6–D5 literature and repository evidence map](docs/Q6_D5_LITERATURE_AND_REPOSITORY_EVIDENCE_MAP.md).
- [Implementation roadmap](docs/IMPLEMENTATION_ROADMAP.md).
- [RGS–STEM interoperability profile](docs/STEM_INTEROPERABILITY_PROFILE.md) — proposal for discussion, not an approved Ontos extension or Stem Framework endorsement.
- [CI failure audit](docs/CI_FAILURE_AUDIT_2026-11-09.md).

The precise modules and callable functions used for a benchmark must be identified by repository commit SHA and recorded in its provenance manifest.
This repository contains source code, focused mathematical specifications, tests and reproducible finite audits. The programme-wide theorem register and cross-project proof architecture are maintained in [THEOREM](https://github.com/fabianleonaressi95/theorem).

> **Evidence policy:** exact proofs, finite computational verification, conditional constructions, candidates and open physical questions are distinct. A passing finite test does not establish a continuum limit or a physical theory.

## Canonical mathematical suite

The 51-step suite is now formalized in **[Unified Mathematical Suite: Formal Statements, Evidence and Dependencies](docs/UNIFIED_MATHEMATICAL_SUITE_FORMALIZATION.md)**.

That document is the canonical index for:
- quaternion and representation foundations;
- Q6 graph, square complex, homology and Hodge theory;
- recovered RGS seed, induced geometry and graph FEM;
- flat and curved Mosco problems;
- qutrit/SU(3), Q3–Q6 and network constructions;
- Kathara, Krystal Spiral, barycentric quadrature and cubical refinement;
- Q6–PRISM obstructions, path-valued maps, \(\Lambda^2\) and GEC-6;
- A5 spectrum, De Aquino candidate normalization;
- shear response, Gibbs benchmark and finite OU Green–Kubo;
- the Pellis two-sheet and \(D_6\) correspondence.

Each record has an evidence label and explicit dependencies. The detailed derivations, implementations, tests and audit outputs remain in their focused files; they are not duplicated here.

## Current open problem for independent mathematical review

**Establish, or identify an obstruction to, Mosco convergence of the curved Ellis–RGS P1 energies on the refinement family specified in this repository.**

The problem is deliberately limited to a scalar Riemannian Dirichlet form. It is not a request to validate the entire RGS–PSU–QLG programme.

- Domain: \(M_L=[-L,L]\times S^2\), with \(L=3r_0\), \(r_0=0.8292\).
- Metric: \(g=du^2+(u^2+r_0^2)d\Omega^2\).
- Common Hilbert space: \(H=L^2(M_L,d\mu_g)\), with \(d\mu_g=(u^2+r_0^2)\,du\,d\Omega\).
- Discrete spaces: nested curved tetrahedral P1 meshes generated from the recovered RGS 64B seed.
- Target form: \(q(u)=\int_{M_L}|\nabla_g u|_g^2\,d\mu_g\), with domain \(H^1(M_L)\) and natural Neumann boundary condition.

For the reconstruction \(J_h:V_h\to H\), prove both the weak liminf inequality and existence of a strong recovery sequence. A valid contribution may instead give a counterexample or identify a missing hypothesis.

**Available starting points**
- [Curved refinement G1](docs/RGS_ELLIS_CURVED_REFINEMENT_G1.md)
- [Common-\(L^2\) embedding G4](docs/RGS_ELLIS_CURVED_G4_COMMON_L2.md)
- [Open proof programme G5–G6](docs/RGS_MOSCO_CURVED_REFINEMENT_PROGRAM.md)
- [Flat periodic P1 Mosco theorem](docs/FLAT_TORUS_G6_P1_MOSCO_PROOF.md) — a separate benchmark, not a proof of the curved result.

Run the curved audit:
```bash
PYTHONPATH=. python examples/rgs_ellis_curved_refinement_audit.py
PYTHONPATH=. python -m pytest tests/test_rgs_ellis_refinement.py -v
```

The missing step is uniform all-level control of element maps and inverses, metric/Jacobian bounds, coercivity and compactness, followed by liminf and recovery. Finite-level audits are not substitutes for these estimates. No Lorentzian, spinorial, teleparallel, Yang–Mills or QCD convergence claim is included.

## Principal implementation areas

| Area | Entry point |
|---|---|
| Quaternion, SU(2), SO(3), binary icosahedral group | `src/`, [THEOREM records](https://github.com/fabianleonaressi95/theorem/tree/main/theorems) |
| Q6 graph and square complex | [Q6 formalization](docs/Q6_COMPLEX_FORMALIZATION.md) |
| Recovered RGS seed and \(O_h\) action | [RGS/Q6 symmetry audit](docs/RGS_Q6_OH_EQUIVARIANCE.md) |
| Combinatorial Hodge theory | [Q6 Hodge specification](docs/Q6_COMBINATORIAL_HODGE_LAPLACIANS.md) |
| Curved Ellis geometry and spin connection | [Ellis refinement](docs/RGS_ELLIS_CURVED_REFINEMENT_G1.md), [spin connection](docs/ELLIS_RGS_LEVI_CIVITA_SPIN_CONNECTION.md) |
| Qutrit chain and gauge-covariant Q6 embedding | [Qutrit specification](docs/QUTRIT_Q6_GAUGE_COVARIANT_EMBEDDING.md) |
| Kathara and Krystal Spiral | [Kathara graph audit](docs/KATHARA_15_LINE_GRAPH_FACE_AUDIT.md), [square geometry](docs/KATHARA_KRYSTAL_SPIRAL_SQUARE_GEOMETRY.md) |
| Q6 dyadic refinement and barycentric quadrature | [Q6 refinement](docs/Q6_DYADIC_SUBDIVISION_CHAIN_MAP.md), [Gate E.1](docs/GATE_E1_BARYCENTRIC_COVARIANCE_QUADRATURE.md) |
| Q6–PRISM and GEC-6 | `docs/research/Q6_PRISM_*.md`, `docs/research/GEC6_*.md` |
| Shear response and Green–Kubo limits | [Q6 action-derived shear](docs/Q6_ACTION_DERIVED_SHEAR_RESPONSE.md), [finite OU/GK record](docs/Q6_HYPERCUBE_GREEN_KUBO_SUPERGOLDEN.md) |
| Pellis correspondence | [Dossier index](docs/PELLIS_D_RGS_PSU_PRISM_INTEGRATION.md) |

## Repository map and provenance

- `src/` — mathematical and computational implementations.
- `tests/` — focused unit, structural and regression tests.
- `examples/` — reproducible audit runners.
- `results/` — generated finite audit outputs.
- `docs/` — mathematical specifications, evidence ledgers and open proof programmes.
- [Experiment inventory](experiments/README.md) — experiment navigation and safe variant-management rules.
- [Unified theorem dependency graph](docs/UNIFIED_THEOREM_DEPENDENCY_GRAPH.md) — typed links to the central theorem records.
- [Master catalogue](docs/MASTER_CATALOGUE_LAWS_THEOREMS.md) — programme-wide laws, theorems, models and open hypotheses.
- [Implementation roadmap](docs/IMPLEMENTATION_ROADMAP.md).

For every reported benchmark, identify the exact repository revision, source module, test/audit command, expected invariant, tolerances and limitations. Do not paste raw test logs into this README; retain them as CI artifacts or dedicated audit records.

## Collaboration

Independent mathematical review, alternative implementations and reproducibility audits are welcome. The most useful contributions state exact hypotheses, identify the relevant finite or analytic object, and distinguish a proof from a numerical check.

## License

See the repository license file. If no license file is present, reuse permissions should not be inferred.
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

This repository is the implementation and computational-audit workspace for the quaternionic, Q6 and QLG constructions. The central **THEOREM** repository organizes the corresponding mathematical statements, proof obligations, evidence levels and cross-project results.

- **Central research programme:** [THEOREM — RGS–PSU–PRISM–Q6 Mathematical Research Programme](https://github.com/fabianleonaressi95/theorem)
- **Central README and theorem index:** [theorem/README.md](https://github.com/fabianleonaressi95/theorem/blob/main/README.md)
- **Lorentz-covariant QLG–Dirac chirality result:** [T049D](https://github.com/fabianleonaressi95/theorem/blob/main/theorems/T049D-lorentz-chiral-qlg-embedding.md) (merged into the central repository's `main` branch).

The two repositories have complementary roles: this repository contains implementation details, numerical audits and model-specific documentation; the central repository records the mathematical claims and their epistemic status. A result documented in either repository does not automatically certify every physical interpretation or establish an equivalence between models. Such claims require their own derivations and tests.
