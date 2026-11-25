# RGS–QLG Quaternion Bridge

**Implementation and reproducibility repository** for the quaternionic, Q6, RGS, lattice-gauge and finite-operator components of the RGS–PSU–QLG programme.

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

Run: \`PYTHONPATH=. python -m pytest tests/test_pellis_d6_operator_gate.py tests/test_fibonacci_d6_channel_intertwiner.py -q\`.


Pellis D6 commutant refinement: src/pellis_d6_commutant.py and tests/test_pellis_d6_commutant.py implement the general singlet-plus-doublet radial transfer allowed by D6 symmetry. The equal three-mode lift is a special case, not a consequence of D6 alone.

Run: PYTHONPATH=. python -m pytest tests/test_pellis_d6_commutant.py tests/test_pellis_d6_operator_gate.py tests/test_fibonacci_d6_channel_intertwiner.py -q.
