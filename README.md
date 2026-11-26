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

This repository is the implementation and computational-audit workspace for the quaternionic, Q6 and QLG constructions. The central **THEOREM** repository organizes the corresponding mathematical statements, proof obligations, evidence levels and cross-project results.

- **Central research programme:** [THEOREM — RGS–PSU–PRISM–Q6 Mathematical Research Programme](https://github.com/fabianleonaressi95/theorem)
- **Central README and theorem index:** [theorem/README.md](https://github.com/fabianleonaressi95/theorem/blob/main/README.md)
- **Lorentz-covariant QLG–Dirac chirality result:** [T049D](https://github.com/fabianleonaressi95/theorem/blob/main/theorems/T049D-lorentz-chiral-qlg-embedding.md) (merged into the central repository's `main` branch).

The two repositories have complementary roles: this repository contains implementation details, numerical audits and model-specific documentation; the central repository records the mathematical claims and their epistemic status. A result documented in either repository does not automatically certify every physical interpretation or establish an equivalence between models. Such claims require their own derivations and tests.
