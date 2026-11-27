# RGS–QLG Quaternion Bridge

**Implementation and reproducibility repository** for the quaternionic, Q6, RGS, lattice-gauge and finite-operator components of the RGS–PSU–QLG programme.

This repository develops finite mathematical constructions and computational audits within the RGS–PSU–QLG programme. It covers quaternionic groups, Q6 cubical geometry, RGS-induced metrics, finite-element operators, lattice gauge covariance, spectral representation theory and selected continuum-limit questions.

> **Evidence policy:** PROVED, FINITE-EXACT, NUMERICAL, CONDITIONAL, CANDIDATE and OPEN are distinct statuses. A finite test is not a continuum theorem or a physical derivation.

## Canonical suite formalization

The 51-step sequence formerly expanded in this README is now maintained as a dedicated mathematical register:

- **[Integrated mathematical suite — 51 formal records](docs/INTEGRATED_MATHEMATICAL_SUITE_FORMALIZATION_2026-11-25.md)** — definitions, precise claims, evidence status, dependencies and proof obligations for every step from quaternionic foundations through Q6, RGS/FEM, qutrit networks, Kathara, Q6–PRISM, A5/De Aquino, shear/Green–Kubo and Pellis.
- [51-step evidence and implementation matrix](docs/51_STEP_EVIDENCE_MATRIX_2026-11-26.md) — source/test/result traceability and outstanding proof obligations. Post-snapshot supplements for merged PRs #233–#234 are recorded at the end of that matrix.
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
```bash
PYTHONPATH=. python examples/rgs_ellis_curved_refinement_audit.py
PYTHONPATH=. python -m pytest tests/test_rgs_ellis_refinement.py -v
```

This is a scalar Riemannian convergence question, not a claim about Lorentzian dynamics, spinorial convergence, teleparallel field equations, Yang–Mills or QCD.

## GEC-6 trigram block decomposition

The exterior-square pair module splits relative to the fixed I-Ching bit blocks (lower 0,1,2; upper 3,4,5) as 3 lower-internal + 9 mixed + 3 upper-internal directions. Its block-preserving wreath symmetry is (S3 x S3) semidirect C2. The mixed sector alone surjects under contraction onto the rank-five A5 quotient; the antipodal complement acts trivially on that quotient. The traditional I-Ching interpretation of the mixed sector and any physical coupling law remain open.

- [GEC-6 trigram block decomposition](docs/research/GEC6_TRIGRAM_BLOCK_DECOMPOSITION.md)
- [Exact implementation](src/gec6_trigram_decomposition.py)
- [Focused tests](tests/test_gec6_trigram_decomposition.py)
- [Audit example](examples/gec6_trigram_decomposition_audit.py)

Reproduce:
~~~bash
PYTHONPATH=. python -m pytest tests/test_gec6_trigram_decomposition.py -q
PYTHONPATH=. python examples/gec6_trigram_decomposition_audit.py
~~~

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


### RGS–PSU carrier connection comparison (27 November 2026)

- [Archived PSU connection versus carrier-derived coframe](docs/RGS_PSU_CARRIER_CONNECTION_COMPARISON_2026-11-27.md) — the A13 shortest-radial-rotation link and the PR #476 carrier link are the same \(SO(3)\) connection in ambient versus local radial frames. The local gauge fixes the radial axis, reducing links to the \(SO(2)\) stabilizer and spin lifts to diagonal \(U(1)\).
- The archived normalized four-component \((\Delta X,\Delta r)\) feature is **not** the carrier's symmetrized three-component spatial coframe. The audit compares both and separately checks the raw spatial Q6 torsion against the previous A13 report.
- Reproduce with: PYTHONPATH=. python scripts/audit_rgs_psu_carrier_connection.py; python -m pytest tests/test_rgs_psu_carrier_connection.py -q.
- Provenance boundary: this is generator-level reconstruction; the original A13 NPZ payload is not bytewise reloaded or compared.


### Icosahedral matching orbits in GEC-6

The six-axis icosahedral action splits the 15 perfect matchings into a distinguished orbit of five and a generic orbit of ten. The diagonal trigram matching is special, but each lower/upper 3+3 split has two special cross-matchings, not one. The full signed A5 action does not descend to the GEC-6 contraction quotient; only a six-element subgroup does. [Dossier and reproduction](docs/research/GEC6_ICOSAHEDRAL_MATCHING_ORBITS.md).
