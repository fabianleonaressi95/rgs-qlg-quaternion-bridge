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
