# Master Research Roadmap — RGS–PSU–PRISM–QLG

**Status date:** 2026-11-06  
**Role:** programme-level decision and navigation document for this repository.  
**Scope:** research prioritisation, dependency order, evidence gates and repository hygiene.  
**Not a theorem:** this roadmap organises work; it does not upgrade any result's scientific status.

## 1. Executive decision

The programme has many implemented finite constructions, audits and exploratory branches. The immediate bottleneck is not a lack of models or numerical experiments: it is the absence of a single, enforced dependency chain from geometry to analysis to continuum limits.

**The next scientific objective is the scalar curved Ellis–RGS Mosco problem.** It is the narrowest well-defined continuum question already backed by a concrete refinement family, a common-L2 reconstruction and curved P1 mass/stiffness forms. Spinorial, gauge, teleparallel, QCD, PRISM and quantum-light extensions must not be used as substitutes for completing this scalar proof gate.

For the next work cycle, do not start a new physical model or a new parallel discretisation. Close one proof obligation at a time, record the result, then advance.

## 2. Programme architecture and repository roles

| Layer | Role | Canonical location |
|---|---|---|
| Theorem statements, assumptions, proof records and programme-wide evidence register | Mathematical source of truth | [THEOREM repository](https://github.com/fabianleonaressi95/theorem) |
| Reusable implementations, finite geometric constructions, numerical audits and regression tests | Computational source of truth | This repository: `src/`, `tests/`, `examples/`, `scripts/` |
| Scientific derivations, conventions, proof programmes and architecture boundaries | Documentation source of truth | `docs/` |
| Exploratory variants and historical computations | Research archive, not production dependencies | `experiments/`, `archive/` |
| Inputs and generated evidence | Immutable inputs plus provenance-tracked outputs | `data/`, `reports/`, `results/` |

A result should have one canonical implementation and one canonical theorem/proof record. Other files should link to it rather than restating it with subtly different conventions.

## 3. Priority roadmap

### P0 — Stabilise and de-duplicate the research surface

**Goal:** make it possible to tell what is active, canonical, validated, duplicated or merely exploratory.

Tasks:
1. Maintain this document as the programme-level index; keep `docs/IMPLEMENTATION_ROADMAP.md` for implementation history and detailed phase checklists.
2. Reconcile the repository map, experiment catalogue, dependency audit and open pull requests. Record the current branch, commit, CI result and evidence status for every active deliverable.
3. Preserve the 311-file experiment inventory. Do not bulk-delete, rename or merge variants based on filenames or apparent similarity.
4. Triage the 15 reported Python parse blockers individually. First decide whether each file is executable code, a notebook/export, a mathematical note or imported material; repair only after recovering its intended meaning.
5. Keep generated reports out of source modules and require provenance (commit, command, inputs, environment, units, tolerances).

**Exit gate:** every active result has a canonical path, owner/reviewer, evidence label, reproducible command and next action; no unresolved duplicate is silently promoted or deleted.

### P1 — Finish the curved Ellis mesh geometry gate (G1–G5 prerequisites)

**Goal:** establish refinement-uniform geometry, not just good finite-level plots.

Canonical references:
- `docs/RGS_ELLIS_CURVED_REFINEMENT_G1.md`
- `docs/ELLIS_RGS_MESH_SHAPE_REGULARITY.md`
- `examples/run_ellis_mesh_shape_regularity_audit.py`
- `src/rgs_ellis_refinement.py`

Tasks, in order:
1. Write the exact recursive definition of the angular midpoint-normalisation map and the fixed tetrahedral staircase split.
2. Prove a level-independent lower bound for spherical face quality and edge comparability under the actual midpoint-normalised subdivision.
3. Prove radial/angular scale comparability for the chosen refinement rule.
4. Transfer those bounds through the Ellis embedding using explicit singular-value bounds for its differential.
5. Derive a uniform element Jacobian/volume comparison, in particular constants (0<c_J\le J_K/(6V_K^{\rm ch})\le C_J<\infty) independent of level.
6. Use the geometric bounds to prove uniform shape regularity and curved-versus-chordal P1 mass equivalence.

**Exit gate:** a written proof with constants and hypotheses that match the generator exactly. Finite-level diagnostics remain supporting evidence, not the proof.

### P2 — Complete scalar common-(L^2) stability (G4–G5)

**Goal:** obtain the functional-analytic estimates needed for Mosco convergence.

Canonical references:
- `docs/RGS_ELLIS_CURVED_G4_COMMON_L2.md`
- `docs/RGS_MOSCO_CURVED_REFINEMENT_PROGRAM.md`
- `examples/run_ellis_lumped_mass_equivalence_audit.py`
- `examples/run_ellis_neumann_poincare_audit.py`

Tasks:
1. Freeze the common Hilbert space (H=L^2(M_L,d\mu_g)), reconstruction (J_h), discrete domain and exact curved P1 forms. No silent switching between consistent and lumped masses.
2. Prove uniform equivalence of the curved mass, reference P1 mass and any auxiliary lumped norm actually used.
3. Prove a level-independent Neumann Poincare inequality on the mean-zero subspace.
4. Prove discrete Rellich compactness for bounded reconstructed (L^2) norm and energy.
5. Verify the boundary convention and show that the discrete forms converge to the intended Neumann form, with no hidden multiplicative constant.

**Exit gate:** uniform analytic estimates, not just a finite list of positive eigenvalues.

### P3 — Prove the scalar curved Mosco theorem (G6)

**Target:** for the explicitly defined extended forms on the common (H), prove both:

- **Liminf:** weak (L^2) convergence implies the lower energy bound.
- **Recovery:** for every (u\in H^1(M_L)), construct (u_h) with strong (L^2) convergence and the upper energy bound.

Use smooth functions first, then density and diagonalisation. State the exact hypotheses inherited from P1 and P2. If a required hypothesis fails, produce a counterexample or revise the mesh/form definition transparently rather than weakening the theorem invisibly.

**Exit gate:** a complete proof with no reliance on spectral-fit extrapolation. Only after this gate may the scalar strong-resolvent and spectral consequences be claimed.

### P4 — Extend to differential forms and spinors (separate theorems)

**Goal:** extend the scalar result only through separately specified complexes and bundles.

Order:
1. Curvilinear cubical/polytopal de Rham complex: define pullbacks, cell maps, metric-dependent inner products and commuting projections.
2. Prove discrete differential identities and uniform stability for the actual non-affine cell maps.
3. Establish the relevant Hodge/Mosco result degree by degree.
4. Define the spin structure and the Levi-Civita spin connection on the same geometric family.
5. Prove convergence of discrete parallel transport and spinorial recovery sequences.
6. Treat Wilson doubler control, overlap locality and Dirac strong-resolvent convergence as distinct gates.

Periodic Fourier Wilson estimates are a reference model only; their Brillouin-corner scaling window must not be transferred to curved graphs without a separate proof.

### P5 — Physical extensions (research tracks, not current blockers)

These tracks remain active as model-development or source-review work, but they do not outrank P1–P3:

- **RGS–PSU / 188D vierbein / Einstein–Cartan / teleparallel:** keep Lorentzian and Riemannian problems separate; derive torsion and field equations from an explicit action before physical interpretation.
- **Q6 / SU(2) / SU(3) / QCD / Yang–Mills:** distinguish finite lattice identities and gauge-covariance tests from dynamical continuum limits.
- **PRISM / Krist-All / QVSM:** preserve the exact finite combinatorics and independent evidence ledger; do not infer universal scaling from finite ratios.
- **2I / (S^3/2I) / Rumin / spectral geometry:** separate standard quotient facts, finite representation computations and uncomputed analytic spectra.
- **PSU–Anti-Twin / QAT / quantum light / emergent spacetime:** treat as conjectural model architecture until equations, observables, limits and empirical discriminants are specified.
- **External source audits (including non-standard ether/magnetic-rhythm claims):** reconstruct definitions and equations critically; do not import claims as established physics.

**Entry rule:** a physical extension enters implementation priority only with a precise mathematical question, a defined object/action, a reproducible test and a stated falsification or proof criterion.

## 4. Active work queue (ordered)

1. **PR #193:** review the finite-level Ellis shape-diagnostic changes and their CI; do not describe them as a uniform shape-regularity proof.
2. **Geometry lemma:** turn the midpoint-normalised spherical subdivision and radial/angular scaling into explicit uniform inequalities.
3. **Jacobian comparison:** prove the uniform curved/chordal element measure bounds.
4. **Mass equivalence:** close the common-(L^2) norm comparison using the geometric lemma.
5. **Uniform Poincare and Rellich:** replace finite-level evidence with level-independent estimates.
6. **Mosco liminf and recovery:** write and independently audit the full proof.
7. **Only then:** begin the curvilinear de Rham and spinorial convergence extensions.

No new exploratory model is to be added to this queue unless it blocks one of these steps or has a separate externally required deadline.

## 5. Evidence and completion policy

Every task and claim must carry exactly one primary status:

- **PROVED:** complete mathematical proof; hypotheses and conventions stated.
- **NUMERICALLY VERIFIED:** reproducible finite computation with commit, command, inputs, residuals and tolerances.
- **CANDIDATE:** construction exists but at least one required proof or validation gate is open.
- **OPEN:** precise unresolved question.
- **UNVERIFIED / IMPORTED:** provenance, executable meaning or independent validation is incomplete.

A passing CI run verifies the configured software checks, not every scientific statement. A plot, regression test, finite spectrum or small residual does not prove a refinement-uniform estimate.

A task is complete only when its source, tests, documentation, evidence status and theorem-level record agree.

## 6. Pull-request and change-control rules

1. One PR should address one coherent gate and list the exact files, mathematical claim, tests and remaining limitations.
2. Work on a named branch; never write directly to `main`.
3. Before editing, inspect the current branch head and file blob SHA. Avoid overwriting intervening changes.
4. Do not merge, close, retarget or rewrite another PR without explicit authorisation.
5. Do not create a second roadmap or duplicate implementation where a canonical file already exists.
6. After each PR, update this active queue and the relevant theorem record; do not update unrelated physical claims.

## 7. Canonical navigation

- Programme-level priorities: this document.
- Detailed implementation history: `docs/IMPLEMENTATION_ROADMAP.md`.
- Repository structure and safe maintenance: `docs/REPOSITORY_MAP.md`.
- Experiment inventory: `experiments/README.md` and `experiments/EXPERIMENTS_CATALOG.json`.
- Dependency findings and parse blockers: `docs/EXPERIMENT_DEPENDENCY_GRAPH_AUDIT.md`.
- Curved scalar Mosco proof programme: `docs/RGS_MOSCO_CURVED_REFINEMENT_PROGRAM.md`.
- Curved common-(L^2) construction: `docs/RGS_ELLIS_CURVED_G4_COMMON_L2.md`.
- Evidence catalogue: `docs/MASTER_CATALOGUE_LAWS_THEOREMS.md`.
- Programme-wide theorem and proof records: [fabianleonaressi95/theorem](https://github.com/fabianleonaressi95/theorem).

---

**Current decision:** concentrate the next technical cycle on P1–P3. The scalar curved Ellis–RGS Mosco theorem is **OPEN**. All spinorial and physical continuum extensions remain separate, downstream research questions.
