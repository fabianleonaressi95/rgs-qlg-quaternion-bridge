# Repository Roadmap — RGS–PSU Gauge–Clifford QLG

This roadmap is an implementation plan for the research proposal in
[docs/RGS_PSU_GAUGE_CLIFFORD_QLG_RESEARCH_PROPOSAL.md](RGS_PSU_GAUGE_CLIFFORD_QLG_RESEARCH_PROPOSAL.md).

## Baseline repository audit

Inspected on the current default branch:

- Existing scalar-first Hamilton quaternion implementation: \`src/quaternions.py\`.
- Active right-handed SO(3) rotations of column vectors: \`src/rotations.py\`.
- Quaternion/SU(2) representation and lift diagnostics: \`src/su2_lifts.py\`.
- SU(2) lattice links, path transport, gauge transformations and Wilson holonomies: \`src/holonomies.py\`.
- Existing tests: \`tests/test_quaternion_algebra.py\`, \`tests/test_su2_so3.py\`, \`tests/test_holonomies.py\`.
- Existing CI workflow compiles \`src\`, \`tests\`, \`examples\` and runs the full pytest suite on pushes and pull requests to \`main\`.
- Dependencies currently declared in \`requirements.txt\`: NumPy and pytest.

Important convention already encoded in \`src/holonomies.py\`: for an oriented link from source g to target h, the stored link transforms as
\`U_gh -> G_h U_gh G_g†\`. Path transport multiplies later links on the left. New code must preserve this convention; do not silently switch to the opposite endpoint convention.

The repository already contains gauge transformation and Wilson-loop tests. The roadmap therefore extends those foundations instead of duplicating or replacing them.

## Implemented in this roadmap commit

### Phase 0 — Baseline and conventions
- [x] Read the existing quaternion, rotation, SU(2)-lift and lattice-holonomy modules.
- [x] Read the current tests and CI workflow.
- [x] Record the existing gauge-transformation and path-order convention.
- [x] Preserve all existing modules and tests.

### Phase 1A — Q6 incidence foundation
- [x] Add a canonical Q6 complex builder with big-endian six-bit indexing.
- [x] Orient edges by increasing integer vertex label.
- [x] Orient each square by its ordered pair of coordinate axes.
- [x] Build integer incidence matrices B1 (64 x 192) and B2 (192 x 240).
- [x] Add tests for cardinalities, degree, shapes and the exact identity B1 B2 = 0.
- [x] Run GitHub Actions: 65 tests passed; lattice holonomy demonstration passed.
- [x] Repair pre-existing 2I audit lookup bug (NumPy rows were used as dictionary keys instead of rounded quaternion tuples).
- [x] Correct the Q6 axis test to reflect that all six single-bit flips are edges.

## Next phases

### Phase 1B — Möbius complex audit
- [ ] Encode the agreed 12-vertex, 18-edge, 6-face subcomplex without conflating it with all Q6 plaquettes.
- [ ] Verify B1 B2 = 0, boundary components, Euler characteristic and Betti numbers.
- [ ] Verify non-orientability by orientation propagation.

### Phase 2 — Gauge-covariance audit
- [ ] Add randomized vertex gauge transformations using the existing SU2Lattice convention.
- [ ] Check open-path covariance and closed-loop conjugation.
- [ ] Check Wilson-trace invariance and inverse-edge consistency.
- [ ] Export residuals, tolerances, random seed and input metadata as JSON.

### Phase 3 — Local frame and feature representations
- [ ] Specify the local frame group separately from the internal SU(2) group.
- [ ] Define scalar, sign-changing and regular two-channel representations for Möbius reflection transitions.
- [ ] Verify transition cocycles and seam-crossing transport.
- [ ] Do not infer torsion or a Pin structure from non-orientability alone.

### Phase 4 — 2I representations and intertwiners
- [ ] Audit the 120 unit-quaternion elements, closure, inverses and multiplication convention.
- [ ] Construct explicit representation matrices.
- [ ] Compute intertwiner spaces and verify their defining equations.
- [ ] Compare dimensions against representation-theoretic predictions.

### Phase 5 — Clifford-valued operators
- [ ] Fix the Clifford relation, basis ordering and metric signature.
- [ ] Implement grade projections and group actions.
- [ ] Build steerable kernel bases as intertwiners.
- [ ] Test equivariance independently of any learning/training code.

### Phase 6 — Discrete differential operators and continuum limit
- [ ] Define discrete inner products and adjoints.
- [ ] Test consistency and stability on controlled smooth fields.
- [ ] State compactness hypotheses and formulate the Mosco convergence problem.
- [ ] Do not claim a continuum theorem from equivariance or spectral plots alone.

## Continuous integration gate

The existing workflow already runs \`python -m compileall -q src tests examples\` and \`python -m pytest tests/ -v\`. Every new module must be covered by this gate. A passing local test suite is not claimed until the workflow result is observed.

## Evidence labels

- **PROVED** — mathematical derivation with assumptions.
- **NUMERICALLY VERIFIED** — reproducible computation with inputs and tolerances.
- **CANDIDATE** — specified but incompletely verified construction.
- **OPEN** — unresolved mathematical or physical question.
