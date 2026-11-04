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

### Phase 1B — Möbius complex audit
- [x] Encode the agreed six face patterns as a 12-vertex, 18-edge, 6-face subcomplex.
- [x] Preserve the cyclic free-bit order 00, 10, 11, 01; lexicographic order would incorrectly introduce crossed diagonals.
- [x] Build integer incidence matrices B1 (12 x 18) and B2 (18 x 6).
- [x] Add an exact chain-complex test B1 B2 = 0.
- [x] Audit 12 boundary edges, 6 interior edges and one boundary component.
- [x] Compute Euler characteristic 0 and Betti numbers (b0,b1,b2) = (1,1,0).
- [x] Detect non-orientability by face-orientation propagation.

### Phase 2 — Gauge-covariance audit
- [x] Generate deterministic random SU(2) link variables on all 18 Möbius edges.
- [x] Generate independent SU(2) gauge transformations at all 12 vertices.
- [x] Check open-path covariance on each of the six face paths.
- [x] Check closed-loop holonomy conjugation on all six square boundaries.
- [x] Check normalized Wilson-trace invariance.
- [x] Check reverse-link consistency before and after gauge transformation.
- [x] Record seed, tolerance, graph counts and maximum residuals in a JSON audit; the module CLI writes to results/json/.

## Next phases

### Phase 3 — Local frame and feature representations
- [x] Derive the six dual-face Z2 transition signs from shared-edge incidence and the fixed cyclic face orientations.
- [x] Keep the local-frame Z2 transition system separate from the independent internal SU(2) lattice gauge field.
- [x] Implement scalar, sign-changing and regular two-channel representations.
- [x] Verify the representation homomorphism laws and local-frame change covariance.
- [x] Verify the dual-face cycle (0,1,2,5,4,3,0) has Z2 holonomy -1, the discrete orientation obstruction.
- [x] Verify regular two-channel transport swaps components across a reversing transition.
- [x] Test invalid signs and feature dimensions; 85 tests passed on the Phase 3 implementation commit.
- [x] Do not infer torsion or a Pin structure from non-orientability alone.


### Phase 4 — 2I representations and intertwiners
- [x] Reuse the exhaustive 120-element 2I audit: unit norm, all 14,400 products closed, inverses present and 60 SO(3) classes.
- [x] Preserve scalar-first Hamilton multiplication and active column-vector rotation conventions.
- [x] Implement the trivial 1D, defining fundamental 2D complex, and vector 3D representations.
- [x] Verify representation homomorphism and unitarity/orthogonality on all 120 elements (and deterministic product samples).
- [x] Verify the central element -1 acts as +1 in trivial/vector representations and -I in the fundamental representation.
- [x] Compute numerical Hom-space bases from the defining equations rho_out(g)K = K rho_in(g).
- [x] Verify selected Schur-lemma dimensions: one-dimensional commutants and zero intertwiners between inequivalent selected irreducibles.
- [x] Keep this explicitly a selected representation set, not the complete irreducible character table of 2I.
- [x] Keep orientation-reversing Möbius transitions outside 2I; no reflection or Pin lift is inferred.


### Phase 5 — Clifford-valued operators
- [x] Fix signature eta=diag(-1,+1,+1,+1) and the anticommutator convention {gamma^mu,gamma^nu}=2 eta^{mu nu} I_4.
- [x] Implement an explicit complex 4x4 Dirac representation and the 16 ordered Clifford blades.
- [x] Verify the 16 anticommutator identities and linear independence of the 16-blade basis.
- [x] Implement grade 0 through grade 4 projections and reconstruction.
- [x] Restrict the 2I defining SU(2) lift to spatial Spin(3) acting on Dirac spinors by diag(U,U).
- [x] Verify the spinor representation law, unitarity, central -1 action, and covariance of the three spatial gamma matrices.
- [x] Compute the 4-dimensional space of 2I-equivariant Dirac-spinor endomorphisms as intertwiners.
- [x] Test Clifford and group equivariance independently of any learning/training code.
- [x] Explicitly limit this phase to spatial rotations: boosts, Möbius reflections, Pin lifts and a full Lorentz-group action remain open.


### Phase 6 — Weiler-compatible geometry and discrete Dirac candidate
- [x] Add the static 3+1 tetrad adapter and pointwise curved gamma fields from the RGS spatial coframe.
- [x] Verify all 16 curved Clifford anticommutators at every vertex, including a non-symmetric coframe test.
- [x] Add compact-support graph kernel-field transform with explicit source-to-target transport.
- [x] Define a polar-frame SU(2) transport candidate, keeping it distinct from the internal gauge field and Möbius Z2 transition system.
- [x] Assemble a volume-weighted Hermitian nearest-neighbour spatial Dirac candidate and test Hermiticity and coordinate-dilation scaling.
- [x] Map the supplied Weiler results (Chs. 3–5, 12–13, 15, 17; Appendices B and J) to implementation obligations and explicit unimplemented hypotheses in docs/WEILER_PHASE6_THEOREM_MAP.md.
- [ ] Identify the subgroup of Q6 automorphisms preserving the reconstructed coordinates and vierbein; then test exact kernel-field invariance/equivariance.
- [ ] Complete the 2I/Clebsch–Gordan kernel basis and invariant bias/pooling spaces.
- [ ] Test smooth-field consistency, spectral stability, fermion doubling and compactness; formulate Mosco convergence.
- [ ] Extend from the static spatial operator to a 4D Lorentzian/time-dependent Dirac operator only after defining lapse, shift and temporal connection.
- [ ] Do not claim a continuum theorem from equivariance or spectral plots alone.

## Continuous integration gate

The existing workflow already runs \`python -m compileall -q src tests examples\` and \`python -m pytest tests/ -v\`. Every new module must be covered by this gate. A passing local test suite is not claimed until the workflow result is observed.

## Evidence labels

- **PROVED** — mathematical derivation with assumptions.
- **NUMERICALLY VERIFIED** — reproducible computation with inputs and tolerances.
- **CANDIDATE** — specified but incompletely verified construction.
- **OPEN** — unresolved mathematical or physical question.


### Phase 9 — Wilson refinement, overlap/Ginsparg–Wilson, and Dirac–Kähler
- [x] Derive the low-momentum Wilson correction and Brillouin-corner masses for the existing centered operator.
- [x] Determine the simultaneous scaling window for a running Wilson parameter: if r(a)=r0(a/L)^alpha, then low-mode correction vanishes and doubler masses diverge iff -1<alpha<1.
- [x] Audit all eight 3D Brillouin corners, including the seven lifted doubler corners and multiplicities by the number of pi-components.
- [x] Add a multi-resolution Fourier-symbol audit at fixed physical momentum (N=16,32,64,128,256).
- [x] Add a finite overlap construction from a gamma5-Hermitian Euclidean Wilson kernel and test the Ginsparg–Wilson identity.
- [x] Add the unweighted combinatorial Dirac–Kähler operator on the Q6 2-skeleton (cochain dimensions 64,192,240), verifying d^2=0, Hermiticity, D_K^2=Hodge Laplacian and degree-parity anticommutation.
- [x] Document the distinctions between the spinor Wilson operator, overlap construction, and cochain Dirac–Kähler operator.
- [x] Build the overlap operator from the archived RGS-Q6 256x256 Wilson/Dirac pair; verify the kernel gap and Ginsparg–Wilson/gamma5-Hermiticity residuals, and record finite graph-distance decay. This is not a uniform locality theorem.
- [ ] Construct a vierbein-weighted DEC Hodge star and metric Dirac–Kähler operator.
- [ ] Specify a refining family of RGS complexes, common Hilbert-space embeddings and quadratic forms; prove Mosco liminf and recovery-sequence conditions.
- [ ] Prove spectral convergence and a continuum index/anomaly statement under explicit hypotheses.
- [ ] No Q6 physical/doubler labels are inferred from a finite spectrum without a momentum or continuum identification.

- [x] Audit the full dispersion of the Hermitian scalar sum D0+W; identify its additional one-dimensional zero at pa=2 arctan(2/r), and distinguish it from the Euclidean Wilson operator.
- [x] Add `assemble_euclidean_wilson_dirac` implementing D_E=iD0+W and verify gamma5-Hermiticity. The overlap construction uses this Euclidean convention.

- [x] Prove a Fourier-embedded Mosco-convergence prototype for the Euclidean Wilson quadratic form on a flat periodic torus, distinguishing Mosco liminf/recovery from the additional UV compactness needed for spectral convergence.
- [x] Add numerical sampling of the fixed-mode quadratic form and the doubler mass along the periodic-torus refinement sequence.
- [ ] Transfer the prototype proof to a nonuniform RGS refinement family; the fixed Q6 graph does not provide such a family.


### Phase 10 — Physics architecture compatibility and local spin connection
- [x] Add a central sector/evidence registry and explicit coupling gates.
- [x] Cross-check the Clifford and EC gamma representations, metric signature and gamma5 convention.
- [x] Re-audit the canonical Q6 incidence dimensions and exact chain identity.
- [x] Derive a local torsion-free Levi-Civita spin connection from a coframe and its first derivatives.
- [x] Check Lorentz antisymmetry, tetrad postulate and vanishing torsion on flat and spherical-coordinate Minkowski benchmarks.
- [x] Convert the spin connection to the Dirac representation and construct coordinate and orthonormal covariant derivatives.
- [x] Feed the orthonormal derivative into the action-derived EC–Dirac kinetic Euler operator in a bridge test.
- [ ] Derive consistent coframe first derivatives on the curved RGS refinement family; the local routine does not infer derivatives from the 64-site point cloud.
- [ ] Construct the teleparallel inertial spin connection separately and verify the TEGR boundary identity.
- [ ] Add internal U(1), SU(2) and SU(3) connections as separate representation-valued terms.
- [ ] Derive the full common action, tetrad/connection equations, Noether identities and total stress conservation.
- [ ] Establish curved discrete-to-continuum convergence and quantum anomaly/doubling gates.

The compatibility registry and local spin-connection scope are documented in
`docs/PHYSICS_ARCHITECTURE_COMPATIBILITY.md` and
`docs/LEVI_CIVITA_SPIN_CONNECTION_CORE.md`. A passing finite benchmark does
not close the RGS refinement or teleparallel gates.
