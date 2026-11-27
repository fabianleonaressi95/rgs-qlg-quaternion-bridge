# 51-step suite — evidence and implementation matrix

**Snapshot:** `main` at the commit used to create this matrix (2026-11-26).
**Purpose:** trace each formal record to implementation, test, retained result and remaining proof obligation.

> **Interpretation rule.** A path association means that a file exists in the repository tree; it does not certify that the implementation proves the mathematical claim. Tests establish only their asserted software properties. The status column preserves the scope and caveats of the formalization register.

## Coverage summary

- Records: **51**.
- Every record has an explicit status and a recorded next obligation.
- A test path is listed only when that exact file exists on `main`.
- A result/certificate is listed only when the exact artifact exists on `main`.
- `—` means no dedicated artifact was located in this repository snapshot; it is a gap to close, not evidence that the claim is false.

## Matrix

| ID | Mathematical item | Status | Implementation | Test | Result / certificate | Remaining obligation |
|---|---|---|---|---|---|---|
| S01 | Hamilton quaternions | PROVED | `src/quaternions.py` | `tests/test_quaternion_algebra.py` | — | Add symbolic/algebraic certificate only if formalized beyond unit tests. |
| S02 | Double cover SU(2) → SO(3) | PROVED | `src/su2_lifts.py`<br>`src/rotations.py` | `tests/test_su2_so3.py` | — | Proof is standard; test verifies implementation conventions. |
| S03 | Binary icosahedral group 2I | PROVED (group); FINITE-EXACT/NUMERICAL (explicit list audit) | `src/binary_icosahedral.py`<br>`src/representations_2i.py` | `tests/test_binary_icosahedral.py`<br>`tests/test_representations_2i.py` | `results/json/quaternion_bridge_2i_audit.json` | Distinguish exact arithmetic closure from floating-point audit. |
| S04 | Intertwiners and Schur | PROVED under irreducibility | `src/q6_2i_intertwiners.py`<br>`src/rgs_icosahedral_q6_intertwiner.py` | `tests/test_q6_2i_intertwiners.py`<br>`tests/test_rgs_icosahedral_q6_intertwiner.py` | — | Record irreducibility evidence for each applied representation. |
| S05 | Six-cube graph | PROVED | `src/q6_complex.py` | `tests/test_q6_complex.py` | — | Graph counts follow directly from definition. |
| S06 | Graph Laplacian spectrum | PROVED | `src/q6_spectral_audit.py` | `tests/test_q6_spectral_audit.py` | — | Retain Walsh-character derivation; numerical diagonalization is cross-check. |
| S07 | Oriented square complex | PROVED; implementation FINITE-EXACT | `src/q6_complex.py` | `tests/test_q6_complex.py` | — | Retain integer incidence matrices and B1 B2 = 0 certificate. |
| S08 | Graph cycle rank vs filled-complex homology | PROVED conditional on integer rank certificate | `src/q6_complex.py` | `tests/test_q6_complex.py` | — | Attach exact integer ranks of B1 and B2; graph rank 129 is not filled-complex beta1. |
| S09 | Weighted combinatorial Hodge Laplacians | PROVED for positive weights | `src/q6_hodge_laplacian.py`<br>`src/q6_geometric_hodge.py` | `tests/test_q6_hodge_laplacian.py`<br>`tests/test_q6_geometric_hodge.py` | — | Keep arbitrary positive cochain weights distinct from geometric DEC stars. |
| S10 | Antipodal involution | FINITE-EXACT audit | `src/q6_complex.py` | `tests/test_q6_complex.py` | — | Explicitly retain signed cell permutations and chain-map identities. |
| S11 | Recovered RGS seed | REPRODUCIBLE NUMERICAL | `src/rgs_q6_induced_metric.py` | `tests/test_rgs_q6_induced_metric.py` | `docs/recovered-rgs-psu-artifacts.md` | Pin source data hash, state ordering, units and provenance manifest. |
| S12 | Octahedral symmetry Oh | FINITE-EXACT (group action); NUMERICAL (coordinates) | `src/rgs_q6_oh_action.py`<br>`src/rgs_icosahedral_q6_intertwiner.py` | `tests/test_rgs_q6_oh_action.py`<br>`tests/test_rgs_icosahedral_q6_intertwiner.py` | — | Separate exact signed-permutation action from floating-coordinate residuals. |
| S13 | Edge and face metrics | NUMERICAL | `src/rgs_q6_induced_metric.py` | `tests/test_rgs_q6_induced_metric.py` | — | Report degeneracies, units, precision and tolerance. |
| S14 | Rank of induced pullback metric | PROVED | `src/rgs_q6_induced_metric.py` | `tests/test_rgs_q6_induced_metric.py` | — | Linear-algebra proof rank((DX)^T DX) ≤ 3; do not call it a 6D nondegenerate metric. |
| S15 | Graph P1 FEM | PROVED (1D element); assembled data FINITE-EXACT/NUMERICAL | `src/q6_barycentric_fem.py` | `tests/test_q6_barycentric_fem.py` | — | Mass matrix requires declared measure/quadrature; not 3D FEM by implication. |
| S16 | Ellis–RGS curved mesh family | FINITE CONSTRUCTION | `src/rgs_ellis_refinement.py` | `tests/test_rgs_ellis_refinement.py` | — | Finite construction does not prove uniform shape regularity. |
| S17 | Flat periodic P1 Mosco theorem | PROVED under stated hypotheses | `src/flat_torus_refinement.py`<br>`src/mosco_torus.py` | `tests/test_flat_torus_mosco_p1.py`<br>`tests/test_flat_torus_refinement.py`<br>`tests/test_mosco_torus.py` | `docs/FLAT_TORUS_G6_P1_MOSCO_PROOF.md` | Keep hypotheses and common-Hilbert-space reconstruction attached to theorem. |
| S18 | Curved Ellis–RGS Mosco problem | OPEN | `src/rgs_ellis_refinement.py` | `tests/test_rgs_ellis_refinement.py` | `docs/RGS_MOSCO_CURVED_REFINEMENT_PROGRAM.md`<br>`docs/RGS_ELLIS_CURVED_G4_COMMON_L2.md` | Missing proof: uniform element-map/Jacobian/metric bounds, coercivity, compactness, liminf and recovery. |
| S19 | Aperiodic qutrit chain | FINITE MODEL; NUMERICAL diagonalization | `src/qutrit_aperiodic.py` | `tests/test_qutrit_aperiodic.py` | `docs/QUTRIT_APERIODIC_ENTANGLEMENT.md` | No QCD derivation; record boundary conditions and symmetry sector. |
| S20 | Gauge-covariant qutrit–Q6 embedding | FINITE-EXACT/NUMERICAL | `src/qutrit_q6_coupling.py` | `tests/test_qutrit_q6_coupling.py` | `docs/QUTRIT_Q6_GAUGE_COVARIANT_EMBEDDING.md` | Gauge covariance of finite model is not dynamical Yang–Mills. |
| S21 | Q3, Omega and Beneš networks | FINITE-EXACT if exhaustive enumeration archived | `src/interconnection_networks_8x8.py`<br>`src/q3_q6_network_bridge.py` | `tests/test_interconnection_networks_8x8.py`<br>`tests/test_q3_q6_network_bridge.py` | `results/interconnection_networks_8x8_audit.json` | Retain wiring, stage order and exhaustive reachable-permutation list. |
| S22 | Trigram–hexagram encoding | PROVED (encoding); instance FINITE-EXACT | `src/trigram_q3_theorems.py`<br>`src/q3_q6_network_bridge.py` | `tests/test_trigram_q3_theorems.py`<br>`tests/test_q3_q6_network_bridge.py` | `results/q3_q6_network_prism_audit.json` | A bijective encoding does not imply identification with PRISM. |
| S23 | Kathara graph transcription | FINITE MODEL / CONDITIONAL | `src/kathara_15_line_graph.py` | `tests/test_kathara_15_line_graph.py` | `docs/KATHARA_15_LINE_GRAPH_FACE_AUDIT.md` | Counts depend on declared source diagram and transcription. |
| S24 | Kathara cellular chain complex | FINITE-EXACT | `src/kathara_cellular_complex.py` | `tests/test_kathara_cellular_complex.py` | `docs/KATHARA_CELLULAR_COMPLEX_AUDIT.md` | Retain integer boundary matrices, ranks and B1 B2 = 0. |
| S25 | Krystal Spiral square similarity | PROVED CONDITIONAL | `src/kathara_square_similarity.py` | `tests/test_kathara_square_similarity.py` | `docs/KATHARA_KRYSTAL_SPIRAL_SQUARE_GEOMETRY.md` | Similarity follows from stated transform; no tribonacci/physical derivation. |
| S26 | Kathara–Q6 chain-map gates | CANDIDATE / FINITE-EXACT; metric canonicity OPEN | `src/kathara_q6_chain_map.py`<br>`src/kathara_q6_nondegenerate_chain_map.py` | `tests/test_kathara_q6_chain_map.py`<br>`tests/test_kathara_q6_nondegenerate_chain_map.py` | `docs/KATHARA_Q6_CHAIN_MAP_GATE_C.md`<br>`docs/KATHARA_Q6_NONDEGENERATE_GATE_D.md` | Separate chain-map identities from injectivity, cellular embedding and metric canonicity. |
| S27 | Barycentric triangle quadrature | PROVED moment identities; implementation audit FINITE-EXACT/NUMERICAL | `src/barycentric_triangle_quadrature.py`<br>`src/barycentric_geometry.py` | `tests/test_barycentric_triangle_quadrature.py` | `docs/GATE_E1_BARYCENTRIC_COVARIANCE_QUADRATURE.md` | Not the 6D covariance of Q6. |
| S28 | Dyadic Q6 2-skeleton subdivision | FINITE-EXACT | `src/q6_dyadic_refinement.py`<br>`src/rgs_q6_cubical_refinement.py` | `tests/test_q6_dyadic_refinement.py`<br>`tests/test_rgs_q6_cubical_refinement.py` | `docs/Q6_DYADIC_SUBDIVISION_CHAIN_MAP.md` | Retain chain-map matrices and integer counts at each level. |
| S29 | Five-factor directional embedding obstruction | PROVED under restricted definition | `src/q6_prism_embedding_obstruction.py` | `tests/test_q6_prism_embedding_obstruction.py` | `docs/research/Q6_PRISM_EMBEDDING_OBSTRUCTION.md`<br>`results/q6_prism_embedding_audit.json` | Not a no-go for arbitrary graph maps or path-valued maps. |
| S30 | Cubical folding Q6 → Q5 | FINITE-EXACT | `src/q6_prism_chain_map.py` | `tests/test_q6_prism_chain_map.py` | `docs/research/Q6_PRISM_CUBICAL_CHAIN_MAP.md`<br>`results/q6_prism_chain_map_audit.json` | Check ∂F=F∂ on complete oriented source-cell list. |
| S31 | Periodic winding obstruction | PROVED under edge-to-edge and face-compatibility hypotheses | `src/q6_prism_periodic_obstruction.py` | `tests/test_q6_prism_periodic_obstruction.py` | `docs/research/Q6_PRISM_PERIODIC_WINDING_OBSTRUCTION.md`<br>`results/q6_prism_periodic_obstruction_audit.json` | Does not apply to path-valued maps. |
| S32 | Path-valued homology map | PROVED as existence construction under path freedom | `src/q6_prism_path_valued_homology.py` | `tests/test_q6_prism_path_valued_homology.py` | `docs/research/Q6_PRISM_PATH_VALUED_HOMOLOGY.md`<br>`results/q6_prism_path_valued_homology_audit.json` | Does not automatically extend over the 240 filled faces. |
| S33 | A6 equivariance obstruction | PROVED CONDITIONAL | `src/q6_prism_symmetry_obstruction.py` | `tests/test_q6_prism_symmetry_obstruction.py` | `docs/research/Q6_PRISM_SYMMETRY_OBSTRUCTION.md`<br>`results/q6_prism_symmetry_obstruction_audit.json` | Keep exact source/target module actions; not a no-go for another target. |
| S34 | Exterior-square target | FINITE-EXACT / PROVED when formula and equivariance retained | `src/q6_prism_exterior_square_map.py` | `tests/test_q6_prism_exterior_square_map.py` | `docs/research/Q6_PRISM_EXTERIOR_SQUARE_EQUIVARIANT_MAP.md`<br>`results/q6_prism_exterior_square_audit.json` | Retain integral formula, S6-equivariance and rank-114 kernel certificate. |
| S35 | GEC-6 root-Cayley model on A5 lattice | PROVED algebraically; quotient FINITE-EXACT | `src/gec6_a5_root_cayley_torus.py` | `tests/test_gec6_a5_root_cayley_torus.py` | `docs/research/GEC6_A5_ROOT_CAYLEY_PERIODIC_MODEL.md`<br>`results/gec6_a5_root_cayley_audit.json` | A5 denotes root lattice {x∈Z⁶:Σxᵢ=0}, not alternating group. |
| S36 | GEC-6 spectrum for m=3 | PROVED by Fourier classification if complete; cross-check NUMERICAL | `src/gec6_a5_spectral.py` | `tests/test_gec6_a5_spectral.py` | `docs/research/GEC6_A5_LAPLACIAN_SPECTRUM_M3.md`<br>`results/gec6_a5_laplacian_spectrum_m3.json` | Retain full character-classification argument; diagonalization is corroboration. |
| S37 | General-modulus spectral gap | OPEN unless global minimizer and uniform remainder certified | `src/gec6_a5_spectral.py` | `tests/test_gec6_a5_spectral.py` | `results/gec6_a5_gap_scaling_m3_m10.json` | Finite m data cannot prove asymptotic limit; needs global minimum certificate and uniform remainder. |
| S38 | Flat A5 torus spectrum | CONDITIONAL on metric and scale normalization | `src/gec6_a5_continuum_torus.py` | `tests/test_gec6_a5_continuum_torus.py` | `docs/research/GEC6_A5_CONTINUUM_TORUS_SPECTRUM_AND_DEAQUINO_BRIDGE.md`<br>`results/gec6_a5_continuum_spectrum.json` | State lattice metric and dual-lattice normalization explicitly. |
| S39 | Algebraic De Aquino threshold | Threshold PROVED; A5 bridge CANDIDATE/OPEN | `src/rgs_d5_scale_identifiability.py` | `tests/test_rgs_d5_scale_identifiability.py` | `docs/research/DEAQUINO_PENTAGON_A5_MOMENTUM_NORMALIZATION.md` | Derive factor 1/4, physical momentum map and dynamical intertwiner; algebraic zero alone does not. |
| S40 | Pentagon frame normalization | PROVED CONDITIONAL | `src/rgs_d5_scale_identifiability.py` | `tests/test_rgs_d5_scale_identifiability.py` | `docs/research/DEAQUINO_PENTAGON_A5_MOMENTUM_NORMALIZATION.md` | Balanced five-amplitude ansatz is not shown to be dynamically selected. |
| S41 | Shear derivative of implemented Q6 action | PROVED for declared action derivative; finite difference NUMERICAL | `src/q6_action_stress.py`<br>`src/unified_discrete_action.py` | `tests/test_q6_action_stress.py`<br>`tests/test_unified_discrete_action.py` | `docs/Q6_ACTION_DERIVED_SHEAR_RESPONSE.md` | Do not identify automatically with full physical stress tensor. |
| S42 | Inverse-square link-weight shear | PROVED CONDITIONAL | `src/q6_metric_wilson_response.py` | `tests/test_q6_metric_wilson_response.py` | `docs/Q6_METRIC_WILSON_RESPONSE.md` | Separate link-weight response from Wilson plaquette derivative. |
| S43 | Metric-weighted plaquette response | CONDITIONAL + NUMERICAL | `src/q6_metric_wilson_response.py` | `tests/test_q6_metric_wilson_response.py` | — | Need compatible dual measures for a full DEC Hodge star. |
| S44 | Two-sheet covariant matter energy | PROVED discrete form; model choice CONDITIONAL | `src/psu_two_sheet_dirac.py` | `tests/test_psu_two_sheet_dirac.py` | `docs/PSU_TWO_SHEET_DIRAC_SPECTRAL_BRIDGE.md` | Not a Dirac action absent spinor, Clifford, measure and limit structures. |
| S45 | Gibbs stress benchmark | Detailed balance PROVED; correlations NUMERICAL | `src/q6_gibbs_dynamics.py` | `tests/test_q6_gibbs_dynamics.py` | `docs/Q6_GIBBS_STRESS_CORRELATIONS.md` | Monte Carlo step is algorithmic time, not physical time. |
| S46 | Finite OU Green–Kubo integral | PROVED under stability; implementation FINITE-EXACT/NUMERICAL | `src/q6_hypercube_gk_supergolden.py`<br>`src/thermodynamics.py` | `tests/test_q6_hypercube_gk_supergolden.py`<br>`tests/test_thermodynamics.py` | `docs/Q6_HYPERCUBE_GREEN_KUBO_SUPERGOLDEN.md` | Physical viscosity needs calibrated time, stress, FDT, equilibrium and limits. |
| S47 | Scalar spectrum on S³/2I | PROVED with action conventions; multiplicities FINITE-EXACT | `src/binary_icosahedral.py`<br>`src/representations_2i.py` | `tests/test_binary_icosahedral.py`<br>`tests/test_representations_2i.py` | `docs/PELLIS_A_2I_SPECTRAL_GEOMETRY.md`<br>`results/json/quaternion_bridge_2i_audit.json` | Retain action conventions and exact Reynolds/Molien multiplicities. |
| S48 | Two-sheet self-adjoint operators | SU(2) intertwiner PROVED; edge exchange CONDITIONAL; operator domain-specific | `src/pellis_su2_fibre.py`<br>`src/psu_two_sheet_dirac.py` | `tests/test_pellis_su2_fibre.py`<br>`tests/test_psu_two_sheet_dirac.py` | `docs/PELLIS_B_TWO_SHEET_SELF_ADJOINT_OPERATORS.md` | Specify dense domain, adjoint blocks, invariant domain and weight symmetry. |
| S49 | D6 channel intertwiner | Commutant FINITE-EXACT; dipole ℓ=1 CANDIDATE; physical selection OPEN | `src/pellis_d6_commutant.py`<br>`src/fibonacci_d6_channel_intertwiner.py` | `tests/test_pellis_d6_commutant.py`<br>`tests/test_fibonacci_d6_channel_intertwiner.py` | `docs/PELLIS_C_D6_CHANNEL_INTERTWINER.md` | D6 allows independent 1 and 2 blocks; does not force triplet degeneracy. |
| S50 | RGS–PSU–PRISM–Pellis bridge discipline | OPEN physical maps; formal separation established | `src/rgs_physics_architecture.py` | `tests/test_rgs_physics_architecture.py` | `docs/PELLIS_D_RGS_PSU_PRISM_INTEGRATION.md` | For every bridge give Hs,Ht,T and preserved property; dimensions/spectra are not maps. |
| S51 | D6 splitting and possible triplet degeneracy | FINITE-EXACT commutant; enhanced-symmetry degeneracy CONDITIONAL | `src/pellis_d6_operator_gate.py`<br>`src/pellis_d6_commutant.py` | `tests/test_pellis_d6_operator_gate.py`<br>`tests/test_pellis_d6_commutant.py` | `docs/PELLIS_C_D6_CHANNEL_INTERTWINER.md` | Triplet degeneracy needs a larger irreducible symmetry plus invariant operator/domain. |

## Coverage interpretation

1. **Formal statements:** the canonical statements, hypotheses and dependency descriptions remain in [the integrated formalization register](INTEGRATED_MATHEMATICAL_SUITE_FORMALIZATION_2026-11-25.md). This matrix is an audit layer, not a second theorem source.
2. **Tests:** listed test modules are present in the repository tree. No test execution is claimed by this document. Run the repository's normal Python test suite on the exact commit before treating software validation as current.
3. **Certificates:** JSON audit outputs and matrix/rank artifacts are computational evidence. For claims marked PROVED, the proof must remain in the mathematical text; a floating-point residual is not a proof.
4. **Open bridges:** S18, S26, S37, S39, S43, S48–S51 retain explicit open or conditional parts. Do not promote them based solely on passing tests.
5. **Cross-model discipline:** Q6 (64 vertices), RGS operator (188D), PRISM (69,120 vertices), Q-ball BVP and Pellis/GSR resonance operator remain separate until a concrete map and its preserved structure are proved.

## Recommended completion order

1. Close the exact finite certificates for S07–S10 and S24, including integer boundary matrices/ranks and signed chain-map identities.
2. Pin RGS seed provenance and exact-vs-floating symmetry evidence for S11–S15.
3. Finish the curved common-`L²` reconstruction and Mosco liminf/recovery proof obligations for S16–S18.
4. Preserve explicit wiring and exhaustive enumeration for S21–S22; separate chain-map, injectivity and metric-canonicality gates for S23–S28.
5. Attach integral formulas and rank/equivariance certificates for S29–S36; treat S37's asymptotic limit as open until a global minimizer and uniform remainder proof exist.
6. Keep the De Aquino bridge and physical momentum normalization open for S39–S40.
7. Separate algorithmic Monte Carlo time from physical time and complete the equilibrium/FDT/limit assumptions before any viscosity claim in S45–S46.
8. For Pellis S47–S51, archive exact representation conventions, self-adjoint domains, commutant blocks and the additional symmetry needed for triplet degeneracy.

## Reproducibility contract

For every future status change, update this matrix and the formalization register together, and provide: (i) exact commit SHA; (ii) source and test paths; (iii) exact command and output; (iv) input/result hashes; (v) proof or certificate; (vi) tolerance where floating-point arithmetic is used; (vii) explicit scope limitations.

Do not change a status from OPEN/CANDIDATE/CONDITIONAL to PROVED merely because a test passes.


## Post-snapshot supplements — PRs #233 and #234 (merged 2026-11-27)

The 51-row table above is preserved as the original 2026-11-26 snapshot. The following merged companion-repository work supplements existing records; it does not create duplicate theorem records or renumber S01–S51.

### S35 — GEC-6 lower/upper trigram decomposition (PR #233)

- Canonical dossier: `docs/research/GEC6_TRIGRAM_BLOCK_DECOMPOSITION.md`.
- Implementation/test: `src/gec6_trigram_decomposition.py`; `tests/test_gec6_trigram_decomposition.py`.
- Exact mathematical content: `Λ²(Z⁶)=Λ²L ⊕ (L⊗U) ⊕ Λ²U`, ranks 3+9+3; block symmetry `(S3×S3)⋊C2` of order 72; the mixed block maps onto `A5` with rank-four kernel.
- Review correction: a general bit translation acts on the oriented pair generator by `(-1)^(t_i+t_j)`. It need not preserve `ker(c)` and does not generally descend to `A5`. The global antipode `t=63` has sign +1 on every pair and therefore acts trivially on the quotient. The dossier and tests retain this distinction.
- Status: exact algebraic statements under the stated orientation convention; the finite test is an implementation audit, not a substitute for the contraction/kernel proof. I-Ching semantics and physical coupling remain OPEN.

### S11 — deterministic RGS–PSU vierbein/tetrad artifacts (PR #234)

- Generator: `examples/reconstruct_rgs_vierbein_tetrads.py`.
- Outputs: `results/vierbein_reconstruction/` (six `.npy` arrays plus JSON report and SHA-256 manifest).
- Test: `tests/test_rgs_vierbein_reconstruction.py`.
- Report: `results/vierbein_reconstruction/vierbein_reconstruction_report.json`; hashes: `results/vierbein_reconstruction/SHA256SUMS.txt`.
- The 64 spatial coframes are deterministic reconstructions from the SHA-256-checked seed using `E_i=principal_sqrt(C_i/c0)`. The report's maximum condition number for `E_i` is 1.733712…; its square, 3.00576…, is the corresponding condition number for `C_i/c0`, reconciling the two conventions in the historical status note.
- The RGS 4D array is only the static lift `diag(1,E_i)`. Ellis and Kerr–Schild arrays are separate analytic benchmark coframes, not solutions of the coupled RGS field equations. The report explicitly leaves TEGR/torsion `NOT_CERTIFIED` pending an inertial spin connection and convention audit.
- Status: deterministic reconstruction and finite numerical consistency checks; not byte-for-byte recovery of historical arrays, not a dynamical tetrad derivation, and not a TEGR proof.

### CI evidence for PRs #233–#235

No pull-request-triggered workflow runs were returned for the three merge commits checked in this audit. This is recorded as **no CI run evidence retrieved**, not as a test failure or a passing test result. The exact reproduction commands remain in the companion dossiers.
