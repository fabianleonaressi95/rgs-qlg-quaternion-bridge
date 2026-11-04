# KILAW source audit: `13` and `naressilaw`

**Audit date:** 4 November 2026  
**Source repository:** [fabianleonaressi95/kilaw](https://github.com/fabianleonaressi95/kilaw)  
**Target:** `rgs-qlg-quaternion-bridge`  
**Scope:** static source inspection and review of outputs already stored in the Kaggle notebooks. This is not a clean-room rerun of the notebooks.

## 1. Immutable source provenance

| Notebook path | Git blob SHA | Format | Cells |
|---|---|---|---:|
| [13](https://github.com/fabianleonaressi95/kilaw/blob/667ccfb54969c3cf0d371e27de4f3817b33be3c7/13) | `667ccfb54969c3cf0d371e27de4f3817b33be3c7` | nbformat 4 JSON, extensionless | 6 |
| [naressilaw](https://github.com/fabianleonaressi95/kilaw/blob/e8fde365216fdb2dd2050009573b56aaa9925cc6/naressilaw) | `e8fde365216fdb2dd2050009573b56aaa9925cc6` | nbformat 4 JSON, extensionless | 104 |

The blob SHA values identify the exact Git objects inspected. The notebook cells contain a mixture of executable code, saved outputs, exploratory calculations and repeated revisions. Saved outputs are historical evidence only; they do not establish that the current notebook executes top-to-bottom in a clean environment.

## 2. Notebook `13`: 64B charges and neutralising cubic jellium

### 2.1 What is actually implemented

The notebook defines a (4\times4\times4) point set (64 charges), pairwise Coulomb energy/force/Hessian, a uniform neutralising cubic background, volume integrals, finite-difference derivatives, Hessian spectra and a constrained Coulomb relaxation. Units are dimensionless with (a=q=(4\pi\epsilon_0)^{-1}=1); physical energy and Hessian scales are stated in the source.

### 2.2 Critical source revisions and findings

1. **Cell 1 is not a reliable baseline.** It uses an analytic rectangular-prism corner formula and finite differences. Its saved output reports enormous background forces (order (10^6)) and Hessian eigenvalues (order (10^9)), with a huge breathing residual. The following cell explicitly calls out an incorrect signed-corner implementation. These results must not be reused.

2. **Cell 2 replaces the potential/gradient calculation with Duffy quadrature**, but explicitly discovers a geometry issue: for (n=4), the coordinates are (pm1.5) and the background cube is ([-1.5,1.5]^3). Thus 56 of the 64 sites lie on the cube boundary (the 8 corners are vertices; the rest include edges/faces). The boundary value of a volume potential may be defined, but its derivative requires a declared one-sided/boundary convention. The notebook itself warns that the original setup needs a boundary convention and that the regulated family (L=L_{sites}+\eta), \(\eta>0\), must be extrapolated.

3. **Cell 3 contains a useful six-sector Duffy map** for the integrals (\int d^3r/|r|) and (\int r_i,d^3r/|r|^3). It is an integration primitive, not by itself a validation of the full model or its boundary limit.

4. **Cell 4 is the most careful jellium revision.** It fixes a background cube larger than the site extent, (L_{background}=L_{sites}+\eta), and resets (\rho_b=-N/(2L_{background})^3) for neutrality. Saved outputs show converged background self-energy quadrature at fixed \(\eta=0.03\) (orders 16–64 agree to roughly (10^{-11})); the representative corner-site gradient also stabilises with quadrature order. However, this is quadrature convergence at fixed regulator, not proof of the \(\eta\to0^+\) limit.

5. **The stored (n=4) regulated results do not establish equilibrium or stability.** For \(\eta=0.1,0.03,0.01,0.003\), the maximum total force remains nonzero (about 22.0, 26.5, 28.0, 28.5 in the notebook's units). The Hessian has negative modes (8 at \(\eta=0.1\), 32 at the smaller listed regulators), and its minimum eigenvalue becomes more negative. Therefore the unrelaxed 64B lattice in this fixed-background model is not shown to be a stable equilibrium. A zero-point-energy calculation based on a purported stable spectrum is not justified for these outputs.

6. **The exact scaling statement is conditional and limited.** If all lengths scale together at fixed dimensionless regulator, the classical Coulomb terms scale as (a^{-1}). That dimensional scaling does not prove mechanical equilibrium, a cosmological equation of state for the physical system, or a quantum vacuum equation of state.

7. **Cell 5 is a separate constrained relaxation** of the 64 like-charge system with centre of mass and moment-of-inertia constraints. It omits the neutralising background. It is therefore a different variational problem and cannot be used as a solution of the jellium model without a separate derivation connecting the functionals.

### 2.3 RGS integration decision

- **REUSABLE after isolated tests:** point-charge pair energy, force and Hessian; Duffy volume-integral primitive; explicit units and scaling checks.
- **CANDIDATE:** regulated jellium energy/force/Hessian as a separate external-potential model.
- **OPEN:** rigorous boundary prescription; convergence/extrapolation \(\eta\to0^+\); stationarity; constrained Hessian after removing the correct symmetry/constraint directions; stability; comparison with the actual RGS 64-site coordinates and 188D tangent Hessian.
- **DO NOT IMPORT:** cell-1 spectrum/forces, any stable-vacuum or ZPE claim from the unstable spectrum, or cell-5 relaxation as if it solved the same jellium functional.

The current RGS coordinates must be compared by an explicit coordinate map and ordering audit. Both objects having 64 sites is not enough: the KILAW point set is a Cartesian (4^3) grid, whereas the RGS (O_h) reconstruction has its own coordinate data and tangent constraints.

## 3. Notebook `naressilaw`: Clifford-like matrices, links, coframes and Q6

### 3.1 Algebraic and lattice-gauge sections

- Early cells construct the 16 tensor-Pauli matrices in (M_4(\mathbb C)). Cell 1 assigns a coarse grading that gives every non-identity Pauli factor the same grade; cell 2 revises this with a (\mathbb Z_2^2) degree map. These are exploratory graded matrix bases, not yet the repository's Lorentzian Clifford convention.
- The reported closure test projects every graded commutator back onto all 16 matrices. Since those matrices span all of (M_4(\mathbb C)), closure in that ambient vector space is expected. It does not establish a new Lie superalgebra, a physical spin representation, or compatibility with the RGS gamma matrices. A proper adapter must compare explicit multiplication tables, grading, adjoint, metric signature and gamma-five convention.
- The link/plaquette test implements (U_{xy}\mapsto\Omega_xU_{xy}\Omega_y^{-1}). Its saved trace difference is (1.36\times10^{-15}). This is a valid finite-dimensional telescoping/cyclic-trace sanity check for the implemented closed loop.
- The separate (3\times3) periodic (U(1)) phase model has 18 link variables and 9 plaquette terms. Its Hessian has 10 numerical zero modes. For the periodic incidence/curl map, rank 8 and nullity 10 are expected: local gauge directions plus two global torus holonomies. The output's description of all ten as simply gauge modes is too coarse; report rank, gauge orbit and harmonic/global modes separately.
- The synthetic metric reconstruction uses randomly generated link vectors and lengths computed from a chosen metric, then solves the same linear system. The near-machine-precision residual is an algebraic round-trip test, not empirical reconstruction of a spacetime metric.

### 3.2 Coframe and torsion prototype

- Cell 6 constructs an aperiodic 64-vertex triangulation, induces local coframes and reports nondegenerate local metric determinants. Its discrete Cartan torsion is near machine zero.
- Cell 7 sets a flat reference spin connection (\(\omega=0\)) but still reports torsion at roughly (10^{-17}). This is consistent with the specific discrete coframe construction being effectively closed/exact at the tested scale; it is not evidence for a nonzero teleparallel torsion field.
- Cells 8 onward explore a modulated/anholonomic coframe and additional discrete torsion constructions. They need a separate orientation/sign/edge-ordering audit and comparison against the continuum identity (T^a=de^a+\omega^a{}_b\wedge e^b). A discrete boundary operator must satisfy the relevant chain identity, and the discrete connection/coframe must be declared before calling the result Weitzenböck/TEGR.
- **Critical distinction:** an internal SU(2) link connection and a spacetime Lorentz/Weitzenböck connection are different bundles and transformation laws. A plaquette Wilson loop is not spacetime torsion.

### 3.3 Q6 and 2I section

- The notebook attempts an embedding of the 64 Q6 states into (S^3\) and comparison with a generated binary-icosahedral target. Its own output says **216 elements generated where 120 were expected**. This fails the basic order check for (2I), so the reported angular-distance table is not an audit of the binary icosahedral group.
- The 64-state embedding and edge-angle statistics are exploratory geometric data. They do not define a group homomorphism, an equivariant map, a quotient (S^3/2I), or an operator intertwiner.
- Any future comparison must first construct exactly 120 distinct unit quaternions, verify closure/inverses/identity and the group multiplication table, then define the Q6-to-orbit map and quantify equivariance and spectral intertwining. The VFD 600-cell and RGS 2I objects remain distinct until this is done.

### 3.4 Q6 spectra and phenomenological outputs

Several later cells report a central spectral gap around (0.01159616), BdG-like gaps, PMNS-like angles, phonon bands, TOV integrations and frequency/music mappings. These are model-dependent outputs from the notebook's own chosen matrices, parameters and post-processing. The displayed PMNS-like angles (approximately (60.79^\circ,56.27^\circ,29.37^\circ)) do not match the comparison values printed by the notebook itself. The TOV examples include failed/zero-radius configurations in earlier revisions. None of these outputs is promoted here to a physical prediction.

### 3.5 CUORE data and falsification record

The notebook includes a CUORE pulse-shape ML pipeline and later tests of predicted decay-time distributions. One structural audit reports 9,000 waveforms and only one class transition in the original ordering, which is a serious block-ordering warning; a later CNN cell globally shuffles samples and reports high cross-validation AUC. A later audit explicitly reports that label-permutation controls average near 0.502, but this does not resolve possible event/run-level leakage or establish external generalisation.

More importantly, the notebook's own later blind/revised tests report:
- original blind test: KS statistic 0.478, (p\approx9.01\times10^{-104}), **falsified**;
- revised clean subset: KS statistic 0.4659, (p\approx1.27\times10^{-17}), **falsified**;
- DOS model: KS statistic 0.961, (p\) underflowing to 0, **falsified**;
- revised macroscopic thermal model: KS statistic 0.426, (p\approx9.01\times10^{-82}), **falsified**.

These are strong internal negative results against those specific decay-time predictions. They must remain visible in any summary; the high AUC classifier result does not rescue a rejected physical timing model.

## 4. Integration contract for the main RGS repository

No notebook code is copied into production modules by this audit. Migration is staged behind explicit interfaces:

1. **Electrostatics adapter:** `sites -> pair energy/force/Hessian`; separate background potential/gradient/Hessian; explicit boundary regulator and unit metadata.
2. **Gauge adapter:** `oriented links -> gauge transform -> plaquette holonomy`; test covariance of the plaquette matrix and invariance of its trace. Keep the group representation explicit.
3. **Coframe adapter:** `coframe + discrete derivative + connection -> torsion`; enforce index order, orientation and chain identities; expose LC and Weitzenböck choices separately.
4. **Group adapter:** exact (2I) constructor and group audit before any Q6 embedding.
5. **Q6 adapter:** derive operators from the declared graph/complex and compare spectra under explicit maps; do not copy fitted gaps or phenomenological labels.
6. **CUORE adapter:** retain preregistered predictions and all failed tests; use grouped/event-level splits, leakage checks and uncertainty-calibrated held-out inference.

## 5. Status ledger

| Item | Status | Reason |
|---|---|---|
| Coulomb pair-force/Hessian formula in `13` | ACTION_DERIVED / NUMERICALLY_TESTABLE | Standard pair potential; implementation still needs independent unit tests |
| Duffy integral primitive | CANDIDATE | Source inspected; convergence and boundary handling must be tested independently |
| 64B jellium equilibrium/stability | OPEN / NOT ESTABLISHED | Nonzero forces and negative Hessian modes in stored regulated outputs |
| Tensor-Pauli basis spans (M_4(\mathbb C)) | PROVED (finite-dimensional) | 16 linearly independent tensor products |
| Gauge covariance of tested plaquette trace | NUMERICALLY VERIFIED | Saved residual (1.36\times10^{-15}) for one seeded test |
| Ten zero modes in 3x3 periodic U(1) Hessian | NUMERICALLY VERIFIED / ACTION_DERIVED | Rank-8 curl map gives nullity 10; decompose gauge vs harmonic modes |
| KILAW discrete teleparallel torsion | OPEN | Current outputs near zero; modulated coframe requires exact discrete audit |
| KILAW generated (2I) | FAILED BASIC CHECK | Generated count 216, expected 120 |
| Q6-to-(S^3/2I) map | HYPOTHESIS / OPEN | No validated equivariant map or intertwiner |
| CUORE decay-time predictions in cells 69–72 | FALSIFIED IN NOTEBOOK'S OWN TESTS | All four stored KS tests reject their corresponding prediction |

## 6. Next reproducible tasks

- Port the pair Coulomb kernel and Duffy primitive only after creating small analytic unit tests (two-charge force/Hessian, symmetry, finite-difference gradient, quadrature order convergence).
- Reproduce the (\eta)-regulated jellium family and fit/extrapolate forces and Hessian eigenvalues with uncertainty; do not assume the limit exists.
- Rebuild the SU(2) link test using the main repository's declared gamma/Clifford conventions only where the representation is genuinely shared.
- Replace the notebook's 216-element generator with a separately tested exact (2I) constructor; then attempt the Q6 map.
- Preserve the CUORE negative tests as regression/falsification artifacts.

**Audit conclusion:** KILAW contains useful mathematical prototypes and several valuable negative results. The cleanest immediate reusable kernels are pair Coulomb interactions and finite lattice-gauge plaquette covariance. The 64B jellium stability claim, the teleparallel torsion interpretation and the Q6/(2I) identification are not established by these notebooks. The CUORE decay-time models recorded in the later cells are explicitly falsified by their own stored tests.
