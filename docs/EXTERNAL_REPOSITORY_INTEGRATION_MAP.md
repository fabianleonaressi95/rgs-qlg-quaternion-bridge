# External repository integration map — RGS–PSU–PRISM–QLG

**Audit date:** 4 November 2026  
**Primary integration target:** [rgs-qlg-quaternion-bridge](https://github.com/fabianleonaressi95/rgs-qlg-quaternion-bridge)

This map records repositories actually discovered and inspected through the linked GitHub account and public GitHub repository metadata. It distinguishes code-bearing projects from empty repositories, research notes, demonstrations and external collaborator work. A thematic overlap is not a mathematical identification.

## 1. Executive map

| Priority | Repository | Concrete overlap | Recommended treatment |
|---|---|---|---|
| P0 | [rgs-qlg-quaternion-bridge](https://github.com/fabianleonaressi95/rgs-qlg-quaternion-bridge) | Canonical RGS/Q6/PRISM/2I/Clifford/vierbein/EC/TEGR/QLG implementation | Keep as integration hub and canonical conventions |
| P1 | [kilaw](https://github.com/fabianleonaressi95/kilaw) | 64B electrostatics, teleparallel torsion prototype, EMS wormhole spectra, QLG, QFI, quantum memory, proton/deuteron double wells | Extract audited modules and datasets into versioned adapters |
| P1 | [vfd-crystallisation](https://github.com/vfd-org/vfd-crystallisation) | 600-cell, binary icosahedral group 2I, spectral closure, Lindblad baseline, quantum recovery and discrete event-order gravity | Independent mathematical comparison; build explicit 2I/600-cell bridge |
| P1 | [torsion-lambda-toolkit](https://github.com/dannyeldridge54/torsion-lambda-toolkit) | Einstein–Cartan-inspired cosmological background, torsion parameter, energy-condition and limiting-case checks | Compare conventions and background equations; do not import as EC derivation |
| P1 | [the-twist-torsion-cosmology](https://github.com/dannyeldridge54/the-twist-torsion-cosmology) | Barbero–Immirzi-motivated scalar, Holst/Nieh–Yan/supergravity context, preregistered cosmological tests | Use as external phenomenology and falsification interface |
| P2 | [AEGIS](https://github.com/dannyeldridge54/AEGIS) | Working cosmology fits, AxiCLASS patches and optimizer | Use the cosmology records and corrected-code audit; keep optimizer separate |
| P2 | [naressi-golden-resonance](https://github.com/fabianleonaressi95/naressi-golden-resonance) | Golden-ratio modulation, aperiodic photon–phonon toy model | Rebuild from explicit Hamiltonian/action and independent data before scientific use |
| P2 | [cosmodex](https://github.com/fabianleonaressi95/cosmodex) | Fibonacci, cut-and-project and Ammann–Beenker-like point-set generators | Reuse only as geometric dataset prototypes after exact geometry audit |
| P2 | [ASTRAQ](https://github.com/fabianleonaressi95/ASTRAQ) | Perturbative nonlinear-response fit scripts; separate orbital SSA demonstrator | Reuse generic fit/audit utilities only; no direct physical identification |
| P2 | [CosmoExplorer](https://github.com/dannyeldridge54/CosmoExplorer) | Cosmological model interface, data likelihood and parameter scan tooling | Historical comparison harness; reconcile against corrected AEGIS/Twist results |
| P3 | [naressi-zuliani-law-quantum-morphodynamics](https://github.com/fabianleonaressi95/naressi-zuliani-law-quantum-morphodynamics) | Toy VBS-vs-QCD jet classifier and information observables | ML/phenomenology sandbox only; not a QCD field model |
| P3 | [q-proof](https://github.com/fabianleonaressi95/q-proof) | Golden-angle graph/consensus toy simulation | Keep outside physics core; benchmark algorithmic claims separately |
| P3 | [quantum-lindblad](https://github.com/fabianleonaressi95/quantum-lindblad) | Repository name suggests open quantum systems, but current main simulation file is empty and downloader targets E. coli CDS data | No usable Lindblad implementation identified |
| P3 | [biverse](https://github.com/fabianleonaressi95/biverse) | Name/README assertion about PSU–Anti-Twin and Yang–Mills gap | Empty repository; recover source and proofs before any integration |
| P3 | [cosmodex-pwa](https://github.com/fabianleonaressi95/cosmodex-pwa) | Companion app name | Empty repository |
| P3 | [astra-q-core](https://github.com/fabianleonaressi95/astra-q-core) | Companion core name | Empty repository |

P0/P1 mean high-priority technical bridges, not validation of the physical claims. P2 are conditional tools or model candidates. P3 are currently non-integrable, unrelated, empty or insufficiently specified.

## 2. Detailed internal repository map

### 2.1 KILAW — highest-value internal source archive

Repository: [fabianleonaressi95/kilaw](https://github.com/fabianleonaressi95/kilaw)

The repository is a collection of extensionless Kaggle notebook blobs and CSV outputs rather than a conventional package. The following source blobs were inspected:

| Path | Existing content | Bridge to RGS–PSU–QLG | Gate before reuse |
|---|---|---|---|
| [13](https://github.com/fabianleonaressi95/kilaw/blob/main/13) | 64B point charges, neutralising cubic background, Duffy quadrature, Hessian and relaxation tests | Compare the 64-site electrostatic equilibrium/Hessian with the RGS 64-site coordinates and 188D tangent Hessian | Verify site ordering, units, boundary conditions, neutrality and whether coordinates are identical |
| [naressilaw](https://github.com/fabianleonaressi95/kilaw/blob/main/naressilaw) | Pauli/Clifford-like 4×4 generators, plaquette gauge-covariance tests, gauge Hessian, aperiodic geometry and a Weitzenböck-style torsion prototype | Compare gamma conventions, SU(2) link/plaquette transformations and discrete torsion with current `clifford_qlg.py`, `holonomies.py`, `rgs_spin_connection.py` and future TEGR module | Re-derive each transformation law; distinguish internal gauge curvature from Lorentz/spacetime torsion |
| [qlgloop](https://github.com/fabianleonaressi95/kilaw/blob/main/qlgloop) | EMS wormhole background, radial coupled Sturm–Liouville spectra and axial-sector diagnostics | Compare with the RGS Ellis/EMS wormhole and PSU–Anti-Twin throat programme | Match the action, charges, radial gauge, asymptotics and self-adjoint boundary conditions |
| [entanglement](https://github.com/fabianleonaressi95/kilaw/blob/main/entanglement) | Gaussian quantum-quench vs classical matched-kernel memory audits | Quantum-information diagnostics for the two-sheet PSU sector | Compare covariance conventions, matched kernels and classical null controls |
| [qfilangevin](https://github.com/fabianleonaressi95/kilaw/blob/main/qfilangevin) | Quantum Fisher information via symmetric logarithmic derivative and Koopman-density toy models | State distinguishability and sensitivity diagnostics for PSU/Anti-Twin dynamics | Audit positivity, trace, SLD residual and identifiability; avoid hard-coded target exponents |
| [quantum](https://github.com/fabianleonaressi95/kilaw/blob/main/quantum) | Three-qubit strong-collision/QME dataset and full/partial/hidden QFI observables | A controlled finite-dimensional quantum-information testbed | Check exact Hamiltonian, subsystem ordering, state preparation and data-generation provenance |
| [morphodynamic](https://github.com/fabianleonaressi95/kilaw/blob/main/morphodynamic) | Proton/deuteron double-well propagation, hydrogen-bond deformation and TSDE/double-cone matching experiments | Matter-sector phenomenology and transition-rate tests | Keep the effective potential distinct from a derived RGS/PSU microscopic action |
| [algebra](https://github.com/fabianleonaressi95/kilaw/blob/main/algebra) | Exact finite matrices, modular factorisation, Fibonacci words and cocycle scans | Algebraic/spectral audit tools | Extract exact statements; separate exploratory scans from proved identities |
| [kiki](https://github.com/fabianleonaressi95/kilaw/blob/main/kiki) | Exact rational block matrices and characteristic-polynomial work | Candidate finite-operator comparisons | Identify the intended operator and invariant bilinear form before mapping |
| [new](https://github.com/fabianleonaressi95/kilaw/blob/main/new) | 576 mutation combinatorics and chemical/tautomer datasets | Possible biological/chemical application of the morphodynamic programme | Separate biological chemistry from fundamental gauge/geometry claims |
| [v1](https://github.com/fabianleonaressi95/kilaw/blob/main/v1) | Falsifiable scale-model and ablation experiments | Generic model-validation methodology | Reproduce from raw data; do not import fitted targets |
| [13 plus CSV outputs](https://github.com/fabianleonaressi95/kilaw/tree/main/figures) | Physical-constrained fits and residual data for the Krasilnikov analysis | Provenance and independent H/D model comparison | Check source datasets, model assumptions and held-out predictions |

**Priority extraction order:** `13` → `naressilaw` → `qlgloop` → `entanglement` → `qfilangevin` → `quantum` → `morphodynamic`.

First create a provenance manifest for every notebook: Git blob SHA, notebook cell hashes, dependencies, input files, units, output artifacts and status. Then migrate one audited model at a time into importable modules under the main repository. Do not bulk-copy notebooks.

### 2.2 VFD Crystallisation — strongest independent geometric comparison

Repository: [vfd-org/vfd-crystallisation](https://github.com/vfd-org/vfd-crystallisation)

Relevant inspected paths:
- `src/vfd/crystallisation/operator.py`: finite-dimensional closure functional and gradient-flow implementation.
- `src/vfd/crystallisation/models/lindblad_model.py`: standard Lindblad baseline.
- `papers/paper-xxii/scripts/run_600cell_spectral_verification.py`: construction of 120 600-cell vertices, degree-12 adjacency, graph Laplacian and 2I/H4 spectral checks.
- `papers/paper-xxviii/paper-xxviii.tex`: quantum/gravity structural bridge hypothesis.

**Exact common mathematical object:** the binary icosahedral group (2I), order 120, and the 600-cell geometry. This is the most direct group-theoretic comparison with the RGS programme's 2I/S³/2I sector.

**Do not identify the complexes:** the VFD 600-cell vertex graph has 120 vertices and degree 12; the Q6 hypercube has 64 vertices, 192 edges and 240 square faces. A bridge requires an explicit equivariant map, common representation basis or shared operator construction. No natural bijection is implied.

**Important scope note:** VFD Paper XXVIII explicitly describes its quantum/gravity correspondences as analogies on different spaces and says it does not prove a common Laplacian, common closure functional or unification theorem. That restraint should be preserved. The repository's claims about mass ratios, couplings, charge radii and Yang–Mills should be independently audited; repository presence and passing local tests are not external validation.

### 2.3 Danny Eldridge — torsion and cosmology repositories

#### TorsionΛ toolkit
Repository: [dannyeldridge54/torsion-lambda-toolkit](https://github.com/dannyeldridge54/torsion-lambda-toolkit)

Useful paths:
- `src/models/torsion-lambda.ts`
- `src/physics/einstein-cartan-nonlinear.ts`
- `src/physics/running-torsion-extension.ts`
- `src/physics/torsion-boltzmann-sources.ts`
- `src/physics/torsion-radiative-transfer.ts`
- `examples/config-einstein-cartan.yaml`

The README describes it as a **pre-simulation consistency toolkit**, not a full Boltzmann solver or observationally validated theory. Its phenomenological `beta` correction to background density and its effective-fluid energy-condition checks must not be equated with the algebraic Einstein–Cartan relation
(T_{abc}=-(\kappa/2)\epsilon_{abcd}J_5^d).
The model's `beta\,\rho_s` term needs an independent derivation from an action and a specified spin ensemble before it can be coupled to the RGS EC–Dirac sector.

#### The Twist — preregistered BI-field cosmology
Repository: [dannyeldridge54/the-twist-torsion-cosmology](https://github.com/dannyeldridge54/the-twist-torsion-cosmology)

Useful records:
- `CLAIMS_INDEX.md`
- `REPRODUCE.md`
- `preregistrations/`
- `results_notes/`
- `code/axiclass/`
- `theory/` (BI-supergravity and Spiral Twist solvers)

The 2 November 2026 claims index reports the corrected-code verdicts: the late field ties ΛCDM on Pantheon+, has a weak DES Y5 preference that weakens on alternative supernova samples and is worse than (w_0-w_a); the coupled model is also worse than (w_0-w_a), and the torsion-natural coupling value (eta=0.816) is excluded. The minimal field does not explain the measured CMB birefringence. These are useful falsification records, not support for a completed unified theory.

**Critical theoretical bridge:** the BI scalar/supergravity sector is dynamical only under its specified extended action and couplings. Minimal Einstein–Cartan torsion remains algebraic. Derive the Holst/Nieh–Yan/scalar action and its connection equation before merging this with the EC–Dirac contact sector. Do not promote a phenomenological cosmological scalar to propagating EC torsion by notation alone.

#### AEGIS
Repository: [dannyeldridge54/AEGIS](https://github.com/dannyeldridge54/AEGIS)

This is a working repository containing the cosmology fit pipeline and a separate generic optimizer. Use the cosmology code, corrected AxiCLASS patch, preregistrations and results as the scientific audit trail. Keep the optimizer as generic tooling, not as a physics result.

#### NEXUS
Repository: [dannyeldridge54/NEXUS](https://github.com/dannyeldridge54/NEXUS)

Generic parameter-space exploration, anomaly detection and model scanning. Potentially reusable as an external driver for RGS parameter audits, provided objective functions, constraints and uncertainty treatment are implemented in the main project. It contributes no physical equations.

#### CosmoExplorer
Repository: [dannyeldridge54/CosmoExplorer](https://github.com/dannyeldridge54/CosmoExplorer)

Cosmological model comparison and sampling harness with UFE torsion models. Its README contains older headline significance claims that do not match the later corrected-code/pre-registered verdicts in AEGIS and The Twist. Treat it as a historical software baseline; reconcile equations, likelihoods, data versions and fit outputs before quoting any result.

#### Other Danny repositories
- [CubeMind-Simple](https://github.com/dannyeldridge54/CubeMind-Simple) — empty repository at inspection.
- [SOUL](https://github.com/dannyeldridge54/SOUL) — empty repository at inspection.

### 2.4 Other repositories under the user's account

#### Naressi Golden Resonance
[naressi-golden-resonance](https://github.com/fabianleonaressi95/naressi-golden-resonance)

`naressi_law.py` implements a golden-angle photon–phonon toy signal, with (arphi)-dependent frequencies and an entropy signal constructed around a preset scale. The frequencies and damping are inserted into the synthetic generator, so an FFT peak ratio from that generator is not an independent discovery of golden selection. Potential bridge: formulate an explicit Hamiltonian/open-system model, derive the couplings, and test against independently generated or measured data.

#### COSMO-DEX
[cosmodex](https://github.com/fabianleonaressi95/cosmodex)

`quasicrystal_data.py` generates Fibonacci-word, cut-and-project and approximate Ammann–Beenker point sets. This can supply geometric test datasets for aperiodic graph/refinement experiments. Its current 2D point sets are not the 3D icosahedral RGS seed or the 4D 600-cell; construct exact coordinate and symmetry audits before identifying them.

#### ASTRAQ
[ASTRAQ](https://github.com/fabianleonaressi95/ASTRAQ)

Contains a separate orbital space-situational-awareness demonstrator and a universal nonlinear scaling-law analysis. The scaling scripts fit
(I(\epsilon,\omega,\phi_0)=\epsilon^2 K_2+O(\epsilon^3)).
This is a reusable perturbative-response analysis pattern, not a derivation that the same exponent or kernel governs RGS, Q-ball, gauge, or gravitational sectors. Require raw data, independent uncertainty estimation, residual analysis and a controlled (epsilon\to0) limit.

#### Naressi–Zuliani quantum morphodynamics
[naressi-zuliani-law-quantum-morphodynamics](https://github.com/fabianleonaressi95/naressi-zuliani-law-quantum-morphodynamics)

The inspected `toy_monte_carlo.py` generates synthetic VBS-like signal and QCD-like background jet kinematics, computes Bessel-weighted modes, IPR, von Neumann entropy and ROC/AUC. This is an ML/phenomenology toy pipeline, not lattice QCD, a QCD action or a derivation of hadronic dynamics. It may be used for classifier stress tests after removing any claim that synthetic separation validates the physics.

#### Quantum Lindblad
[quantum-lindblad](https://github.com/fabianleonaressi95/quantum-lindblad)

At inspection, `main_simulation.py` was empty; `download.py` downloads E. coli CDS data and the repository contains an E. coli sequence file. No working Lindblad solver was identified. Do not treat this repository as the quantum open-system implementation; use the explicit VFD Lindblad baseline or the main project's audited quantum modules instead.

#### Q-PROOF
[q-proof](https://github.com/fabianleonaressi95/q-proof)

The inspected simulation constructs a golden-angle spatial network but assigns headline capture probability, TPS, finality and energy metrics as fixed protocol constants. It is a protocol/consensus toy, not a physics substrate. Keep separate from RGS physical claims.

#### Biverse and empty companion repositories
- [biverse](https://github.com/fabianleonaressi95/biverse) — README states a PSU/Anti-Twin Yang–Mills-gap result, but the repository has no code or proof files at inspection. Recover the actual derivation and reproducible artifacts before linking it as evidence.
- [cosmodex-pwa](https://github.com/fabianleonaressi95/cosmodex-pwa) — empty.
- [astra-q-core](https://github.com/fabianleonaressi95/astra-q-core) — empty.

## 3. Scientific integration graph

The arrows below indicate candidate data/derivation interfaces, not established physical equivalences.

```text
                         VFD: 600-cell / 2I / closure
                                      |
                         common group / representation audit
                                      |
RGS coordinates -> Q6/PRISM -> coframe E -> LC / EC / TEGR connection
       |              |                    |
       |              |                    +--> spinor D_mu, EC-Dirac
       |              +--> SU(2) links / plaquettes / Wilson action
       |                                   |
       +--> 64B electrostatics (KILAW)     +--> U(1), SU(2), SU(3) sectors
       |                                        |
       +--> Ellis refinement / EMS (KILAW) -----+--> matter / stress / Noether
                                                |
PSU–Anti-Twin <--> QFI / Gaussian memory / QME (KILAW)
       |                    |
       +--> Q-ball / QAT ---+--> open-system and response tests
                                  |
BI/Holst/Nieh–Yan scalar (The Twist) -- action-level bridge required --> EC/TEGR
                                  |
VFD crystallisation / Lindblad ------> independent selection/decoherence baseline

Aperiodic phi generators (Golden Resonance, COSMO-DEX)
       --> geometric test fixtures only until a physical action is derived

ASTRAQ K2, NEXUS, CosmoExplorer
       --> analysis/scan/likelihood tooling, not physical sectors
```

## 4. Canonical integration contracts

Every migrated module must implement or document these interfaces:

1. **Geometry:** coordinate arrays, vertex ordering, units, orientation, metric signature, frame gauge and refinement level.
2. **Cochains:** exact integer (B_1,B_2,B_3), dimensions and (B_1B_2=B_2B_3=0).
3. **Frames:** (e^a{}_{\mu}), inverse (e_a{}^\mu), first jet (\partial_\nu e^a{}_{\mu}), Lorentz gauge transition and tetrad-postulate residual.
4. **Connections:** separate records for LC spin connection, EC independent connection/contorsion, TEGR inertial spin connection and internal gauge connections.
5. **Matter:** explicit representation, gamma convention, covariant derivative, action, Euler operator, currents and stress tensor.
6. **Quantum states:** Hilbert-space basis/order, density-matrix positivity, trace, Hamiltonian, dissipators, QFI/SLD residuals and boundary conditions.
7. **Spectral data:** operator definition, inner product, self-adjointness/Hermiticity, spectrum, multiplicities, zero modes, convergence and error bars.
8. **Cosmology:** action-derived background equations, units, parameter dictionary, likelihood/data versions, priors, preregistration hashes and corrected-code provenance.
9. **Provenance:** source repository + commit/blob SHA, exact source path, licence, migration commit, tests and evidence label.

## 5. Immediate work programme

### A — Recover and audit KILAW source material
- [ ] Materialise and hash the extensionless notebooks; create a machine-readable inventory.
- [ ] Reconstruct the 64B coordinate/order mapping and compare with the canonical RGS seed.
- [ ] Port the 64B Coulomb/Hessian audit as an independent cross-check, not as replacement data.
- [ ] Rebuild the `naressilaw` plaquette and torsion examples with explicit index and gauge conventions.
- [ ] Compare the EMS wormhole background and spectra with the Ellis-RGS implementation.
- [ ] Extract QFI, Gaussian memory, QME and double-well experiments into separate test modules.

### B — Build the 2I/600-cell comparison with VFD
- [ ] Independently regenerate the 120 600-cell vertices and 2I action.
- [ ] Compare group multiplication, central element, irreducible representations and character projectors.
- [ ] Identify whether any RGS 64-site/Q6 object admits a natural 2I-equivariant map; do not assume one.
- [ ] Compare VFD closure functional with Krist-All and PRISM operators on a common, explicitly defined space.
- [ ] Reproduce VFD's finite spectral claims independently before citing them.
- [ ] Preserve the VFD paper's own distinction between structural analogy and theorem.

### C — Build the EC cosmology bridge with Danny's work
- [ ] Compare TorsionΛ background equations with the exact EC spin-current convention used in the main repository.
- [ ] Derive the BI/Holst/Nieh–Yan scalar action and connection equation from the same declared conventions.
- [ ] Check the minimal EC algebraic torsion limit and the extended propagating-scalar limit separately.
- [ ] Re-run the preregistered corrected-code checks in The Twist before using cosmological fit claims.
- [ ] Reconcile CosmoExplorer historical UFE claims with AEGIS/The Twist results.

### D — Use auxiliary tooling without confusing it with physics
- [ ] Audit the Golden Resonance synthetic generator and replace imposed golden peaks with blind predictions.
- [ ] Test COSMO-DEX point sets against exact cut-and-project and symmetry criteria.
- [ ] Validate ASTRAQ scaling fits with raw data and independent uncertainty estimates.
- [ ] Keep NEXUS/CosmoExplorer as scan/likelihood drivers with a strict model interface.

## 6. Licensing and collaboration constraints

- The linked GitHub account has read access but no push permission to `vfd-org/vfd-crystallisation` and Danny Eldridge's public repositories. External changes require the maintainers' approval.
- `torsion-lambda-toolkit` declares LGPL-3.0; preserve its licence and notices if code is ever reused.
- `the-twist-torsion-cosmology` declares Hippocratic License 3.0 with an additional restrictive addendum. Do not copy its code into the main repository without explicit licence compatibility review. Reproduction and citation can be performed without code migration.
- VFD repository licensing must be checked at the exact files/release before code reuse.
- For the user's repositories, inspect each repository's licence before copying code; repository ownership does not itself specify the licence of every artifact.
- Prefer independent reimplementation of mathematical equations with source attribution, and keep external code/data as separately versioned dependencies when licensing is unclear.

## 7. Evidence policy

- **PROVED:** mathematical derivation with stated hypotheses and conventions.
- **NUMERICALLY VERIFIED:** reproducible calculation with independent tests, residuals and tolerances.
- **CANDIDATE:** explicit construction not yet fully validated.
- **HYPOTHESIS:** proposed physical interpretation or coupling.
- **OPEN:** missing derivation, data, continuum limit or falsification test.

A repository's README, paper, passing self-tests or numerical headline is not independent evidence by itself. No cross-repository link becomes a physical coupling until a common action/operator, map between spaces, compatible conventions and the corresponding field equations or convergence statement have been supplied.
