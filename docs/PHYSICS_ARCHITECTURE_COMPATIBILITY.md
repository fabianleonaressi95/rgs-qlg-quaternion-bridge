# RGS–PSU–PRISM–QLG: physics architecture and compatibility register

**Status:** repository architecture/control document, not a claim of a completed unified physical theory.  
**Convention baseline:** metric signature \(\eta=(-,+,+,+)\), natural units only where a model explicitly declares them.

## 1. Purpose

The repository contains independently developed mechanics: discrete geometry, quaternionic and Clifford representations, lattice gauge transport, Dirac and Einstein–Cartan operators, teleparallel geometry, gauged scalar Q-balls, QLG/Dirac–Kähler operators, QAT state models, and effective continuum-mechanics targets.

This document and the module src/rgs_physics_architecture.py provide a central compatibility register. They:

1. cross-check conventions that must agree before modules can be coupled;
2. classify each physical sector by what its implementation currently supports;
3. list mathematical gates that must be passed before claiming a common theory.

They do **not** assemble an action by merely summing unrelated models.

## 2. Core convention audit

The audit compares two independently maintained Dirac representations:

- src.clifford_qlg.gamma_matrices;
- src.rgs_ec_algebraic_torsion.dirac_gamma_matrices_minus_plus.

It checks the metric, all four gamma matrices, gamma5, Clifford anticommutators, and the canonical Q6 chain identity B1 B2 = 0.

The expected finite Q6 structure is 64 vertices, 192 edges, 240 square faces, with B1 in Z^(64x192) and B2 in Z^(192x240). This combinatorial identity is not a continuum or gauge-theory result.

Run: python -m pytest tests/test_rgs_physics_architecture.py -q

Or run the audit from Python with audit_physics_architecture().

## 3. Sector registry and exact distinctions

| Sector | Current repository evidence | Not yet established |
|---|---|---|
| RGS geometry | finite coordinate, tangent and coframe audits | refining family and continuum spacetime |
| Q6 / PRISM | finite incidence and graph constructions | common geometric realization and refinement maps |
| Quaternions / 2I / Clifford | finite representation and algebra audits | full Lorentz spin connection and global bundle |
| Dirac–Einstein–Cartan | action-derived local kinetic/contact Euler operators; algebraic axial torsion | solved D_a psi, full tetrad/connection equations and total conservation |
| TEGR | separate teleparallel construction | matching to EC through an explicitly derived connection map |
| U(1) / Maxwell / Q-ball | specified two-complex-scalar gauged candidate and radial equations | coupled gravity, stability and full matter backreaction |
| SU(2) lattice gauge | oriented links, path transport and Wilson holonomies | derivation from the same connection used by the spinor |
| SU(3) / QCD | selected flavour/representation algebra | dynamical SU(3) links, chiral fermions and continuum QCD |
| QLG / Dirac–Kähler | finite spatial, Wilson/overlap and cochain operators | curved 4D continuum convergence |
| MHD / fluids / superfluids / elasticity | effective target theories and separate component audits | constitutive closure and common-action derivation |
| QAT / thermodynamics | finite-state and entropy calculations | physical time, partition measure and gravitational backreaction |

### Non-negotiable distinctions

- Minimal Einstein–Cartan torsion is algebraic in the current construction. It is not a propagating torsion field.
- TEGR torsion belongs to a teleparallel connection and is not the Einstein–Cartan spin-sourced torsion by notation.
- The polar-frame SU(2) transporter in src/discrete_dirac.py is a candidate spatial frame transporter. It is not automatically the independent internal SU(2) gauge connection in src/holonomies.py.
- The Q-ball diagonal U(1) charge is not automatically the Dirac axial current, electric current of every matter species, or an MHD current.
- The SU(3) flavour representation module is not a lattice QCD implementation.
- Spinor Dirac, Wilson/overlap and cochain Dirac–Kähler operators are distinct constructions with different Hilbert spaces and convergence obligations.

## 4. Common action: target, not assumption

The intended continuum target may be written schematically as

\[
S= S_{\rm grav}[e,\omega]
+S_{\rm Dirac}[e,\omega,A,\psi]
+S_{\rm YM}[A]
+S_{\rm Higgs}[e,A,H]
+S_{\rm PSU}[e,\Phi_+,\Phi_-,A]
+S_{\rm matter}[e,\chi].
\]

This is a list of sectors, not yet a valid single action. Before it can be called the RGS unified action, specify:

1. exact field content and representations under SU(3)_C x SU(2)_L x U(1)_Y;
2. whether omega is independent, constrained, or inertial, and which gravitational formulation is used;
3. the common spacetime/discrete measure and inner products;
4. all gauge-covariant derivatives and their signs;
5. boundary terms and variational variables;
6. couplings and dimensions;
7. symmetry transformations and Noether identities.

Do not include both Einstein–Hilbert/Einstein–Cartan and TEGR terms as independent gravitational sectors unless a non-equivalent extension is deliberately intended. In equivalent formulations, adding both would double-count gravitational dynamics.

## 5. Required coupling gates

| Gate | Required derivation | Status |
|---|---|---|
| RGS coframe → spin connection | construct omega_mu and its local Lorentz transformation | **OPEN** |
| spin connection → Dirac derivative | implement D_mu psi = partial_mu psi + Omega_mu psi + A_mu psi with explicit representation | **OPEN** |
| EC ↔ TEGR | derive connection map and boundary identity | **OPEN** |
| U(1), SU(2), SU(3) | independent generators, links, couplings and matter representations | **OPEN** |
| matter → Hilbert tensor | vary the complete action before imposing equations or fluid closure | **OPEN** |
| Noether/Bianchi | prove total gauge, Lorentz and diffeomorphism identities | **OPEN** |
| discrete → continuum | refining complexes, common embeddings, Mosco liminf/recovery and operator convergence | **OPEN** |
| quantum consistency | doubling, chirality, anomalies, positivity and measure | **OPEN** |

## 6. Immediate implementation order

1. Construct the RGS spin connection from the actual coframe and a declared connection choice.
2. Implement a single covariant-derivative interface with separate spin and internal gauge representations.
3. Feed that derivative into the merged EC–Dirac kinetic/contact Euler operator.
4. Derive tetrad and connection equations from the same action.
5. Add Maxwell and scalar Q-ball stress/current to the same variational system.
6. Add non-Abelian gauge groups only after the U(1)/SU(2) interface and representation contracts are stable.
7. Derive fluid, MHD and superfluid limits from matter equations and stress balance.
8. Quantize and test anomalies only after the classical action and field content are fixed.

## 7. Evidence policy

- **PROVED:** derivation under stated assumptions.
- **ACTION_DERIVED:** Euler equation or identity obtained from the declared action.
- **NUMERICALLY_VERIFIED:** reproducible finite computation with inputs and tolerances.
- **CANDIDATE:** specified construction with unresolved mathematical gates.
- **HYPOTHESIS:** physical interpretation not yet derived.
- **OPEN:** missing proof, equation, or independent validation.

A passing convention audit establishes compatibility of selected algebraic interfaces only. It does not establish a common action, solved field configuration, conservation, continuum limit, or agreement with experiment.
