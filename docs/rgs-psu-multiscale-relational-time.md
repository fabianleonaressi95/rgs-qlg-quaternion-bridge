# RGS–PSU–Q6: multiscale geometry, relational clock and thermodynamic arrows

**Status:** reproducible finite-model audit and research programme. It is not a unified theory, a derivation of quantum gravity, or a proof that physical time is unreal.

## Scientific inputs

1. **Critical gravitational collapse.** C. Ecker, F. Ecker and D. Grumiller, “Analytic Discrete Self-Similar Solutions of Einstein-Klein-Gordon at Large D,” *Physical Review Letters* **136**, 191401 (2026), DOI [10.1103/qgl5-5l3t](https://doi.org/10.1103/qgl5-5l3t). The paper constructs analytic discretely self-similar solutions in a large-spacetime-dimension limit for Einstein gravity coupled to a massless scalar. This is not a demonstration of a four-dimensional microscopic black-hole population or a derivation of microscopic fluctuating time.
2. **Superfluid analogue gravity.** P. Švančara et al., “Rotating curved spacetime signatures from a giant quantum vortex,” *Nature* **628**, 66–70 (2024), DOI [10.1038/s41586-024-07176-8](https://doi.org/10.1038/s41586-024-07176-8). The experiment stabilizes a giant quantum vortex in superfluid helium-4 and measures wave–vortex effects described by an effective acoustic metric. It is an analogue curved spacetime, not an Einstein gravitational field.
3. **Fluctuation thermodynamics.** The Markov module is a controlled stochastic-thermodynamics example. It does not identify thermodynamic fluctuations with gravitational fluctuations.

## Implemented audits

### A. Exact Q6 combinatorics

Vertices are six-bit strings, with one edge for each single-bit flip: |V|=64, |E|=192, degree 6, graph cycle rank E−V+1=129, and combinatorial Laplacian eigenvalues 2j with multiplicity C(6,j). The adjacency matrix is constructed canonically from bit flips and hashed in the JSON report.

### B. Finite relational-clock example

Split Q6=Q2 □ Q4. Set H_C=A(Q2) and H_S=−A(Q4). Their spectra admit paired energies and a state satisfying

    (H_C ⊗ I + I ⊗ H_S)|Psi> = 0.

Conditioning on phase states of the clock factor produces distinguishable conditional states. This is a finite Page–Wootters-style example. The clock has only three distinct energies (−2,0,+2) and period π in the chosen units; it is not a continuous ideal clock and theta is not calibrated to proper time. This demonstrates relational correlations, not emergence of physical spacetime.

### C. Reversible and driven Markov controls

A symmetric-rate Q6 chain with uniform stationary distribution has detailed balance and zero stationary entropy-production rate. A second chain adds a nonzero affinity around one square face while preserving the uniform stationary distribution. Its stationary entropy-production rate is positive. Seeded Gillespie trajectories are recorded as diagnostics. The affinity is an explicit model input, not derived from PSU geometry.

### D. Acoustic metric control

For a 2+1 dimensional effective acoustic metric,

    g = [[−c_s² + v², −v_x, −v_y], [−v_x,1,0], [−v_y,0,1]],

the audit verifies det(g)=−c_s² and Lorentzian signature. A pure circulating flow (v_r=0) has an ergoregion but no radial acoustic event horizon; a draining component is required for a horizon. This avoids conflating an ergosurface with an event horizon.

### E. Discrete-scale toy identity

A log-periodic test function is checked under one discrete rescaling. This is only a unit test of the mathematical definition of discrete scale invariance; it is **not** a numerical reproduction of the Ecker–Ecker–Grumiller solution. The physical echoing period must be extracted from the paper's equations or supplemental material and tested in a faithful Einstein–Klein–Gordon implementation.

### F. Dimensionless self-gravity diagnostic

The report computes alpha_g=Gm²/(hbar c) and compactness at the reduced Compton length, where C=2 alpha_g. This illustrates why elementary-particle self-gravity is tiny while the Planck-mass coupling is order unity. It does **not** imply that gravity causes all microscopic thermodynamic or temporal fluctuations: those also occur in systems where gravity is negligible.

## Reproduction

~~~bash
python src/rgs_multiscale_relational_audit.py --output reports/rgs_multiscale_relational_audit.json --seed 20261104
pytest -q tests/test_rgs_multiscale_relational_audit.py
~~~

The JSON records seed, tolerances, Q6 adjacency SHA-256, section residuals and pass/fail statuses.

## What remains open

1. Replace the chosen Q2 clock split with a physically motivated relational-clock observable and test robustness under alternative subsystem splits.
2. Derive transition rates and cycle affinities from an RGS–PSU action or experimental data, rather than prescribing them.
3. Derive a continuum vierbein and metric from the discrete model, retaining existing RGS tangent-space conventions and auditing torsion separately.
4. Implement the actual large-D Einstein–Klein–Gordon solution and compare finite-D numerical solutions; do not substitute the log-periodic toy.
5. Reproduce the superfluid acoustic metric using measured flow profiles and dispersion, preserving the distinction between analogue and real gravity.
6. Establish whether any scale-dependent coupling between gravity, clock correlations and entropy follows from the equations. It is not assumed.

## Interpretation rule

**Proved by this code:** finite Q6 spectral identities, the selected global constraint, Markov detailed-balance controls, and algebraic acoustic-metric identities.

**Not proved:** emergent physical time, a new theory of gravity, microscopic black-hole production in four dimensions, or a gravitational origin of thermodynamic fluctuations.
