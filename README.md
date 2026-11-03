# RGS–QLG Quaternion Bridge

**Research codebase for quaternionic geometry, discrete lattice structures, spectral operators and reproducible computational validation.**

This repository is part of the RGS–PSU–QLG research programme. It develops explicit mathematical constructions and computational tests involving unit quaternions, SU(2)/SO(3), discrete incidence geometry, lattice holonomies, symmetry representations and Dirac-type operators.

> **Evidence policy:** mathematical identities, numerical verification, candidate constructions and open physical hypotheses are reported separately. Passing a finite numerical test does not establish a continuum limit or a physical theory.

## Research scope

### Quaternionic and group geometry
- Hamilton quaternion conventions and composition order.
- SU(2) and its double covering of SO(3).
- Binary icosahedral group 2I, order 120, and selected representations.
- Representation intertwiners and symmetry-sector checks.

### Discrete geometry
- Q6: 64 vertices, 192 edges and 240 square faces.
- Integer incidence operators and chain-complex identities.
- Möbius-complex test constructions and local-frame conventions.
- Graph Laplacians, cochain operators and topology.

### Spectral and spinorial operators
- Clifford algebra conventions.
- Discrete Dirac-type and Dirac–Kähler constructions.
- Wilson-kernel and overlap-operator experiments.
- Spectral residuals, multiplicities and symmetry checks.

### Geometric field structures
- Lattice holonomies and gauge covariance.
- Discrete vierbein and connection constructions.
- Torsion and teleparallel research directions.

## Mathematical starting point: SU(2) to SO(3)

Unit quaternions form a group isomorphic to SU(2). The standard map

    rho: SU(2) -> SO(3)

is a double covering:

    rho(q) = rho(-q).

Thus q and -q represent the same spatial rotation, while remaining distinct elements of SU(2). This distinction is essential in spinorial calculations and in the treatment of the central element of 2I.

## Q6 reference structure

The six-dimensional hypercube graph has

    |V| = 64, |E| = 192, degree(v) = 6.

With the six-bit vertex convention and the oriented incidence matrices B1 and B2, the chain-complex condition is

    B1 B2 = 0.

Its graph Laplacian has exact eigenvalues

    lambda_j = 2j, multiplicity = binomial(6,j), j = 0,...,6.

These identities provide exact finite-dimensional reference tests. They do not, by themselves, establish a continuum geometric limit.

## Repository map

- src/ — mathematical and computational implementations.
- tests/ — unit and structural tests.
- docs/RGS_PSU_GAUGE_CLIFFORD_QLG_RESEARCH_PROPOSAL.md — research direction.
- docs/IMPLEMENTATION_ROADMAP.md — implementation roadmap.
- docs/STEM_INTEROPERABILITY_PROFILE.md — proposed interoperability profile.
- [Issue #3](https://github.com/fabianleonaressi95/rgs-qlg-quaternion-bridge/issues/3) — tracked implementation discussion.

The precise modules and callable functions used for a benchmark must always be identified by repository commit SHA and recorded in its provenance manifest.

## RGS–STEM interoperability

The [RGS–STEM Interoperability Profile](docs/STEM_INTEROPERABILITY_PROFILE.md) proposes a conservative sidecar-based exchange format for finite algebraic, geometric and spectral data.

It is a proposal for technical discussion, **not an approved Ontos extension or an endorsement by Stem Framework**. The current Ontos stream specification remains unchanged unless its maintainers agree to a formal extension.

## Evidence and reproducibility standards

Every reported result should state:

1. the mathematical object and conventions;
2. the implementation and exact repository revision;
3. the test method and expected invariant;
4. numerical residuals and tolerances, where relevant;
5. whether verification uses the same code, a separate module or an independent implementation;
6. limitations and unresolved assumptions.

Claims are classified as:

- **PROVED** — established by a mathematical proof.
- **NUMERICALLY VERIFIED** — supported by reproducible computation with stated tolerances.
- **CANDIDATE** — explicitly constructed but not fully established.
- **OPEN** — requiring further mathematical, numerical or empirical work.

In particular, finite Q6 or 2I computations do not automatically prove continuum quantum field theory, general relativity, Einstein–Cartan dynamics, antimatter phenomenology or cosmological claims.

## Collaboration

Independent mathematical review, alternative implementations, code review and reproducibility audits are welcome. The most useful contributions are precise, testable and accompanied by explicit conventions.

For the broader research overview, see the [research profile](docs/GITHUB_PROFILE_README.md).

## License

See the repository license file. If no license file is present, reuse permissions should not be inferred from this README.

## Master catalogue of laws and theorems

The consolidated [RGS–PSU–PRISM–QLG Master Catalogue](docs/MASTER_CATALOGUE_LAWS_THEOREMS.md) records mathematical theorems, exact finite identities, numerical results, model definitions, conjectures and open physical hypotheses with explicit evidence labels. It is versioned and deliberately does not present unproved physical interpretations as established laws.

The [Phase 6–10 stable-ID register](docs/MASTER_CATALOGUE_PHASE6_10_ADDENDUM.md) expands the catalogue with theorem statements, hypotheses, numerical residuals, source files and explicit limits for Q6 incidence, curved Clifford fields, discrete Dirac, Wilson, overlap/Ginsparg–Wilson, Dirac–Kähler, the flat-torus Mosco prototype and Möbius finite-spectrum thermodynamics.

Machine-readable index: [theorem registry JSON](results/catalogue/laws_theorems_registry.json). The next analytical gate is documented in the [curved RGS Mosco refinement program](docs/RGS_MOSCO_CURVED_REFINEMENT_PROGRAM.md); its status remains OPEN until a curved refinement family and liminf/recovery proof are supplied.

The [flat-torus G0–G3 refinement benchmark](docs/FLAT_TORUS_G0_G3_AUDIT.md) now provides a periodic Freudenthal mesh generator, incidence matrices, P1 stiffness, positive lumped 0-form Hodge star, refinement diagnostics and dedicated CI tests. It is a flat reference benchmark, not a curved RGS Mosco result.
