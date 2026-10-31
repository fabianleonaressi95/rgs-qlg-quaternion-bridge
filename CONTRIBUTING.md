# Contributing to RGS–QLG Quaternion Bridge

Thank you for your interest in contributing to this project.

The purpose of this repository is to develop a rigorous, reproducible mathematical and computational framework for unit quaternions, SU(2)/SO(3) rotations and lattice holonomies.

We welcome mathematical corrections, independent implementations, numerical tests, documentation improvements and constructive scientific criticism.

## 1. Contribution principles

All contributions should follow these principles:

- Mathematical definitions must be explicit.
- Quaternion conventions must be documented.
- Matrix multiplication and composition order must be stated.
- Numerical results must be reproducible.
- Tests must include tolerances and explain their numerical meaning.
- Mathematical identities must be distinguished from numerical observations.
- Physical interpretations must be identified as hypotheses unless independently established.

## 2. Development environment

The initial implementation uses Python 3.

Recommended dependencies:

- NumPy
- SciPy
- pytest

Install the dependencies with:

```bash
python -m pip install numpy scipy pytest
```

Run the test suite with:

```bash
python -m pytest -v
```

## 3. Code conventions

All quaternion implementations must explicitly document:

1. Component ordering.
2. Scalar-first or scalar-last representation.
3. Hamilton or alternative multiplication convention.
4. Active or passive rotation convention.
5. Coordinate-system handedness.
6. Composition order.

The initial project convention is:

\[
q=(w,x,y,z)
\]

with scalar-first storage and Hamilton multiplication.

Any alternative convention must be explicitly documented and tested.

## 4. Numerical validation

Contributions involving rotations or holonomies should include tests for:

- Quaternion normalization.
- Quaternion conjugation and inversion.
- Associativity of quaternion multiplication.
- SU(2) unitarity.
- SU(2) determinant.
- SO(3) orthogonality.
- SO(3) determinant.
- SU(2) to SO(3) consistency.
- Composition and inverse transport.
- Gauge covariance where applicable.

Numerical tolerances must be justified and should not conceal systematic errors.

## 5. Scientific status labels

Use the following labels in documentation:

- **PROVED**: established by a mathematical argument.
- **NUMERICALLY VERIFIED**: supported by reproducible computational tests.
- **CANDIDATE**: a proposed construction requiring further validation.
- **OPEN**: an unresolved mathematical or computational question.
- **PHYSICAL HYPOTHESIS**: a proposed physical interpretation not yet established.

A passing numerical test is not, by itself, a mathematical proof.

## 6. Pull requests

Please submit changes through pull requests whenever possible.

A pull request should contain:

- A concise description of the change.
- The mathematical or computational motivation.
- Relevant equations and conventions.
- Tests covering the new functionality.
- Expected numerical tolerances.
- Any known limitations.

Keep pull requests focused on one mathematical or computational topic.

## 7. Discussion and collaboration

Before introducing a substantial new mathematical construction, open a discussion describing:

- The problem being addressed.
- The proposed definitions.
- The assumptions being made.
- The expected mathematical properties.
- The tests required for independent verification.

Alternative formulations and independent implementations are particularly welcome.

## 8. Scope

This repository focuses initially on quaternion algebra, SU(2)/SO(3) representations and lattice holonomies.

Connections to the broader RGS–PSU and QLG research programme are welcome, provided that established mathematical results are clearly separated from conjectures and physical hypotheses.

Thank you for helping make this project mathematically precise, transparent and reproducible.
