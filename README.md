# rgs-qlg-quaternion-bridge
Quaternionic rotations, SU(2)/SO(3), lattice holonomies and reproducible computational validation
# RGS–QLG Quaternion Bridge

**An open mathematical and computational collaboration on unit quaternions, SU(2), SO(3), lattice rotations and holonomy structures.**

[![Status](https://img.shields.io/badge/status-research%20proposal-blue)]()
[![License](https://img.shields.io/badge/license-MIT-green)]()

## Purpose

This repository develops a reproducible computational framework for quaternionic representations of three-dimensional rotations and their application to discrete lattice connections and holonomies.

The project focuses on a precise mathematical question:

**How can quaternionic rotation conventions, SU(2) lifts, and lattice holonomies be implemented and independently validated without conflating quaternion equality with equality of the represented SO(3) rotations?**

The repository is an invitation to mathematical and computational collaboration. It does not assume that the physical interpretations of the broader RGS–PSU programme have been established.

## Mathematical starting point

A unit quaternion is written as

\[
q=(\mathbf v,w),\qquad
\|\mathbf v\|^2+w^2=1.
\]

Its quaternionic conjugate and inverse are

\[
q^*=(-\mathbf v,w),\qquad q^{-1}=q^*.
\]

The unit-quaternion group is isomorphic to \(SU(2)\), and the standard rotation map

\[
\rho:SU(2)\longrightarrow SO(3)
\]

is a double covering:

\[
\rho(q)=\rho(-q).
\]

Consequently,

\[
-q^*=(\mathbf v,-w)
\]

represents the same inverse rotation in \(SO(3)\) as \(q^*\), but it is not the quaternionic inverse of \(q\):

\[
q(-q^*)=-1,\qquad qq^*=1.
\]

This distinction is central to the conventions and tests developed here.

## Research objectives

1. Establish unambiguous quaternion conventions and composition order.
2. Validate the SU(2) to SO(3) map and its double-cover property.
3. Compare direct matrix rotations with quaternionic rotations.
4. Implement and test lattice parallel transport and plaquette holonomies.
5. Examine complex conjugation and orientation reversal in a two-sheet construction.
6. Publish reproducible numerical tests and document all limitations.

## Proposed computational programme

The initial implementation will include:

- quaternion multiplication, conjugation, inverse and normalization;
- quaternion-to-rotation-matrix conversion;
- SO(3) orthogonality and determinant tests;
- SU(2) representation and lift consistency;
- composition-order tests;
- closed-loop holonomy calculations;
- gauge covariance checks;
- numerical stability and edge-case tests.

## Connection to the RGS–PSU / QLG programme

The broader research programme contains a discrete geometric model with an SU(2) connection, 192 edge rotations and 240 plaquette holonomies.

These structures provide a test case for the mathematical framework developed here. Their physical interpretation remains a separate research question and will not be treated as a consequence of the quaternionic identities alone.

## Collaboration proposal

We welcome independent review, corrections, alternative implementations and discussion of quaternion conventions, rotation representations and numerical validation.

A particularly valuable contribution would be an independent review of the relationship between quaternionic inversion, the SU(2) double cover and the composition of lattice holonomies.

## Evidence policy

Results will be labelled explicitly as:

- **Mathematical identity** — established by proof.
- **Numerically verified** — supported by reproducible tests.
- **Candidate construction** — mathematically specified but not fully established.
- **Open question** — requiring further analysis.

No numerical experiment will be presented as a proof of a physical theorem.

## License

MIT License — subject to confirmation before publication.
