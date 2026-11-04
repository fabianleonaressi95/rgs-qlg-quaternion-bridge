# Positive C6 constitutive coupling on Q6 directional differences

## Motivation and correction

The diagonal rule previously tested, w_(v,d)=1+|C6[d,d]|, is degenerate:
the constructed C6 has zero diagonal. This module does not replace it with
another scalar edge weight. It uses the full off-diagonal C6 matrix as a
constitutive tensor on the six directional differences at each Q6 vertex.

## Definition

For the canonical big-endian Q6 axis masks m_d=2^(5-d), define
\[
(Df)_{v d}=f_{v\oplus m_d}-f_v,\qquad
D:\mathbb R^{64}\to\mathbb R^{384}.
\]
The icosahedral signed axis coupling C6 is pulled back to the Q6 direction
frame by a declared permutation a in S6:
\[
(C_Q)_{de}=(C_6)_{a(d),a(e)}.
\]
This is a frame identification, not a canonical embedding and not a claim of
SO(3) or O_h equivariance. Relabelling the frame and its components together
is a coordinate covariance statement only.

Because C6 has eigenvalues +/-sqrt(5), it is indefinite and cannot by itself
be a positive kinetic metric. Define
\[
H=\mu I_6+\beta C_Q,\qquad
\mu=\sqrt5|\beta|+\epsilon,\quad \epsilon>0.
\]
Then H is strictly positive, with eigenvalues epsilon and
2 sqrt(5)|beta|+epsilon (up to ordering). The quadratic energy and operator are
\[
E[f]=\frac12\sum_{v=0}^{63}(Df)_v^T H(Df)_v
=\frac12 f^T Kf,\qquad
K=D^T(I_{64}\otimes H)D.
\]

## What this establishes

- The full off-diagonal C6 coupling enters a well-defined Q6 directional
  constitutive form; the former diagonal-weight defect is avoided.
- H is strictly positive by construction.
- K is symmetric positive semidefinite and annihilates constants. Since the
  Q6 graph is connected and H is positive definite, its kernel is exactly
  the constant mode.
- The six-channel frame identification is explicit and can be relabelled
  covariantly.

## What remains open

This construction is a discrete constitutive coupling, not yet a geometric
identification of icosahedral axes with physical Q6 bit directions. The
permutation a is not selected by a proved intertwiner. No teleparallel
connection, torsion, Weyl/Dirac operator, SU(2) gauge covariance, or
continuum/Mosco limit is derived here. The parameter beta is a model coupling,
not a measured physical constant. The strict shift mu is a stability
regularization, not a fitted mass.

## Evidence status

- Algebraic construction and positivity: **PROVED** from the definitions.
- Numerical tests of the finite Q6 operator: **NUMERICALLY VERIFIED** by CI.
- Physical icosahedral-to-Q6 frame selection: **OPEN**.
- Teleparallel fermion and fine-structure constant derivations: **OPEN**.
