# Weyl phase closure, icosahedral C6 and Q6: algebraic bridge audit

## Scope

This audit connects (1) massless two-component Weyl helicity, (2) an oriented
antipodal-incidence matrix built from the 12 vertices of an icosahedron, and
(3) the 64-vertex, 192-edge hypercube Q6. It is a reproducible algebraic
bridge, **not** yet a teleparallel fermion model and not a derivation of alpha.

## 1. Weyl helicity

For nonzero spatial momentum, define
\[
\Sigma(\hat p)=\boldsymbol\sigma\cdot\hat p,\qquad
\Pi_\pm=\frac12(I_2\pm\Sigma).
\]
The tests verify idempotence, orthogonality, completeness and the Weyl equations
under the convention
\[
(EI_2-\boldsymbol\sigma\cdot\mathbf p)u_R=0,\qquad
(EI_2+\boldsymbol\sigma\cdot\mathbf p)u_L=0,\quad E=|\mathbf p|>0.
\]
Thus u_R has helicity +1/2 and u_L helicity -1/2 in this convention.
This identification is for massless positive-energy modes; it does not hold
unchanged for massive fermions.

## 2. Icosahedral C6 construction

The source reproduces the 12 normalized icosahedral vertices, nearest-neighbour
adjacency A, six oriented antipodal axes and incidence matrix B, then sets
\[
C_6=\frac12 B^T A B.
\]
For this ordering, the audit finds C6 symmetric, C6 squared equal to 5 I6,
and diag(C6)=0. Its eigenvalues are therefore +sqrt(5) and -sqrt(5), each
with multiplicity three. C6 is a signed off-diagonal coupling among six axes;
its square identity alone does not make sqrt(5) a fermion mass or a
propagator normalization.

## 3. Critical defect in the supplied Q6 toy script

The supplied script defines every Q6 edge weight by
\[
w_{(u,v,d)}=|C_{6,dd}|+1.
\]
But C6[d,d]=0 for all six axes. Consequently all 192 edge weights are exactly
one, their standard deviation is zero, and the Q6 operator contains **no C6
anisotropy at all**. The C6 off-diagonal entries are never used.

A mapping from the six icosahedral axes to the six Q6 bit directions must be
specified and tested for equivariance before off-diagonal C6 entries can
define a geometrically meaningful Q6 coupling. An arbitrary ordering is not
such a proof.

The random SU(2) link matrices in the toy script are internal gauge links.
They are not a spacetime Lorentz/Weitzenbock connection and do not define
torsion. Its Hermitian 128x128 matrix is a massive covariant graph
adjacency/Laplacian-type operator, not a two-component Weyl Dirac operator.
The inverse is a resolvent of that finite matrix; a small linear-solver
residual verifies only the inversion.

## 4. Correct integration path in the RGS–PSU–QLG repository

The repository already distinguishes:
- src/discrete_dirac: four-component spatial Dirac candidates and Wilson terms;
- src/holonomies: internal SU(2) link transport and Wilson loops;
- src/rgs_spin_connection: local torsion-free Levi-Civita spin connection;
- src/dirac_kahler: a 304-cochain finite Q6 Dirac–Kahler operator;
- src/q6_complex: canonical Q6 incidence and chain identity.

The next physical step is not to relabel random SU(2) links as torsion. It is
to supply a discrete coframe on Q6, define an independent Lorentz/teleparallel
connection and its edge transport, construct the Clifford-valued covariant
Dirac operator, and audit its gauge covariance, hermiticity/gamma-hermiticity,
chirality, spectral flow and continuum/refinement behaviour. The internal
electromagnetic U(1) connection must remain a separate bundle connection.

## 5. The Wyler / alpha bridge

Wyler's number is a geometric candidate, not a consequence of the Weyl
projectors or of C6 squared = 5 I6. A genuine bridge requires an explicit map
from the relevant bounded symmetric domain and its invariant measure to a
normalized phase-space/holonomy functional on the chosen compact geometry,
followed by a derivation of the coefficient of F_munu F^munu in the same
action that generates the fermion operator. Only after continuum matching and
renormalization to the low-energy Thomson coupling could that coefficient be
compared with alpha.

No such map or action derivation is established by this audit. Status:
Weyl projector algebra **PROVED/NUMERICALLY VERIFIED**; displayed C6 identity
**NUMERICALLY VERIFIED**; proposed C6-to-Q6 anisotropic weights **FALSIFIED
for the supplied diagonal-weight rule**; teleparallel fermion coupling and
alpha derivation **OPEN**.
