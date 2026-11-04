# Q6 to relative sphere: A5-equivariant finite direction map

## Construction

The six unoriented icosahedral axes are the columns of a tight frame
A in R^(3x6), with A A^T = 2 I_3.

A Q6 vertex v=(b_1,...,b_6) in {0,1}^6 determines signs
s_j=(-1)^(b_j). Define x(v)=A s(v) and n(v)=x(v)/||x(v)|| in S^2.

The implementation uses the same signed icosahedral-axis representation
already acting on Q6. For each proper icosahedral rotation R, its signed
permutation S_R obeys A S_R=R A. Since the bit action is
s(gv)=S_R s(v), it follows that x(gv)=R x(v) and n(gv)=R n(v).

This is an explicit A5-equivariant map from the 64 Q6 states to a finite
subset of the relative-direction sphere.

## Finite image and angular quadrature

The computed image has 32 distinct unit directions, each with multiplicity
two among the 64 Q6 vertices. The 32-point image is antipodal (16 antipodal
pairs). All 64 signed sums are nonzero and occur on four radial shells
before normalization.

With equal weights on the 32 distinct directions, the moment tensors through
degree five agree with the uniform sphere moments:

    (1/32) sum_a n_a = 0
    (1/32) sum_a n_a n_a^T = I_3/3

    (1/32) sum_a n_a tensor4 =
      (delta_ij delta_kl + delta_ik delta_jl + delta_il delta_jk)/15

Odd moments through degree five vanish by antipodality. Thus the finite
direction set is numerically verified as a spherical 5-design. The sixth
moment differs from the uniform-sphere tensor; the finite set is not an
exact all-orders angular measure. Tests explicitly require this nonzero
sixth-order defect rather than hiding it.

## What this does and does not establish

NUMERICALLY VERIFIED: the signed-axis map, rank-three linear projection,
A5 equivariance, 32-point image, two-to-one fibres, antipodality, and
spherical moment identities through degree five.

CANDIDATE: use of this finite 5-design as an angular quadrature rule for
the relative S^2 coordinate in the two-body shell. A quadrature rule is
not a continuum parametrization; interpolation and convergence must still
be supplied.

OPEN: a continuous Q6-to-S^2 extension; a measure-preserving identification
with Smilga's on-shell phase space; a map from this angular set to the
five-dimensional coexact Hodge eigenspace; a global embedding into the D5
Lie ball; and any derivation of the fine-structure constant.

The 5-dimensional A5 coexact eigenspace in
rgs_q6_d5_phase_space.py and the 3-dimensional standard A5 direction
representation used here are different irreducible representations. This
module does not identify them. Their tensor-product/intertwiner structure
must be computed before a claim that the relative-direction modes are
contained in the D5 spectral sector.
