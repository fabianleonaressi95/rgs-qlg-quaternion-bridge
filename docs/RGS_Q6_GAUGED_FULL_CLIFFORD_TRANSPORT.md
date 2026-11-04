# Gauge-covariant full cubical Q6 Dirac-Kaehler transport

## Construction

The full cubical 6-cube has cell counts (64, 192, 240, 160, 60, 12, 1),
for 729 cells in degrees 0 through 6.

Put an SU(2) fibre at the lower-vertex anchor of each cell. For a
(k+1)-cell c=(b; i_0<...<i_k), the lower face normal to i_p has anchor b,
while the upper face has anchor b+e_i_p. The link U_(b,b+e_i_p) transports
the upper-face fibre back to b. The covariant coboundary is

    (D_U omega)_c = sum_p (-1)^p [
        U_(b,b+e_i_p) omega_(upper face p) - omega_(lower face p)
    ].

For a vertex gauge transformation G_b in SU(2), every cell cochain
transforms with the matrix at its anchor. Thus D_U' = P_(k+1) D_U P_k^dagger.
The Hermitian operator is D_U + D_U^dagger. The implementation stores the
full 1458 by 1458 operator sparsely and evaluates all 240 square holonomies.

## Curvature and nilpotency

For identity links, D_U squared is exactly zero, recovering the ordinary
cubical cochain complex. For a curved connection, D_U squared is generally
nonzero: on each square it contains the difference between the two parallel
transports from the upper corner to the lower anchor. This is the expected
plaquette-curvature obstruction, not a failure of the ordinary incidence
complex.

## Relation to the exterior Clifford module

The six Hermitian generators gamma_i = epsilon_i + contraction_i act on
the 64-dimensional exterior algebra and satisfy the Euclidean Clifford
relations. The SU(2) link acts on a separate internal factor, so
U tensor I_64 commutes with I_2 tensor gamma_i. This is a consistent product
internal-gauge/Clifford module. It does not construct a Spin(6) connection:
a spin connection must act through a specified Spin(6) lift and must be tied
to a discrete frame/coframe geometry.

## Evidence boundary

Exact by construction: cell-anchor gauge representation, covariant incidence
formula, Hermiticity of D_U + D_U^dagger.

Numerically audited: local SU(2) covariance, 240 plaquette conjugations,
Wilson trace invariance, flat nilpotency, nonzero curvature obstruction,
Clifford anticommutation and commutation of the independent tensor factors.

Still open: deriving the connection from the RGS coframe/action, a Spin(6)
or Lorentzian Spin(1,3) lift, the continuum limit, and fermion doubling.
The operator is Dirac-Kaehler, not yet a spin-Dirac operator.

## Reproduction

Run:
    PYTHONPATH=. python -m pytest tests/test_rgs_q6_gauged_full_clifford.py -q
