# Discrete frame alignment and Spin(4)/Clifford transport

## Implemented convention

Coframes are real (4\times4) matrices with internal index first and
Euclidean metric (g_v=e_v^T e_v). For an oriented edge (u\to v), define

[
A_{v\leftarrow u}=e_v e_u^{-1}=R_{v\leftarrow u}P_{v\leftarrow u},
]

where (R\in SO(4)) is the orientation-preserving orthogonal polar factor.
The module lifts (R) to (S\in\mathrm{Spin}(4)) using

[
\Omega=\log R,\qquad
S=\exp\!\left(\frac14\sum_{a,b}\Omega_{ab}\gamma_a\gamma_b\right).
]

The convention is checked numerically by

[
S\gamma_cS^\dagger=\sum_bR_{bc}\gamma_b,
\qquad
\{\gamma_a,\gamma_b\}=2\delta_{ab}I.
]

The four gamma matrices are explicit Hermitian (4\times4) Euclidean
representatives. The module returns all 192 Q6 edge transports and reports
metric mismatch and coframe-alignment residuals rather than suppressing them.

## Geometric status and limitations

This is a **polar frame-alignment prescription**, not automatically a
Levi-Civita connection. It is an exact orthogonal frame change when the two
endpoint metrics coincide and (e_v=R e_u). For generic endpoint metrics,
the polar factor is a chosen nearest-frame alignment; it does not encode
the missing metric variation or affine connection data.

The lift is Euclidean Spin(4), not Lorentzian Spin(1,3). No discrete curvature,
torsion tensor, Cartan structure equation, TEGR action, or Dirac kinetic
operator is claimed here. Those require an independently defined discrete
coframe difference, connection/holonomy convention, and a compatible
Clifford-valued Dirac operator.

## Tests

Tests check the Clifford anticommutators, nontrivial SO(4) spin intertwining,
inverse-edge transport, all 192 Q6 edges, orientation/non-degeneracy checks,
and explicit reporting of metric mismatch.
