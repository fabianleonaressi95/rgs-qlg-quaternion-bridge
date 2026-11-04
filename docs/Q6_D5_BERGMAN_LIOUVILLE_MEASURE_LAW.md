# Exact Q6 Maxwell–D5 Bergman/Liouville measure law

## Theorem

Let (q,p) be canonical coordinates on the five-mode Q6 Maxwell sector, with
H = 1/2 [ beta lambda q.q + beta^(-1) p.p ], beta,lambda>0.
Put omega=sqrt(beta lambda) and define

    Phi(q,p) = z = (sqrt(omega) q + i p/sqrt(omega))/sqrt(20).

For the type-IV domain D5, use generic norm

    N(z)=1-2 z* z + |z^T z|^2

and Bergman kernel K=C N^(-5). With Bergman Kähler form
omega_B=i partial bar-partial log K and volume dV_B=omega_B^5/5!,
and canonical Liouville measure dmu_L=product dq_j dp_j, the exact law is

    Phi^*(dV_B) = N(Phi(q,p))^(-5) dmu_L,

wherever Phi(q,p) belongs to D5. At the origin N=1, so the two volume
forms agree exactly.

## Proof

Differentiation gives

    partial_i N = -2 bar(z_i) + 2 z_i overline(z^T z)
    partial_i partial_bar_j N = -2 delta_ij + 4 z_i bar(z_j).

Since log K=log C-5 log N, the Bergman metric is

    g_i_barj = -5 [ (partial_i partial_bar_j N)/N
                   - (partial_i N)(partial_bar_j N)/N^2 ].

The type-IV determinant identity (genus five) is

    det(g)=10^5 N^(-5).

The real Jacobian of Phi from (q,p) to (Re z, Im z) is independent of
beta and lambda:

    |J_Phi|=(1/20)^5=1/(32*10^5).

In complex coordinates the Kähler volume is
dV_B=2^5 det(g) d^10(Re z,Im z). Hence

    Phi^*dV_B = 2^5 (10^5 N^(-5))/(32*10^5) dmu_L
              = N^(-5) dmu_L.

At z=0, g(0)=10 I_5 and direct pullback gives

    Phi^*omega_B at 0 = sum_j dq_j wedge dp_j = omega_can.

This proves the volume-density identity and the infinitesimal symplectic
normalization.

## Interpretation and limits

- Proved: the Bergman volume pulls back to the five-mode Liouville volume
  with the exact, non-adjustable density N^(-5).
- Proved at the tangent origin: the Bergman Kähler form equals the canonical
  symplectic form.
- Not proved: equality of the Kähler and canonical symplectic forms away
  from the origin; a global nonlinear embedding of the compact Q6 gauge
  quotient; equality with the finite-beta_t Bessel/compact Wilson measure;
  or a physical normalization/derivation of alpha.

The factor N^(-5) is the precise measure correction any proposed global
Q6-to-D5 completion must reproduce.

## Reproducibility

    pytest -q tests/test_rgs_q6_d5_liouville_measure.py

The audit function is q6_d5_bergman_liouville_report in
src/rgs_q6_d5_liouville_measure.py.
