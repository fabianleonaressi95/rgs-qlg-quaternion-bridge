# Exact Q6 Maxwell–D5 Bergman/Liouville measure law

## Theorem

Let (q,p) be canonical coordinates on the five-mode Q6 Maxwell sector, with
H = 1/2 [ beta lambda q.q + beta^(-1) p.p ], beta,lambda>0.
Put Omega=sqrt(lambda), kappa=beta*Omega and define

    Phi(q,p) = z = (sqrt(kappa) q + i p/sqrt(kappa))/sqrt(20).

This is the flow-adapted canonical scale for this Maxwell Hamiltonian.
The older scale sqrt(beta*lambda) gives the same constant volume Jacobian,
but it does not produce the Hamiltonian common-phase flow unless beta=1.

For the type-IV domain D5, use generic norm

    N(z)=1-2 z* z + |z^T z|^2

and Bergman kernel K=C N^(-5). With Bergman Kähler form
omega_B=i partial bar-partial log K and volume dV_B=omega_B^5/5!,
and canonical Liouville measure dmu_L=product dq_j dp_j, the exact law is

    Phi^*(dV_B) = N(Phi(q,p))^(-5) dmu_L,

wherever Phi(q,p) belongs to D5. At the origin N=1, so the two volume
forms agree exactly. Hamilton's equations have frequency Omega=sqrt(lambda)
and the chart obeys Phi(q(t),p(t))=exp(-i Omega t)Phi(q(0),p(0)).

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
Q6-to-D5 completion must reproduce. The weighted measure is preserved by
this isotropic oscillator because N is constant along its common-phase flow;
for a general Hamiltonian the exact criterion is X_H(N)=0.

The PR #145 volume identity is independent of the positive reciprocal scale
used in the linear chart. The dynamical statement is not: it requires the
flow-adapted scale kappa=beta*sqrt(lambda). This distinction reconciles the
older chart convention with the Maxwell-normalized Hamiltonian in
RGS_Q6_D5_MAXWELL_HAMILTONIAN.md.

## Reproducibility

    pytest -q tests/test_rgs_q6_d5_liouville_measure.py

The audit function is q6_d5_bergman_liouville_report in
src/rgs_q6_d5_liouville_measure.py.


## Literature and scope

The type-IV domain, its generic norm and Bergman-kernel exponent are standard Cartan-domain results; see De Fabritiis (1990), https://www.bdim.eu/item?id=RLIN_1990_9_1_2_131_0, and Roos–Yin (arXiv:math/0301322), https://arxiv.org/abs/math/0301322. The preservation of canonical Liouville volume by Hamiltonian flow is the standard symplectic Liouville theorem; see the Hamiltonian-systems chapter in Cieliebak's symplectic geometry notes, https://ammann.app.uni-regensburg.de/lehre/2017w_symplectic/cieliebak_symplectic-geometry-a.pdf.

The additional statement that the weighted density N^(-5) is preserved requires N to be a first integral. In this implementation that is proved for the isotropic five-mode oscillator because the flow is a common phase rotation and N is phase invariant. It is not asserted for a general Hamiltonian.

The identity proved here is local to the explicit linear chart on the preimage of D5. It is not a pushforward theorem for the complete compact Wilson quotient, nor does it identify a physical cutoff.
