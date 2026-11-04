# Exact Q6 integer-fiber quotient onto the D5 five-mode sector

## 1. Rational spectral map

Let C be the 129 x 192 integral fundamental-cycle matrix and
H = D1^T D1 the integral edge Hodge matrix. The selected D5 tangent
sector is the five-dimensional eigenspace ker_Q(H - 12 I). Choose
any rational basis V (192 x 5) of that eigenspace. The exact mode map is

    Pi_5 : Z^129 -> Q^5,   c |-> c C V.

The choice of rational basis changes the target coordinates by an
invertible rational 5 x 5 matrix but does not change the kernel. Exact
linear algebra gives rank(Pi_5)=5 and

    K_Z = Z^129 intersect ker_Q(Pi_5),   rank(K_Z)=124.

Thus Z^129/K_Z is torsion-free and isomorphic to the rank-five image
lattice Pi_5(Z^129) inside Q^5. This is an exact arithmetic quotient,
not merely a real rank count.

## 2. Why the pointwise transfer weight cannot be quotiented

For n=c C, the normalized temporal Wilson weight is

    w(c) = product_e I_|n_e|(beta_t) / I_0(beta_t).

It is not invariant under c -> c+k for k in K_Z. The audit constructs
an explicit nonzero integer k in the exact kernel and compares w(0)=1
with w(k)<1. Therefore the five-mode quotient alone does not inherit
the original pointwise Wilson weight.

## 3. The correct fiber-summed weight

The canonical positive marginal on a quotient class is

    W_eff([c]) = sum_{k in K_Z} w(c+k).

It is independent of the representative c because translation by an
element of K_Z permutes the summation set. It is strictly positive.

For fixed beta_t > 0,

    I_m(beta_t)/I_0(beta_t) <= (beta_t/2)^m/m!.

The integer cycle map C is injective, so ||k C||_1 grows at least
linearly with ||k|| on the finite-dimensional kernel lattice. The
factorial Bessel bound then gives absolute convergence of the sum over
K_Z. This proves existence of the fiber-summed quotient weight.

The implementation includes a finite-subset evaluator and checks its
translation consistency. It does not mislabel a finite subset as the
complete 124-dimensional fiber sum.

## 4. What remains before the Hua cutoff

This quotient is the electric five-dimensional lattice sector only.
It does not yet construct the spatial Wilson potential on the quotient,
the canonical coordinate q sector, or the full ten-real-dimensional
D5 phase-space measure. In particular it does not determine the Hua
generic-norm regulator epsilon. That requires a compatible q-sector and
a derived pushforward measure, not an identification by notation.

Run:

    python scripts/run_q6_d5_integer_fiber_quotient.py

Tests:

    pytest -q tests/test_rgs_q6_d5_integer_fiber_quotient.py

Exact results: rational top eigenspace, rank-five image, rank-124
integer kernel, non-invariance of the pointwise transfer weight, and
absolute convergence of the complete fiber-summed weight.

Open: controlled evaluation of the complete fiber sum, spatial sector,
full D5 phase-space measure, global map, Hua cutoff identification, and
any physical 3+1D reduction.
