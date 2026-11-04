# Exact finite-cutoff normalization of the Q6-D5 Bergman pullback

## Distribution of the type-IV generic norm

Let D5 be the Cartan domain of type IV in complex dimension five, with
normalized Euclidean volume dnu and generic norm
N(z)=1-2 z*z+|z^T z|^2. The Hua identity is

    integral_D5 N(z)^s dnu(z)
      = 60/[(s+1)(s+2)(s+3)(s+4)(s+5/2)], Re(s)>-1.

Writing t=N(z), the pushforward of normalized Euclidean volume under N has the exact probability density (obtained from the Hua Mellin transform and its Beta-product factorization):

    f_N(t) = (20/3)(1-sqrt(t))^4(1+4sqrt(t)+t).

Equivalently,

    f_N(t) = 20/3 - 60t - 60t^2 + (320/3)t^(3/2) + (20/3)t^3.

It is nonnegative, integrates to one, and its Mellin transform reproduces
the Hua polynomial exactly. In particular f_N(0)=20/3, the residue at the
sharp Hua convergence boundary.

## Finite cutoff

For 0<epsilon<1, restrict the Bergman density to N(z)>=epsilon. Its
normalizer relative to normalized Euclidean volume is exactly

    Z(epsilon) = integral_epsilon^1 t^-5 f_N(t) dt

    = (5/3)(epsilon^-4-1) - 20(epsilon^-3-1)
      - 30(epsilon^-2-1) + (128/3)(epsilon^-5/2-1)
      + (20/3)(epsilon^-1-1).

Therefore Z(epsilon) is finite and positive for every fixed positive
cutoff, but

    Z(epsilon) ~ (5/3) epsilon^-4, as epsilon tends to zero.

The normalized cutoff distribution of the generic norm is

    p_epsilon(t) = 1_[epsilon,1](t) t^-5 f_N(t) / Z(epsilon).

This is a well-defined probability density for each fixed cutoff. It does
not have a finite uncut Bergman-volume limit: its normalizer diverges as
epsilon^-4, and the probability mass concentrates toward N=0.

## Connection to the Q6 five-mode chart

The exact local pullback law from the preceding result is

    Phi^*(dV_B) = N(Phi(q,p))^-5 dmu_L.

The present calculation specifies precisely what a generic-norm regulator
does to this density. It does not identify epsilon with the Q6 lattice
spacing, transfer-matrix coupling, PSU scale, or a physical ultraviolet or
infrared cutoff. Such an identification must follow from an independent
derivation from the compact Wilson/transfer sector.

## Status

- Exact analytic result: factorized generic-norm density, cutoff normalizer,
  normalized cutoff distribution, and 5/3 leading divergence coefficient.
- Numerically cross-checked: cutoff integral against independent
  Gauss-Legendre quadrature and the small-cutoff asymptotic.
- Still open: derive a regulator from Q6 Wilson dynamics; prove a global
  compact Q6 gauge-quotient map into D5; compare normalized finite-cutoff
  measures with the finite-beta_t Bessel transfer operator; derive any
  physical scale or fine-structure constant from this construction.

## Reproducibility

    pytest -q tests/test_rgs_q6_d5_hua_cutoff.py

Executable API:
- d5_hua_generic_norm_density(t)
- d5_hua_cutoff_normalizer(epsilon)
- d5_hua_cutoff_probability_density(t, epsilon)
- d5_hua_cutoff_report()


## Literature and logical dependency

The Hua integral and Hua polynomial for irreducible bounded symmetric domains are given by Roos and Yin, arXiv:math/0301322, Proposition 2.1 and Theorem 2.3: https://arxiv.org/abs/math/0301322. For type IV, the generic norm and Bergman kernel form are also described by De Fabritiis (1990): https://www.bdim.eu/item?id=RLIN_1990_9_1_2_131_0. The moment uniqueness step is the Hausdorff moment theorem; see Diaconis and Freedman (2004), https://doi.org/10.1007/s00209-003-0633-9.

Terminology is important: f_N is the radial pushforward density of normalized Euclidean volume. The Bergman radial density is t^(-5) f_N(t), before cutoff normalization. The literature supplies the Hua integral and moment theorem; the explicit factorization, density and cutoff primitive are derived in this repository.
