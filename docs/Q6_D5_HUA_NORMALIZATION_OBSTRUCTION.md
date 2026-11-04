# D5 Hua integral: exact normalization threshold and Bergman-volume obstruction

## Theorem

For the Cartan domain of type IV_5, the rank is r=2, the characteristic
multiplicities are a=3 and b=0, the genus is g=5, and the complex dimension
is 5. The Hua polynomial is

    chi_5(s) = (s+1)_4 (s+5/2)
             = (s+1)(s+2)(s+3)(s+4)(s+5/2),

with chi_5(0)=60.

Let dnu be the normalized Euclidean volume measure on D5 and
N(z)=1-2 z* z+|z^T z|^2 its generic norm. The Hua integral formula gives

    I(s) = integral_D5 N(z)^s dnu(z)
         = chi_5(0)/chi_5(s),       Re(s)>-1.

The convergence threshold is sharp. At s=-1+epsilon,

    I(-1+epsilon) ~ 20/(3 epsilon).

The D5 Bergman volume density relative to Euclidean volume is proportional
to N(z)^(-g)=N(z)^(-5). Since -5 is outside the Hua convergence half-line,
the total Bergman volume of D5 is infinite.

## Consequence for the Q6-D5 pullback

The previously proved local-chart identity

    Phi^*(dVol_B) = N(Phi(q,p))^(-5) dmu_L

is exact wherever the Maxwell chart maps into D5. This does not turn dVol_B
into a finite probability measure on the whole noncompact D5 domain.

In particular, it cannot be identified directly with a normalized finite
Haar measure on compact Wilson link-configuration space. Any such comparison
must specify an additional ingredient: a finite-domain cutoff, a quotient
with finite volume, or a compensating weight. That ingredient must be derived,
not fitted.

This is a genuine obstruction, not a numerical mismatch.

## Proof of the polynomial

For type IV_n, r=2, a=n-2, b=0. The Hua polynomial formula is

    chi_n(s) = (s+1)_(n-1) (s+n/2).

Setting n=5 yields the displayed degree-five polynomial. The standard Hua
integral identity gives I(s)/I(0)=chi(0)/chi(s) for Re(s)>-1. Since
chi_5(s) has a simple zero at s=-1 and

    chi_5(-1+epsilon) = 9 epsilon + O(epsilon^2),

the normalized integral has residue 60/9=20/3. For s=-5 the integral lies
outside the convergence range and diverges.

## Reproducibility

    pytest -q tests/test_rgs_q6_d5_hua_normalization.py

The executable report is
d5_hua_wilson_normalization_report() in
src/rgs_q6_d5_hua_normalization.py.

## Status separation

- Proved from the Hua integral theorem: the convergence threshold, polynomial,
  residue and infinite total Bergman volume.
- Numerically checked: executable evaluation of the polynomial and boundary
  asymptotic.
- Open: deriving a physically justified cutoff/quotient/weight from the
  compact Wilson transfer operator, and any resulting Wyler/alpha formula.


## Derived cutoff asymptotic (no fitted regulator)

The Hua Mellin transform admits an exact probabilistic factorization. If
independent random variables have laws

    X_1 ~ Beta(1,1), X_2 ~ Beta(2,1), X_3 ~ Beta(3,1),
    X_4 ~ Beta(4,1), X_5 ~ Beta(5/2,1),

then their product X=X_1 X_2 X_3 X_4 X_5 has Mellin transform

    E[X^s] = product_j a_j/(a_j+s)
           = 60/[(s+1)(s+2)(s+3)(s+4)(s+5/2)]
           = I(s)/I(0).

Because X and the generic norm N under normalized Euclidean volume both
take values in [0,1], equality of all nonnegative integer moments implies
equality of their distributions (Hausdorff moment determinacy).

Put Y=X_2 X_3 X_4 X_5, so X=X_1 Y and X_1 is uniform on [0,1].
Conditioning on Y gives the density

    f_X(t) = E[Y^-1 1_{Y>t}],  0<t<1.

Monotone convergence and independence give

    lim_{t downarrow 0} f_X(t)
      = E[Y^-1]
      = product_{a in {2,3,4,5/2}} a/(a-1)
      = 20/3.

Consequently, for the norm cutoff 0<delta<1, define the dimensionless
truncated Bergman integral relative to normalized Euclidean volume by

    Z_5(delta) = E[N^-5 1_{N>=delta}]
               = integral_delta^1 t^-5 f_X(t) dt.

The boundary-density limit implies the sharp leading divergence

    lim_{delta downarrow 0} delta^4 Z_5(delta) = (20/3)/4 = 5/3,

or

    Z_5(delta) ~ (5/3) delta^-4.

This coefficient is derived from the Hua law; it is not fitted. The code
exposes the leading term as an asymptotic approximation only, not as an
exact finite-delta integral.

### Interpretation and limits

- The cutoff N>=delta gives a mathematically explicit family of finite
  integrals and its divergence rate.
- It does not select a physical value of delta.
- It does not show that the Q6 transfer operator induces this norm cutoff.
- It does not identify the compact Wilson Haar measure with the D5
  Bergman measure; those remain different measures on different spaces.
- A physical regulator must still be derived from a specified discrete
  transfer kernel or pushforward, not chosen to match alpha.
