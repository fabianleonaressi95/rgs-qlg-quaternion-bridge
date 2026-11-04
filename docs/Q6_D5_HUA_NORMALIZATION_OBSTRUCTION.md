# D5 Hua integral: exact normalization threshold and Bergman-volume obstruction

## Theorem

For the Cartan domain of type IV_5, the rank is r=2, the characteristic
multiplicities are a=3 and b=0, the genus is g=5, and the complex dimension
is 5. The Hua polynomial is

    chi_5(s) = (s+1)_4 (s+5/2)
             = (s+1)(s+2)(s+3)(s+4)(s+5/2),

with chi_5(0)=60.

Let dnu be the normalized invariant Euclidean volume measure on D5 and
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
