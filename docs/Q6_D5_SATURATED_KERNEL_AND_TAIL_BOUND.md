# Saturated Q6-D5 integer kernel and certified Wilson tail

This extension strengthens the exact quotient introduced in PR #153.

## 1. Saturated integer kernel

Write the rational mode map as M: Z^129 -> Q^5. Clear a common denominator
to obtain an integer matrix A in Z^(5x129), without changing its kernel.

The Smith decomposition supplies unimodular matrices U,V such that U A V = S.
If rank(A)=5, the last 124 columns of V form a basis of the full integer
kernel, not merely 124 independent integer vectors in a finite-index
sublattice. The implementation verifies the Smith identity and exact
annihilation, then applies row Hermite normal form to canonicalize the row lattice.

## 2. A rigorous global tail majorant

For n in Z^192, put s=||n||_1. The normalized Bessel bound gives

    w_beta(n) <= (beta/2)^s / product_e |n_e|!.

By the multinomial theorem, s!/product_e |n_e|! <= 192^s, hence

    w_beta(n) <= (96 beta)^s / s!.

The number of integer vectors in Z^192 with l1 norm s is at most
2^192 * binomial(s+191,191). Therefore the omitted mass of any subset,
including one Wilson fibre, outside the edge-flux l1 ball of radius S-1
is bounded by

    T_beta(S) = sum_{s=S}^infinity 2^192 binomial(s+191,191) (96 beta)^s/s!.

The ratio of consecutive summands is

    r_s = 96 beta (s+192)/(s+1)^2.

Once r_s<1, it decreases thereafter; the implementation sums the finite
pre-asymptotic segment and bounds the remaining tail by a geometric series.
The result is reported in logarithmic form to avoid overflow.

## 3. Scope

- Exact: saturated integer kernel basis via Smith normal form over ZZ.
- Exact: rank and annihilation checks.
- Rigorous bound: Wilson mass omitted beyond a specified edge-flux l1 cutoff.
- Not claimed: complete numerical evaluation of the 124-dimensional fibre sum.
  Enumerate all fibre points inside the edge-flux ball, then add the certified tail.
- Not claimed: derivation of the Hua generic-norm cutoff. The electric marginal
  still needs the spatial/canonical coordinate sector for a full D5 measure.

## Reproduction

    python scripts/run_q6_d5_saturated_kernel_bounds.py
    pytest -q tests/test_rgs_q6_d5_saturated_kernel_bounds.py
