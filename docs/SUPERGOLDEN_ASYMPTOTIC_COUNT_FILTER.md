# Stage-1 supergolden asymptotic count filter

This stage uses only the number of objects in each channel at each refinement
level. It requires no element IDs, child incidence, geometry, or orientation.

## Input

CSV, one row per consecutive level:

    n,A,B,C
    0,1,0,0
    1,1,1,0
    2,2,1,1
    3,3,2,1

Counts must be nonnegative integers; level labels must be unique and contiguous.

## Exact recurrence

For
\[
M=\begin{pmatrix}1&0&1\\1&0&0\\0&1&0\end{pmatrix},
\]
Cayley-Hamilton gives \(M^3=M^2+I\). Hence the count vector and its total obey
\[
v_{n+3}=v_{n+2}+v_n,\qquad N_{n+3}=N_{n+2}+N_n.
\]
Both checks use exact integer arithmetic. The vector test is stronger than
checking totals alone.

## Perron frequencies

The positive right Perron vector is proportional to
\((1,s^{-1},s^{-2})\), where \(s^3-s^2-1=0\). Normalizing gives
\[
(A,B,C)\longrightarrow(s^{-2},s^{-3},s^{-4}),
\]
whose entries sum to one. The script computes s using Decimal Newton
iteration and reports observed fractions and absolute errors at every level.
The optional --frequency-tolerance argument checks the last available level;
it is a finite-data diagnostic, not a proof of convergence.

## Run

    python -m pytest tests/test_supergolden_level_count_filter.py -q
    python scripts/filter_supergolden_level_counts.py data/level_counts.csv --output count_audit.json
    python scripts/filter_supergolden_level_counts.py data/level_counts.csv --frequency-tolerance 0.01

Exit status is 0 when the exact vector and total recurrences pass, 2 otherwise.
The optional frequency tolerance does not affect that exact gate.

## Mathematical caveats

- Since s has irreducible minimal polynomial \(x^3-x^2-1\) over Q, any
  integer matrix having s as an eigenvalue has characteristic polynomial
  divisible by this degree-three polynomial. Its dimension is therefore
  at least three (indeed a multiple of three).
- The displayed M has characteristic polynomial \(x^3-x^2-1\), trace 1,
  sum of principal 2x2 minors 0, and determinant 1; it is unimodular.
- The other two eigenvalues have modulus \(s^{-1/2}\), so the normalized
  Perron convergence has subdominant spectral ratio \(s^{-3/2}\). Their
  arguments yield an oscillatory correction.
- These are necessary signatures, not sufficient evidence of geometric
  realization. Different matrices can share a characteristic polynomial.
  Stage 2 must still inspect actual incidence, multiplicities, orientation,
  and Anti-Twin action.
- The existing midpoint-subdivision benchmark \(D=2I_3\) does not pass this
  supergolden recurrence filter.
