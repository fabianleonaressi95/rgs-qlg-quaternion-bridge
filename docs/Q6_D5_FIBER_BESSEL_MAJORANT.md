# Q6–D5 complete-fiber convergence: explicit Bessel majorant

## Purpose and status

The exact integer quotient in PR #153 has a rank-124 kernel. The complete
positive marginal is

\[
W_{\rm eff}([c])=\sum_{k\in K_{\mathbb Z}}
\prod_{e=1}^{192}\frac{I_{|(c+k)C)_e|}(\beta_t)}{I_0(\beta_t)}.
\]

This note strengthens the convergence argument with an explicit
one-dimensional majorant. It does **not** claim that the 124-dimensional
sum has been numerically evaluated.

## Factorial bound

For every integer \(m\geq0\) and \(\beta>0\), termwise comparison of
the power series gives

\[
I_m(\beta)\leq \frac{(\beta/2)^m}{m!}I_0(\beta),
\qquad
\frac{I_m(\beta)}{I_0(\beta)}
\leq\frac{(\beta/2)^m}{m!}.
\]

For an edge-flux vector \(n\in\mathbb Z^{192}\), put
\(N=\|n\|_1\). By the multinomial theorem,

\[
\prod_e\frac{I_{|n_e|}(\beta)}{I_0(\beta)}
\leq
\frac{(\beta/2)^N}{\prod_e|n_e|!}
\leq \frac{(96\beta)^N}{N!}.
\]

## Exact coefficient-count bound

The fundamental-cycle matrix \(C\in\mathbb Z^{129\times192}\) has full
row rank. Select 129 independent columns, forming an invertible integer
minor \(A\). Write \(M=\|A^{-1}\|_\infty\), an exactly computed rational
number. For \(n=(c+k)C\), the selected edge coordinates obey
\((c+k)=n_P A^{-1}\), hence

\[
\|k\|_\infty
\leq M\|n\|_1+\|c\|_\infty.
\]

Therefore the number of kernel translations with flux norm exactly \(N\)
is bounded above by the number of all integer coefficient vectors in the
corresponding box:

\[
\#\{k\in K_{\mathbb Z}:\|(c+k)C\|_1=N\}
\leq
\left(2(MN+\|c\|_\infty)+1\right)^{129}.
\]

The kernel restriction can only reduce this count. Combining the two
inequalities yields the convergent series

\[
W_{\rm eff}([c])\leq
\sum_{N=0}^{\infty}
\left(2(MN+\|c\|_\infty)+1\right)^{129}
\frac{(96\beta_t)^N}{N!}.
\]

The ratio of successive smooth majorant terms is bounded by

\[
q_N=
\left(
\frac{2(M(N+1)+\|c\|_\infty)+1}
     {2(MN+\|c\|_\infty)+1}
\right)^{129}
\frac{96\beta_t}{N+1},
\]

which tends to zero. Once \(q_N<1\), monotonicity gives a geometric
bound on the remaining tail. The implementation evaluates this bound in
log space and computes \(M\) from an exact rational inverse.

## Computational contract

Run:

    python -m pytest -q tests/test_rgs_q6_d5_fiber_bessel_majorant.py

The audit checks the Bessel inequality over orders 0–80, selects the exact
cycle minor, and returns a finite logarithmic upper bound.

**Important limitation:** this count bound is deliberately coarse. It may
be far too large to provide a useful numerical interval. It is a rigorous
analytic domination strategy, not an enumeration of the complete
rank-124 fiber and not a numerical value for \(W_{\rm eff}\).

## Bergman connection

The Hua cutoff normalizer from PR #149 is a separate D5 radial integral.
The present result bounds an electric lattice marginal only. A justified
identification of a Wilson-derived radial density with the D5 generic-norm
density still requires a constructed spatial/canonical sector and a
derived pushforward measure. No equality between \(\beta_t\), a lattice
spacing, and the Hua cutoff \(\epsilon\) is asserted.
