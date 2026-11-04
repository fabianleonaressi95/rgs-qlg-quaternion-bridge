# RGS–PSU rank-three bundle and multifractal audit

## Finite vector bundle

The eight base labels form the vertices of the base cube $Q_3$ and the
fibre has rank three. The archived RGS seed provides positive orbit
representatives $p_b\in\mathbb R^3_{>0}$. Define local frames and fibre
metrics by

\[
D_b=\operatorname{diag}(p_b),\qquad h_b=D_b^T D_b.
\]

For an oriented base edge $a\to b$, the frame-induced transport is

\[
U_{b\leftarrow a}=D_b^{-1}D_a.
\]

It is metric compatible:
$U_{b\leftarrow a}^T h_b U_{b\leftarrow a}=h_a$.
For every base square, the ordered product telescopes to the identity.
Thus this canonical connection is **flat and pure gauge**. This is an
important no-go result: the archived coordinates alone do not generate
nonzero Yang–Mills curvature. A curved connection requires additional
independent link variables or a variational principle that determines them.

The gauge transformation at each base vertex is $G_b\in GL(3,\mathbb R)$:
$U'_{b\leftarrow a}=G_b U_{b\leftarrow a}G_a^{-1}$ and
$h'_b=G_b^{-T}h_bG_b^{-1}$. Face holonomy transforms by conjugation.

## Multifractal scaling

For normalized positive masses $\mu_i(\ell_n)$ on compatible nested
partitions, compute

\[
Z_n(q)=\sum_i\mu_i(\ell_n)^q,\qquad
\tau(q)=\lim_{n\to\infty}\frac{\log Z_n(q)}{\log\ell_n},
\qquad D_q=\frac{\tau(q)}{q-1}.
\]

At $q=1$, $D_1$ is obtained from the slope of Shannon entropy
$H_n=-\sum_i\mu_i\log\mu_i$ against $\log(1/\ell_n)$.
The implementation estimates slopes by least squares over at least three
scales. It requires strictly positive masses; zero-mass cells must be
removed consistently at every level or handled with a documented limiting
convention.

A nonlinear $D_q$ is evidence of multifractal scaling only if the partitions
are nested, the scaling range is stable, and finite-size/crossover effects
are controlled. The routine is a diagnostic; it does not infer a physical
measure from the 64-point seed.

## Reproducible tests

The tests check metric compatibility, flatness, gauge covariance, and two
synthetic scaling controls: a uniform one-dimensional measure ($D_q=1$)
and a non-uniform binomial cascade (nonconstant $D_q$). The latter is a
method validation, not evidence that RGS–PSU itself is multifractal.
