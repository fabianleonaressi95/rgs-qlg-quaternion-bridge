# Mosco prototype for the Euclidean Wilson form on a flat periodic torus

## Scope

This note proves a reference result for a **regular flat periodic lattice family**. It is not a Mosco theorem for the fixed curved RGS-Q6 graph. The proof fixes the Hilbert-space embeddings explicitly and separates Mosco convergence from the stronger compactness needed for eigenvalue convergence.

## 1. Discrete Fourier family

Let \(\mathbb T_L^d=(\mathbb R/L\mathbb Z)^d\), \(a=L/N\), and let
\(\mathcal H_a\) be the finite-dimensional subspace of
\(L^2(\mathbb T_L^d;\mathbb C^4)\) spanned by the Fourier modes
\[
e_k(x)=L^{-d/2}e^{2\pi i k\cdot x/L},
\qquad
k\in\Lambda_N,
\]
where \(\Lambda_N\) is a centered Brillouin index set. This trigonometric
interpolation is an isometric embedding of the discrete Fourier space into
one fixed Hilbert space. Forms are extended by \(+\infty\) outside
\(\mathcal H_a\).

For Hermitian spatial Clifford matrices satisfying
\(\{\gamma^\mu,\gamma^\nu\}=2\delta^{\mu\nu}I\), define the Euclidean
Wilson symbol
\[
D_{E,a}(p)=
i\frac{2}{a}\sum_{\mu=1}^d\gamma^\mu\sin(p_\mu a)
+\frac{r(a)}{a}\sum_{\mu=1}^d[1-\cos(p_\mu a)]I.
\]
Since the scalar Wilson term commutes with the momentum-space naive term,
\[
D_{E,a}(p)^\dagger D_{E,a}(p)
=\sigma_a(p)^2I,
\]
where
\[
\sigma_a(p)^2=
\frac{4}{a^2}\sum_\mu\sin^2(p_\mu a)
+\frac{r(a)^2}{a^2}
\left(\sum_\mu[1-\cos(p_\mu a)]\right)^2.
\]

The discrete quadratic form is
\[
q_a(u)=\sum_{k\in\Lambda_N}
\sigma_a(2\pi k/L)^2\,\|\widehat u(k)\|_{\mathbb C^4}^2,
\qquad u\in\mathcal H_a,
\]
and \(q_a(u)=+\infty\) for \(u\notin\mathcal H_a\).

The candidate limit is
\[
q(u)=4\int_{\mathbb T_L^d}|\nabla u(x)|^2\,dx,
\qquad
\mathcal D(q)=H^1(\mathbb T_L^d;\mathbb C^4).
\]

## 2. Scaling and pointwise symbol convergence

Take
\[
r(a)=r_0(a/L)^\alpha,\qquad r_0>0.
\]
For each fixed physical momentum \(p\), Taylor expansion gives
\[
\frac{4}{a^2}\sum_\mu\sin^2(p_\mu a)
\longrightarrow 4|p|^2,
\]
and
\[
\frac{r(a)^2}{a^2}
\left(\sum_\mu[1-\cos(p_\mu a)]\right)^2
=O(r(a)^2a^2|p|^4).
\]
Thus the Wilson contribution to the fixed-mode form vanishes when
\(r(a)a\to0\), i.e. \(\alpha>-1\), and
\[
\sigma_a(p)^2\longrightarrow4|p|^2.
\]

At a Brillouin corner with k components equal to \(\pi/a\), the naive
symbol vanishes and the Wilson mass is \(2kr(a)/a\). It diverges iff
\(\alpha<1\). The joint window
\[
\boxed{-1<\alpha<1}
\]
therefore gives fixed-mode consistency and ultraviolet doubler decoupling.

## 3. Mosco liminf inequality

Suppose \(u_a\rightharpoonup u\) weakly in the fixed \(L^2\) Hilbert space.
For any finite set of continuum Fourier indices F, all its modes belong to
\(\Lambda_N\) for sufficiently large N. Positivity of every symbol gives
\[
q_a(u_a)\ge
\sum_{k\in F}\sigma_a(2\pi k/L)^2
\|\widehat u_a(k)\|^2.
\]
Weak convergence implies convergence of each fixed Fourier coefficient.
Using the pointwise symbol limit,
\[
\liminf_{a\to0}q_a(u_a)
\ge
\sum_{k\in F}4|2\pi k/L|^2\|\widehat u(k)\|^2.
\]
Taking the supremum over finite F yields
\[
q(u)\le\liminf_{a\to0}q_a(u_a).
\]
If the right side is finite, the Fourier characterization also implies
\(u\in H^1\).

## 4. Recovery sequence

First take a trigonometric polynomial u. For sufficiently large N all its
Fourier modes are represented exactly in \(\mathcal H_a\). Set \(u_a=u\).
The finite sum and pointwise symbol convergence imply
\[
q_a(u_a)\longrightarrow q(u).
\]

Trigonometric polynomials are dense in \(H^1(\mathbb T_L^d;\mathbb C^4)\).
For general u in the form domain, approximate by trigonometric polynomials
and use a diagonal sequence. This gives a strong recovery sequence. If
u is outside H1, the recovery upper bound is vacuous because q(u)=+infty.

Hence, for \(\alpha>-1\), the Fourier-embedded forms \(q_a\) Mosco-converge
to q. The lower endpoint is the consistency threshold for this Wilson
scaling.

## 5. Spectral compactness and the upper endpoint

Mosco convergence alone should not be confused with convergence of every
eigenvalue of the finite matrices. If \(\alpha\ge1\), corner masses need
not diverge; high-frequency lattice states can remain at bounded energy.
This can obstruct compactness of form sublevel sets and leave spurious
low-energy spectral branches.

For \(\alpha<1\), \(r(a)/a\to\infty\). If a sequence of Brillouin
momenta has physical magnitude tending to infinity, either its dimensionless
momentum approaches zero, in which case the centered derivative term grows,
or it approaches a nonzero point of the Brillouin torus. At such a point,
either a sine component stays nonzero or at least one component approaches
a pi-corner and the Wilson term grows like r(a)/a. Thus the energy escapes
to infinity along ultraviolet sequences. This supplies the additional
high-frequency compactness needed for spectral convergence in this flat
periodic model.

Therefore \(-1<\alpha<1\) is the natural joint window for Mosco form
convergence plus decoupling/compactness of the doubler sector in this model.

## 6. What is not proved for RGS-Q6

The proof above uses translation invariance, an exact Brillouin torus,
trivial links, a flat coframe, and a Fourier interpolation into a fixed
Hilbert space. None of these structures is currently available as a
refinement family for the reconstructed 64-vertex RGS-Q6 graph.

To transfer the argument to RGS–PSU, the project still needs:
1. a sequence of Q6-derived or related complexes with mesh parameter h tending
   to zero and controlled geometric distortion;
2. compatible discrete Hilbert spaces and explicit embeddings into a common
   continuum space;
3. a vierbein-weighted discrete Hodge/volume measure with uniform ellipticity;
4. convergence and uniform control of polar spin links to the intended spin
   connection, including holonomy consistency;
5. uniform lower bounds and ultraviolet coercivity for the Wilson forms;
6. liminf and recovery-sequence proofs on the actual nonuniform graph family.

Thus this is a **proved periodic-torus prototype**, not a Mosco theorem for
RGS-Q6. The numerical sampler is `src/mosco_torus.py`, with tests in
`tests/test_mosco_torus.py`.
