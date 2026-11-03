# Phase 9 — Wilson refinement, overlap/Ginsparg–Wilson, and Dirac–Kähler

## 1. Refinement scaling of the Wilson parameter

For the centered spatial symbol already used in Phase 8,

\[
D_0(p)=\frac{2}{a}\sum_{\mu=1}^{d}\gamma^\mu\sin(p_\mu a),
\qquad
W(p)=\frac{r(a)}{a}\sum_{\mu=1}^{d}[1-\cos(p_\mu a)]I.
\]

At fixed physical momentum p with ap small,

\[
W(p)=\frac{r(a)a}{2}|p|^2+O(r(a)a^3|p|^4).
\]

At a Brillouin corner with k components equal to pi/a, the naive symbol vanishes and

\[
M_k(a)=\frac{2kr(a)}{a}.
\]

Set \(r(a)=r_0(a/L)^\alpha\), with fixed physical box length L. Then

- the Wilson correction on fixed physical modes vanishes iff \(\alpha>-1\);
- the corner masses diverge iff \(\alpha<1\).

Therefore the simultaneous scaling window is

\[
\boxed{-1<\alpha<1.}
\]

The conventional choice \(r(a)=r_0\) corresponds to \(\alpha=0\): the low-momentum correction is O(a), while the doubler masses grow as O(1/a). The endpoints matter: alpha=1 leaves finite corner masses; alpha=-1 leaves a non-vanishing low-mode Wilson correction. This scaling statement concerns the symbol and does not depend on identifying Fourier momenta on Q6.

The implementation is `src/wilson_refinement.py`; the reproducible multi-resolution audit is `examples/wilson_refinement_audit.py`. It samples N=16,32,64,128,256 at fixed physical box length and records the naive dispersion error, the low-mode Wilson correction, and all 2^3 corner masses.

## 2. Three-dimensional corner accounting

There are 2^3=8 corners. One has no pi-components and is the physical low-momentum corner. The other seven have k=1,2,3 pi-components and acquire masses

\[
2r/a,\quad 4r/a,\quad 6r/a
\]

with multiplicities 3,3,1 respectively. Each corner retains the spinor multiplicity of the chosen representation. This is the exact Fourier-symbol account for the periodic cubic lattice. It is not a mode classification on the finite curved Q6 graph.

## 3. Finite overlap / Ginsparg–Wilson algebra

The Phase 8 matrix is Hermitian Hamiltonian-like and anticommutes with Gamma5 before adding Wilson. To use the standard overlap construction, the implementation explicitly converts this convention to a Euclidean gamma5-Hermitian kernel:

\[
D_E=iD_{0}+W,\qquad
H_W=\Gamma_5(aD_E-\rho I),
\]

where \(D_0=D_0^\dagger\), \(\{D_0,\Gamma_5\}=0\), \(W=W^\dagger\), and \([W,\Gamma_5]=0\). Thus \(H_W=H_W^\dagger\). Provided H_W has no zero eigenvalue, define

\[
D_{ov}=\frac{\rho}{a}
\left[I+\Gamma_5\operatorname{sign}(H_W)\right].
\]

The finite matrix audit checks

\[
\Gamma_5D_{ov}+D_{ov}\Gamma_5
=\frac{a}{\rho}D_{ov}\Gamma_5D_{ov}
\]

and gamma5-hermiticity. The function is `src/overlap_dirac.py`.

**Scope:** the current CI test verifies the exact finite-matrix Ginsparg–Wilson identity on a transparent 2x2 test kernel. It does not yet construct or diagonalize the overlap operator from the archived RGS-Q6 256x256 Wilson/Dirac pair. It also does not prove exponential locality; that requires a uniform spectral gap for H_W and suitable locality bounds on the underlying graph family. The continuum index theorem and anomaly statement remain open.

## 4. Dirac–Kähler comparison on the Q6 two-skeleton

The canonical Q6 complex provides boundary matrices B1 (64x192) and B2 (192x240). On cochains,

\[
d_0=B_1^T:C^0\to C^1,\qquad
d_1=B_2^T:C^1\to C^2,\qquad d_1d_0=0.
\]

With the unweighted Euclidean cochain inner products, define

\[
D_{DK}=d+d^\dagger
\]

on \(C^0\oplus C^1\oplus C^2\). Its dimension is 64+192+240=496. Its square is the block-diagonal Hodge Laplacian:

\[
D_{DK}^2=\operatorname{diag}
(d_0^\dagger d_0,
d_0d_0^\dagger+d_1^\dagger d_1,
d_1d_1^\dagger).
\]

The degree-parity operator \(\Gamma_{DK}=\operatorname{diag}(+I_{64},-I_{192},+I_{240})\) anticommutes with D_DK. Implementation: `src/dirac_kahler.py`.

**Scope:** this is a valid finite combinatorial Dirac–Kähler operator on the Q6 2-skeleton. It is not yet the metric-weighted DEC operator for the RGS vierbein. A physically geometric DEC version requires primal/dual cell volumes or a discrete Hodge star, compatible inner products, and a refinement family. Dirac–Kähler cochains and the four-component spinor operator are different representations; no equivalence is asserted here.

## 5. What this establishes, and what remains open

### PROVED algebraically / NUMERICALLY VERIFIED in finite tests
- Wilson symbol: fixed-momentum correction is O(r(a)a); Brillouin corner mass is exactly 2kr(a)/a.
- For r(a)=r0(a/L)^alpha, the simultaneous scaling window is -1<alpha<1.
- In d=3, seven nonphysical Brillouin corners are lifted with the expected k-dependent masses.
- The finite Q6 Dirac–Kähler matrix is Hermitian, squares to the unweighted Hodge Laplacian, and anticommutes with degree parity.
- The finite overlap construction satisfies the Ginsparg–Wilson relation when the kernel is gapped and its gamma5-hermiticity hypotheses hold.

### CANDIDATE / NOT YET ESTABLISHED
- Applying the overlap construction to the archived RGS-Q6 Wilson operator with the actual spin transporters.
- A vierbein-weighted DEC Hodge star and RGS Dirac–Kähler operator.
- Uniform overlap-kernel spectral gap and exponential locality on a sequence of refined RGS graphs.

### OPEN: Mosco and continuum limit
A Mosco proof cannot be inferred from a few Fourier eigenvalues or a finite-dimensional residual. It requires:
1. a specified sequence of geometric complexes with mesh h tending to zero;
2. Hilbert spaces and explicit isometric/controlled embeddings into a common continuum Hilbert space;
3. discrete quadratic forms with uniform lower bounds and equicoercivity;
4. the Mosco liminf inequality for every weakly convergent sequence;
5. a recovery sequence for every continuum form-domain element;
6. consistency of discrete Hodge stars, vierbein, link transport and spin connection;
7. treatment of boundary conditions, topology, and possible zero modes.

The Phase 9 refinement audit is a Fourier-symbol consistency and scaling study, **not a Mosco-convergence proof**. The finite Q6 graph is one fixed complex and does not itself supply a mesh-refinement sequence.

## Reproducibility

Run:

```bash
PYTHONPATH=. python -m pytest tests/ -v
PYTHONPATH=. python examples/wilson_refinement_audit.py
```

The JSON report is written to `results/json/wilson_refinement_audit.json`.
