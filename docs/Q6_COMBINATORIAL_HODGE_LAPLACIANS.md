# Weighted combinatorial Hodge Laplacians on Q6

**Status:** finite-dimensional construction; algebraic and numerical tests are reproducible. A metric-derived discrete Hodge star, continuum limit and physical interpretation remain OPEN.

## 1. Cochain complex

The canonical Q6 square complex has
\[
C^0\xrightarrow{d_0}C^1\xrightarrow{d_1}C^2,
\qquad d_0=B_1^T,\quad d_1=B_2^T,
\]
with cell counts
\[
(\dim C^0,\dim C^1,\dim C^2)=(64,192,240).
\]
The incidence convention is inherited unchanged from `src/q6_complex.py`. The chain identity (B_1B_2=0) is checked exactly over integers, hence (d_1d_0=0).

## 2. Positive cochain inner products

Choose strictly positive diagonal weights
\[
H_k=\operatorname{diag}(w_k),\qquad w_k(\sigma)>0,
\]
for (k=0,1,2). They define finite-dimensional inner products
\[
\langle u,v\rangle_k=u^T H_kv.
\]
These weights are explicit input data. Unit weights give the standard unweighted combinatorial inner products. Positive nonuniform weights provide a controlled weighted family.

The implementation returns symmetric representatives in orthonormalized cochain coordinates:
\[
\widehat\Delta_0=H_0^{-1/2}d_0^TH_1d_0H_0^{-1/2},
\]
\[
\widehat\Delta_1=H_1^{-1/2}
(d_0H_0d_0^T+d_1^TH_2d_1)H_1^{-1/2},
\]
\[
\widehat\Delta_2=H_2^{-1/2}d_1H_1d_1^TH_2^{-1/2}.
\]
Each is a Gram sum and is therefore positive semidefinite in exact arithmetic.

For unit weights these reduce to
\[
\Delta_0=B_1B_1^T,\quad
\Delta_1=B_1^TB_1+B_2B_2^T,\quad
\Delta_2=B_2^TB_2.
\]

## 3. Harmonic dimensions and cellular topology

For any strictly positive inner products, finite-dimensional Hodge theory identifies the kernel of the degree-(k) Laplacian with the (k)-th cohomology. The Q6 incidence ranks are
\[
\operatorname{rank}B_1=63,\qquad \operatorname{rank}B_2=129.
\]
Thus
\[
(b_0,b_1,b_2)=(64-63,192-63-129,240-129)
=(1,0,111).
\]
The graph alone has cycle rank 129; after attaching the 240 square faces, the cellular first Betti number is zero. The 111 harmonic 2-cochains are not 111 independent graph cycles.

## 4. Reproduction

From the repository root:
```bash
python -m pytest tests/test_q6_hodge_laplacian.py -q
python examples/q6_hodge_laplacian_audit.py
```

The tests verify the incidence formulas for unit weights, exact boundary-of-boundary, positive semidefiniteness, and harmonic dimensions for both unit and positive nonuniform weights.

## 5. Scope and next proof gates

**Established for this finite construction:** the incidence identities, the weighted Gram formulas, and the expected harmonic dimensions under positive diagonal inner products.

**Not established:** a geometric Hodge star from a six-dimensional primal/dual cell decomposition; a circumcentric or barycentric dual metric; convergence under refinement; Mosco convergence; a Yang–Mills continuum action; a TEGR vierbein interpretation; or any physical identification of the 111 harmonic modes.

The weights are not silently inferred from the Q6 bit embedding. A future metric-dependent construction must state its cell volumes, dual volumes, orientation and consistency conditions, and must compare against independent analytic identities.
