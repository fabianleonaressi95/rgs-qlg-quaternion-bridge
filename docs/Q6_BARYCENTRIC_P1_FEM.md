# Q6 barycentric P1 finite elements on the edge skeleton

**Status:** NUMERICALLY VERIFIED finite-dimensional construction; continuum and physical interpretation OPEN.

## Geometric embedding

Use the canonical six-bit embedding
\[
X_v=(b_0(v),\ldots,b_5(v))\in\{0,1\}^6\subset\mathbb R^6,
\]
with the same big-endian bit convention as `src/q6_complex.py`. The 64 vertices and 192 Q6 edges are retained in the repository's canonical order. Every edge has Euclidean length one. The uniform vertex barycentre and covariance are
\[
X_B=(1/2,\ldots,1/2),\qquad
\Sigma_B=\frac14 I_6.
\]

## Element matrices

Each Q6 edge is treated as a one-dimensional affine simplex. For an oriented edge (e=(u,v)), length (L_e=|X_v-X_u|>0), the P1 element matrices are
\[
M_e=\frac{L_e}{6}\begin{pmatrix}2&1\\1&2\end{pmatrix},\qquad
K_e=\frac1{L_e}\begin{pmatrix}1&-1\\-1&1\end{pmatrix}.
\]
Assembly gives global (M) and (K). With (B_1) the canonical node-edge incidence matrix,
\[
K=B_1\,\operatorname{diag}(L_e^{-1})B_1^T.
\]
For the canonical hypercube all (L_e=1), so (K=B_1B_1^T), the unweighted Q6 graph Laplacian. Since Q6 is connected, (K\succeq0) and its kernel is the constant mode; the assembled consistent edge mass is positive definite because every vertex is incident to edges.

The edgewise ambient gradient of a scalar nodal field is
\[
(\nabla_e f_h)=\frac{f_v-f_u}{L_e^2}(X_v-X_u).
\]
This is the tangential gradient on each edge, not a full-dimensional gradient in \(\mathbb R^6\).

## Scope and limits

This module connects the existing affine barycentric P1 element to the exact Q6 incidence ordering. It uses only the one-skeleton. The 240 square faces remain available in the chain complex, but this module does not construct a 2-form Hodge star or a 6D volume finite element. It does not identify the covariance with a physical metric, nor the stiffness with a Lorentzian/gravitational operator. It makes no continuum, Mosco, TEGR, QCD or TOE claim.

## Reproducible audit

Run:

```bash
python -m pytest tests/test_q6_barycentric_fem.py -q
```

The tests check canonical counts and edge lengths, barycentre/covariance, positive definiteness of (M), the constant kernel and graph factorisation of (K), exact (B_1B_2=0), coordinate-field gradients/energies, translation invariance, anisotropic edge lengths, and rejection of collapsed embeddings.
