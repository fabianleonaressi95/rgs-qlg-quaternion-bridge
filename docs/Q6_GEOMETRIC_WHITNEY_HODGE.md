# Metric-derived cubical Hodge operators for Q6

**Status:** exact finite cubical incidence construction plus metric-derived finite-dimensional L2 inner products. A primal-to-dual Hodge-star map, curved RGS geometry and continuum convergence remain open.

## Why the 2-skeleton is insufficient

The existing Q6 square complex contains only degrees 0, 1 and 2:
`(64,192,240)`. It is not a six-dimensional cellulation. In particular, the 111-dimensional second cohomology of this 2-skeleton is not a topological invariant of the full six-cube; it is filled by higher-dimensional cells.

The full cubical decomposition of `[0,1]^6` has
[
N_k=\binom{6}{k}2^{6-k},\qquad
(N_0,\ldots,N_6)=(64,192,240,160,60,12,1).
]
Its cellular boundary maps satisfy `B_{k-1} B_k = 0` exactly. Their ranks are
[
(63,129,91,35,7,1),
]
so the full cube has Betti numbers `(1,0,0,0,0,0,0)`, as expected for a contractible space.

## Metric-derived inner products

Use the Euclidean metric on the unit cube. For each active axis set (I), the lowest-order tensor-product cubical differential-form basis has the form
[
\omega_{I,\epsilon}(x)=
\prod_{j\notin I}\lambda_{\epsilon_j}(x_j)\,dx_I,
\quad \lambda_0(x)=1-x,\quad\lambda_1(x)=x.
]
Different active-axis sets are pointwise orthogonal. For each fixed coordinate,
[
\int_0^1\lambda_a\lambda_b\,dx=
\begin{pmatrix}1/3&1/6\\1/6&1/3\end{pmatrix}_{ab}.
]
The degree-(k) mass matrix (H_k) is assembled as the tensor product of these one-dimensional Gram factors over fixed axes, with unit factors over active axes. It is therefore derived from the Euclidean L2 inner product, not selected as arbitrary positive weights.

With (d_k=B_{k+1}^T), the cochain Laplacian is
[
\Delta_k=H_k^{-1}d_{k-1}^TH_{k-1}d_{k-1}
+d_k^TH_{k+1}d_k,
]
omitting terms outside degrees 0 and 6. The code returns its symmetric representative in L2-orthonormal coordinates.

## Key distinction

The former 2-skeleton has (b_2=111). The full six-cube has (b_2=0): its 3-cells and higher-dimensional cells fill those 2-cycles. This is a useful consistency result, not a contradiction.

## Reproduction

```bash
python -m pytest tests/test_q6_geometric_hodge.py -q
python examples/q6_geometric_hodge_audit.py
```

## Scope

**Verified finite construction:** cell counts, orientations, exact boundary-of-boundary identities, positive-definite Euclidean Whitney mass matrices, positive semidefinite Hodge Laplacians, and harmonic dimensions equal to the full cube's Betti numbers.

**Not established:** a primal-to-dual discrete Hodge-star map, a dual cellulation, a non-flat RGS metric, vierbein/TEGR coupling, a refinement family, Mosco convergence, a continuum Laplace–Beltrami limit, Yang–Mills dynamics, or a physical interpretation of Q6.
