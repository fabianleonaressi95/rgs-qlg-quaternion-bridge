# Q6 geometric Hodge metrics from cubical Whitney forms

## Construction

The canonical Q6 complex is the 2-skeleton of the unit six-cube
\([0,1]^6\subset\mathbb R^6\). Every p-cell has a canonical set of active
coordinate axes and fixed endpoint coordinates.

For each one-dimensional fixed coordinate, the endpoint hat functions are
\(\lambda_0(x)=1-x\) and \(\lambda_1(x)=x\), with Gram matrix

\[
M^{(0)}_{ij}=\int_0^1\lambda_i(x)\lambda_j(x)\,dx
=\begin{pmatrix}1/3&1/6\\1/6&1/3\end{pmatrix}.
\]

In an active coordinate the cubical Whitney factor is \(dx_i\), whose
squared norm integrates to one. A p-cell basis form is the tensor product
of these factors. Distinct active-axis orientations are orthogonal in the
Euclidean metric. Consequently each \(H_p\) is block diagonal by
orientation, with each block a tensor product of copies of \(M^{(0)}\)
over the fixed coordinates.

The canonical ordering is inherited from the Q6 incidence constructor:
64 vertices, 192 increasingly oriented edges, and 240 faces oriented by
axis pairs \(a<b\). No cell reordering is introduced.

## Metric adjoints and Hodge Laplacians

With coboundaries \(d_0=B_1^T\), \(d_1=B_2^T\) and positive Gram matrices
\(H_p\), the adjoint is

\[
d_p^*=H_p^{-1}d_p^T H_{p+1}.
\]

The implementation uses Cholesky factors \(H_p=L_pL_p^T\) to transform
cochains to orthonormal coordinates. The transformed coboundary is
\(A_p=L_{p+1}^T d_p L_p^{-T}\), and the symmetric Hodge operators are

\[
\Delta_0=A_0^TA_0,\quad
\Delta_1=A_0A_0^T+A_1^TA_1,\quad
\Delta_2=A_1A_1^T.
\]

This makes symmetry and positive semidefiniteness structural, up to floating
point roundoff. Harmonic dimensions are invariant under any positive
definite cochain metrics and equal the cellular Betti numbers \((1,0,111)\).

## Evidence and limits

- **PROVED algebraically:** the tensor-product Gram construction is positive
  definite; the adjoints follow from the stated inner products; finite
  dimensional Hodge decomposition gives harmonic dimensions equal to
  cellular cohomology.
- **NUMERICALLY VERIFIED:** the finite 64/192/240 matrices, symmetry,
  positivity, and harmonic dimensions are checked by the audit and tests.
- **OPEN:** a curved RGS-induced metric, a geometrically defined dual-cell
  construction, refinement compatibility, Mosco convergence, and coupling
  to a dynamical Yang-Mills action.
- **Not claimed:** this is not a circumcentric dual Hodge star, curved
  spacetime metric, TEGR construction, or physical QCD model.

Run:

```bash
PYTHONPATH=. python -m pytest tests/test_q6_whitney_hodge.py -q
PYTHONPATH=. python examples/q6_whitney_hodge_audit.py
```
