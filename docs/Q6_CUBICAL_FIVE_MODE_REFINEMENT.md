# Q6 cubical refinement: exact five-mode block

## Statement

Take the cubical product complex formed by six copies of an interval
subdivided into m edges, with the standard tensor-product cochain inner
product and absolute (Neumann) boundary convention. Put h=1/m.

For one interval, the nonzero spectrum of the 0-form and 1-form Hodge
Laplacians is
\[
\mu_j(m)=2-2\cos\left(\frac{j\pi}{m+1}\right),
\qquad j=1,\ldots,m.
\]
The scalar 0-form spectrum also has its zero mode.

In the six-dimensional tensor-product complex, consider the **joint**
fundamental index (1,1,1,1,1,1). In the 1-form space there are six
component vectors, one for each differential direction, all with total
eigenvalue
\[
\lambda_{\rm fund}(m)=6\mu_1(m).
\]
The gradient of the scalar product eigenfunction spans the diagonal vector
(1,1,1,1,1,1)/sqrt(6) in this six-dimensional multiplicity space.
Its orthogonal complement is coexact and has dimension
\[
6-1=5.
\]
Thus the fundamental **joint** block has an exact rank-five coexact sector
for every m>=1. A fixed orthonormal basis is supplied by the 6-by-5
Helmert matrix in src/rgs_q6_cubical_refinement.py.

The rescaled eigenvalue obeys
\[
\lim_{m\to\infty}h^{-2}\lambda_{\rm fund}(m)
=\lim_{m\to\infty}6m^2
\left[2-2\cos\left(\frac{\pi}{m+1}\right)\right]
=6\pi^2.
\]

At m=1, the cell counts are 64 vertices, 192 edges and 240 square
faces, and the eigenvalue is 12. These recover the Q6 cubical data used by
the existing five-mode Maxwell matching.

## What is proved, and what is not

**Proved in this discrete tensor-product model:** the cell counts; the
six-dimensional fundamental joint block; the exact rank-one/coexact
rank-five splitting; and the rescaled eigenvalue limit.

**Not proved here:** that the total-eigenvalue-12 eigenspace is globally
isolated from all other joint blocks (degeneracies can occur); convergence
of Whitney finite-element Hodge Laplacians with their mass matrices;
Mosco convergence of the full quadratic forms; convergence of Wilson
measures or of the Wilson-marginalized effective dynamics; and any
identification of the limiting measure with the Bergman measure on D5.

The joint spectral labels are important: selecting only by the *sum*
eigenvalue can mix distinct tensor-product blocks. This is therefore a
rigorous finite-complex refinement statement, not yet the Wilson
continuum-limit theorem.

## Reproduce

Run:

    pytest -q tests/test_rgs_q6_cubical_refinement.py

The tests check the Q6 base counts and eigenvalue, the rank-five splitting,
orthonormality, and convergence to 6*pi^2. No workflow or generated
artifact is required.
