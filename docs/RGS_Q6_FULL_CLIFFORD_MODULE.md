# Full cubical Q6 complex and six-dimensional Clifford module

## Construction

Q6 is the 1-skeleton of the 6-cube, but the full cubical complex contains
cells in every degree:
\[
(n_0,\ldots,n_6)=(64,192,240,160,60,12,1),\qquad
\sum_k n_k=729.
\]
Cells are indexed by a base vertex and an ordered subset of free coordinate
directions. Orientations follow increasing coordinate order. The boundary
operator uses the cubical alternating-face formula, and the coboundary is its
transpose. The implementation verifies \(d^2=0\).

## Clifford module

The fibre is the exterior algebra \(\Lambda^*\mathbb C^6\), of dimension 64.
For each coordinate, the Hermitian generator is
\[
\gamma_j=\varepsilon_j+\iota_j,
\qquad
\gamma_i\gamma_j+\gamma_j\gamma_i=2\delta_{ij}I.
\]
The grading \(\Gamma=(-1)^{\mathrm{degree}}\) anticommutes with every
\(\gamma_j\). This is a concrete complex Clifford module, not merely a
name for the graph incidence operator.

## Relation to the graph Dirac–Kähler operator

The full cubical incidence operator gives \(D_K=d+d^\dagger\) on all 729
cochains. It is distinct from the 512-dimensional vertex-plus-edge truncation
introduced in the SU(2) link PR. This first Clifford milestone is algebraic:
it does not yet define a gauge-covariant spin Dirac operator on Q6, because
parallel transport on higher-dimensional cells and a compatible spin
connection have not yet been specified.

## Fermion doubling and scope

No claim of eliminating fermion doubling is made. The finite 6-cube has
boundaries and no periodic momentum torus, so the standard Brillouin-zone
doubling count cannot be imported unchanged. A physical doubling audit needs
a specified kinetic operator, boundary conditions, continuum scaling and
low-energy pole/zero count. Chiral grading alone does not establish a
single chiral fermion or evade Nielsen–Ninomiya assumptions.

## Audit

The test suite checks cell counts, \(d^2=0\), the Clifford anticommutator,
Hermiticity and chirality anticommutation. Boundary ranks are expected to be
\((63,129,111,49,11,1)\), consistent with a contractible 6-cube.
