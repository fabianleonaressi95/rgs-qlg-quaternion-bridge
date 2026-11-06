# Intrinsic Ellis-RGS spinorial Wilson-Dirac operator

**Status:** finite-mesh operator assembled from the explicit Ellis Levi-Civita SU(2) connection. This is not a continuum convergence theorem.

## Definition

At each vertex choose the oriented orthonormal tangent frame of the Ellis hypersurface. For each canonical edge i<j use the same affine-u / shortest-great-circle path as the spin connection, its intrinsic Ellis length, the unit path tangents at both endpoints expressed in their local frames, and the SU(2) Levi-Civita transporter U_ji from i to j.

Let c_i and c_j be the endpoint Clifford contractions with spatial gamma matrices, and let S(U_ji)=diag(U_ji,U_ji) be the four-component spinor lift. The Hermitian naive spatial matrix has edge block

H_ji = i (c_j S(U_ji) + S(U_ji)c_i)/(2 ell_ij),

with the reverse block its adjoint.

The Wilson graph term uses diagonal weights r/(2 ell_ij), off-diagonal blocks -r S(U_ji)/(2 ell_ij), and the adjoint reverse block. Both matrices are normalized by the lumped P1 vertex measure, one quarter of each incident tetrahedron volume. The Euclidean operator is D_E=iH+W, hence gamma5-Hermitian.

## Implementation

- src/ellis_rgs_spin_dirac.py
- tests/test_ellis_rgs_spin_dirac.py

Run focused tests:
PYTHONPATH=. python -m pytest tests/test_ellis_rgs_spin_dirac.py -q

Run the complete level-0 refinement graph:
PYTHONPATH=. python examples/run_ellis_rgs_spin_wilson_dirac_audit.py --level 0 --radial-base 2 --steps 12 --output ellis-spin-wilson-audit.json

The assembler uses the actual embedded R4 hypersurface, local tangent frames, intrinsic path lengths and geometrically integrated SU(2) links. It does not reinterpret the archived 64x3 edge coframe as an Ellis vierbein.

## Audits and limits

Tests cover gamma5 hermiticity, positivity of the Wilson graph term, local SU(2) frame-gauge covariance and rejection of incomplete link data. These are finite-dimensional algebraic checks.

Open: prove uniform shape regularity and measure equivalence, derive the consistency order of the path/integrator scheme, establish discrete Poincare/Rellich compactness, construct recovery sequences and prove Mosco/strong-resolvent convergence. The finite graph Wilson term has not yet been shown to remove all curved-lattice doublers.
