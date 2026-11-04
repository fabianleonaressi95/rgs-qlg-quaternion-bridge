# Exact integer five-mode sector, A5 symmetry, and the D5 measure

## 1. Rationality and integer fluxes

Let D0 be the integer vertex-edge incidence matrix of Q6 and let
L1=D1^T D1 be its integer 1-cochain Hodge operator. The selected five-mode
sector is the solution space of the simultaneous equations

D0^T x = 0,
(L1 - 12 I_192) x = 0.

Both matrices have integer entries. Computing the nullspace over Q gives
dimension five. Every rational null vector can be multiplied by a common
denominator to produce an integer vector. Therefore

dim_Q ker_Q([D0^T; L1-12I]) = 5
and
rank_Z(ker_Z(D0^T) intersect ker_Z(L1-12I)) = 5.

The implementation uses exact SymPy rational row reduction, clears
denominators, and verifies the resulting integer vectors with exact
integer matrix products. This is stronger than rounding the floating-point
eigenvectors from the SVD.

Important lattice qualification: the five displayed integer vectors are
independent generators of a rank-five sublattice. They are not claimed to
be a saturated basis of the complete integer intersection. A Smith/Hermite
normal-form saturation audit remains open.

## 2. A5 invariance

The 60 proper icosahedral rotations act on Q6 edges by signed permutation
matrices P1 with integer entries. Since these are cubical cochain
automorphisms,

D0^T(P1 x)=0 and L1(P1 x)=12(P1 x)

whenever x satisfies the two defining equations. Thus the full integer
intersection is A5-invariant, not merely invariant up to floating-point
tolerance.

The orthonormal real five-mode frame realizes the irreducible five-
dimensional A5 representation. Its matrices rho obey rho^T rho=I.
Consequently, on the complexified five-mode coordinates,

z_bar.z and z.z

are invariant under z -> rho z.

## 3. D5 type-IV Bergman measure

For the Cartan type-IV domain D5, define the generic norm

N(z)=1-2 z_bar.z+|z.z|^2,  z in C^5.

The Bergman kernel has the form K(z,z)=C N(z)^(-5). With the standard
Bergman metric g_i_barj=partial_i partial_barj log K, its determinant is

det(g_i_barj)=10^5 N(z)^(-6).

Therefore the Bergman volume density is A5-invariant because A5 preserves
both z_bar.z and z.z. The implementation audits the generic norm directly
under all 60 group elements and reports the density 10^5 N^-6.

At the origin, g=10 I_5, so det(g)=10^5. Under the tangent normalization
z_D=z/sqrt(10), the real ten-dimensional coordinate Jacobian is 10^-5.
The product is

det(g(0)) * |Jac_R(z -> z_D)| = 10^5 * 10^-5 = 1.

This is an exact local normalization identity. It does not determine a
global nonlinear map from the compact Q6 quotient to D5, nor the global
volume of D5.

## 4. Epistemic status

**Exact algebraic result:** the Hodge-12 Gauss sector is defined over Q and
contains five independent integer electric-flux vectors.

**Exact symmetry result:** the full A5 action preserves the integer
intersection because it acts by integer cochain automorphisms preserving
D0 and L1.

**Finite-dimensional geometric result:** the A5 five-dimensional
representation preserves the D5 type-IV generic norm and hence its Bergman
density. The origin tangent-coordinate Jacobian cancels the Bergman
density normalization.

**Still open:** saturation of the rank-five integer lattice; the relation
between the compact integer-flux Hilbert space and the real Gaussian
oscillator sector; a global Q6-to-D5 map; derivation of the Bergman measure
from the microscopic RGS-PSU action; Wyler's global volume prescription;
and any derivation of the fine-structure constant.

Run:
pytest -q tests/test_rgs_q6_integer_a5_d5_measure.py
