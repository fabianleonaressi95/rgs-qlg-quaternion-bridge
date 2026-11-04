# Minimal Einstein–Cartan–Dirac: effective four-fermion term

## Convention map

The upstream EC module defines the real current
    J5^a = i bar(psi) gamma5 gamma^a psi.
Define the raw bilinear A^a=bar(psi) gamma5 gamma^a psi. Thus J5=iA and
A=-iJ5. In the frequently used raw-bilinear convention, eliminating the
algebraic contorsion of minimally coupled Dirac matter gives
    L_4f = -(3 kappa/16) A_a A^a.
Since A_a A^a=-J5_a J5^a, the repository's real-current convention is
    **L_4f = +(3 kappa/16) J5_a J5^a.**
The sign conversion is essential: importing the raw coefficient without
the factor of i reverses the interaction.

Here kappa=8 pi G in c=1 units. The coefficient has mass dimension -2;
the current-current operator has dimension 6. Minimal EC has no torsion
kinetic term, so the contact interaction is local and torsion adds no
propagating initial data.

## Auxiliary-field elimination

Split the independent connection as omega=omega_LC+K. The curvature
contains terms quadratic in K plus a divergence; the symmetrized minimally
coupled Dirac action contributes the axial source linear in K. Varying K
gives the algebraic Cartan equation. Substituting its solution into the
quadratic-plus-linear auxiliary sector (completing the square) yields the
contact operator above.

The code implements the coefficient, density, functional derivative,
raw/real current conversion and multi-flavour sum-before-square. Summing
species currents before squaring retains cross-flavour contact terms.

## Evidence

- PROVED within minimal EC and the declared conventions: algebraic torsion
  elimination gives a local axial-axial interaction.
- NUMERICALLY VERIFIED: current-convention equivalence, finite-difference
  functional derivative, spinor-current-to-torsion/contact pipeline and
  multi-flavour cross terms.
- OPEN: explicit component-by-component substitution of full contorsion into
  the repository's first-order action with boundary terms; quantum loops/RG;
  Holst/Nieh-Yan/nonminimal couplings; cosmological finite-density averaging.
  Q-ball U(1) charge is not fermionic spin density.

## References

Bonder & Corral, Phys. Rev. D 97, 084001 (2018), doi:10.1103/PhysRevD.97.084001:
first-order EC and algebraic spin-sourced torsion.

*Torsion at Different Scales: From Materials to the Universe*, Universe 9
(2023) 516: explicit integration of algebraic torsion and axial contact term.
For nonminimal Holst/Immirzi extensions, vector-axial terms and coefficients
change; do not reuse this minimal-EC coefficient.
