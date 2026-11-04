# PSU two-sheet Dirac spectral bridge

## Existing repository components reused

This extension builds on src/discrete_dirac.py, src/overlap_dirac.py,
src/dirac_kahler.py, src/discrete_spin_geometry.py,
src/spinor_kirchhoff.py, and src/unified_discrete_action.py. It does not
replace or relabel any of them. In particular, the existing unified Q6
matter hopping is explicitly a covariant graph energy, not a Dirac kinetic
term; the existing frame lift is Euclidean Spin(4) polar alignment, not
automatically a Levi-Civita or TEGR connection.

## Definition

For Hermitian one-sheet operators D+, D-, unitary W, and real g,m+,m-:

H_PSU = [[D+ + m+ I, g W], [g W†, D- + m- I]].

The construction is Hermitian by block adjunction. W acts on the complete
site-spinor space. It is not identified with spatial SU(2) spin transport,
electromagnetic U(1), colour SU(3), or the Mobius orientation cocycle.

For identical sheets D+=D-=D, W=I and m+=m-=0, symmetric/antisymmetric
sheet combinations diagonalise the sheet factor, giving the exact identity

spec(H_PSU) = spec(D+gI) union spec(D-gI).

This is a finite-dimensional spectral identity, not a mass prediction.

## Chirality

The optional Gamma5 audit reports each sheet anticommutator, [Gamma5,W],
and the full block anticommutator. A scalar inter-sheet coupling commuting
with Gamma5 is mass-like and breaks naive anticommutation for nonzero g.
A chiral inter-sheet term requires a separately derived Clifford-odd block;
it must not be introduced as an arbitrary sign convention.

## Scope

- Fixed finite graph; no continuum or Mosco theorem.
- Hermitian spatial operator; not Lorentzian 4D Dirac evolution.
- No physical fermion-mass identification or curved-Q6 no-doubling theorem.
- W and g are not yet derived from the RGS/PSU action.
- U(1), SU(2), SU(3), Spin(4)/Spin(1,3), and sheet transport remain distinct.

## Reproduction

Run: python -m pytest tests/test_psu_two_sheet_dirac.py -q
