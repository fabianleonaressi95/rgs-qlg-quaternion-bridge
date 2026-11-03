# Chiral and Wilson audit — Phase 7

Date: 2026-11-03

## Scope

This phase implements a finite-graph covariant Wilson-Laplacian candidate on
the existing Q6/RGS edge set. It does **not** claim a Nielsen–Ninomiya theorem
for Q6, an exact count of eight Q6 fermion species, a continuum limit, or a
physical chiral theory.

## Implemented operator

For an oriented edge i -> j, let S_ji be a unitary fibre transporter and
let l_ij = ||x_j-x_i||. The unweighted graph Wilson form is

  W_ji = - r/(2 l_ij) S_ji,
  W_ij = - r/(2 l_ij) S_ji^dagger,
  W_ii += r/(2 l_ij) I,
  W_jj += r/(2 l_ij) I.

The volume-weighted operator is W_V = V^(-1/2) W V^(-1/2), with
V_i = det(E_i) > 0 repeated over the four spinor components.

For each edge its quadratic form contribution is proportional to

  || psi_j - S_ji psi_i ||^2,

so W is positive semidefinite for unitary links and r >= 0. It is Hermitian.
Under independent unitary changes of fibre frame G_i, with
S_ji -> G_j S_ji G_i^dagger, it transforms by block conjugation
W -> G W G^dagger. These are finite-dimensional algebraic statements.

The module now exposes:
- assemble_wilson_laplacian(...)
- assemble_dirac_wilson(...)

The latter sums the already implemented Hermitian spatial Dirac candidate
and the Wilson graph Laplacian. The construction uses the same geometric
link convention as the existing discrete Dirac module.

## What this does not establish

1. **Nielsen–Ninomiya applicability.** The classical no-go theorem assumes a
   translation-invariant periodic lattice and a momentum-space symbol with
   specified locality, hermiticity, and chiral properties. The finite
   reconstructed Q6 graph with nonuniform coordinates and vertex-dependent
   coframes is not such a lattice. Its finite spectrum cannot be labelled by
   a three-dimensional Brillouin zone without constructing an additional
   periodic family and proving its hypotheses.

2. **Eight species on Q6.** The formula
   D(p) proportional to sum_mu gamma_mu sin(p_mu a), and its 2^3 corner
   counting, describe a particular regular-lattice symbol. They are not a
   theorem about the 64-vertex Q6 operator. No 160^3 momentum-grid result or
   DOS ratio is produced by this repository phase.

3. **Wilson decoupling.** On a regular lattice, the standard Wilson symbol
   gives corner masses proportional to 2 r k/a. To claim the same mechanism
   here requires a periodic/refinement family, a symbol or controlled
   low-energy expansion, and spectral convergence as a -> 0. The finite
   graph Laplacian positivity test is not that derivation.

4. **Chirality.** The Wilson term is a scalar covariant Laplacian in spinor
   space. In the usual Dirac representation it commutes with gamma_5, while
   the naive kinetic term anticommutes with gamma_5; consequently the sum
   explicitly breaks the naive anticommutation relation. A modified
   Ginsparg–Wilson relation is not implemented.

5. **Lorentzian gauge group.** The stored links in the current implementation
   are SU(2) lifts of spatial polar rotations. They are not generic
   SL(2,C) Spin(3,1) links. No boost sector, noncompact Haar measure, or
   Lorentzian Palatini-Holst path integral is implemented.

6. **Möbius axial flux.** Non-orientability (w_1 != 0) does not alone imply a
   nonzero axial anomaly or a non-cancelling face flux. Such a statement
   requires a globally well-defined current on the relevant bundle, an
   orientation/local-system convention for each face, and a derived discrete
   Ward identity including regulator terms and boundary contributions.
   No such anomaly calculation is included in this phase.

## Correct status of the supplied report

| Claim in the supplied report | Status after repository audit |
|---|---|
| Pointwise curved Clifford relation | PROVED algebraically from invertible coframe; numerically tested |
| 2I finite representation checks | NUMERICALLY VERIFIED for selected representations; not a full irrep catalogue |
| Discrete spatial Dirac hermiticity | NUMERICALLY VERIFIED by construction/tests |
| Wilson graph Laplacian hermiticity and positivity | PROVED edgewise for unitary links, tested numerically |
| Wilson link gauge covariance | PROVED algebraically for fibre-frame transformations, tested numerically |
| Eight doublers on Q6 | OPEN / not established |
| DOS ratio 8.0000 from a 160^3 grid | NOT PRESENT in the repository evidence inspected |
| Chiral staggering suppression on Möbius | NOT ESTABLISHED |
| Net axial flux -0.041880 and its topological interpretation | NOT ESTABLISHED |
| Full Weiler theorem compliance | NOT ESTABLISHED; the theorem map records hypotheses and unfinished obligations |
| Spin(3,1) Palatini-Holst dynamics | OPEN; no such implementation currently exists |

## Next mathematical gates

- Specify a sequence of periodic 3D/4D complexes with mesh scale a -> 0,
  bounded geometry, and a controlled coframe/connection family.
- Derive the free and gauge-covariant Fourier symbols on the periodic
  reference lattice; count zeros and chirality indices explicitly.
- Compare naive, Wilson, and (if feasible) overlap/Ginsparg–Wilson operators.
- Prove or numerically substantiate spectral convergence with mesh refinement,
  rather than infer a continuum result from one finite graph.
- Define a discrete axial current from the action by local chiral variation
  and derive the exact lattice Ward identity, including the Wilson breaking
  term and any boundary/local-system contribution.
- Keep the Möbius orientation cocycle, spatial SU(2) frame transport, and any
  internal Yang–Mills gauge field as distinct structures unless a coupling
  and its transformation law are explicitly defined.
