# T031 — 2I-equivariant Rumin bridge: geometry, representations and proof boundary

**Purpose.** This document is the implementation-side companion to `theroem/theorems/T031-rumin-2i-equivariance.md`. It connects the existing 120-quaternion 2I implementation to the geometric descent of the three-dimensional Rumin complex, without duplicating the analytic torsion proofs T021/T024.

## Canonical conventions

- Unit quaternions, scalar-first `(w,x,y,z)`.
- Hamilton product.
- Deck action on `S^3=SU(2)`: **left** multiplication `q -> g*q`.
- Hopf/Reeb circle: **right** multiplication `q -> q*h(t)`.
- The two actions commute by associativity. The binary icosahedral group is **not** a subgroup of the Hopf circle.
- Flat twist convention: lifted coefficient forms obey `u(g*q) = rho(g) u(q)`. The associated-bundle quotient must use the matching convention.

## Existing implementation (do not duplicate)

- `src/binary_icosahedral.py`: deterministic 120-element table and exhaustive floating-point closure/inverse/unit-norm/SO(3)-quotient audit.
- `src/representations_2i.py`: trivial, fundamental 2D SU(2), and vector 3D SO(3) representations; homomorphism/unitarity checks; selected intertwiner spaces.
- `tests/test_binary_icosahedral.py` and `tests/test_representations_2i.py`: current finite tests.

These are selected representations only. They are **not** the complete nine-irrep unitary catalogue used by the Molien/torsion calculation in T021.

## Geometric theorem imported from the analytic side

For a contactomorphism preserving the contact data, pullback commutes with the Rumin differentials. In particular, for the middle differential, the unique corrected representative `beta_alpha` is preserved by pullback, so
`D_R(F^* alpha) = F^*(D_R alpha)`.

The free left 2I action therefore makes the Rumin complex on the cover a 2I-complex. Its invariant part is the untwisted quotient complex; rho-equivariant coefficient forms descend to the flat bundle complex for a unitary representation rho.

The right Hopf circle commutes with the left deck action and descends to the quotient Reeb flow. Do not write that 2I lies inside the Hopf circle.

## Isotypic projector contract

For a unitary 2I-module U and an irreducible character chi_rho of dimension d_rho:
`P_rho = (d_rho/120) sum_g conjugate(chi_rho(g)) U(g)`.

Required identities:
- `P_rho^2 = P_rho = P_rho^*`;
- `P_rho P_sigma = 0` for inequivalent irreducibles;
- `sum_rho P_rho = I` when the complete irreducible list is used;
- an equivariant Rumin Laplacian commutes with every `P_rho`.

A character projector is not a matrix-unit basis and does not generate the missing irreducible matrices by itself.

## Work split

### This repository (implementation)
1. Preserve the canonical 2I element order and quaternion conventions.
2. Add a deterministic left/right action audit on unit quaternions.
3. Add a small, explicit test that left multiplication by every 2I element commutes with right Hopf multiplication for fixed sample unit quaternions and circle parameters.
4. Record selected-representation scope; do not label the three implemented representations as all nine irreducibles.
5. Later: export the full nine-irrep data only from a source with audited character ordering and exact algebraic entries.

### `theroem` (analysis)
1. Prove naturality/descent of the Rumin differential from its normalization definition.
2. Keep the analytic hypotheses and determinant-line/torsion conventions in T021/T024.
3. Do not claim a Q6/PRISM continuum limit from finite group equivariance.

## Reproducible acceptance checks

- Existing finite group tests: `pytest -q tests/test_binary_icosahedral.py tests/test_representations_2i.py`.
- New left/right action test: verify `g*(q*h) == (g*q)*h` and unit norm to tolerance; this is a finite arithmetic check of the implemented convention, not a proof of CR analytic properties.
- Report maximum residual and tolerance; no rounded output may be used to assert exact closure.
- The analytic statement is proved in the theorem document, not by this test.

## Evidence classification

- **PROVED (mathematical):** associativity gives commutation of left/right multiplication; the free left 2I action gives the regular quotient; naturality gives Rumin descent, under the stated contact-preservation hypotheses.
- **NUMERICALLY VERIFIED:** existing 120-element group and selected representation audits.
- **OPEN:** full nine-irrep model export, exact arithmetic projector/matrix-unit checks, and any Q6/PRISM-to-Rumin convergence theorem.
