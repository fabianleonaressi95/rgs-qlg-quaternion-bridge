# RGS–PSU–PRISM–QLG: recovered original artifacts

This integration records original artifacts recovered from the user's Library bundles. It distinguishes exact imported data from independently recomputed checks and from still-unavailable objects. No synthetic replacement matrices are introduced.

## Recovered sources

### A14 Jacobian bundle
Source archive: `A14_RGS_PSU_QLG_JACOBIAN_BUNDLE.zip` (SHA-256 `80442238349cd6ddc720dabae93d083075dea88decba8acf0a78c3ec03e0abcc`).

Contains:
- `A14_RGS_PSU_QLG_JACOBIAN.py`
- `A14_jacobian_matrices.npz`
- `A14_jacobian_report.json`
- `A14_jacobian_derivation.md`
- `RGS_PSU_all_matrices.npz`

The A14 matrices archive has:
- `X`: (64,3), float64
- `J_constraint`: (4,192), float64
- `Q_archived`: (192,188), float64
- `normal_basis_svd`: (192,4), float64
- `normal_basis_labelled`: (192,4), float64
- `Q_tangent_svd`: (192,188), float64
- `P_tangent_archived`, `P_tangent_svd`: (192,192), float64
- singular values and centroid.

A14 defines three translation constraints and one fixed-radius constraint:
`F=(sum x_i, sum y_i, sum z_i, (sum ||X_i||^2-R_ref^2)/2)`.
Thus the last Jacobian row is the flattened embedding and the first three rows are uniform translations.

Reported A14 results:
- rank(J)=4; nullity=188
- singular values: 15.49193338482967, 8, 8, 8
- ||J Q_archived||_F = 1.2385465290863038e-14
- ||Q_archived^T Q_archived-I||_F = 1.3145467710370518e-14
- archived-vs-SVD tangent projector difference = 0
- projector complement residual = 9.013175987547113e-15
- finite-difference directional derivative relative residual = 4.110283102810646e-09
- centroid residual max = 2.220446049250313e-16.

These are the reported results in the recovered A14 report; the archive and report are preserved as source artifacts. They do not by themselves establish the Hessian spectrum.

### A13 QLG bundle
Source archive: `A13_RGS_PSU_QLG_KAGGLE_BUNDLE.zip` (SHA-256 `7ca0cdfeadbe1d693d4f776a5e22d957630b206b842f4b822f3afb578c488922`).

Its NPZ files include:
- `Q`: (192,188), float64
- `K0`: (188,188), float64
- `K_geometry_tangent`: (188,188), float64
- `X`: (64,3), float64
- `edges`: (192,3)
- `R_SO3`: (192,3,3)
- `U_SU2`: (192,2,2)
- 240 SU(2) holonomies on each of two sheets.
- Other fibre, SEA and QLG operators.

Important: `K0` and `K_geometry_tangent` are not automatically the same object as the original constrained energy Hessian (K_{tangent}=Q^T H Q). Their derivations and definitions must be compared before assigning that name.

### Exact Q6 incidence bundle
Recovered files:
- `Q6_exact_incidence_and_H241.npz`
- `RGS_Q6_incidence_matrices.npz)

The first archive contains:
- edges (192,3)
- squares (240,6)
- B1 (64,192), int8
- B2 (192,240), int8
- plaquette adjacency (240,240)
- Q6 adjacency and Laplacian (64,64)
- H241 (241,241) and its eigenvalues.

Independent checks performed on the recovered arrays:
- B1 @ B2 = 0 exactly (Frobenius norm 0)
- rank(B1)=63
- rank(B2)=129
- nonzero counts: B1=384, B2=960
- entries belong to {-1,0,1}.

The second archive contains a higher-dimensional incidence chain D1,...,D6 with shapes (64,192), (192,240), (240,160), (160,60), (60,12), (12,1). These are distinct from B1/B2 and must retain their own cell-complex interpretation.

## Recovery status matrix

| Requested object | Recovered status | Notes |
|---|---|---|
| J (4x192) | FOUND / REPORTED VALIDATED | Exact A14 source |
| Q (192x188) | FOUND | A14 archived basis and A13 copies |
| K_tangent (188x188) | PARTIAL | A13 has K0 and K_geometry_tangent; exact relation to constrained energy Hessian not established here |
| Original H (192x192) | NOT FOUND IN THESE BUNDLES | Do not infer it from K |
| Spectral projectors | PARTIAL | A14 has P_T/P_N; A13 has P_G; full negative/zero/positive spectral projectors not identified here |
| Q6 incidence | FOUND | B1/B2 and edges/squares recovered |
| Explicit 2I (120 quaternions) | NOT YET RECOVERED | Existing quaternion/SU2 code is not itself a finite 2I element table |

## Required next validation
1. Compare A14 Q, A13 Q and all archived tangent bases by subspace projectors.
2. Determine the exact definitions and provenance of K0 and K_geometry_tangent.
3. Locate the original energy functional and H (192x192), including any Lagrange-multiplier correction.
4. Recover complete spectral projectors and verify idempotence, orthogonality, ranks and completeness.
5. Recover the explicit 120-element binary icosahedral quaternion set and test closure, inverses and the projection to A5.
6. Preserve all NPZ bundles byte-for-byte and attach their SHA-256 hashes when they are added to the Git repository.

## Provenance
The A14/A13 archives and Q6 NPZ files were recovered from the user's Library artifacts during this integration. Binary source archives are not silently regenerated. This document is an index and validation record, not a substitute for the original binary arrays.
