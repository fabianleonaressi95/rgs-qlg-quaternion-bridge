# Extended recovered matrix inventory (A13/A14)

This file records array-level contents inspected in the recovered original archives. Shapes and dtypes below are read from the NPZ members, not inferred from filenames.

## A14: Jacobian and common matrix source

Archive: `A14_RGS_PSU_QLG_JACOBIAN_BUNDLE.zip`
SHA-256: `80442238349cd6ddc720dabae93d083075dea88decba8acf0a78c3ec03e0abcc`

`A14_jacobian_matrices.npz`:
- X (64,3) float64
- J_constraint (4,192) float64
- Q_archived, Q_tangent_svd (192,188) float64
- normal_basis_svd, normal_basis_labelled (192,4) float64
- P_tangent_archived, P_tangent_svd (192,192) float64
- singular_values (4,) float64; centroid (3,) float64

The companion `RGS_PSU_all_matrices.npz` contains:
- X_Oh_reconstructed (64,3)
- Q_tangent (192,188)
- A_geometry_orthonormal, A_raw (188,64)
- P_G, K0, R_geo_tangent (188,188)
- Pi_AT (64,64), B_PG_A (188,64), U_polar (188,64)
- T_transport (188,188)
- Q6_states (64,6), Q6_adjacency and Q6_laplacian (64,64)
- Q6_antipodal (64,64)

## A13: QLG, fibre and SEA arrays

Archive: `A13_RGS_PSU_QLG_KAGGLE_BUNDLE.zip`
SHA-256: `7ca0cdfeadbe1d693d4f776a5e22d957630b206b842f4b822f3afb578c488922`

`psu_fibre_connection_matrices.npz`:
- X (64,3), Q (192,188), K0 (188,188)
- edges (192,3), R_SO3 (192,3,3), U_SU2 (192,2,2)
- D_fibre (384,128), D_fibre_biverso (768,128)
- K_fibre (128,128), K_fibre_biverso (256,256)
- holonomy_sheet0 and holonomy_sheet1 (240,2,2)
- wilson_angles (240,), D_tangent (576,188)
- K_geometry_tangent (188,188)

`sea_operators.npz`:
- D_kin (512,512) complex128
- V_fermion (512,512), K_B0 (128,128), V_boson and H_cross (128,128)
- edge_coframe (192,4), tau (64,), Q (192,188), K0 (188,188)
- fermion_singular_values (512,), boson_eigenvalues (128,)

`A13_RGS_PSU_QLG_matrices.npz`:
- Q (192,188), normal_basis (192,4)
- P_T and P_N (192,192)
- edges (192,3), edge_coframe and raw_coframe (192,4)
- R_SO3 (192,3,3), U_SU2 (192,2,2), tau (64,)
- correlation_twin and correlation_antitwin (4,3)
- exchange (8,8), correlation_observables (3,8,8)

The A13 archive also includes `A13_report.json`, `A13_tau_spectra.csv`, `A13_multifractal_partition.csv`, source Python, SEA and fibre audit bundles, and a Giacosa twin/antitwin bridge NPZ.

## Hopf/SU(2) two-seas bundle

Archive: `RGS_PRISM_two_seas_Hopf_SU2_bundle.zip`
SHA-256: `a2dfa00e6a8ece2a35e96cb5683b5bfb8482fb2661338c446dcfa68c02d4b10c`

Contains a report, README, builder script and `two_seas_hopf_su2_bundle_operator.npz` (94,930 bytes). This is an additional SU(2)/Hopf construction. Its existence does not establish that it contains a complete enumerated 120-element binary icosahedral subgroup. That claim requires explicit elements plus closure, inverse and group-order tests.

## Naming and mathematical safeguards

- `K0`, `R_geo_tangent`, and `K_geometry_tangent` retain their source names until their defining derivations establish equivalence to a constrained energy Hessian.
- (P_T) and (P_N) are tangent/normal projectors, not the complete negative/zero/positive spectral projectors of a Hessian.
- (P_G) is a geometric-sector projector, not automatically a Hessian spectral projector.
- SU(2) link matrices and 240 plaquette holonomies are not the same thing as a canonical explicit 2I element table.
- No arrays are regenerated or reordered in this inventory.
