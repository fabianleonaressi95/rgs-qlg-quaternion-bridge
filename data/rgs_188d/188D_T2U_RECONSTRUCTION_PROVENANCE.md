# 188D T2u Reconstruction Provenance

## Executive conclusion

The original `J_A` / `J_B` definitions are still not recoverable from the verified chain; they remain `NOT_RECOVERABLE`.

The mathematically explicit object that is recoverable from the verified local source is the T2u tangent triplet basis:

```python
B_T2u = T2U_DATA__B_T2u         # shape (188, 39)
K_T2u = T2U_DATA__K_T2u         # shape (39, 39)
first_idx = T2U_DATA__first_triplet_indices   # [0, 1, 2]
V_triplet = T2U_DATA__V_first_triplet         # shape (188, 3)
B_T2U_TAN = T2U_DATA__eigenvectors_T2u_tan[:, first_idx]
```

This is source-backed by the archived T2u diagnostic bundle. The direct check gives:

```python
np.linalg.norm(B_T2U_TAN - V_triplet) / np.linalg.norm(V_triplet)
```

which is effectively machine precision (≈ 8.5e-16), so the triplet basis is consistent with the archived T2u construction.

## Formula and dimensions

- `B_T2u`: `(188, 39)`
- `K_T2u`: `(39, 39)`
- `first_triplet_indices`: `(3,)` with values `[0, 1, 2]`
- `V_first_triplet`: `(188, 3)`
- `B_T2U_TAN` / `J0`: `(188, 3)`

The 3D triplet is selected by diagonalizing the projected T2u block and taking the first degenerate triplet:

```python
K_T2u = B_T2u.T @ K_coulomb @ B_T2u
K_T2u = 0.5 * (K_T2u + K_T2u.T)
vals, vecs = np.linalg.eigh(K_T2u)
B_T2U_TAN = B_T2u @ vecs[:, :3]
```

or equivalently, in the archived source, the same triplet is already stored as:

```python
B_T2U_TAN = T2U_DATA__V_first_triplet
```

## Provenance verdict

- `J_A` / `J_B` original provenance: `NOT_RECOVERABLE`
- `B_T2U_TAN` / `J0` provenance: `RECONSTRUCTED_FROM_EXPLICIT_FORMULA`
- `M_geom = B_T2U_TAN^T K^{-1} B_T2U_TAN`: mathematically valid downstream object if `K` is the tangent-space stiffness matrix.
