# 188D J_A / J_B first-principles reconstruction

## Executive conclusion

The original mathematical definition of `J_A` and `J_B` is not recoverable from the verified local sources.

The only explicit formulas found in the recovered notebook are synthetic mock constructions, not a verified physical definition:

```python
J_A = J + 0.01 * rng.normal(size=J.shape)
J_B = J + 0.01 * rng.normal(size=J.shape)
J_A, _ = la.qr(J_A, mode='economic')
J_B, _ = la.qr(J_B, mode='economic')
```

These lines live in the recovered notebook and are not provenance of the original `J_A/J_B` bundle. In the same context, `J` is itself either a random QR construction or a synthetic `B_T2U_TAN` proxy; both are not verified original sources. The bundle elements `J_A.npy`, `J_B.npy`, `K0.npy`, and `Kinv.npy` remain UNVERIFIED and therefore cannot be used as truth data.

## Verified upstream chain

Verified objects in local replay:

- `X_Oh` from `X_Oh_reconstructed.npy`: shape `(64, 3)`
- `Q` from `Q_rebuilt.npy`: shape `(192, 188)`
- `K_tangent` from `K_tangent_rebuilt.npy`: shape `(188, 188)`

These are the only locally verified inputs. No verified derivation of an original `J_A` or `J_B` was found that starts from these objects and ends at a unique, documented formula.

## Evidence from source files

The closest source text found is in the recovered notebook snippet:

- `tutti_gli_output/kaggle/working/188D_FORENSIC_REPLAY/notebook_source_recovered.py`
- `tutti_gli_output/kaggle/working/188D_RECOVERY/notebook_source_recovered.py`

The relevant block is:

```python
# 4. AUDIT UNIVERSALITÀ (ACCOPPIAMENTO DUE SPECIE DI MATERIA PHI_A, PHI_B)
# Specie A e Specie B accoppiate a blocchi differenti del settore 185D
J_A = J + 0.01 * rng.normal(size=J.shape)
J_B = J + 0.01 * rng.normal(size=J.shape)

J_A, _ = la.qr(J_A, mode='economic')
J_B, _ = la.qr(J_B, mode='economic')

M_A = J_A.T @ G0 @ J_A
M_B = J_B.T @ G0 @ J_B
```

This is a synthetic audit construction, not original physics. The same file also contains a fallback mock generator:

```python
J_A = np.random.randn(N_COLL, N_DIM)
J_B = np.random.randn(N_COLL, N_DIM)
```

Therefore the local corpus explicitly demonstrates that these matrices are not uniquely defined by the verified objects but are mock/synthetic artifacts.

## A/B semantics

The code comments support only a partial semantic role:

- `A` and `B` are two matter-species coupling channels (`phi_A`, `phi_B`) in the 185D internal sector.
- They are not two verified geometric subspaces or two verified basis vectors of a mathematically fixed original model.

This is therefore `PARTIALLY_RECOVERED` only; it does not justify a unique original `J_A/J_B` formula.

## Shape audit

The verified objects impose only the following dimensional logic:

- `X_Oh`: `(64, 3)`
- `Q`: `(192, 188)`
- `K_tangent`: `(188, 188)`

From the recovered code, the nominal `J` used in the synthetic audit is a triplet basis of shape `(188, 3)`, but this `J` is not proven to be the original source object. Because the provenance is unresolved, the shape `(188, 3)` cannot be promoted to an original definition.

## Final verdict

`J_A` / `J_B` original mathematical definition = `NOT_RECOVERABLE`.

No formula may be invented as if it were original. The correct classification is:

- Original provenance: `NOT_RECOVERABLE`
- Mathematical reproducibility: `NOT_RECOVERABLE`
- Numerical reproducibility: only for synthetic mock constructions, not original physics
- Historical bundle compatibility: `MISMATCH`

## Required downstream status

Because the original formula is not recoverable, all downstream objects depending on `J_A` / `J_B` must remain flagged as non-original until a verified external source is found.
