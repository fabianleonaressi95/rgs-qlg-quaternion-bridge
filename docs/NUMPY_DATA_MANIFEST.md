# NumPy data manifest

Run from the repository root:

```bash
python tools/build_npy_manifest.py
```

The command scans `data/**/*.npy` and writes
`data/manifest/npy_manifest.json`. It records the actual NPY header (format
version, shape, dtype, Fortran order), byte counts, SHA-256, Git blob SHA-1,
payload-size consistency, exact-content duplicate groups, and Python files
that mention the same basename.

## Interpretation limits

- A matching SHA-256/Git blob ID establishes byte identity, not semantic
  equivalence. Keep differently named duplicates until provenance is resolved.
- `producer_candidates` are source files that mention a basename. They are
  leads, not verified generators.
- `units`, `coordinate_convention`, and `semantic_role` are deliberately
  left unresolved. Do not infer them from filenames or matrix dimensions.
- Header/payload checks validate container consistency only; they do not
  validate an array's mathematical or physical interpretation.
- Object-dtype arrays are identified, but payload size is not inferred because
  their data are pickle-encoded.

The manifest is generated from the checked-out files, so its SHA-256 values
refer to the exact bytes present at generation time. Re-run after any data
change. The generator requires NumPy.
