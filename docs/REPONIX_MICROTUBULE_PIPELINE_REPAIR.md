# Reponix follow-up: microtubule pipeline hardening

Date: 2026-11-05  
Scope: `experiments/microtubule_pipeline` numbered analysis chain.

## Changes in this patch

- The phase 1E weighted pipelines (`mvp_middleware49.py` and
  `mvp_middleware50.py`) now use a shared downloader that streams to a
  temporary sibling file, applies a finite timeout, validates gzip streams,
  optionally verifies a caller-pinned SHA-256, and publishes with an atomic
  rename only after validation.
- Their output directory, and the output directory in middleware phases 1F–11
  (`mvp_middleware51.py`–`mvp_middleware69.py`), is now configurable through
  `RGS_QIP_OUTDIR`. The default is
  `~/Scaricati/QIP_QC_EMD25156`, not a machine-specific absolute home path.
- Regression tests cover atomic publication, checksum mismatch preservation,
  gzip validation, finite timeout handling, and the coordinate-function call
  signatures.

## Findings reviewed rather than blindly changed

- The reported `density_coordinates` signature mismatch is not reproduced in
  the inspected phase-1 variants: 47/48 use the two-argument form and their
  callers pass two arguments; 49/50 use the weighted three-argument form and
  their callers pass three. The AST regression test locks this contract.
- `matplotlib.use("Agg")` occurs before importing `matplotlib.pyplot` in the
  inspected scripts. That ordering is valid for headless CI and is not itself
  a defect, so it was not changed.
- Exact duplicate experiment paths are retained. The repository's experiment
  archive policy explicitly requires reference and entrypoint analysis before
  deleting or renaming historical paths.

## Integrity boundary and remaining work

The downloader computes and reports SHA-256 for every transfer, and verifies
it when an independently obtained expected digest is supplied through
`RGS_QIP_MAP_SHA256` or `RGS_QIP_PDB_SHA256`. No trusted digest was present
in the repository for these upstream files, so this patch does **not** claim
cryptographic source verification when those variables are unset. In that
case the code reports the limitation; gzip CRC/trailer validation still runs
for the compressed map.

This is a focused repair of the microtubule pipeline family, not a claim that
all static-analysis findings across the full repository have been resolved.
Other downloaders and experiment families need separate provenance-aware
patches.
