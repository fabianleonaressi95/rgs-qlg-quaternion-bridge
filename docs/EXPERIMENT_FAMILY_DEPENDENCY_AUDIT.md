# Experiment family dependency and consolidation audit

Date: 2026-11-04  
Scope: targeted static audit of selected numbered experiment families after the 311-file catalog. This is **not** a complete execution or scientific validation of all 311 scripts.

## Decision

**Do not bulk-rename, move, delete, or auto-select the highest-numbered script.** Numeric suffixes are historical names, not validated semantic version numbers. Preserve original paths until provenance and downstream usage are mapped. The inventory's five identical-file pairs are candidates for aliasing, not immediate deletion.

## Findings

### 1. Exact duplicate family members

The catalog identified byte-identical pairs (same Git blob SHA):

- `analysis_sequence/scriptanal48.py` = `scriptanal49.py`
- `astronomy/hst_stis_cos/02_HST_STIS_COS_SPECTRUM_EXTRACTOR3.py` = `...EXTRACTOR4.py`
- `microtubule_pipeline/mvp_middleware48.py` = `qip_qc_phase1_microtubule6.py`
- `microtubule_pipeline/mvp_middleware64.py` = `mvp_middleware65.py`
- `prism_q6/ks_Q3_genuine.py` = `ks_Q3_genuine_1.py`

Interpretation: exact content equivalence is established at the byte level; scientific equivalence, path usage, and deletion safety are not. Keep both path names until the reference audit and provenance mapping are complete.

### 2. `analysis_sequence/scriptanal48.py` / `scriptanal49.py`

- Confirmed exact byte duplicate.
- These are standalone separation/non-normality audits; not a shared package/module.
- Recommended action: after reference scan, select one canonical implementation and retain the other as a compatibility shim or historical alias. Do not remove either in this pass.

### 3. `aperiodic_topology/Dockerfile*.py` naming family

- These are Python source files, not Docker build recipes (e.g. the base `Dockerfile.py` imports NumPy/SciPy and defines an Aubry–André chain and modal entropy).
- The family contains gaps in suffixes and large jumps (e.g. `Dockerfile57.py`), so the numbering cannot be interpreted as a complete linear release history.
- `Dockerfile57.py` identifies itself in its docstring as “Dockerfile56.py” and describes a Q-PROOF L4.6 nested Monte Carlo experiment. This filename/docstring mismatch is a provenance warning, not proof of a bug.
- Recommended action: eventually migrate under descriptive names (e.g. `qproof_l4_6_nested_mc.py`) with old-path wrappers. Do this only after checking path references and capturing source SHA/provenance in a manifest.

### 4. `peierls/peierls45.py`

Static inspection shows this script invokes names that are not defined or imported in the file:

- `generate_amorphous_lattice`
- `build_hamiltonian`
- `find_local_optimal_gap`

The script calls `scan_magnetic_phase_diagram()` under `__main__`, so direct execution will fail with `NameError` unless those functions are injected externally (not a normal Python execution model). The comments say to include previously validated functions, but those functions are absent from this file.

Classification: **execution-incomplete / external dependency unresolved**. This is not evidence that the physical model or result is invalid; it means the artifact is not self-contained as committed. Recommended action: recover the exact referenced implementation from history/adjacent files, add explicit imports or a local module, then add a smoke test before calling it runnable.

### 5. `microtubule_pipeline/mvp_middleware47.py`

- This is a substantial executable pipeline for EMDB-25156 / PDB 7SJ7, with download paths, generated output directories, compressed arrays and plots.
- It performs filesystem writes and external data retrieval as part of its workflow; it should be treated as a data pipeline, not a pure numerical function.
- Recommended action: isolate I/O in a CLI entrypoint, support an offline/cache mode, record source data checksums and versions, and add a small synthetic-data smoke test. Do not execute the full pipeline in CI by default.

### 6. `microtubule_pipeline/mvp_middleware77.py`

- This is a distinct screw-symmetry/RMT experiment for a 13-protofilament, 64-ring lattice (832 sites), not simply a later rewrite of file 47.
- It should be classified by scientific objective, not numeric suffix.
- Recommended action: add metadata (model, data provenance, assumptions, expected invariants, runtime) and a small deterministic test for the screw operator/commutator before any relocation.

### 7. `peierls/peierls.py` vs `peierls45.py` / `peierls50.py`

- The base file builds a Peierls-phase tight-binding Hamiltonian from supplied coordinates.
- `peierls45.py` and `peierls50.py` are separate higher-level Chern/phase-diagram workflows, not established drop-in successors.
- `peierls45.py` has unresolved missing helper functions (above).
- `peierls50.py` includes its own amorphous-lattice generator and Delaunay-based workflow; do not merge it into the base builder without a shared API contract and numerical tests.

## Reference scan status

GitHub code search did not surface direct import/call sites for the targeted numbered filenames in the indexed repository. This is **not** a proof that references are absent: search indexing may omit files, notebooks, shell commands, workflow strings, external scripts, or historical users.

Before source relocation, perform a complete local checkout scan over:
- Python imports and literal path strings;
- notebooks, shell scripts, YAML workflows, docs, and JSON manifests;
- Git history and external user-facing paths;
- generated data/output dependencies.

## Suggested next safe sequence

1. Add a per-file manifest schema: `path`, `blob_sha`, `topic`, `entrypoint`, `imports`, `external_data`, `writes_files`, `status`, `provenance`, `scientific_claim_level`.
2. Build a dependency/reference graph from a local checkout (AST + text scan), not only GitHub search.
3. Recover missing Peierls helpers and add a smoke test.
4. Add reproducibility wrappers and metadata to high-value experiments.
5. Only then create focused migration PRs with compatibility shims and tests.

## Evidence labels

- **Confirmed:** exact blob equality; observed imports/definitions/call sites in inspected source; filename/docstring mismatch.
- **Static concern:** missing helper functions in standalone file; side-effecting pipeline; ambiguous numbering.
- **Not assessed:** scientific validity, full runtime, data provenance for all experiments, whether all external references exist.
