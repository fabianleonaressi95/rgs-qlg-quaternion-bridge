# Experiments archive and navigation

Snapshot: 2026-11-04, `main` tree `c4995f693a50a4b478cfb087efefb44401a1196e`.

This folder contains **311 tracked files** across topic directories. It is a research/workbench archive—not a single validated software package. The authoritative machine-readable file map, including Git blob SHA and byte size, is [EXPERIMENTS_CATALOG.json](EXPERIMENTS_CATALOG.json).

## Directory map

| Directory | Files | Topic / organization |
|---|---:|---|
| `analysis_sequence/` | 59 | Iterative `scriptanal*` research sequence; version lineage is not inferred from numeric suffix alone |
| `analysis_utilities/` | 11 | Matrix, topology, proof and audit utilities |
| `aperiodic_topology/` | 28 | Aperiodic topology exploratory variants (many `DockerfileN.py` historical filenames) |
| `astronomy/` | 21 | HST/STIS/COS spectral processing and MAST UV downloads |
| `ecoli_mutation/` | 2 | E. coli mutation data/build experiments |
| `engineering/` | 10 | Coil, thermal FEM and turbine iterations |
| `fea/` | 3 | Finite-element simulations |
| `legacy_imports/` | 3 | Imported legacy tests/summary |
| `legacy_proofs/` | 17 | Historical proof iterations |
| `microtubule_pipeline/` | 30 | MVP middleware and QIP/QC versions |
| `misc/` | 16 | One-off audits, recovery and external-data scripts |
| `peierls/` | 14 | Peierls model iterations |
| `prism_q6/` | 20 | PRISM/Q6, KS/Yang–Mills, 69120 and spectral experiments |
| `psu_bounce/` | 2 | PSU bounce/BVP experiments |
| `qproof/` | 2 | Q-proof matrix variants |
| `quantum_models/` | 26 | Quantum models, QNN, Koopman, QFIM and jellium experiments |
| `quasicrystal/` | 7 | Quasicrystal geometry and audit |
| `quaternion_kristall/` | 12 | Hurwitz/quaternion/Krist-All investigations |
| `qvsm_rnr/` | 11 | QVSM/RNR sources, reweighting and universality tests |
| `rgs_188d/` | 5 | 188D and O_h geometry experiments |
| `rgs_psu_core/` | 12 | RGS–PSU analytical/numerical bridge, Mosco, GH and spectral gap |

## Exact byte duplicates found

The following pairs have the same Git blob SHA in this snapshot:

1. `analysis_sequence/scriptanal48.py` and `scriptanal49.py`
2. `astronomy/hst_stis_cos/02_HST_STIS_COS_SPECTRUM_EXTRACTOR3.py` and `...EXTRACTOR4.py`
3. `microtubule_pipeline/mvp_middleware48.py` and `qip_qc_phase1_microtubule6.py`
4. `microtubule_pipeline/mvp_middleware64.py` and `mvp_middleware65.py`
5. `prism_q6/ks_Q3_genuine.py` and `ks_Q3_genuine_1.py`

**No duplicate file is deleted or renamed in this pass.** These entries are tagged in the JSON catalog for later reference/path checks. Byte-identical files may still be intentional aliases for workflows.

## Rules for handling experiment variants

- A suffix such as `1`, `2`, `fixed`, `corrected`, `final` or `v2` is not sufficient to establish chronology or correctness.
- Preserve all variants until their code diff, outputs, references, dependencies and scientific role are inspected.
- Do not promote an experiment into `src/` just because it is newer or passes a local run. Require a clean API, tests, provenance and reproducible validation.
- Keep raw/generated data beside the experiment or in `data/` with a manifest; do not silently overwrite data.
- For a new experiment, use a descriptive snake_case filename and a topic directory. Put short context, command, inputs, outputs and evidence status in its docstring/README.
- Treat physical interpretations as hypotheses unless derived and independently validated. A numerical audit passing establishes only the stated computational check.

## Next cleanup pass

The inventory is a first step. A safe deduplication/migration requires:
1. search all imports, shell commands, notebooks, workflow references and external entry points;
2. inspect diff/headers for numbered series;
3. identify canonical maintained implementation versus historical variants;
4. add compatibility wrappers or redirects if paths are referenced;
5. run relevant tests and full audit;
6. only then move or remove files in separate reviewable PRs.
