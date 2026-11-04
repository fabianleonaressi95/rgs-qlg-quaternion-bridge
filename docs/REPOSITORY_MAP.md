# Repository map and maintenance rules

Snapshot: 2026-11-04, default branch `main`, tree `68c4e671ca19f6aad12d756b16f8373bbc9644b0`.

## Maintained code

| Path | Contents | Policy |
|---|---|---|
| `src/` | 70 importable Python modules: quaternion/SU(2), 2I, Q6, Clifford/Dirac, Wilson/overlap, RGS–PSU, Einstein–Cartan, vierbein and spin connections | Production/reusable implementations; avoid parallel duplicate modules |
| `tests/` | 71 pytest modules | Canonical regression and structural tests |
| `examples/` | 30 runnable audit/demo scripts | Reproducible entry points; document command and output |
| `scripts/` | 20 Python tools plus one NPZ artifact | Maintained automation and audit utilities |
| `.github/workflows/` | CI workflows | Required validation gates |

## Research and evidence

| Path | Contents | Policy |
|---|---|---|
| `docs/` | 91 Markdown documents, including architecture, derivations, conventions, roadmap and imported notes | Source of truth for scientific scope and evidence labels |
| `experiments/` | 311 files, predominantly exploratory Python variants | Preserve historical work; do not import into `src/` without review |
| `data/` | 190 files, primarily CSV/NPY/NPZ numerical artifacts | Preserve bytes; require provenance, shape, dtype, units and generator metadata |
| `reports/` | 27 files, mainly logs and imported reports | Results and audit outputs; distinguish source reports from regenerated outputs |
| `results/` | 7 JSON result/catalogue files | Machine-readable outputs |
| `papers/` | Draft manuscripts | Draft status must remain explicit |
| `references/` | Bibliographic/imported references | Keep source and licence metadata |
| `archive/` | Bundles and transcripts | Immutable historical material; no scientific status implied |

The recursive Git tree contains 986 entries (files and directories). There are 77 top-level entries, including 70 root-level files. The root-level legacy artifact inventory is documented in [ROOT_ARTIFACT_TRIAGE.md](ROOT_ARTIFACT_TRIAGE.md).

## Scientific architecture map

1. **Algebra and symmetry:** `quaternions.py`, `rotations.py`, `su2_lifts.py`, `binary_icosahedral.py`, `representations_2i.py`.
2. **Discrete geometry and topology:** `q6_complex.py`, `mobius_complex.py`, `discrete_differential.py`, `q6_geometric_hodge.py`, `q6_whitney_hodge.py`.
3. **Spinors and gauge fields:** `clifford_qlg.py`, `discrete_dirac.py`, `holonomies.py`, `rgs_q6_su2_wilson.py`, `su3_flavour.py`.
4. **RGS–PSU geometry and refinement:** `rgs_q6_oh_action.py`, `rgs_psu_cut_project.py`, `rgs_psu_intertwiner.py`, `rgs_psu_fibred_intertwiner.py`, `rgs_psu_multifractal_bundle.py`, `rgs_ellis_refinement.py`.
5. **Gravity and matter coupling:** `rgs_spin_connection.py`, `rgs_ellis_spin_connection.py`, `rgs_ec_*`, `rgs_psu_hilbert_qhd_mhd.py`, `spinor_kirchhoff.py`, `rgs_psu_qball_*`.
6. **Multiscale and statistical mechanics:** `rgs_multiscale.py`, `rgs_unified_multiscale_operator.py`, `koopman_rgs_psu.py`, `spin_glass` and thermodynamics modules.

This is a navigation map, not a claim that all sectors are derived from one action or physically unified.

## Required evidence labels

- **PROVED:** proof with explicit assumptions and conventions.
- **NUMERICALLY VERIFIED:** reproducible calculation, tolerances and independent checks identified.
- **CANDIDATE:** specified construction with outstanding verification gates.
- **OPEN:** unresolved mathematical, numerical or physical question.
- **UNVERIFIED / IMPORTED:** provenance or validation is incomplete.

Do not promote a result based only on a passing unit test. In particular, finite graph identities do not establish continuum convergence, and compatibility adapters do not constitute a unified field action.

## File placement rules

- New reusable library code: `src/`; one canonical implementation per mathematical object.
- Regression tests: `tests/`; names must track the module or invariant under test.
- Reproducible runnable audits: `examples/`.
- Maintained command-line and batch utilities: `scripts/`.
- Numerical source inputs: `data/<domain>/`; never silently overwrite a source artifact.
- Generated reports: `reports/<domain>/`; machine-readable result payloads may go in `results/`.
- Exploratory alternatives: `experiments/<domain>/`; preserve variant suffixes until a semantic and byte-level audit is complete.
- Imported bundles and transcripts: `archive/`.

## Safe cleanup protocol

1. Inventory paths, blob SHA, size, file type and references.
2. Classify each file as source, test, example, experiment, data, report, reference or archive.
3. Detect byte-identical blobs separately from semantically equivalent files.
4. Check code, workflow and Markdown references before moving a path.
5. Move by reusing the original Git blob SHA; do not rewrite binary artifacts.
6. Keep a path-migration map and provenance note.
7. Run compile, full pytest, examples and relevant domain audits before merge.
8. Do not delete or merge variants solely because their names look duplicated.

## CI baseline at this snapshot

The latest `python-package.yml` run for this `main` tree completed successfully (run 402, 2026-11-04). This confirms the configured software checks only; it does not execute or certify every file under `experiments/`, every imported dataset, or every scientific claim.
