# Full static experiment dependency-graph audit

Date: 2026-11-04  
Source: GitHub Actions run [41836984040](https://github.com/fabianleonaressi95/rgs-qlg-quaternion-bridge/actions/runs/41836984040), artifact `experiment-dependency-graph`  
Commit scanned: `a23b6125f4c704c1f76427ff5daf1c377dacbe79` (PR #108 head; subsequently merged)

## Scope and counts

The workflow scanned repository text files and built a static graph. It did not execute experiments.

| Measure | Count |
|---|---:|
| Text files scanned | 711 |
| Python files scanned | 510 |
| Catalogued experiment files | 311 |
| Python experiment files | 307 |
| Binary `.npz` experiment data files excluded from AST scan | 4 |
| Dependency/reference edges emitted | 1,176 |
| Python import edges (whole repository) | 265 |
| Literal path references | 622 |
| Literal basename references | 289 |
| Python parse errors in experiment sources | 15 |
| Exact duplicate groups across scanned text files | 6 |

The graph's `experiment_files=309` counter is the number of text files under `experiments/` (including the README and catalog), not the 311-item scientific-artifact inventory. The inventory has 307 Python files and four binary NPZ datasets. The NPZ files are catalogued but intentionally excluded from text/AST parsing.

## Key structural result: weakly connected experiment scripts

Only **one direct Python import edge between experiment source files** was resolved:

`experiments/peierls/peierls45.py` → `experiments/peierls/peierls46.py`

This is consistent with the explicit helper recovery in PR #108. The graph therefore does **not** support treating the numbered experiment folders as a coherent imported software package. Most files appear to be standalone scripts or notebooks/CLI artifacts linked by textual references and shared concepts rather than Python module imports.

This is not proof that experiments are independent: dynamic imports, runtime path construction, subprocesses, notebook execution, data-file loading and external scripts are not resolved by this static scan.

## Important interpretation warning: textual references are not runtime dependencies

The 1,176 edges must not be read as 1,176 executable dependencies.

- The experiment catalog alone accounts for 308 literal references because it lists artifact paths.
- Documentation, README files, workflow YAML and reports also contain filename mentions.
- Basename matches can be false positives.
- The useful executable-import signal is the AST import subset, not the raw literal-reference total.

A future graph revision should keep separate edge classes for Python imports, explicit data loads, subprocess/CLI calls, documentation references and catalog membership. It should also label uncertain basename matches rather than treating them as dependency edges.

## Python parse blockers

The AST pass found 15 experiment Python files that do not parse under the workflow's Python 3.12 parser:

| File | Parser diagnostic |
|---|---|
| `experiments/rgs_psu_core/RGS_PSU_Equidistribution_Proof.py` | Unicode subscript character U+2083 in source at line 74 |
| `experiments/rgs_psu_core/RGS_PSU_closure.py` | unexpected indent, line 874 |
| `experiments/misc/rgs_supergolden_a6.py` | invalid decimal literal, line 318 |
| `experiments/analysis_sequence/scriptanal39.py` | mismatched closing bracket, line 61 |
| `experiments/analysis_sequence/scriptanal40.py` | mismatched closing bracket, line 61 |
| `experiments/analysis_sequence/scriptanal42.py` | mismatched closing bracket, line 82 |
| `experiments/analysis_sequence/scriptanal44.py` | mismatched closing bracket, line 70 |
| `experiments/aperiodic_topology/Dockerfile19.py` | invalid syntax, line 105 |
| `experiments/peierls/peierls4.py` | invalid syntax, line 88 |
| `experiments/quantum_models/orchestrate.py` | unexpected indent, line 12 |
| `experiments/qvsm_rnr/qvsm_universality_audit.py` | unclosed parenthesis, line 841 |
| `experiments/qvsm_rnr/build_qvsm_rnr_neo_reduction.py` | mismatched bracket/brace, lines 386–397 |
| `experiments/legacy_proofs/proof2.py` | invalid syntax, line 198 |
| `experiments/microtubule_pipeline/qip_qc_phase1_microtubule4.py` | unexpected character after line continuation, line 341 |
| `experiments/astronomy/hst_stis_cos/02_HST_STIS_COS_SPECTRUM_EXTRACTOR.py` | invalid syntax, line 102 |

These are **syntax/encoding parse blockers**, not a verdict on the mathematical content. The Unicode case may be a prose/formula-in-code issue rather than an intended executable module. Each file needs manual inspection before repair: preserve the mathematical expression, identify whether the file is meant to be executable Python, and add a minimal parse/import test only after its role is clear.

### Priority repair order

1. **Core RGS–PSU:** `RGS_PSU_closure.py` and `RGS_PSU_Equidistribution_Proof.py`. These names imply central closure/proof artifacts; recover their intended derivations and make executable code distinct from mathematical prose.
2. **QvSM/RNR:** `qvsm_universality_audit.py` and `build_qvsm_rnr_neo_reduction.py`. Repair only after checking whether the malformed sections are incomplete experiments or accidental truncations.
3. **Supergolden / Peierls:** `rgs_supergolden_a6.py` and `peierls4.py`. Compare with adjacent versions and Git history; do not copy the newest numbered file blindly.
4. **Sequence family:** `scriptanal39.py`, `scriptanal40.py`, `scriptanal42.py`, `scriptanal44.py`. Determine whether each is an abandoned draft or a required result before patching.
5. **Other isolated files:** the astronomy extractor, Dockerfile19.py, orchestrate.py, proof2.py and qip_qc_phase1_microtubule4.py. Establish intended entrypoint and dependencies first.

## Exact duplicate groups

Six byte-identical groups were detected across scanned text files. Five are experiment pairs already listed in the catalog; one is a duplicated report artifact.

- `experiments/analysis_sequence/scriptanal48.py` = `scriptanal49.py`
- `experiments/astronomy/hst_stis_cos/02_HST_STIS_COS_SPECTRUM_EXTRACTOR3.py` = `...EXTRACTOR4.py`
- `experiments/microtubule_pipeline/mvp_middleware48.py` = `qip_qc_phase1_microtubule6.py`
- `experiments/microtubule_pipeline/mvp_middleware64.py` = `mvp_middleware65.py`
- `experiments/prism_q6/ks_Q3_genuine.py` = `ks_Q3_genuine_1.py`
- `reports/rgs_188d/testB_Oh_constrained_report (1).json` = `reports/rgs_188d/testB_Oh_constrained_report.json`

Exact byte equality establishes identical contents only. It does not establish which path is canonical, whether external consumers use both paths, or whether deletion is safe. No files are deleted in this audit.

## Architecture decision

Do not flatten all experiments into one package and do not promote a numbered variant to canonical solely by suffix. The repository currently contains a large **research archive plus a smaller set of reusable modules**. Separate these roles explicitly:

- **Core library:** stable mathematical operators, geometry, gauge/Dirac structures and shared numerical primitives. These belong under `src/` with versioned APIs and tests.
- **Reproducible experiment drivers:** scripts with explicit inputs, seeds, outputs, provenance and command-line entrypoints. Keep them under `experiments/`.
- **Data pipelines:** scripts with downloads, filesystem side effects or large external datasets. Give them cache/offline modes and checksums; never run them as default unit tests.
- **Historical/proof artifacts:** preserve original source and derivations; mark non-executable or incomplete material clearly rather than forcing it into imports.
- **Generated reports/data:** keep under `reports/` and `results/`, with source commit and provenance manifest.

## Next implementation sequence

1. Inspect the 15 parse blockers individually, starting with the two RGS–PSU core files and two QvSM/RNR files.
2. Extend the scanner to distinguish executable imports, data reads, subprocess calls, docs references and catalog membership; reduce false-positive basename edges.
3. Add a manifest per experiment: topic, status, entrypoint, imports, input data, outputs/side effects, seed, provenance SHA, tests and evidence class (PROVED / NUMERICALLY VERIFIED / CANDIDATE / OPEN).
4. Promote only genuinely shared and tested functions into `src/`; keep old experiment paths as compatibility wrappers when needed.
5. Do not delete exact duplicates until all in-repository and external path consumers have been considered.

## Evidence boundary

- **Confirmed:** file counts, parser diagnostics, exact SHA duplicate groups and one experiment-to-experiment AST import edge in this snapshot.
- **Static indication:** scripts are mostly loosely coupled by imports; documentation and catalog references dominate literal edges.
- **Not established:** runtime dependency completeness, scientific validity, numerical reproducibility, external consumer references or whether any malformed file is obsolete.
