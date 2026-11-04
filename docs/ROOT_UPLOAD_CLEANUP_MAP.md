# Root cleanup — organization map

This cleanup moves the remaining loose research files from repository root into topic folders. The following root-level project files remain intentionally: `README.md`, `requirements.txt`, `CONTRIBUTING.md`, and `license`.

## Folder map
- `experiments/analysis_sequence/`: scriptanal numbered series
- `experiments/peierls/`: Peierls model scripts
- `experiments/microtubule_pipeline/`: MVP middleware and QIP microtubule experiments
- `experiments/prism_q6/`: PRISM, Q6, 69120, PSU, Yang–Mills and related audits
- `experiments/quantum_models/`: QNN, QFI, quantum pipeline, simulator and related models
- `experiments/quasicrystal/`: quasicrystal geometry/audits
- `experiments/misc/`: other Python research prototypes
- `scripts/audit/`: general audit, scan and solver scripts
- `data/rgs_188d/`: 188D and O_h reconstruction arrays/provenance
- `data/raw_artifacts/`: other NPZ/NPY artifacts
- `papers/drafts/`: TeX manuscript and computational appendix files
- `archive/bundles/`: uploaded ZIP archives
- `reports/imported/` and `docs/imported/`: loose reports/documentation

**Scope:** organization only. Contents are not rewritten, duplicate data are not removed, and relocation does not establish correctness. Some scripts may use root-relative paths; check/update those references in a separate compatibility pass.
