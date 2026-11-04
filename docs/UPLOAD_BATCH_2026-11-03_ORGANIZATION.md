# Uploaded research batch — organized

This directory layout groups the files added in commit `281ad75995baafd32e2a601df7902e05702bc419`.

## Important provenance
- This commit is a **file relocation / repository hygiene pass**.
- Original contents are preserved byte-for-byte; filenames may be normalized only for imported bibliography/report paths where spaces/parentheses existed.
- Relocation is **not scientific validation** and does not imply the scripts run or the claims are established.
- Original upload commit: [281ad75995baafd32e2a601df7902e05702bc419](https://github.com/fabianleonaressi95/rgs-qlg-quaternion-bridge/commit/281ad75995baafd32e2a601df7902e05702bc419).

## Layout
- `experiments/aperiodic_topology/`: Dockerfile-named Python Aubry–André/topology experiments.
- `experiments/microtubule_pipeline/`: sequential MVP middleware scripts.
- `experiments/fea/`: finite-element/impact prototypes.
- `experiments/quantum_models/`: QNN, Koopman, jellium and related model scripts.
- `experiments/quasicrystal/`: QC/truss experiments.
- `scripts/`: reusable audit/recovery scripts.
- `reports/`: logs and numerical reports.
- `papers/drafts/`: Typst drafts.
- `references/imported/`: imported bibliography files.
- `docs/imported/`: imported markdown proposals/audits.
- `archive/imported/`: miscellaneous notes and files not confidently classified.

## Next steps
1. Fix syntax/runtime errors and add tests after relocation.
2. Convert hard-coded paths to configurable paths.
3. Mark each result as reproduced / exploratory / unverified.
4. Consolidate duplicate or successive script versions only after comparing them.
