# Upload batch organization and triage — 2026-11-04

This change relocates the newly uploaded root-level scripts and NumPy artifacts
without changing their bytes. Git blob IDs are reused, including for binary
.npy files; duplicate suffixes such as "(1)" are retained rather than guessed
to be identical or overwritten.

## Directory map

- `experiments/astronomy/hst_stis_cos/`: HST/STIS/COS spectrum and extractor variants.
- `experiments/astronomy/mast_uv/`: MAST UV retrieval utilities.
- `experiments/rgs_psu_core/`: analytical derivations, numerical checks, spectral-gap, GH, Mosco and multifractal scripts.
- `experiments/rgs_188d/` and `experiments/rgs_188d/oh64_analysis/`: 188D / O_h / 64B analysis.
- `experiments/qvsm_rnr/`: QVSM/RNR reduction, bounce and universality variants.
- `experiments/quaternion_kristall/`: Hopf, Hurwitz, quaternion and Krist-All exploratory scripts.
- `experiments/qproof/`: Q-Proof versions.
- `experiments/psu_bounce/`: PSU bounce and physical-test scripts.
- `experiments/engineering/`: turbine, coil and thermal FEM scripts.
- `experiments/ecoli_mutation/`: E. coli mutation HDF5 builders.
- `experiments/legacy_proofs/`: numbered proof scripts; retained as historical variants, not certified proofs.
- `experiments/analysis_utilities/`: miscellaneous analysis helpers.
- `experiments/legacy_imports/`: imported summary/test scripts requiring independent validation.
- `archive/transcripts/`: transcript text-processing scripts.
- `data/rgs_188d/`, `data/q6/`, `data/qnn_controls/`, `data/quaternion_geometry/`: numerical artifacts grouped by subject.

## Initial technical triage

1. Filename families with numbered suffixes are version/import batches, not
   established sequential improvements. Do not assume the largest suffix is best.
2. HST/STIS/COS extractors need common input provenance, wavelength/flux/error
   schema, calibration metadata, masks, and independent reference spectra.
3. RGS/PSU and 188D arrays require a manifest with shape, dtype, units,
   coordinate convention, generating script, and SHA-256; similarly named
   files are not interchangeable.
4. QVSM/RNR and PSU bounce scripts require dimensional analysis, explicit
   equations, parameter provenance, solver tolerances and convergence tests
   before physical interpretation.
5. Quaternion/Hopf/Hurwitz/Krist-All scripts need exact group/action definitions,
   orbit-stabilizer checks, normalization conventions, and independent algebraic
   identities. Numerical agreement alone is not a proof.
6. Turbine/FEM and coil scripts need mesh/time-step convergence, boundary
   conditions, material parameters and conservation checks.
7. The legacy proof series and Q-Proof files are exploratory imports until each
   theorem statement, hypotheses and proof/test coverage are reviewed.
8. No numerical artifact was deleted, renamed internally, or deduplicated.

This is an organizational pass and a triage map, not a claim that every imported
script has been fully executed or scientifically validated.
