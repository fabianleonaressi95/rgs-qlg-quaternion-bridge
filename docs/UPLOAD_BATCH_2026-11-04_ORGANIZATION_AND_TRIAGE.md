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


## Focused mathematical audit of the RGS–PSU imports

The following are code-level findings from reading the imported source, not
results of executing every script.

### Critical: equidistribution claim is false as stated

`RGS_PSU_Equidistribution_Proof.py` calls the finite 2I orbit (120 points on
S³, or 60 points after antipodal quotient) uniformly distributed with respect
to continuous Haar measure. A finite orbit measure is atomic and cannot equal
or converge to continuous Haar measure at a single fixed level. A finite
subgroup action being transitive on its own orbit, or preserving the set under
symmetry, does not imply equidistribution on the whole sphere. The script's
purported Weyl theorem is therefore misapplied. The 600-cell is highly
symmetric and is a spherical design of finite strength, but that is not the
same as continuous equidistribution. A valid route is to define a sequence of
refinements and prove convergence of empirical measures against a determining
class of test functions, with quantitative discrepancy bounds if possible.

### GH and measure convergence remain conditional

`RGS_PSU_GH_Measure_Convergence.py` explicitly says equidistribution is
plausible and requires proof. Its later proof sketch treats mesh decay and
covering of the target as consequences of a refinement construction, but the
code does not establish those facts for every point of the target at all
levels. The finite level-0 set of 60 quotient vertices is not a GH
approximation of the full continuous Poincaré sphere without a controlled
refinement sequence and metric embeddings.

### Critical normalization inconsistency in the Mosco script

`RGS_PSU_Mosco_Renormalization.py` defines in one section
(L_n^{ren}=L_n/N_n), while another section invokes an (N_n^{-2/3})
renormalization for a three-dimensional Weyl scaling. These are different
normalizations and cannot both define the same operator without an explicit
change of graph metric, edge weights, and measure convention. In addition,
checking a single level-0 eigenvalue cannot establish a uniform Poincaré
bound or spectral/resolvent convergence over the refinement sequence.
Mosco liminf and recovery-sequence conditions still need independent
construction and verification.

### Fluid–gauge bridge is currently formal/phenomenological

`RGS_PSU_Analytical_Derivations.py`,
`RGS_PSU_Shear_Deformation_Bridge.py` and
`RGS_PSU_Multifractal_Shear_Bridge_Integrated.py` state the desired
Navier–Stokes/Wilson-stress correspondence, but a common shear parameter
alone does not derive a gauge-to-fluid equivalence. The gauge action's
dependence on the deformed coframe/metric, the fluid configuration space,
measure, dissipative/forcing terms, and the Green–Kubo response must be
specified from one action or a controlled statistical-mechanical limit.
A computed cross-correlator is evidence for a chosen model, not proof of
dynamical equivalence.

### Teleparallel closure records unresolved contradictions

`RGS_PSU_closure.py` is valuable because its header explicitly records
inconsistent tetrad, NLED field/energy-density, and Misner–Sharp mass
formulas. Its concluding printout nevertheless announces a teleparallel
identity verified at machine precision, while the header says the corrected
identity leaves a nonzero (-2R''/R) term for generic profiles. This must be
resolved by deriving the exact convention/signature, tetrad, inertial spin
connection, torsion scalar, and boundary divergence from the same geometry;
a small residual for a special profile is not a general identity proof.

### Exact duplicate blobs

The relocated tree contains 12 groups of byte-identical NumPy blobs. Examples:
the three `qnn_6_6_6_tensor_paths` variants are identical; multiple
`X_Oh_reconstructed` and `64B_Oh_reconstructed` variants share one blob;
and `B_T2U_TAN_rebuilt.npy` and `J0_rebuilt.npy` are byte-identical despite
different names. These were retained, not deleted. The differing semantic
names need provenance checks before deciding whether they are aliases,
mislabelled copies, or intentionally duplicated artifacts.

## Validation status

- The move operation reuses original Git blob SHAs; file bytes are preserved.
- Root-level Python and NumPy files from this upload batch have been moved into
  topic folders.
- The current change does not execute all imported research scripts. Existing
  repository CI checks the maintained package/tests, not a scientific
  validation of every legacy upload.
- Next audit should produce a machine-readable manifest for each `.npy`
  (shape, dtype, byte count, SHA-256, semantic role, generating script) and
  add syntax/runtime smoke tests for imported Python files in isolated groups.
