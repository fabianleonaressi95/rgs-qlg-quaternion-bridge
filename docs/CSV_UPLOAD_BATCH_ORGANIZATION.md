# CSV upload organization

Files from upload commit `9ca584b13b439e256246f3748b68aed5fca677e0` have been relocated from repository root into thematic `data/` directories. Original CSV contents are preserved byte-for-byte; names containing spaces/parentheses are normalized for portability. The upload commit remains in Git history.

Directories:
- `astronomy_forecasting/`: CR3BP, conjunctions and forecasts
- `prism_fibonacci/`: Q6/PRISM/Fibonacci/69120 audits
- `krasilnikov/`: critical fits, residuals and theta scans
- `qnn_controls/`: QNN, entanglement, equivariance, null and control datasets
- `morphodynamics_information/`: morphodynamic, HOM, Fisher and descriptor data
- `model_scans/`: parameter scans and D2/D3 variants
- `koopman_edmd/`: Koopman/EDMD data
- `imported_misc/`: not safely classifiable from filename alone

This is file organization only, not scientific validation. CSV schemas and provenance should be documented separately.
