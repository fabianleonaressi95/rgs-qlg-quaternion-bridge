# RGS–PSU Vierbein reconstruction — status and recovery map

**Status:** archival reconstruction specification. This document does not claim that the original binary arrays have been recovered or regenerated, and it does not certify a TEGR identity.

## 1. Target artifacts

The historical work refers to these files:

- `RGS_PRISM_FINAL_vierbein_E.npy` — spatial RGS–PRISM vierbein, historically described as a collection of 64 local 3D frames.
- `RGS_PRISM_tetrad_4d.npy` — 4D tetrad candidate.

The filenames and historical numerical reports are known, but the original `.npy` payloads and the complete source matrices needed to reproduce them are not currently present in the recovered repository context. They must not be represented as available artifacts until retrieved and checksum-verified.

## 2. Historical construction data (reported, not re-run here)

The archived generation note reports:

- `ell2 = 7.035514834154695`
- `c0 = 43.96495130838031`
- principal-square-root prescription: `E_i = principal_sqrt(C_i / c0)`
- finite-difference steps: coarse `1e-5`, fine `5e-6`
- relative Jacobian step stability: `2.9658161229881756e-09`
- `J_ambient` shape `(576, 192)`
- `J_tangent` shape `(576, 188)`
- `T2u_response` shape `(576, 3)`
- tangent-frame orthogonality residual: `1.3121453838922012e-14`
- reported ranks: `rank(J_ambient)=192`, `rank(JQ)=188`

These are historical reported outputs, not fresh computations in this commit.

## 3. Historical geometric diagnostics

Previously reported diagnostics:

| Quantity | Historical value |
|---|---:|
| 3D determinant range | 0.793963 – 1.010983 |
| minimum singular value | 0.7692575 |
| condition number | 3.00576 |
| 4D metric-signature convention | `eta = diag(-1,+1,+1,+1)` |
| minimum eigenvalue of reported `C/c0` | 0.591757 |
| reported 4D condition number | 1.7337 |

The 4D eigenvalue/condition diagnostics alone do not establish a Lorentzian tetrad: the tetrad must reproduce a metric of the stated signature through `g = e^T eta e` (with index convention documented), and the determinant/orientation must be checked.

## 4. Required inputs for exact regeneration

### Spatial 3D frames

Recover the original 64 matrices `C_i` (or the exact upstream arrays from which they are built), their vertex ordering, normalization, and the precise matrix-square-root convention. Then regenerate each `E_i` and test:

- symmetry and positive definiteness of each `C_i/c0`;
- `E_i^T E_i = C_i/c0` (or the documented alternative convention);
- determinant, singular values, condition number;
- equivariance under the stated `O_h` action;
- sensitivity to the archived finite-difference steps.

### 4D tetrad

Recover the actual 4D metric/profile data and the coordinate and internal-index conventions. Construct a coframe `e^a_mu` satisfying

`g_mu_nu = eta_ab e^a_mu e^b_nu`, with `eta = diag(-1,1,1,1)`.

Record the tetrad gauge, orientation, determinant, inverse tetrad, and residual `||g - e^T eta e||`. Do not identify a fourth SO(3)-singlet mode with physical time without an independent Lorentzian construction.

### Teleparallel / TEGR audit

A tetrad alone is not a complete teleparallel specification. Fix the inertial spin connection `omega^a_{b mu}` consistently with the tetrad gauge. Then compute torsion, the torsion scalar, boundary term, and Levi-Civita Ricci scalar, and test the TEGR identity with explicit sign and index conventions. The Weitzenböck gauge `omega=0` is admissible only when the chosen tetrad/gauge makes that choice consistent.

## 5. Scope and verification labels

- **Recovered from historical notes:** filenames, construction prescription, constants, reported dimensions and residuals above.
- **Not recovered:** original `C_i` matrices, original `.npy` vierbein/tetrad payloads, complete 4D metric/profile inputs.
- **Not revalidated in this commit:** all historical numerical diagnostics.
- **Open:** exact 3D regeneration, Lorentzian 4D tetrad reconstruction, inertial spin-connection selection, and independent TEGR identity audit.

This status note is deliberately conservative: it preserves the research trail without converting an old report into a new computational PASS.
