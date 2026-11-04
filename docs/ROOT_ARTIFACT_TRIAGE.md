# Root artifact triage

Initial inventory snapshot: `main` tree `68c4e671ca19f6aad12d756b16f8373bbc9644b0` (2026-11-04).
Relocation branch: `chore/root-json-report-relocation-2026-11-04`.

## Applied relocation

The root initially contained 70 files, including legacy JSON reports. **57 JSON artifacts (approximately 3.42 MB) have now been relocated** to domain-specific folders under `reports/` and `data/`, preserving original Git blob SHAs and contents. The relocation is staged in PR #105 and is not yet merged.

Current branch root files:
- `README.md`
- `CONTRIBUTING.md`
- `requirements.txt`
- `license`
- `PROVENANCE_STATUS.json`
- `RGS_PSU_FILE_INVENTORY.json`

See [reports index](../reports/INDEX.md) and [repository map](REPOSITORY_MAP.md).

## Destination families

| Family | Destination |
|---|---|
| `188D_*`, `J_*`, `testB_Oh_*`, PSU basis/spinorial lift | `reports/rgs_188d/` |
| `64B_*` | `reports/rgs_188d/oh64/` |
| `K0_*`, Hodge/Killing audits | `reports/spectral/` |
| DNS conservation variants | `reports/dns/` |
| RGS–PSU master/unified/Green–Kubo/Q6 reports | `reports/rgs_psu/` |
| SU2 beta/L8 | `reports/q6_su2/` |
| generic/imported audit reports | `reports/imported_root/` |
| manifests, model/map parameters | `data/manifests/` |
| OMM station metadata | `data/astronomy_forecasting/` |
| barycentric gradient | `reports/discrete_geometry/` |

## Duplicate and variant policy

- `testB_Oh_constrained_report.json` and `testB_Oh_constrained_report (1).json` share blob SHA `6fce2715e2f99c1553b57d203cf02982a30c0189`. Both paths are retained.
- Green–Kubo report variants have different blob SHAs and sizes; they remain separate.
- `manifest.json`, `manifest (1).json`, and `manifest (2).json` have different blob SHAs; no precedence is inferred from suffixes.
- No semantic deduplication was performed.

## Provenance warning

`PROVENANCE_STATUS.json` marks several 188D arrays (including `K0.npy`, `Kinv.npy`, `J_A.npy`, `J_B.npy`, `J.npy` and `P_G.npy`) as UNVERIFIED. Presence in the repository or a passing test does not certify provenance or physical validity.

## Validation still required before merge

- Run and inspect Python tests and Full Audit on the relocated tree.
- Search for scripts/notebooks with hard-coded former root paths; GitHub code search may not index all branches or files.
- Repair path references only when destination is unambiguous.
- Keep imported reports labeled as imported/generated outputs, not established physical results.
