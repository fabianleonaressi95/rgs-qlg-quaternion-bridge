# Root artifact triage

Snapshot: `main` tree `68c4e671ca19f6aad12d756b16f8373bbc9644b0` (2026-11-04).

## Decision

The 70 root-level files were inventoried. The project essentials remain at root: `README.md`, `CONTRIBUTING.md`, `requirements.txt` and `license`. The remaining root-level JSON files are legacy audit outputs, numerical reports, manifests or data sidecars.

They are **not moved in this pass**. Existing Markdown, scripts and external notebook references may depend on their current paths. Moving them before a reference audit would create broken links and make provenance harder to follow. This document makes the root contents navigable and establishes safe destinations for the next relocation pass.

## Recommended destinations

| Filename family / file | Destination for a future reference-audited move |
|---|---|
| `188D_*`, `J_*`, `testB_Oh_*`, `RGS_PSU_basis_resolution.json`, `RGS_PSU_SPINORIAL_LIFT_REPORT.json` | `reports/rgs_188d/` |
| `64B_*` | `reports/rgs_188d/oh64/` |
| `K0_*` | `reports/spectral/` |
| `DNS_CONSERVATION_AUDIT*.json` | `reports/dns/` |
| `RGS_PSU_MASTER_REPORT.json`, `rgs_psu_*.json`, `rgs_q6_susy_mobius_report.json` | `reports/rgs_psu/` |
| `su2_b*.json` | `reports/q6_su2/` |
| `audit_*.json`, `*_audit.json`, `*_summary.json`, `report.json`, `quantic_report.json`, `results*.json`, `validation_metrics.json` | `reports/imported_root/` or the owning domain report folder after provenance review |
| `manifest*.json` | `data/manifests/` |
| `stations_omm.json` | `data/astronomy_forecasting/` |
| `model_parameters.json`, `map_derivation_report.json`, `rnr_qvsm_collapse_audit_template.json` | `data/manifests/` or domain-specific configuration folders after content inspection |

## Duplicate and version findings

- `testB_Oh_constrained_report.json` and `testB_Oh_constrained_report (1).json` have the same Git blob SHA `6fce2715e2f99c1553b57d203cf02982a30c0189`: they are byte-identical. Both names are retained for now to avoid breaking references.
- `rgs_psu_green_kubo_geometry_audit.json` and `rgs_psu_green_kubo_geometry_audit (1).json` have different blob SHAs and sizes. They are not duplicates and must not be collapsed without a content/provenance comparison.
- `manifest.json`, `manifest (1).json` and `manifest (2).json` have different blob SHAs. Suffixes are not evidence of duplicate content or chronological precedence.
- Similar filenames elsewhere in `data/` may describe different arrays. Use SHA-256, shape, dtype, units, semantic role and generating script before aliasing or deduplication.

## Known provenance caution

`PROVENANCE_STATUS.json` explicitly marks several 188D arrays (including `K0.npy`, `Kinv.npy`, `J_A.npy`, `J_B.npy`, `J.npy` and `P_G.npy`) as UNVERIFIED. They must not be promoted to certified inputs merely because they are present in the repository. See [the provenance record](../PROVENANCE_STATUS.json) and [the upload triage](UPLOAD_BATCH_2026-11-04_ORGANIZATION_AND_TRIAGE.md).

## Follow-up required before physical relocation

- Search all tracked source, test, example, workflow and Markdown references to each proposed path.
- Inspect JSON schemas and provenance to assign ambiguous files to a domain.
- Compare the two Green–Kubo reports and manifest variants semantically.
- Decide whether the byte-identical `testB_Oh` alias is still needed by external consumers.
- Move only after the path migration map and compatibility policy are prepared.
