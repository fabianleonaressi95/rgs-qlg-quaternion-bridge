# Programme-wide theorem dependency reconciliation — 27 November 2026

This report records the first record-level reconciliation of the primary theorem corpus and the companion dependency graph. It is an audit of repository metadata, not a new mathematical theorem.

## Inventory

- Primary `main`: 284 theorem Markdown files.
- Primary working branch `rgs-q6-antitwin-spacetime-bridge`: 289 theorem Markdown files (five branch-only proposals: T146–T149 and T161).
- Generated primary theorem register: 280 valid rows covering 280 unique files.
- Four files on primary `main` are missing from that register: T157 and T159 vertex-fixed Jitterbug metric-obstruction cards; T158 and T160 snub-cube twist-angle cards.
- The 51-step formal suite S01–S51 remains represented separately.

The generated register must be regenerated from the current theorem tree before its status counts are used as corpus-wide totals. Its labels are machine-assisted heuristics.

## Declared dependency extraction

The 280 valid register rows contain 410 dependency tokens across 142 rows:

| Resolution | Count | Graph treatment |
|---|---:|---|
| Unique display-ID match | 306 | Imported as `declared-dependency-unreviewed` |
| Ambiguous display-ID match | 101 | Not converted to an edge; all candidate paths retained in primary audit JSON |
| No matching display ID | 3 | Not converted to an edge; requires repair |

The register has 32 duplicated display-ID groups. Consequently, a bare dependency such as `T021` or `T129` is not a stable key in this corpus. Full paths are required to disambiguate records.

The 306 imported arrows are **not proof dependencies**. They faithfully preserve unique declarations from the generated CSV, with a warning on every edge. A card-level audit must inspect both endpoints and the exact statement before upgrading an edge to `proves`, `uses-definition`, `implements`, `application-of` or another mathematical relation.

## T161 interface

The T161 spine retains six separate reviewed links:
- T129 → T161: reuse of the golden cut-and-project frame;
- T146 → T161: shared Q6–IVM carrier, with 64 projected points distinguished from 32 zonotope vertices;
- T147 → T161: Euclidean Spin(6) input, not a Lorentzian identification;
- T148 → T161: conditional PSU–Anti-Twin timelike/spacelike construction;
- T149 → T161: exact gapped-operator obstruction;
- T150 → T161: comparison-only D6/Fibonacci suborbit.

T161 is not a universal parent theorem. Rumin, D5, FEEC, Morandotti, arithmetic, PRISM and MHD remain connected through their own definitions and application interfaces, not through an invented implication from the triacontahedron.

## Machine-readable graph

- [Unified graph JSON](UNIFIED_THEOREM_EDGES.json): 393 nodes and 734 edges after importing the 306 provisional register declarations.
- [Human-readable dependency graph](UNIFIED_THEOREM_DEPENDENCY_GRAPH.md): preserves the hub and T161 architecture.
- [Primary register audit report](https://github.com/fabianleonaressi95/theorem/blob/rgs-q6-antitwin-spacetime-bridge/docs/T161_RECORD_LEVEL_INTEGRATION_AUDIT_2026-11-27.md).
- [Primary machine-readable audit](https://github.com/fabianleonaressi95/theorem/blob/rgs-q6-antitwin-spacetime-bridge/docs/T161_REGISTER_DEPENDENCY_AUDIT_2026-11-27.json).

## Next mandatory gates

1. Regenerate and validate the primary theorem register; add the four omitted main-tree cards.
2. Replace display-ID dependencies with path-stable identifiers; resolve 101 ambiguous and 3 missing references.
3. Review the 306 unique declarations against the endpoint cards before changing their status.
4. Continue the chapter-by-chapter crosswalk across all 27 chapters and the full 51-step suite.
5. Maintain separate evidence labels: standard result, proved internal theorem, exact finite identity, numerical audit, source-reported result, conditional bridge, candidate, physical hypothesis, open.

No theorem has been renumbered, deleted, or promoted in status by this audit.
