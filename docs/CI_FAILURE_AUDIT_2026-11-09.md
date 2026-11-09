# CI failure audit — 2026-11-09

## Executive result

The supplied failed run is the workflow **RGS-QLG Quaternion Bridge — Full Audit**, job **Quaternion / SU2 / SO3 / QLG tests**. Its first Python validation step is:

```bash
python -m compileall -q src tests examples scripts
```

The observed failure is a syntax error in `tests/test_q6_dyadic_refinement.py`, not a failed quaternion, SU(2), SO(3), or QLG mathematical test. The test suite and all later audit steps are downstream of compilation and therefore were not executed in that run.

## 1. Root cause fixed

The import statement had literal backslash-n sequences inside the source line:

```python
from src.q6_dyadic_refinement import (\n    build_q6_dyadic_refinement, ...
```

Python treats the backslash as explicit line continuation and then encounters `n`, producing:

```
SyntaxError: unexpected character after line continuation character
```

The import has been rewritten using real line breaks. The file was fetched back from GitHub and verified to contain the corrected multi-line import.

## 2. Additional defects found and corrected during the audit

These were not reported as failures in the supplied excerpt; they were latent defects found by scanning nearby changed/relevant files.

1. `examples/kathara_cellular_complex_audit.py`: JSON serialization appended the two literal characters `\\n` rather than a newline. Corrected to append the newline escape `\n`.
2. `scripts/run_q6_d5_saturated_kernel_bounds.py`: same invalid trailing literal `\\n` in the JSON payload. Corrected to append `\n`.
3. `tests/test_quaternion_algebra.py`: duplicate import of `quaternion_dot` removed.
4. `theorem/README.md`: the T060 index row and T061 row were joined by a literal `\\n`, corrupting the Markdown table. Replaced with an actual line break.

These corrections do not establish that every audit passes; they remove the identified source/serialization/index defects.

## 3. Workflow inventory and what the supplied failure means

| Workflow | Trigger | Consequence for the supplied run |
|---|---|---|
| `.github/workflows/python-package.yml` — Full Audit | Push to `main` or manual dispatch | This is the workflow named in the supplied annotation. Compilation is first; later steps did not run after the syntax error. |
| `.github/workflows/tests.yml` — Python tests | Pull request targeting `main` | Separate PR workflow; it is not the push-to-main full-audit workflow. |
| `.github/workflows/wilson-periodic-audit.yml` | Dedicated audit branch push or PR to `main` touching listed paths | Separate focused Wilson/curved-Q6 checks; no failure log for this workflow was supplied. |
| `.github/workflows/experiment-dependency-audit.yml` | Push to `main` touching `experiments/**`, its audit script, or workflow file | Separate dependency graph audit; the files corrected in this pass do not match its path filter. |

The full-audit workflow orders the following after compileall: quaternion algebra tests, SU(2)–SO(3) tests, QLG holonomy tests, full pytest suite, lattice demo, Wilson/refinement, flat and curved Mosco/geometry audits, Möbius/integrated audits, and Q6–D5 numerical scripts/artifact uploads. Since compileall stopped the job, these later steps are **not classified as failed**; they are **not yet tested by that run**.

## 4. Theorem repository audit

The exact finite Q6 midpoint-subdivision theorem was already present as **T055** in `fabianleonaressi95/theorem` and was already linked from its theorem index. No duplicate theorem was created.

A missing foundational theorem-level record for the implementation's quaternion/SU(2)/SO(3) convention was added as **T063**:
- `theorems/T063-quaternion-su2-so3-double-cover.md`
- theorem index row added to `theorem/README.md`

T063 proves the scalar-first unit-quaternion isomorphism with SU(2), the conjugation homomorphism onto SO(3), its kernel \(\{\pm1\}\), and the two-lift property. It formalizes standard mathematics and the implementation conventions; it does not claim a new physical spin law or QLG continuum theorem.

A second missing theorem record was added as **T064** for the finite SU(2) lattice layer: path transport transforms by endpoint gauge factors, closed holonomy by conjugation, and the Wilson trace is gauge invariant and real with normalized value in [-1,1]. It formalizes `src/holonomies.py` and `tests/test_holonomies.py`; it is not a continuum Yang–Mills theorem.

## 5. Remaining verification

- The corrected source files and theorem records have been fetched back from GitHub and checked for the intended text.
- No local Python test execution is claimed in this audit.
- The available GitHub connector exposes retrieval of workflow runs only for pull-request-triggered runs, not a repository-wide listing of push-triggered runs. The supplied Full Audit is a push-triggered workflow. Therefore the post-fix run's conclusion and the exact text of the one warning and one notice cannot be independently read from the available run-log interface here.
- A new successful run must show, in order: compileall; quaternion tests; SU(2)–SO(3) tests; holonomy tests; complete pytest; and all downstream audit/artifact steps.

## 6. Evidence links

- [Full Audit workflow](../.github/workflows/python-package.yml)
- [PR test workflow](../.github/workflows/tests.yml)
- [Wilson/curved-Q6 workflow](../.github/workflows/wilson-periodic-audit.yml)
- [Experiment dependency workflow](../.github/workflows/experiment-dependency-audit.yml)
- [Corrected Q6 dyadic tests](../tests/test_q6_dyadic_refinement.py)
- [Q6 dyadic theorem T055](https://github.com/fabianleonaressi95/theorem/blob/main/theorems/T055-q6-dyadic-subdivision-chain-map.md)
- [Quaternion/SU(2)/SO(3) theorem T063](https://github.com/fabianleonaressi95/theorem/blob/main/theorems/T063-quaternion-su2-so3-double-cover.md)
