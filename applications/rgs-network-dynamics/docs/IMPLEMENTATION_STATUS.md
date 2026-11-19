# Implementation status — research protocol v1.0

## Added

- Research protocol with hypotheses H1–H7, evidence grading and Practices-specific evidence gates.
- Data dictionary for agents, skills, relations, decisions and outcomes.
- BibTeX starter bibliography and source register with source-type/evidence limitations.
- Competence concentration and coverage indicators.
- Explicit knowledge-state update function with validation.
- Task-based competence coverage and agent-removal primitives.
- Unit tests for graph construction, reproducibility and research metrics.

## Verification performed

- Retrieved the source and test files from branch `feat/rgs-network-dynamics-prototype`.
- Executed the 10-test suite in an isolated local Python environment with the project's `conftest.py` import-path setup: **10 passed**.
- Python compilation check: passed.
- Executed the default synthetic experiment (30 agents, 8 skills, 45 edges per architecture, 160 update steps, 80 random-removal repetitions).

The GitHub repository could not be cloned into the local runner because outbound DNS/network access to github.com was unavailable. Therefore this is a local verification of the retrieved branch files, not a GitHub Actions run. No CI status is implied.

## Default synthetic run

| Architecture | Components | Mean knowledge | Coverage > 0.60 | Within-skill dispersion | Largest-component retention |
|---|---:|---:|---:|---:|---:|
| Hierarchical | 1 | 0.25998 | 0.0000 | 0.02902 | 0.96458 |
| Matrix | 3 | 0.24635 | 0.0250 | 0.09726 | 0.89688 |
| Practices-like | 5 | 0.27686 | 0.0375 | 0.15337 | 0.74635 |

These figures are **synthetic outputs of the declared toy model only**. They do not show that any real organizational design is better. In particular, the Practices-like graph has the highest mean knowledge and threshold coverage in this run but the lowest graph-connectivity retention; metrics trade off and must not be collapsed into a ranking.

## Not yet done

- No real organizational dataset has been collected or analyzed.
- No Practices historical claim has been independently verified.
- No cross-case statistical model has been estimated.
- No RGS embedding or comparison with standard graph embeddings has been implemented.
- No empirical advantage or causal effect is claimed.

## Reproduce

From the repository root:

```bash
python -m pip install -r applications/rgs-network-dynamics/requirements.txt
python -m pip install pytest
python -m pytest applications/rgs-network-dynamics/tests -q
python applications/rgs-network-dynamics/src/network_dynamics.py --seed 20261119 --output applications/rgs-network-dynamics/results
```
