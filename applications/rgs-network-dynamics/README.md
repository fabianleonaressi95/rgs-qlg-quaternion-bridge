# RGS Network Dynamics — synthetic prototype v0.2

An isolated research application for comparing organizational network structures, competence distribution, knowledge transfer and task-based resilience.

**Evidence boundary:** all current executable experiments use synthetic data. The code does not reconstruct Sun Microsystems Italia, validate Giovanni Vota's Practices empirically, establish that self-management is superior, or demonstrate an advantage of RGS embeddings.

## Research documents

- [Research protocol](docs/RESEARCH_PROTOCOL.md) — research question, hypotheses H1–H7, comparative cases, mixed-method design, evidence grading and validation gates.
- [Data dictionary](docs/DATA_DICTIONARY.md) — schemas for agents, skills, relations, decisions and outcomes.
- [Source register](bibliography/SOURCE_REGISTER.md) — bibliography quality and evidence limitations.
- [BibTeX bibliography](bibliography/references.bib).
- [Mathematical specification](docs/MODEL_SPEC.md) — v0.1 synthetic graph model.
- [Implementation status](docs/IMPLEMENTATION_STATUS.md).

## Current computational components

- `src/network_dynamics.py`: synthetic hierarchical, matrix and Practices-like graphs; bounded knowledge-state diffusion; connectivity metrics; random vertex-removal graph robustness.
- `src/competence_model.py`: thresholded skill coverage, Herfindahl-style competence concentration and effective-agent count.
- `src/knowledge_dynamics.py`: explicit bounded diffusion/decay/external-input update for exploratory longitudinal models.
- `src/resilience.py`: task-based capability coverage under agent unavailability and ordering by uniquely covered skills.
- `tests/`: finite software checks for graph invariants and research metrics.

The original graph-connectivity robustness metric and the new task-based competence resilience metric answer different questions. They must not be reported as interchangeable.

## Install and test

From the repository root:

```bash
python -m pip install -r applications/rgs-network-dynamics/requirements.txt
python -m pip install pytest
python -m pytest applications/rgs-network-dynamics/tests -q
```

## Run the synthetic baseline

```bash
python applications/rgs-network-dynamics/src/network_dynamics.py \
  --seed 20261119 \
  --output applications/rgs-network-dynamics/results
```

The command writes JSON and CSV summaries locally. Do not commit generated output as verified evidence unless the command, environment, dependency versions and commit SHA are recorded.

## Research priorities

1. Recover independent historical documentation for Practices / Sun Microsystems Italia.
2. Obtain comparable team-level data for Buurtzorg, Morning Star and FAVI.
3. Add standard graph baselines and preregister sensitivity ranges.
4. Separate decision autonomy, knowledge transfer, task allocation, decision latency and output quality.
5. Validate the task-based resilience measures against real operational tasks.
6. Only then define an explicit RGS embedding and compare it against standard graph methods using prespecified distance, spectral, symmetry and robustness criteria.

## Evidence classification

| Item | Status |
|---|---|
| Model definitions and source code | SPECIFIED |
| Software invariants | TESTS WRITTEN; execution status must be checked |
| Generated numerical outputs | SYNTHETIC ONLY |
| Practices rollout and profit claims | UNVERIFIED TESTIMONIAL CLAIMS |
| Comparative organizational effects | OPEN EMPIRICAL HYPOTHESES |
| RGS-specific modelling advantage | OPEN COMPUTATIONAL HYPOTHESIS |

No empirical validation, causal conclusion or commercial-readiness claim is made.
