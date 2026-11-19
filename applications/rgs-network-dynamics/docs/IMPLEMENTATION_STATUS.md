# Implementation status — research protocol v1.0

## Added in this change

- Research protocol with hypotheses H1–H7, evidence grading and Practices-specific evidence gates.
- Data dictionary for agents, skills, relations, decisions and outcomes.
- BibTeX starter bibliography and source register with source-type/evidence limitations.
- Competence concentration and coverage indicators.
- Explicit knowledge-state update function with validation.
- Task-based competence coverage and agent-removal primitives.
- Unit tests for the new mathematical utilities.

## Not yet done

- No real organizational dataset has been collected or analyzed.
- No Practices historical claim has been independently verified.
- No cross-case statistical model has been estimated.
- No RGS embedding or comparison with standard graph embeddings has been implemented.
- No empirical advantage or causal effect is claimed.

## Verification

Run from the repository root:

```bash
python -m pip install -r applications/rgs-network-dynamics/requirements.txt
python -m pip install pytest
PYTHONPATH=applications/rgs-network-dynamics/src python -m pytest applications/rgs-network-dynamics/tests -q
```

Test status must be reported only after executing the command in a clean environment.
