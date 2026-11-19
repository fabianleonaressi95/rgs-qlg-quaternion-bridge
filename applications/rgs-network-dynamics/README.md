# RGS Network Dynamics — prototype v0.1

A small, reproducible Python baseline for comparing three **synthetic** organizational network topologies: hierarchical, matrix / cross-functional, and Practices-like competence clusters.

This prototype tests graph construction, knowledge-state diffusion and structural robustness. It does **not** establish that one organizational model is superior, reproduce Sun Microsystems' historical results, validate the empirical effectiveness of Practices, or prove that RGS geometry improves network modelling.

## Model

There are 30 agents and 8 skills. The same synthetic skill profiles and initial knowledge matrix are used across the three architectures. Each architecture has 45 undirected edges.

For each time step:

```text
K[t+1] = clip(K[t] + eta * (P @ K[t] - K[t]) + rho * (S - K[t]), 0, 1)
```

Here `K` is the agent-by-skill knowledge matrix, `P` is a row-stochastic neighbor-mixing matrix, `S` is the fixed synthetic skill-profile matrix, `eta` controls neighbor diffusion and `rho` controls relaxation toward the profile. This is a modelling choice, not an established law of organizational behaviour.

## Metrics

- skill coverage: fraction of agent-skill entries above 0.60;
- mean knowledge state;
- within-skill dispersion across agents;
- number of connected components and mean degree;
- robustness: mean fraction of retained nodes in the largest connected component after random removal of 20% of nodes.

Metrics are interpreted separately. Lower dispersion is not automatically better, and graph connectivity is not equivalent to organizational performance.

## Run

From the repository root:

```bash
python -m pip install -r applications/rgs-network-dynamics/requirements.txt
python applications/rgs-network-dynamics/src/network_dynamics.py --seed 20261119 --output applications/rgs-network-dynamics/results
```

The command writes `summary.json` and `summary.csv` into the output directory and prints the CSV table. Results are generated locally and are not committed as pre-verified data.

## Test

```bash
python -m pip install pytest
python -m pytest applications/rgs-network-dynamics/tests -q
```

## Next research gates

1. Sensitivity analysis over diffusion, learning, thresholds and graph-generation rules.
2. Add standard graph baselines: Erdős–Rényi, small-world and degree-preserving rewiring.
3. Separate knowledge diffusion from task allocation, decision latency and output quality.
4. Pre-register metrics and objective weights before empirical comparison.
5. Test RGS embeddings against standard graph embeddings using explicit preservation criteria: graph distances, Laplacian spectra, symmetries and robustness.
6. Validate against real organizational data only under an agreed protocol, consent and appropriate anonymization.

## Evidence classification

| Item | Status |
|---|---|
| Graph and update-rule definitions | MODEL DEFINITION |
| Code-level invariants | TESTED only when the test suite passes |
| Numerical output | SYNTHETIC COMPUTATION; generate locally |
| Practices performance advantage | OPEN EMPIRICAL HYPOTHESIS |
| RGS geometric advantage | OPEN MATHEMATICAL / COMPUTATIONAL HYPOTHESIS |
| Commercial demand | NOT ASSESSED |

No claim of empirical validation or startup readiness is made by this prototype.
