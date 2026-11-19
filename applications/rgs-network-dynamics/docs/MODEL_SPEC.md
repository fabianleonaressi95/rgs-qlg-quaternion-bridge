# Mathematical specification — RGS Network Dynamics

## 1. Scope and status

This application contains two related but distinct synthetic model components:

1. **v0.1 comparative graph baseline:** equal-budget graph generation and a fixed-profile knowledge-state update.
2. **v0.2 research utilities:** competence coverage/concentration, an exploratory diffusion-decay-input update, and task-based resilience under agent removal.

They are software and modelling baselines, not validated theories of human learning or organizational performance.

## 2. Comparative graphs

Let (V=\{1,\ldots,n\}) be agents and (E\subseteq\{\{i,j\}:i\ne j\}) be undirected, loop-free links. The adjacency matrix (A\in\{0,1\}^{n\times n}) is symmetric with zero diagonal.

Three stylized graph-generation rules are compared at equal (n) and equal (|E|):

1. **Hierarchical:** binary-tree spanning backbone plus links preferentially connecting nodes at similar tree depths.
2. **Matrix:** edges ranked by complementarity of synthetic skill profiles.
3. **Practices-like:** edges ranked by overlap of synthetic skill profiles.

These are operational toy definitions, not reconstructions of any actual organization.

## 3. v0.1 knowledge-state update

For (s) skills, let (K_t\in[0,1]^{n\times s}) be the agent-by-skill knowledge state and (S\in\{0,1\}^{n\times s}) the synthetic target profile.

For row-stochastic neighbor-mixing matrix (P), the baseline update is

[
K_{t+1}=\Pi_{[0,1]}\left[K_t+\eta(PK_t-K_t)+\rho(S-K_t)\right].
]

Here (eta) controls neighbor mixing and (ho) relaxation toward a fixed synthetic profile. This is a modelling convention, not an empirically established learning law.

## 4. v0.2 exploratory knowledge update

The separate utility in `src/knowledge_dynamics.py` implements

[
X_{t+1}=\Pi_{[0,1]}\left[X_t+\eta(P_KX_t-X_t)-\delta X_t+U_t\right],
]

where (P_K) is row-stochastic, (delta\in[0,1]) is a decay parameter and (U_t) is an external learning input. It is intentionally not conflated with the v0.1 relaxation-to-profile term. Parameters require empirical estimation or transparent sensitivity analysis.

## 5. Competence distribution

Let (C\in[0,1]^{n\times m}) be an assessed proficiency matrix. For threshold (	au),

[
N_k(\tau)=\sum_i\mathbf 1(C_{ik}\geq\tau).
]

For (s_k=\sum_i C_{ik}>0), define (D_{ik}=C_{ik}/s_k) and

[
H_k=\sum_iD_{ik}^2.
]

(H_k) is a Herfindahl-style concentration indicator. It is undefined for an uncovered skill. (1/H_k) is an effective-number interpretation, not proof of practical substitutability.

## 6. Task-based resilience

A task is represented by a nonempty set of required skills. It is feasible if every required skill is available above threshold from at least one available agent.

For a collection of tasks, (F(C,\mathcal T,\tau)) is the fraction feasible. After removing agents (R), report

[
F_{\mathrm{after}}=F(C_{V\setminus R},\mathcal T,\tau).
]

If baseline coverage is positive, the retention ratio is (F_{\mathrm{after}}/F_{\mathrm{before}}). If baseline coverage is zero, the ratio is undefined and represented as NaN. Test random removal, high-centrality removal and removal of agents uniquely covering critical skills as separate scenarios.

This measure is task capability, not merely graph connectivity. It still assumes competence is transferable to the task and must be validated against operational records.

## 7. Existing graph metrics

The original prototype reports mean knowledge state, threshold coverage, within-skill dispersion, connected components, mean degree and largest-component retention after random vertex removal. Lower dispersion is not automatically better; graph connectivity is not equivalent to business continuity or output quality.

## 8. Reproducibility and limitations

- Synthetic profiles are not employee data.
- Graph-generation rules encode assumptions about organizational archetypes.
- The models omit incentives, deadlines, task dependencies, emotional states, communication costs and quality outcomes.
- No composite score is defined.
- No empirical superiority claim follows from simulation.
- No RGS geometric embedding is implemented.

A future RGS embedding must be explicitly defined and compared with standard graph baselines under prespecified preservation criteria, such as graph distances, Laplacian spectra, symmetry actions and perturbation robustness.

## 9. Evidence labels

- **MODEL DEFINITION:** equations and graph rules specified.
- **CODE TEST:** finite software invariants checked by the test suite.
- **SYNTHETIC NUMERICAL RESULT:** reproducible output on generated data.
- **OPEN EMPIRICAL HYPOTHESIS:** claims about real organizations and Practices.
- **OPEN COMPUTATIONAL HYPOTHESIS:** claims about RGS representation advantages.
