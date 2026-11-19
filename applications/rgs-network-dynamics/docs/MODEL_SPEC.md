# Mathematical specification — RGS Network Dynamics v0.1

## 1. Scope

A finite synthetic model of knowledge-state propagation on alternative organizational graphs. It is a baseline for software and mathematical review, not a validated theory of human learning or organizational performance.

## 2. Graphs

Let (V=\{1,\ldots,n\}) be agents and (E\subseteq\{\{i,j\}:i\ne j\}\) be undirected, loop-free links. The adjacency matrix (A\in\{0,1\}^{n\times n}) is symmetric with zero diagonal.

Three graph-generation rules are compared at equal (n) and equal (|E|):

1. **Hierarchical:** a deterministic binary-tree spanning backbone, followed by additional edges preferentially connecting nodes at similar tree depths.
2. **Matrix:** edges ranked by complementarity of synthetic skill profiles, with seeded tie-breaking noise.
3. **Practices-like:** edges ranked by overlap of synthetic skill profiles, with seeded tie-breaking noise.

These are operational toy definitions. They are not a reconstruction of Sun Microsystems' internal organization or Giovanni Vota's historical implementation.

## 3. State space

For (s) skills, let (K_t\in[0,1]^{n\times s}) be the normalized knowledge-state matrix. Let (S\in\{0,1\}^{n\times s}) be a synthetic target profile, with each agent assigned two or three initial skills. The initial state is a seeded random fraction of (S).

Define the row-stochastic neighbor-mixing matrix (P) by

[
P_{ij}=\begin{cases}
A_{ij}/d_i,&d_i>0,\\
1,&d_i=0\text{ and }j=i,\\
0,&\text{otherwise}.
\end{cases}
]

The update is

[
K_{t+1}=\Pi_{[0,1]}\left[K_t+\eta(PK_t-K_t)+\rho(S-K_t)\right],
]

where (Pi_{[0,1]}) clips entries to the unit interval, (0\le\eta,\rho\le1), and (eta,\rho) are declared model parameters.

The first increment represents neighbor mixing; the second represents relaxation toward the fixed synthetic profile. This is a modelling convention, not an empirically established learning law.

## 4. Observables

- Mean knowledge: mean of all entries of (K_T).
- Threshold coverage: fraction of entries of (K_T) strictly greater than 0.60.
- Within-skill dispersion: mean, over skills, of the standard deviation across agents.
- Connectivity: number of connected components of (A).
- Robustness: for uniformly random deletion of 20% of vertices, the expected fraction of retained vertices in the largest connected component, estimated by 80 seeded repetitions.

No scalar composite score is defined. The metrics are not interchangeable and should not be assigned subjective weights after seeing the results.

## 5. Reproducibility

The default seed is 20261119. The code writes JSON and CSV outputs with the run parameters represented by the CLI invocation. Generated output files should be committed only with the command, environment, software versions and commit SHA recorded.

## 6. Limitations and open questions

- Graph-generation rules are stylized and encode assumptions about each organizational archetype.
- The model does not represent tasks, deadlines, incentives, emotional states, communication costs or output quality.
- Synthetic profiles are not employee data.
- The robustness statistic measures graph connectivity only, not business continuity.
- Lower knowledge dispersion is not necessarily desirable.
- No empirical superiority claim follows from a simulation under these assumptions.
- No RGS geometric embedding is implemented in v0.1.

## 7. RGS research gate

A later RGS embedding must be defined as an explicit map from graph vertices or graph features into a specified RGS space. It must then be tested against standard baselines for stated preservation properties, for example graph distances, Laplacian spectra, symmetry actions or perturbation robustness. No such advantage is assumed in this version.

## 8. Evidence labels

- **MODEL DEFINITION:** equations and graph-generation rules specified here.
- **CODE TEST:** finite software invariants checked by the test suite.
- **SYNTHETIC NUMERICAL RESULT:** output of a reproducible run on generated data.
- **OPEN EMPIRICAL HYPOTHESIS:** claims about real organizations and Practices.
- **OPEN COMPUTATIONAL HYPOTHESIS:** claims about advantages of an RGS representation.
