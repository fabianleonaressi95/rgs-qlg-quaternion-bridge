# Unified finite RGS–PSU–Q6 multiscale operator

## Definition

The finite state space is the Cartesian product

\[
\Omega = Q_6 \times K_2 \times P_3,\qquad |\Omega|=64\cdot2\cdot3=384.
\]

Here \(Q_6\) is the 6-cube configuration graph, \(K_2\) is a two-sheet fibre used as an Anti-Twin bookkeeping factor, and \(P_3\) is a three-level scale index. This graph product is a model choice. In particular, \(K_2\) alone does not establish a Möbius bundle, non-orientability, torsion, or a Pin structure.

The undirected conductance graph has 1,600 edges:
- Q6 edges: \(192\cdot2\cdot3=1,152\);
- sheet edges: \(64\cdot3=192\);
- scale edges: \(64\cdot2\cdot2=256\).

The continuous-time row generator is assembled from these symmetric unit conductances plus a divergence-free oriented square current on the first two Q6 bits. On each square \(0\to1\to3\to2\to0\), forward and reverse rates are \(1+a\) and \(1-a\), with \(a=0.2\). The cycle is replicated across all sheet and scale sectors. Since \(0\le a<1\), all transition rates remain nonnegative.

## Single operator and observables

The assembled Markov generator is the common operator acting on the entire product state space. The Q2 relational-clock observable is the adjacency operator of the first two (least significant) Q6 bits, lifted identically over sheet and scale factors. Its commutator with the assembled generator vanishes for this construction. This is an exact finite compatibility identity, not a derivation of proper time.

The uniform distribution is stationary because the symmetric conductance graph is regular in each Q6/sheet sector and the imposed square current is divergence-free. The symmetric part of the generator is negative semidefinite. The stationary Markov entropy-production rate is evaluated by the Schnakenberg expression

\[
\dot S_{\rm path}/k_B =
\frac{1}{2}\sum_{i,j}
(\pi_i q_{ij}-\pi_jq_{ji})
\log\frac{\pi_iq_{ij}}{\pi_jq_{ji}}.
\]

For the uniform stationary distribution, the audit reports:
- reversible control: 0 nats per model-time unit;
- driven square: 0.010136627702704109 nats per model-time unit.

The affinity is imposed externally. It is not derived from PSU geometry or a gravitational action.

## Reproducibility

Run:

~~~bash
python -m src.rgs_unified_multiscale_operator --output reports/rgs_unified_multiscale_operator.json
pytest -q tests/test_rgs_unified_multiscale_operator.py
~~~

The report includes matrix hashes, state/edge counts, generator row-sum residual, stationarity residual, minimum off-diagonal rate, relative clock commutator, symmetric-part spectral bound, and entropy-production controls.

## Audit status and limits

The finite operator tests:
1. state count and product-edge count;
2. Markov generator row sums;
3. uniform stationarity;
4. nonnegative off-diagonal transition rates;
5. Q2-clock commutator;
6. negative-semidefinite symmetric generator part;
7. zero entropy production in the reversible control;
8. positive entropy production in the driven model.

It does **not** yet derive:
- a continuum limit or a vierbein from the graph;
- the existing 188-dimensional RGS tangent-space vierbein;
- a teleparallel torsion tensor or its relation to the two-sheet factor;
- a physical renormalization flow from the three scale labels;
- the Einstein–Klein–Gordon critical-collapse solution;
- a gravitational cause for microscopic thermodynamic fluctuations;
- a physical arrow of time or the nonexistence of fundamental time.

The next mathematically meaningful step is a convergence theorem or numerical Mosco-type test for a specified family of operators, together with an independently defined map from discrete edge data to the existing RGS tangent frame. A finite commutator identity is not a substitute for either.
