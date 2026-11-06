# Recovered RGS/Q6 embedded-graph Wilson baseline

**Status:** finite baseline on the recovered 64-vertex graph; not a curved-tetrad or continuum result.

## Why this baseline is useful

The archived RGS seed provides 64 coordinates and the canonical 192 Q6 edges, with source hashes checked by `load_rgs_seed_data()`. It also provides an edge coframe with four components per edge, not a vertexwise invertible 3x3 tetrad. It would be mathematically incorrect to silently reinterpret that edge array as a vertex tetrad.

The baseline therefore makes the following explicit reference choices:
- recovered RGS coordinates and Q6 edge incidence;
- actual Euclidean lengths of the 192 embedded edges;
- identity 3x3 coframe at every vertex;
- identity SU(2) spin links.

This isolates the effect of nonuniform embedded edge lengths in the Wilson graph kernel. It is a useful reproducible finite reference before introducing a separately justified vertex metric and spin connection.

## Run

```bash
PYTHONPATH=. python examples/rgs_q6_wilson_spectral_baseline.py
```

The JSON report records the source geometry, edge-length range, Wilson-kernel eigenvalue bounds, near-zero and below-cutoff multiplicities, Euclidean Wilson–Dirac singular values and gamma5-Hermiticity residual.

## Evidence boundary

The output is a finite-dimensional baseline only. It does not establish:
- that the embedded graph carries the intended curved RGS metric;
- that the edge coframe determines a unique vertexwise tetrad;
- that identity links approximate Levi-Civita spin parallel transport;
- a curved Brillouin-zone species count;
- a refinement-uniform spectral gap or strong-resolvent limit.

The next geometric input must be a declared vertexwise metric/tetrad and a spin connection with a justified edge integration rule. That data cannot be inferred uniquely from the current four-component edge coframe without additional mathematical choices.
