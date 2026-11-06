# Levi-Civita spin connection on the Ellis-RGS refinement

**Status:** explicit geometric edge-transport construction and finite-mesh tests. Refinement-uniform estimates and continuum Dirac convergence remain OPEN.

## 1. Smooth target geometry

The target is
`g=du^2+R(u)^2d\Omega^2,quad R(u)=\sqrt{u^2+r_0^2}`,
with exact isometric embedding
`F(u,n)=(R(u)n,r_0\,\operatorname{arsinh}(u/r_0))\in\mathbb R^4`.

Its unit normal is
`N(u,n)=(r_0 n/R,-u/R)`.
For a tangent vector `V`, Levi-Civita parallel transport along a curve in the hypersurface is characterized by vanishing tangential component of its ambient derivative. Differentiating `<V,N>=0` gives the exact ODE
`V'=-<V,N'>N`.

This avoids guessing a spin connection from endpoint coframes.

## 2. Discrete edge rule

For each oriented edge `i -> j`:
1. Recover `(u_i,n_i)` from the exact embedded R4 vertices.
2. Fix the path: affine interpolation in `u`, shortest great-circle interpolation in `n\in S^2`.
3. Integrate the ambient parallel-transport ODE by classical RK4.
4. Project/re-orthonormalize the transported tangent frame at the endpoint.
5. Express the transported source frame in the target frame to obtain `R_{ji}\in SO(3)`.
6. Lift `R_{ji}` to the repository's scalar-first SU(2) convention.

The path rule is part of the discretization and is recorded in every audit. Antipodal angular endpoints are rejected because the shortest path is not unique.

## 3. Code and tests

- `src/ellis_spin_connection.py`
- `tests/test_ellis_spin_connection.py`

Run:
```bash
PYTHONPATH=. python -m pytest tests/test_ellis_spin_connection.py -q
```

Tests cover frame orthonormality/tangency, exact radial transport, edge reversal, SU(2) unitarity, and link construction on the actual Ellis refinement chain.

## 4. What is and is not established

The construction gives a well-defined finite-edge SU(2) parallel transporter for the stated smooth Ellis metric and the stated edge-path convention. It is geometrically derived, unlike the previous polar coframe comparison.

It does **not** yet prove:
- convergence of the edge transport to the continuum connection uniformly over all refinement edges;
- gauge-independent curvature estimates or plaquette holonomy convergence;
- a spinorial Dirac operator consistent with this connection;
- Mosco convergence, recovery sequences, or strong-resolvent convergence.

## 5. Next proof obligation

Prove a local consistency estimate for a shape-regular refinement family:
`\|U_{ij}^{(h)}-\mathcal P_{\gamma_{ij}}\|\le C h_{ij}^{p+1}`
for the chosen path and RK quadrature order, then combine it with a uniform discrete Poincare/Rellich estimate and a common-Hilbert-space embedding. The order `p` must be derived from the actual edge path and integrator, not assumed.
