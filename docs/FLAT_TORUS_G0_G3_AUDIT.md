# Flat-torus G0–G3 generator and audit

**Status:** implemented benchmark / finite numerical verification.  
**Not claimed:** curved RGS refinement, Mosco convergence theorem, or physical continuum limit.

## Construction

The generator `src/flat_torus_refinement.py` builds a periodic Freudenthal triangulation of
\[
\mathbb T_L^d=(\mathbb R/L\mathbb Z)^d,\qquad d\in\{2,3\}.
\]
Each periodic cube is split into \(d!\) simplices by axis permutations. Element coordinates are unwrapped locally to compute exact Euclidean simplex volumes and P1 gradients; vertex identifiers are reduced modulo the periodic grid.

For every simplex, the code assembles the continuous piecewise-linear stiffness matrix
\[
K_{ij}=\sum_{\sigma}\operatorname{vol}(\sigma)
\,\nabla\lambda_i^\sigma\cdot\nabla\lambda_j^\sigma
\]
and a positive lumped 0-form mass/Hodge star
\[
(M_0)_{ii}=\sum_{\sigma\ni i}\frac{\operatorname{vol}(\sigma)}{d+1}.
\]
The scalar energy is \(q_h(u)=u^TKu\). No adjustable cubo-sphere weights are used.

The code also builds canonically oriented edge and triangle incidence matrices. Its exact chain condition is \(B_1B_2=0\).

## G0–G3 mapping

- **G0:** fixed flat torus, side length, dimension and coordinate convention.
- **G1:** periodic Freudenthal meshes on at least three increasing resolutions.
- **G2:** edge lengths, simplex volumes/conditioning, mass and stiffness derived from the same Euclidean metric.
- **G3:** incidence matrices, positive lumped \(M_0\), symmetry/row-sum/semidefinite diagnostics.

## Reproduction

Run:

```bash
PYTHONPATH=. python examples/flat_torus_refinement_audit.py
python -m pytest tests/test_flat_torus_refinement.py -q
```

The JSON report records the grid levels, counts, spacing, incidence residual, mass totals, stiffness residuals, simplex conditioning and smooth Fourier-mode energy error. When run in GitHub Actions, the commit SHA is read from `GITHUB_SHA`.

## Expected 2D counts

For an \(N\times N\) periodic square grid triangulated with two Freudenthal triangles per cell:
- vertices: \(N^2\);
- edges: \(3N^2\);
- triangles: \(2N^2\);
- simplices: \(2N^2\).

The total lumped mass equals \(L^2\), up to floating-point roundoff.

## Limits and next mathematical obligations

1. This is a flat reference mesh, not an RGS-derived curved geometry.
2. Only the 0-form Hodge star is implemented. A positive edge/2-form Hodge star and a DEC/Whitney comparison are separate tasks.
3. The finite-element energy audit is a consistency diagnostic, not the liminf/recovery proof.
4. A decreasing error across three grids is not a proof of Mosco convergence.
5. The next step is a curved embedding/coframe family with uniform shape regularity, metric consistency, compactness and analytical liminf/recovery estimates.
