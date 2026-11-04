# Q6 action: analytic directional derivatives and verification gate

This increment adds an analytic *directional derivative evaluator* for the
finite-dimensional Euclidean action in `src/unified_discrete_action.py`.
It is not an optimiser or a claim that a stationary solution has been found.

## What is differentiated

- the fixed-frame coframe edge regulariser;
- the determinant-volume dependence of hopping, twin locking and quartic terms;
- both complex PSU sheets under independent real/imaginary perturbations;
- SU(2) links under left-trivialised variations (U(t)=e^{tX}U);
- Wilson plaquettes with the actual path-product order and both forward and
  inverse occurrences of a stored edge.

The API `directional_action_derivative(...)` returns a derivative for each
action component and their sum. Link generators are keyed by canonical Q6
edges and must be anti-Hermitian and traceless.

## Independent numerical gate

`tests/test_unified_discrete_gradients.py` compares the analytic directional
derivative against centred differences of the *full action evaluator*, with
independently perturbed fields and group-valued links. It covers coframe-only,
plus-sheet-only, minus-sheet-only, link-only and combined directions over a
step-size sweep.

The finite-difference checks validate first derivatives at sampled generic
configurations. They do not prove global correctness, solve the Euler–Lagrange
equations, or establish a continuum limit.

## Spin connection / Clifford boundary

No spin connection is silently inferred from the present SU(2) links. A
subsequent geometric extension must first specify the coframe convention and
metric signature, then define a discrete orthonormal-frame transport and test
its covariance and torsion. A polar factor of neighbouring coframes alone is
only a frame-alignment candidate; it is not automatically the torsion-free
Levi-Civita spin connection, nor the inertial spin connection required for a
covariant teleparallel formulation. Clifford transport and a fermionic
operator must be built only after that structure is fixed.
