# Q6 analytic gradients and independent finite-difference audit

This increment implements the first derivatives of the finite Euclidean action
in `src/unified_discrete_action.py` and independently checks their directional
contractions against centred differences of the action itself.

## Gradient conventions

- Coframe: ordinary Frobenius gradient with respect to each real (4\times4)
  matrix (e_v), on a fixed determinant-sign sector.
- PSU sheets: complex gradients (g_{\sigma,v}) defined by
  (\delta S=2\operatorname{Re}\sum_v g_{\sigma,v}^\dagger
  \delta\psi_{\sigma,v}).
- Links: left-trivialised coordinates (X_a=i\sigma_a), with
  (U_{uv}(t)=\exp(t\sum_a x_aX_a)U_{uv}). The stored array has shape
  ((192,3)) and follows the canonical Q6 edge order.

The link derivative applies the product rule to each ordered face holonomy.
It handles both (U) and (U^\dagger) occurrences, rather than treating a
face as an unordered set of edges.

## Validation

The test draws independent directions for all 64 coframes, both 64-by-2
complex sheets, and all 192 SU(2) links. It compares the analytic directional
contraction with a centred finite difference of the *full action* at three
step sizes. The link perturbations are exponentiated in SU(2), so the test
stays on the group manifold.

This is a numerical derivative audit for the implemented finite model. It is
not a proof of a continuum limit and does not establish TEGR, a spin
connection, a Clifford/Dirac operator, or Lorentzian dynamics.
