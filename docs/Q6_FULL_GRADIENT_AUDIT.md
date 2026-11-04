# Full Q6 gradients: implementation and numerical audit

This extension complements (does not replace) the directional-derivative API
introduced in PR #48.

## Returned arrays

- `coframe`: real Frobenius gradient, shape `(64,4,4)`.
- `psi_plus`, `psi_minus`: complex gradients, shape `(64,2)`, with
  convention `delta S = 2 Re(g^dagger delta psi)`.
- `links`: real array, shape `(192,3)`, in canonical `q6.edges` order.
  Components are derivatives along `U(t)=exp(t i sigma_a) U`.

The coframe gradient includes the derivative of every determinant-volume
factor in hopping, twin-locking and quartic terms. The link gradient includes
both matter hopping and Wilson terms. Wilson derivatives are computed from
the ordered face product, with separate forward and inverse occurrences.

## Independent check

A new test contracts all four gradient blocks with independent random
directions and compares the result with a centred finite difference of the
full action. Every link is perturbed on SU(2), using its exponential map.
The test sweeps three step sizes and checks both a tight best error and a
bounded error envelope. At the smallest step, roundoff can dominate truncation,
so strict monotonic improvement is not required. The pre-existing directional-derivative tests remain intact.

This validates sampled first derivatives of the finite Euclidean action.
It does not solve the stationary equations or establish TEGR, a coframe-derived
spin connection, Clifford/Dirac transport, Lorentzian dynamics, or a continuum
limit. Those require separate geometric definitions and independent audits.
