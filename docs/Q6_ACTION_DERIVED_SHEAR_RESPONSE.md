# Q6 action-derived shear response: exact finite derivative and physical obstruction

## Scope

The repository already contains a finite Euclidean Q6 action in
`src/unified_discrete_action.py` and analytic directional derivatives in
`src/unified_discrete_gradients.py`. This record reuses those objects.
It does not introduce arbitrary edge weights, a fitted force law, or a
new dynamics.

For a declared symmetric shear generator (H), define the coframe
variation
[
e_v(\gamma)=(I+\gamma H)e_v,qquad
\dot e_v(0)=He_v,
]
holding the two PSU matter sheets and SU(2) links fixed. The action response is
[
\boxed{\Pi_H^{\rm action}
=-\frac1V\left.\frac{d}{d\gamma}
S_{Q6}[e(\gamma),\psi_+,\psi_-,U]\right|_{\gamma=0}.}
]
For the default (xy) convention in a 4-by-4 coframe,
(H=(E_{12}+E_{21})/2); index 0 is not strained.

## What is exact now

The derivative is evaluated by the existing analytic
`directional_action_derivative`, with direction
(dot e_v=He_v). The new module
`src/q6_action_stress.py` exposes this response and an independent central
finite-difference audit. The derivative is split by action sector:
coframe, Wilson, hopping, twin and quartic.

Because the strain generator is symmetric and traceless,
(mathrm{tr}(H)=0), the first variation of every local determinant volume
vanishes:
[
\left.\frac{d}{d\gamma}\log|\det e_v(\gamma)|\right|_0
=\mathrm{tr}(e_v^{-1}He_v)=\mathrm{tr}(H)=0.
]
Consequently the volume-only factors in the current hopping, twin and
quartic terms have zero first-order shear response when matter fields are
held fixed. The current Wilson term depends only on link holonomies and
has no coframe/metric dependence, so its response is also exactly zero.
The present action's nonzero response can therefore come from the
fixed-frame coframe-difference regularizer. This is a property of the
implemented action and strain convention, not a universal stress law.

## Decisive obstruction

The finite action currently has no metric-dependent edge kinetic/hopping
coefficient and no metric-dependent Wilson plaquette coupling. It therefore
does not yet provide the full shear stress of a mechanical/gauge medium.
In particular, (Pi_H^{\rm action}) must not be identified with the
microscopic (Pi_{xy}) in Green--Kubo without adding and deriving the
missing metric dependence and specifying the equilibrium ensemble.

The next model-level extension must derive, rather than choose:
1. the edge metric/length dependence in the matter kinetic term;
2. the metric dependence of the gauge/plaquette term (or an equivalent
   discrete Hodge-star formulation);
3. the physical volume measure and its strain transformation;
4. momenta or a justified overdamped stochastic dynamics and its invariant
   Gibbs measure.

Only then can one form a centered stress covariance and calibrate time,
temperature and volume for Green--Kubo. Stokes matching remains a later
hydrodynamic limit.

## Reproduction

```bash
PYTHONPATH=. python -m pytest tests/test_q6_action_stress.py -q
```

## Evidence ledger

- **EXACT BY DIFFERENTIATION:** action response for the explicitly declared
  linear coframe shear and fixed matter/link fields.
- **NUMERICALLY VERIFIED:** analytic derivative against central differences
  when the test suite is run.
- **OPEN:** physical metric stress, equilibrium dynamics, FDT normalization,
  Green--Kubo viscosity and Stokes limit.


## Previously derived inverse-square edge-weight variation

A distinct geometric edge-energy convention used in the earlier Q6 stress
work is (w_e(\gamma)=|d_e(\gamma)|^{-2}), with
(d_e(\gamma)=(I+\gamma H)d_e) and
(H=(E_{xy}+E_{yx})/2). Direct differentiation gives
[
\left.\frac{dw_e}{d\gamma}\right|_0
=-\frac{2d_{e,x}d_{e,y}}{|d_e|^4}.
]
If an edge Hamiltonian is explicitly of the form
(H_{\rm link}(\gamma)=\sum_e w_e(\gamma)E_e), with the edge energies
(E_e) held fixed under this geometric variation, then
[
\Pi_{xy}^{Q6}
=-\frac1V\left.\frac{dH_{\rm link}}{d\gamma}\right|_0
=\frac{2}{V}\sum_e
\frac{d_{e,x}d_{e,y}}{|d_e|^4}E_e.
]
The code exposes this exact conditional derivative through
`inverse_square_edge_weight_shear_derivative` and
`q6_inverse_square_link_shear_response`. It does **not** silently identify
this inverse-square weighted edge Hamiltonian with the separate plaquette
Wilson action in `unified_discrete_action.py`. To combine them, the actual
(E_e) must be derived from the declared link action and its dependence on
the strained geometry must be included.
