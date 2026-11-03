# Recovered RGS seed: induced O_h action and Q6 equivariance

## Purpose

Phase 18 checks the actual coordinate action on the recovered 64-point
RGS seed. It does not synthesize a point cloud or infer a nearest-neighbour
graph. The source is reconstructed by
`src.rgs_ellis_refinement.load_rgs_seed_data()`, which verifies the
archived SHA-256 hashes for the 64x3 coordinate array, the 192 Q6 edges
and the normalized four-component edge coframe.

## The finite group

The full cubic point group is the group of signed coordinate permutations:
[
O_h={Sin{-1,0,1}^{3	imes3}:S^TS=I, 
	ext{one nonzero entry per row and column}}.
]
There are (3!2^3=48) such matrices. For each (S), the audit maps
[
X_imapsto SX_i
]
to a unique coordinate row (X_{pi_S(i)}), using a KD-tree only to
locate the nearest row and then requiring a residual below tolerance,
a bijection, and a strict nearest/second-nearest gap.

It checks:
- all 48 matrices induce unique permutations of the 64 coordinates;
- the action is faithful;
- the induced permutations obey the group multiplication law;
- the central inversion (-I) is a fixed-point-free involution;
- the coordinate antipodal residual is zero to numerical tolerance.

## Q6 action is a separate test

The canonical Q6 2-skeleton has 64 vertices, 192 edges and 240 square
faces. The audit applies each induced coordinate permutation to the
actual Q6 edge and square sets and counts how many group elements preserve
each set.

This separation matters. The fact that the point cloud has an (O_h)
symmetry does **not** logically imply that an independently labelled
Q6 graph is equivariant under the same vertex action. If fewer than 48
elements preserve the edge set, the report marks the Q6 labelling as
non-equivariant; it does not silently relabel edges or claim a symmetry
that the data do not show.

The output fields `full_Oh_graph_automorphism` and
`full_Oh_cellular_automorphism` are therefore explicit measured gates.

## Evidence labels

- **NUMERICALLY VERIFIED:** the 48 signed-permutation matrices act as
  unique bijections on the recovered coordinate array, the action is
  faithful and respects composition, and central inversion is the
  coordinate antipode.
- **INDEPENDENT EQUIVARIANCE GATE:** preservation of the Q6 1- and
  2-cells is checked directly for all 48 induced permutations.
- **OPEN:** SU(2) link variables, plaquette holonomies and Wilson action
  have not yet been coupled equivariantly to this recovered action.

## Reproduction

```bash
PYTHONPATH=. python -m pytest tests/test_rgs_q6_oh_action.py -q
PYTHONPATH=. python -m src.rgs_q6_oh_action --output results/json/rgs_q6_oh_equivariance.json
```
