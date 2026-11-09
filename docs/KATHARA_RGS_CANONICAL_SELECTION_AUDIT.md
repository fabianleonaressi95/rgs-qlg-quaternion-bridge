# Kathara–Q3 canonical-selection audit against the recovered RGS/Q6 geometry

**Date:** 2026-11-09  
**Status:** repository-level structural audit; the proposed RGS selection test is specified, but no numerical selection result is claimed.

## 1. Purpose and correction

The local comparison of the two Kathara-to-Q3 path conventions (bit-0-first and bit-1-first) is not a sufficient basis for a global negative conclusion about canonical selection. In particular, a tie for
\[
E_{\rm local}(S)=\alpha\|S_1\|_1+\beta\|S_2\|_1
\]
only says that this functional does not distinguish those two representatives. It does not test the recovered RGS seed, its tangent constraints, the Q6 incidence geometry, or the existing icosahedral spectral sectors.

This note connects those existing constructions and states exactly what they do—and do not—select. It deliberately does not invent a new coordinate threshold or claim that a numerical result has been obtained from missing arrays.

## 2. Recovered RGS seed and its actual symmetry

The seed is not merely a hypothetical 64-point cloud. It is loaded by
`src.rgs_ellis_refinement.load_rgs_seed_data()`, which checks the archived hashes. The provenance record is
`data/rgs_188d/ellis_curved_seed_manifest.json`; the source archive is A14, array `X_Oh_reconstructed`, shape 64×3, with recorded SHA-256
`6a158cb30a671a5c994d8ab223e9283b57ac703808041d46e67f286126e419cf`.

The six Q6 bits use the integer label
\[
g=\sum_{a=0}^{5}2^a b_a.
\]
The repository's Q6 geometric lookup is an exact 64-state lookup, not an affine assumption. The affine fit in `derive_q6_projection_from_rgs_seed` is diagnostic: a nonzero residual measures the non-affine part of the recovered embedding and must not replace the exact lookup.

The existing `src/rgs_q6_oh_action.py` constructs all 48 signed coordinate permutations of \(O_h\), matches their action on the recovered coordinates, and separately tests whether the induced permutations preserve the Q6 edges and square faces. This distinction is essential:
- the coordinate seed has a recovered \(O_h\) action;
- preservation of the independently labelled Q6 cells is a separate measured gate;
- one must not infer Q6-label equivariance from coordinate symmetry alone.

Thus the symmetry test in the previous local Kathara note must be replaced by an actual comparison with the recovered induced actions. In particular, a reflection of the planar Kathara drawing is not automatically an element of the RGS \(O_h\) action on Q6 labels.

## 3. The 188D tangent is a constraint tangent, not a canonical frame

The A14 construction has
\[
J_{\rm constraint}\in\mathbb R^{4\times192},\qquad
Q\in\mathbb R^{192\times188},\qquad
Q^TQ=I_{188},
\]
and the tangent projector
\[
P_T=QQ^T.
\]
The reported residuals are \(\|JQ\|_F\simeq1.24\times10^{-14}\) and
\(\|Q^TQ-I\|_F\simeq1.31\times10^{-14}\). These validate the finite constraint/tangent construction, not a preferred basis of the 188D space.

The distinction matters for canonicity: replacing \(Q\) by \(QU\), with \(U\in O(188)\), changes the frame but leaves \(P_T\) unchanged. Any proposed selection energy must therefore be written in terms of invariant projectors/operators (or transform covariantly), not depend on a particular SVD basis or column ordering.

## 4. Existing Q6-to-tangent operator

The repository already defines the vertex-major ambient space
\[
\mathcal H_{\rm amb}=\mathbb R^{64}\otimes\mathbb R^3,
\qquad L_{\rm amb}=L_{Q6}\otimes I_3,
\]
and the Galerkin compression
\[
L_{\rm tan}=Q^T L_{\rm amb}Q.
\]
The leakage
\[
\ell_T=\|(I-QQ^T)L_{\rm amb}Q\|_F
\]
distinguishes an invariant tangent subspace (zero leakage) from a mere compression (nonzero leakage). This operator is already implemented in `src/koopman_rgs_psu.py`.

This is a legitimate RGS-derived finite quadratic form on the 188D tangent. It is not, by itself, a functional on Kathara chain maps: a chain map \(S:C_*(K)\to C_*(Q3)\) has not been assigned a linear lift into \(\mathcal H_{\rm amb}\) or into the 188D tangent. Such a lift is required before evaluating \(L_{\rm tan}\) on A/B.

## 5. Existing non-arbitrary Q6 spectral sector

The module `src/rgs_q6_d5_phase_space.py` constructs the Coulomb-gauge Q6 link tangent
\[
\ker D_0^T\subset C^1(Q6),\qquad \dim=192-63=129,
\]
then selects the top coexact Hodge eigenspace at eigenvalue 12. Its multiplicity is five and its spectral projector is basis-independent. The 60 proper icosahedral rotations act on this space with character
\[
(5,1,-1,0,0),
\]
the irreducible real 5-dimensional \(A_5\) representation.

Its cotangent lift is a 10-real-dimensional phase space. The unit electric kinetic metric is an explicit convention, not a derived physical normalization. This construction is a canonical spectral subspace relative to the stated Q6 Hodge operator and symmetry action; it does not identify a Kathara centre or edge with one of those five modes.

A further existing module, `src/rgs_q6_d5_sphere_intertwiner.py`, constructs a unique (up to sign and scale, then normalized and sign-fixed by a deterministic convention) \(A_5\)-intertwiner
\[
\Lambda^2(V_5)\longrightarrow\mathbb R^3.
\]
It sends \(q\wedge p\) to a geometric three-vector, invariant under the oscillator \(SO(2)\) phase action and equivariant under \(A_5\). This is the strongest existing candidate for a symmetry-controlled geometric direction. But it acts on the D5 phase-space sector, not on Kathara's 12 centres or its 21 elementary edges. A map from Kathara chains into this sector is still missing.

## 6. What a valid RGS variational comparison must do

A selection functional should be independent of arbitrary bases and should compare complete chain maps. A future candidate can have the form
\[
\mathcal E_{\rm RGS}(S)=
\lambda_{\rm geom}\,\mathcal D_{\rm geom}(S)
+\lambda_{\rm tan}\,\mathcal D_{\rm tan}(S)
+\lambda_{\rm sym}\,\mathcal D_{\rm sym}(S)
+\lambda_{\rm inc}\,\mathcal D_{\rm inc}(S),
\]
but none of these terms may be filled by an ad hoc choice.

The terms require the following explicit data:
1. **Geometric lift:** a specified map from Kathara vertices/edges/cells to the recovered 64-state Q6 geometry (or a stated finite-dimensional interpolation into \(\mathcal H_{\rm amb}\)). The 12-to-64 assignment must be derived from incidence and geometry, not nearest-neighbour tie-breaking.
2. **Tangent compatibility:** a basis-independent defect using \(P_T\), for example a commutator/intertwining norm once the lift is defined. Do not use a particular \(Q\) frame as physical data.
3. **Symmetry compatibility:** the actual induced group actions on source and target, followed by an equivariance defect. Coordinate \(O_h\) and Q6-cell \(O_h\) must remain separate until the latter passes.
4. **Incidence compatibility:** the chain identities must be imposed exactly over integers; they are constraints, not soft penalties.
5. **Weights:** the \(\lambda\)'s must follow from a previously defined action/inner product or remain symbolic. Choosing them to break the A/B tie would simply reintroduce the arbitrariness under investigation.

The two existing path conventions A/B should be compared only after the same lift, target incidence ordering, metric and group action have been fixed. If the induced functional remains tied, the result is a genuine non-uniqueness for that RGS functional. If it separates them, the proof must show that the separation is invariant under changes of tangent frame, Q6 coordinate relabelling allowed by the group, and equivalent orientations.

## 7. Why the current data do not yet select A or B

The available repository objects have different domains:
- Kathara: a 12-vertex, 21-elementary-edge planar graph with a chosen 10-face filling;
- Q3: an 8-vertex, 12-edge, 6-square 2-skeleton;
- Q6: a 64-vertex, 192-edge, 240-square 2-skeleton;
- RGS seed: a 64-point embedding of Q6 labels in \(\mathbb R^3\);
- A14: a 188D tangent subspace in the flattened 64×3 coordinate space;
- Q6–D5: a 5D coexact eigenspace and its 10D cotangent lift.

The known Kathara–Q3 map is coordinate-defined by the half-plane tests in `src/kathara_q3_chain_map.py`. Its third trigram bit is fixed to zero, so it factors through a square. The two path choices are chain-homotopic as recorded in T061, but the map is not yet derived from the RGS seed or the A14 tangent.

There is no incidence-preserving map in the repository from all 12 Kathara centres/cells into Q6, nor a chain-level map from Kathara into the 5D A5 eigenspace or its cotangent lift. Without such a map, evaluating an RGS tangent or D5 spectral functional on A/B is undefined—not merely numerically unfinished.

The RGS symmetry therefore does **not** currently prove either uniqueness or non-uniqueness of a Kathara representative. It invalidates the earlier overstrong claim that local threshold asymmetry alone establishes a global impossibility theorem.

## 8. Exact next implementation gate

The next implementation must be staged and report failures rather than silently repair them:

**Gate A — provenance and symmetry**
- load the seed through `load_rgs_seed_data()` (hash checks mandatory);
- run `audit_rgs_q6_oh()`;
- record separately the number of group elements preserving Q6 edges and faces;
- recover the canonical Q6 incidence order and compare it with the seed edge order by endpoint sets and orientations.

**Gate B — tangent/spectral invariants**
- load the archived A14 \(Q\) and verify \(JQ=0\), \(Q^TQ=I\);
- compute \(P_T\), \(L_{\rm tan}\), and \(\ell_T\);
- construct the Q6 top-coexact projector and A5 character checks;
- verify the \(\Lambda^2(V_5)\to\mathbb R^3\) intertwiner residual.

**Gate C — missing map**
- define a source incidence map from Kathara's 12 vertices/21 edges/10 filled faces into Q6, with explicit integer matrices \(R_0,R_1,R_2\);
- check \(B_1^{Q6}R_1=R_0B_1^K\) and \(B_2^{Q6}R_2=R_1B_2^K\);
- specify whether the map is injective, a cellular embedding, or a quotient; do not conflate these;
- require equivariance under a stated subgroup only if the source action on Kathara is explicitly defined.

**Gate D — canonical comparison**
- lift the two Q3 chain maps through the same Q3-to-Q6 inclusion (lower and upper factors both audited);
- compose with the same Kathara-to-Q6 map only if the chain-level diagrams commute;
- evaluate a projector-based RGS functional and report all terms separately;
- perform frame-change and group-relabeling invariance checks;
- classify the outcome as unique, tied, or undefined because the required map does not exist.

## 9. Evidence ledger

| Statement | Status |
|---|---|
| 64-point seed and coordinate/coframe provenance | **RECOVERED / HASH-CHECKED BY LOADER** |
| 48-element coordinate \(O_h\) action | **IMPLEMENTED; tests are present** |
| Equivariance of the independently labelled Q6 cells under that action | **SEPARATE GATE; read audit output, do not assume** |
| A14 188D constraint tangent/projector | **REPORTED NUMERICALLY VERIFIED** |
| Q6 top coexact 5D eigenspace and irreducible \(A_5\) character | **IMPLEMENTED FINITE SPECTRAL CONSTRUCTION** |
| \(\Lambda^2(V_5)\to\mathbb R^3\) intertwiner | **IMPLEMENTED CANDIDATE WITH EXPLICIT RESIDUAL AUDIT** |
| Kathara-to-Q3 chain maps A/B | **CANDIDATE; local chain identities are separately audited** |
| RGS-derived map Kathara-to-Q6 / Kathara-to-D5 | **MISSING** |
| RGS variational functional selecting A or B | **UNDEFINED UNTIL THE LIFT EXISTS** |
| Unique canonical Kathara representative | **OPEN** |

## 10. Reproduction entry points

```bash
PYTHONPATH=. python -m pytest tests/test_rgs_q6_oh_action.py -q
PYTHONPATH=. python examples/rgs_q6_oh_equivariance_audit.py
PYTHONPATH=. python -m pytest tests/test_rgs_q6_induced_metric.py -q
PYTHONPATH=. python -m pytest tests/test_rgs_q6_d5_phase_space.py -q
```

For the A5 relative-direction map, inspect `src/rgs_q6_d5_sphere_intertwiner.py` and its dedicated tests. The Kathara chain-map and chain-homotopy records remain T060/T061 in the central theorem repository.

**No test execution is claimed by this documentation commit.** The commands above are reproduction instructions; numerical outputs must be generated in the repository environment and attached before upgrading any status label.
