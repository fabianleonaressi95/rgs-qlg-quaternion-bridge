# Omega-8 and Beneš-8: exhaustive finite routing audit

**Status:** NUMERICALLY VERIFIED by the committed exhaustive enumerator.  
**Model:** 8 wires labelled by three-bit integers; four independent 2×2 switches per stage.

## 1. Explicit topology convention

A stage acts on one address-bit dimension \(d\in\{0,1,2\}\). The four switches pair wires that differ only in bit \(d\). Each switch has one binary setting:

- 0: preserve bit \(d\) (straight);
- 1: flip bit \(d\) (cross).

Switch IDs are formed by packing the two wire bits other than \(d\), from low to high. The stage is therefore a product of four disjoint 2×2 switches.

The network is described in **shuffle-normalized coordinates**: relabelling wires after each perfect-shuffle interstage permutation turns the Omega stages into the dimension sequence below. This convention is fixed in code and avoids an implicit/ambiguous wiring diagram.

- Omega-8: dimensions \((0,1,2)\), 3 stages, 12 switches.
- Beneš-8: dimensions \((0,1,2,1,0)\), 5 stages, 20 switches.

The Beneš sequence is the symmetric butterfly / reverse-butterfly construction.

## 2. Exhaustive enumeration method

The implementation does not count raw switch configurations as if they were distinct permutations. It dynamically builds the set of distinct input-output permutations after each stage:

\[
\mathcal P_0=\{\mathrm{id}\},\qquad
\mathcal P_{s+1}=\{f_{d_s,c}\circ p:
p\in\mathcal P_s,\ c\in\{0,\ldots,15\}\}.
\]

Here \(f_{d,c}\) is the stage permutation for dimension \(d\) and four-bit switch mask \(c\). Duplicate permutations are removed at every stage. A witness (one sequence of stage masks) is retained for every reachable permutation.

This is exhaustive: every one of the 16 possible settings of each stage is applied to every distinct permutation reachable at the preceding stage.

## 3. Results

| Network | Stage dimensions | Switches | Distinct permutations by stage | Final routable |
|---|---|---:|---|---:|
| Omega-8 | 0, 1, 2 | 12 | 1, 16, 256, 4,096 | 4,096 / 40,320 |
| Beneš-8 | 0, 1, 2, 1, 0 | 20 | 1, 16, 256, 4,096, 9,216, 40,320 | 40,320 / 40,320 |

The 8×8 bit-reversal permutation is
\[
(0,4,2,6,1,5,3,7).
\]
It is not in the Omega-8 reachable set and is in the Beneš-8 reachable set.

**Interpretation:** within the declared stage model, Omega-8 is blocking and reaches exactly 4,096 distinct permutations; Beneš-8 is rearrangeable and reaches all \(8!=40,320\) permutations.

## 4. Reproduce

From repository root:

\`\`\`bash
PYTHONPATH=. python -m pytest tests/test_interconnection_networks_8x8.py -q
PYTHONPATH=. python examples/interconnection_networks_8x8_audit.py
\`\`\`

The audit writes \`results/interconnection_networks_8x8_audit.json\`.

## 5. Evidence and scope boundary

- The 4,096 and 40,320 counts are now reproduced from the explicit model, rather than merely copied from the historical session.
- The result counts distinct permutations, not the number of switch-setting strings; multiple settings can realize the same permutation.
- A routing witness is returned for every permutation reachable in the network.
- This is a finite combinatorial result. It does not establish an equivalence between Q6 and a Beneš/Omega switching fabric, nor a Green–Kubo transport theorem, nor a physical implementation performance claim.
- The bridge to Q3/Q6 requires an explicit state/routing map and separate preservation tests.
