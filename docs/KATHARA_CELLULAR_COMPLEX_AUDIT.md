# Kathara cellular-complex audit

This note records the exact finite chain construction in src/kathara_cellular_complex.py.

## Input and convention

The 12-center / 15-primary-line graph and its planar embedding are the T058 object. Expanding the three polyline paths gives 21 elementary edges. We fill its ten bounded planar faces as 2-cells. This filling is a chosen CW-complex model; it is not asserted to be part of the source's physical ontology.

## Exact audit output

| Quantity | Result |
|---|---:|
| dim C0 | 12 |
| dim C1 | 21 |
| dim C2 | 10 |
| B1 shape | 12 x 21 |
| B2 shape | 21 x 10 |
| B1 B2 | exactly zero |
| rank B1 over Q | 11 |
| rank B2 over Q | 10 |
| b0 | 1 |
| b1 | 0 |
| b2 | 0 |
| Euler characteristic | 1 |

The graph-only cycle rank is 21-12+1=10. Attaching the ten bounded faces kills all ten independent graph cycles, yielding the homology of a disk.

## Hodge operators

With standard Euclidean cellular inner products:
- Delta0 = B1 B1^T, size 12 x 12;
- Delta1 = B1^T B1 + B2 B2^T, size 21 x 21;
- Delta2 = B2^T B2, size 10 x 10.

The harmonic kernel dimensions are exactly (1,0,0). All matrices are assembled with integer arithmetic. The audit's exact rank routine works over Q.

## Reproduction

From the repository root:

```bash
PYTHONPATH=. python -m pytest tests/test_kathara_15_line_graph.py tests/test_kathara_cellular_complex.py -q
PYTHONPATH=. python examples/kathara_cellular_complex_audit.py
```

The exporter writes the full integer B1, B2 and three Hodge matrices to results/kathara_cellular_complex/audit.json.

## Scientific boundary

This gives a rigorously defined cellular Hodge system on the filled planar Kathara graph. It does not prove that it equals the user's Krist-All operator, does not derive the supergolden transfer matrix, and does not establish a chain map to the Q3 trigram cube or Q6 cubical complex. Those require explicit incidence-preserving maps and separate proofs.
