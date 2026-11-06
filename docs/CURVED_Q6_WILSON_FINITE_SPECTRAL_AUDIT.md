# Finite curved-graph Wilson spectral audit

**Status:** finite-matrix diagnostic added; refinement-uniform and continuum claims remain OPEN.

## Purpose

The periodic Fourier-symbol audit answers a different question from a curved finite graph. On a curved graph with vertex-dependent coframes and edge spin transports, momentum is generally not a good quantum number. The correct immediate task is therefore to report the spectrum of the actual finite Wilson kernel and the singular spectrum of its Euclidean Wilson–Dirac matrix, without labelling eigenvectors as Brillouin doublers.

## Implementation

- `src/curved_wilson_spectral_audit.py`: calls the existing finite graph assemblers in `src/discrete_dirac.py`.
- `tests/test_curved_wilson_spectral_audit.py`: exercises a nonuniform coframe sample and invalid cutoff handling.

The audit records:
- all eigenvalues of the positive Wilson kernel through its minimum/maximum and multiplicities below a declared cutoff;
- singular-value extrema and count below the same cutoff for the Euclidean Wilson–Dirac matrix;
- Hermiticity and gamma5-Hermiticity residuals;
- unitarity residuals for explicitly supplied spin links;
- a machine-readable finite-check status and explicit limitations.

When spin links are omitted, the existing polar coframe prescription is used by the assembler. That prescription remains a candidate transport rule, not a proven discretization of continuum spin parallel transport.

## Run

```bash
PYTHONPATH=. python -m pytest tests/test_curved_wilson_spectral_audit.py -q
```

## Interpretation boundary

A successful report proves only that the specified finite matrices satisfy the reported algebraic/numerical checks at the chosen tolerance. In particular it does not establish:
- a Brillouin-zone or species count on curved Q6;
- a lower bound uniform over a refinement family;
- compactness of low-energy spectral projectors;
- Mosco convergence of a spinorial form;
- strong-resolvent convergence or convergence of an index.

## Next proof gate

To compare refinement levels, define a family `(K_h, E_h, U_h)` with:
1. explicit curved cells and edge lengths;
2. a spin connection and link transport defined from it, rather than only polar coframe matching;
3. compatible spinor Hilbert spaces and common-space maps;
4. a fixed physical spectral window and normalization convention;
5. mesh-uniform geometric and Wilson-kernel estimates.

Only after those data are fixed does it make sense to test convergence of spectral projectors and formulate a resolvent theorem.
