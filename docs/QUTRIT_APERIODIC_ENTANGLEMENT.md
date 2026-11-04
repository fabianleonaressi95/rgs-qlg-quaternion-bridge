# Qutrit entanglement on an aperiodic RGS–PSU substrate

## Mathematical definition

A qutrit is a local three-dimensional Hilbert space
\[
\mathcal H_i=\mathbb C^3,\qquad
|0\rangle,|1\rangle,|2\rangle.
\]
For two qutrits, the canonical maximally entangled state is
\[
|\Phi_3\rangle=\frac{1}{\sqrt3}
(|00\rangle+|11\rangle+|22\rangle),
\]
with Schmidt probabilities \((1/3,1/3,1/3)\) and entropy
\(S=\log_2 3\) bits.

## Aperiodic coupling sequence

The finite substrate uses the Fibonacci substitution
\[
A\mapsto AB,\qquad B\mapsto A.
\]
For an open chain of \(N\) qutrits, the bond sequence is the prefix of length
\(N-1\). The Hamiltonian is
\[
H_{\rm Fib}=\sum_{i=0}^{N-2}J_iP_{i,i+1},\qquad
J_i=\begin{cases}J_A&A,\\J_B&B,\end{cases}
\]
where \(P_{i,i+1}\) swaps adjacent \(\mathbb C^3\) factors. This is an
SU(3)-invariant permutation-exchange model; it is a well-defined finite
quantum-spin Hamiltonian, not a QCD colour Hamiltonian.

The code diagonalises the finite matrix and computes the Schmidt spectrum and
von Neumann entropy across a selected bipartition. It separately checks
finite-prefix periods. A finite word can never certify infinite aperiodicity;
the infinite Fibonacci fixed point is aperiodic by a separate combinatorial
argument.

## Relation to the existing RGS–PSU–Q6 programme

- The qutrit local space is distinct from the SU(3) flavour decomposition
  \(3\otimes3\otimes3=10\oplus8\oplus8\oplus1\): the former is one
  three-level degree of freedom; the latter is a three-particle representation.
- The SU(3) permutation chain is a controlled bridge between three-level
  quantum states and nonuniform/aperiodic couplings.
- Coupling this chain to Q6 cochains, barycentric geometry, spin connections,
  or the PSU twin sheets requires an explicitly defined tensor-product map
  and interaction term. No such identification is assumed here.
- SU(3) flavour \(SU(3)_F\) and QCD colour \(SU(3)_C\) must remain distinct.

## Evidence status and limits

- **PROVED (finite algebra):** \(|\Phi_3\rangle\) has three equal Schmidt
  coefficients and entropy \(\log_2 3\).
- **NUMERICALLY VERIFIED (finite model):** Hermiticity, spectrum, Schmidt
  probabilities and entanglement entropy for each specified finite chain.
- **CANDIDATE:** use of Fibonacci-modulated qutrit exchange as an aperiodic
  quantum sector in the broader RGS–PSU model.
- **OPEN:** coupling to Q6/PRISM incidence, PSU–Anti-Twin, spinorial currents,
  vierbein and a continuum limit; robustness of entanglement under size and
  boundary changes; spectral/transport signatures distinguishing the
  aperiodic model from periodic controls.

## Reproducibility

Run:
```bash
PYTHONPATH=. python -m pytest tests/test_qutrit_aperiodic.py -q
PYTHONPATH=. python examples/qutrit_aperiodic_audit.py
```

The dense exact diagonalisation is capped at a small number of sites. Larger
systems require sparse/Krylov or tensor-network methods and convergence tests.
