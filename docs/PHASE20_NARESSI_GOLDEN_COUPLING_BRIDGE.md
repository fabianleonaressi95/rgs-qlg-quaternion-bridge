# Phase 20 insert — Naressi golden coupling bridge

This is an extension of the existing Phase 20 insert in the canonical 27-chapter monograph. It does not create a new chapter or assert completed physical unification.

## Golden coupling

Define
\[
q_\phi=\phi^{-1}e^{-1/\phi^2},\quad
J_n=\alpha q_\phi^n,\quad
\theta_\phi=2\pi/\phi^2,\quad
\chi_n=J_n\cos\theta_\phi.
\]
At the golden ratio and alpha=1, n=1, q_phi=J_1=0.4218194769 and chi_1 is approximately -0.3110. The constant ratio J_(n+1)/J_n=q_phi follows exactly from the definition.

## QAT hermiticity

The capsule expression with b but without its adjoint is not Hermitian. The audit therefore implements the explicitly corrected finite-Fock Hamiltonian
\[
H_{QAT}=\omega_p n_a+\omega_f n_b+
\chi(n_an_b+b+b^\dagger).
\]
It checks the Hermiticity residual, propagator unitarity and state norm. This is a corrected model for audit, not a claim that the original capsule used this correction. Fock-cutoff convergence remains to be checked.

## TFIM

The open-chain sparse model is
\[
H_N=-\sum_{n=1}^{N-1}\alpha q_\phi^n
\sigma_n^x\sigma_{n+1}^x-g\sum_{n=1}^N\sigma_n^z.
\]
The CSR implementation uses big-endian site ordering and is restricted to N<=24. It returns E0, E1, gap and half-chain entropy. MPS/DMRG is not implemented or verified by this module.

## Explicit Q6 spectral dictionary

The Q6 heat-kernel reduced state has entropy
\[
S_{Q6}(\tau)=6h\left(\frac{1}{1+e^{2\tau}}\right).
\]
Identifying the Naressi decay ratio with the Q6 heat-kernel level-amplitude ratio gives
\[
q_\phi=e^{-\tau_\phi}
\quad\Rightarrow\quad
\boxed{\tau_\phi=-\log q_\phi=\log\phi+\phi^{-2}.}
\]
At phi golden, tau=0.8631778363 and S_Q6=2.5472014005 nats.

This is a parameter-free spectral dictionary, but remains a model ansatz. It is not yet derived from the QAT tetrad/action, and it does not prove that the TFIM ground-state reduced density matrix equals the Q6 heat-kernel reduced state. An intertwiner between these states/operators is still required.

## Evidence ledger

**Exact from definitions:** golden decay ratio, phase identity, chi projection, Q6 entropy evaluated at tau_phi.

**Numerically audited by the added test suite:** corrected QAT finite-Fock Hermiticity/unitarity; CSR TFIM symmetry, low-energy spectrum and half-chain entropy for tested N. The test run for these new commits is pending confirmation.

**Model choices:** displacement quadrature b+b†, open-chain boundaries, spectral dictionary q_phi=exp(-tau).

**Open:** derive tau from the QAT coframe/action; construct a TFIM-to-Q6 state/operator map; derive a geometric area law; implement and converge MPS/DMRG; compare against uniform, random and non-golden exponential controls.

## Reproduction

Run the focused tests:
PYTHONPATH=. python -m pytest tests/test_naressi_golden_bridge.py -q

Generate the JSON audit (default sizes N=4,8,12,16):
PYTHONPATH=. python examples/naressi_golden_bridge_audit.py
