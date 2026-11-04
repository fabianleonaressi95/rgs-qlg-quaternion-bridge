# Qutrit–Q6 gauge-covariant embedding

## Embedding

Use the first vertices of the six-cube binary-reflected Gray path
\[
v_k=k\oplus(k\!\gg\!1),\qquad k=0,1,\ldots,N-1.
\]
Consecutive vertices differ in exactly one bit, hence define a genuine Q6
edge. The active qutrit Hilbert space is
\[
\mathcal H_N=\bigotimes_{k=0}^{N-1}\mathbb C^3.
\]
This is an active finite path in Q6, not a claim that all 64 vertices have
been diagonalised as qutrit sites (that Hilbert space has dimension 3^64).

## Aperiodic exchange and SU(3) links

The Fibonacci word supplies bond strengths J_k. Store a link
U_{v<-u} in SU(3) for each canonical edge u<v, with local gauge rule
\[
U_{v\leftarrow u}\mapsto g_vU_{v\leftarrow u}g_u^\dagger.
\]
For a path edge oriented from source u to target v, use its oriented
transporter U_{v<-u} (or its adjoint if the path traverses the canonical
edge backwards). The transported swap on the two endpoint qutrits is
\[
P_U=(U^\dagger\otimes U)P,
\]
where P exchanges the two tensor factors. It is Hermitian and transforms
covariantly under independent SU(3) transformations at the endpoints.
The Hamiltonian is
\[
H_{Q6,qutrit}=\sum_{k=0}^{N-2}J_kP_{U_k}.
\]
Consequently
\[
H[U^g]=G H[U]G^\dagger,\qquad
G=\bigotimes_k g_{v_k}.
\]
The implementation includes a numerical covariance residual for this
identity and a finite-size ground-state Schmidt audit.

## Relation to the RGS–PSU programme

- Q6 supplies the actual incidence graph and a concrete embedding of the
  qutrit sites.
- Fibonacci modulation supplies deterministic aperiodic bond data.
- SU(3) links make the exchange interaction locally gauge covariant.
- This SU(3) is a finite lattice internal symmetry construction. It is not
  yet the Yang–Mills plaquette action, QCD with dynamical quarks, or the
  SU(3)_F baryon octet/decuplet calculation.
- Coupling the qutrit path to the PSU two sheets, 188D vierbein, discrete
  spin transport and barycentric FEM requires additional explicit
  interaction terms. Those remain OPEN.

## Evidence labels

- **PROVED:** Gray-code consecutive vertices are adjacent in Q6; the
  transported-swap covariance follows algebraically from the link rule.
- **NUMERICALLY VERIFIED:** the implementation tests Hermiticity and
  \(\|H[U^g]-GH[U]G^\dagger\|_F\) for finite examples.
- **CANDIDATE:** the Fibonacci-modulated Q6 qutrit chain as an aperiodic
  quantum sector in RGS–PSU.
- **OPEN:** extension to the full Q6 graph, large-system tensor networks,
  dynamical SU(3) gauge action, coupling to PSU sheets/vierbein/spinors,
  and continuum or physical interpretation.

Run:
`PYTHONPATH=. python -m pytest tests/test_qutrit_q6_coupling.py -q`
