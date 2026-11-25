# Pellis B — Geometric self-adjoint operators and the two-sheeted PSU geometry

**Version:** 1.0 — 25 November 2026  
**Status:** CLOSED for the explicitly defined finite models and their self-adjointness. A canonical identification of all Q6 states with the separate internal Twin/Anti-Twin space, and a two-sheeted continuum limit, remain OPEN.

## 1. Scope and spaces

This note separates three objects which must not be conflated:

1. the 64-vertex Q6 graph with a rank-two \(SU(2)\) fibre and its conjugate on the second sheet;
2. a separate finite-dimensional internal Twin/Anti-Twin exchange algebra on \(\mathbb C^8\), used in a 16-dimensional effective throat Hamiltonian;
3. a proposed smooth continuum double cover, which is not yet derived from the finite data.

Let \(Q_6=(V,E)\) be the six-dimensional hypercube graph, with \(|V|=64\), \(|E|=192\). Its two-sheeted lift is
\[
\widetilde Q_6=Q_6\times\{+,-\}.
\]
The deck involution is \(\tau(g,+)=(g,-)\), \(\tau(g,-)=(g,+)\).

For scalar fields,
\[
\mathcal H_0=\ell^2(V)\oplus\ell^2(V).
\]
For the archived rank-two fibre connection,
\[
\mathcal H_{\rm fib}
=\ell^2(V;\mathbb C^2)\oplus\ell^2(V;\overline{\mathbb C^2}).
\]
The conjugate vector space is essential: it makes the sheet exchange complex-linear on the direct sum.

## 2. Sheet involution and parity sectors

Let \(c:\mathbb C^2\to\overline{\mathbb C^2}\) be the canonical complex-linear identification. Define
\[
J(\psi,\bar\chi)=\big(c^{-1}(\bar\chi),c(\psi)\big).
\]
Then
\[
J^*=J,\qquad J^2=I.
\]
The orthogonal projections
\[
P_\pm=\frac12(I\pm J)
\]
give
\[
\mathcal H_{\rm fib}=\ker(J-I)\oplus\ker(J+I).
\]
The two parity sectors have equal dimension. This is an exact property of the defined finite direct-sum construction.

## 3. Covariant graph form and self-adjoint operator

Let \(U_e\in SU(2)\) be the first-sheet transport and \(\bar U_e\) its conjugate on the second sheet. For an oriented edge \(e=(g,h)\), set
\[
(D_U\Psi)_e=
\left(\psi_h-U_e\psi_g,\;
\bar\chi_h-\bar U_e\bar\chi_g\right).
\]
For positive edge coefficients \(\kappa_e>0\), define
\[
q_U[\Psi]=\sum_{e\in E}\kappa_e
\left(\|\psi_h-U_e\psi_g\|^2+
\|\bar\chi_h-\bar U_e\bar\chi_g\|^2\right).
\]
On a finite graph this form is defined on the entire finite-dimensional Hilbert space, is non-negative, and is closed. Its associated operator is
\[
L_U=D_U^*WD_U,\qquad W=\operatorname{diag}(\kappa_e I).
\]
Therefore
\[
L_U=L_U^*\ge0,\qquad [L_U,J]=0.
\]
The even and odd sheet sectors reduce \(L_U\).

The baseline vertex measure is counting measure. The recovered embedding gives edge lengths \(\ell_e=\|X_h-X_g\|\), but a canonical dual-cell vertex volume has not been fixed. No geometric volume measure is silently inserted into the Hilbert inner product. If positive vertex weights are later selected, the adjoint and operator must be recomputed in the weighted Hilbert space.

The convention \(\kappa_e=\ell_e^{-1}\) is the one-dimensional P1 edge-stiffness convention. It is a declared graph Dirichlet-energy model, not yet a uniquely derived four-dimensional PSU action.

## 4. Separate internal Twin/Anti-Twin exchange model

A distinct recovered finite archive contains \(K\cong\mathbb C^8\), a Hermitian involution \(E\), and Hermitian observables \(O_a\) satisfying
\[
E^*=E,\quad E^2=I_8,\qquad
O_a^*=O_a,\quad O_a^2=I_8,\quad EO_aE=-O_a.
\]
On \(\mathbb C^2_{\rm sheet}\otimes K\), define
\[
H(t,\delta)=-t\,\sigma_x\otimes I_8
+\delta\,\sigma_z\otimes O_a,\qquad t>0,\quad\delta\in\mathbb R,
\]
and
\[
\mathcal J=\sigma_x\otimes E.
\]
Direct multiplication gives
\[
\mathcal J^*=\mathcal J,\quad\mathcal J^2=I_{16},
\quad[H,\mathcal J]=0,
\qquad H^2=(t^2+\delta^2)I_{16}.
\]
Consequently
\[
\operatorname{spec}(H)=
\{-\sqrt{t^2+\delta^2},+\sqrt{t^2+\delta^2}\},
\]
each eigenvalue with multiplicity eight. The lower branch has expansion
\[
E_-=-t-\frac{\delta^2}{2t}+O(\delta^4/t^3).
\]
This is an exact algebraic statement for the declared finite model. It is not a derivation of \(t\), a physical throat equation, or a continuum Dirac operator.

The 16-dimensional model is not the 256-dimensional rank-two Q6 fibre space. In particular, no canonical equivariant map from the 64 Q6 vertices to the 8-dimensional internal space has been established.

## 5. Continuum self-adjointness template

For a compact smooth Riemannian manifold \((M,g)\) without boundary and a Hermitian bundle \(E\to M\) with unitary connection, the connection Laplacian on
\[
L^2(M,E\oplus\bar E;d\mu_g)
\]
is non-negative and self-adjoint on its standard \(H^2\) domain. If a unitary bundle involution preserves the metric, measure, connection and domain, it commutes with the operator and its parity sectors reduce it.

For a manifold with boundary, the boundary conditions and the corresponding domain must be specified; self-adjointness cannot be claimed from the differential expression alone.

This is a standard continuum construction conditional on the existence of the smooth bundle, metric, connection and deck involution. Those data have not yet been derived as a limit of the Q6 seed.

## 6. Continuum-limit gate

The fixed 64-vertex graph has no continuum limit by itself. A theorem requires a refinement sequence \(G_h\), a fixed target Hilbert space, reconstruction maps, compatible measures, and closed forms \(q_h\). One must prove Mosco liminf and recovery, then deduce strong-resolvent convergence under the appropriate hypotheses. The discrete sheet involutions must also intertwine with the continuum involution under reconstruction.

The existing curved Ellis–RGS scalar refinement programme is a separate result and does not supply a spinorial two-sheeted Mosco theorem.

## 7. Status of deliverable B

**Closed and ready to send**
- the finite two-sheet Hilbert spaces and their inner products;
- the complex-linear sheet involution and orthogonal parity projections;
- the finite covariant graph form and its non-negative self-adjoint operator;
- the exact 16-dimensional Twin/Anti-Twin Hamiltonian identities and spectrum;
- the conditional continuum self-adjointness template, with its hypotheses explicit.

**Still open**
- a canonical map from the 64-state Q6 space to the separate 8-dimensional internal space;
- a derived geometric vertex measure and constitutive edge weights;
- a continuum refinement sequence and two-sheeted Mosco/strong-resolvent theorem;
- identification of the finite graph operator with a unique physical PSU or Pellis resonance operator.

Thus B is closed as a finite operator-theoretic specification, not as a completed continuum physical theory.

## 8. Relation to PRISM

The comparison with PRISM must preserve the distinction between the 64-vertex Q6 seed and the periodic product graph
\[
C_{64}\square C_6^3\square C_5,
\]
which has \(69\,120\) vertices and degree \(10\). A map between their Hilbert spaces or operators must be written explicitly and checked for equivariance, inner products, domains and operator intertwining. Shared dimensions, labels or golden-ratio values do not constitute such a map.
