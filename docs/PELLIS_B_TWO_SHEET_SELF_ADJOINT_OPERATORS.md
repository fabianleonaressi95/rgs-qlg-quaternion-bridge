# Pellis B — Geometric self-adjoint operators and the two-sheeted PSU geometry

**Version:** 1.1 — 25 November 2026  
**Status:** CLOSED for the explicitly defined finite Q6 scalar and fibre graph operators. A canonical identification with the Pellis resonance operator and a two-sheeted continuum limit remain OPEN.

## 1. Scope and source objects

This note formalizes the finite two-sheet construction supported by the recovered Q6/QLG archive. It does not combine distinct models into one operator.

Let \(Q_6=(V,E)\) be the six-dimensional hypercube graph, with \(|V|=64\), \(|E|=192\), and canonical incidence data. Its two-sheet lift is
\[
\widetilde Q_6=Q_6\times\{+,-\}.
\]
The deck involution is
\[
\tau(g,+)=(g,-),\qquad \tau(g,-)=(g,+),\qquad \tau^2=I.
\]

For scalar fields, the Hilbert space with counting measure is
\[
\mathcal H_0=\ell^2(V)\oplus\ell^2(V).
\]
For the recovered rank-two fibre connection, use
\[
\mathcal H_{\rm fib}
=\ell^2(V;\mathbb C^2)\oplus
\ell^2(V;\overline{\mathbb C^2}).
\]
The conjugate vector space on the second sheet is part of the definition; it is not an independently chosen second complex structure.

## 2. Scalar graph operator

Let \(L_{Q_6}\) be the combinatorial graph Laplacian. On the double cover define
\[
L_0=L_{Q_6}\oplus L_{Q_6}.
\]
The sheet exchange
\[
\tau(\psi_+,\psi_-)=(\psi_-,\psi_+)
\]
is a self-adjoint unitary involution and satisfies
\[
[L_0,\tau]=0.
\]
The parity projections
\[
P_\pm=\frac12(I\pm\tau)
\]
are orthogonal and reduce \(L_0\). Each parity sector is unitarily equivalent to one copy of \(L_{Q_6}\).

For the unweighted six-cube,
\[
\operatorname{spec}(L_{Q_6})=\{2k:k=0,\ldots,6\},
\qquad
m_k=\binom6k.
\]
Thus the two-sheet operator has the same eigenvalues with doubled total multiplicities, while each sheet-parity sector carries one copy of the single-sheet spectrum. This is an exact finite graph statement.

## 3. Fibre involution

Let \(c:\mathbb C^2\to\overline{\mathbb C^2}\) be the canonical complex-linear identification with the conjugate vector space. Define
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
give the even and odd fibre-sheet sectors. They have equal dimension.

### 3.1 Explicit SU(2)-equivariant realization and uniqueness

Write the defining representation space as \\(V=\\mathbb C^2\\), and let \\(\\bar V\\) carry the conjugate representation \\(\\bar\\rho(U)=\\bar U\\). Set
\\[
\\varepsilon=i\\sigma_2=
\\begin{pmatrix}0&1\\\\-1&0\\end{pmatrix}.
\\]
For
\\[
U=\\begin{pmatrix}a&b\\\\-\\bar b&\\bar a\\end{pmatrix}
\\in SU(2),
\\qquad |a|^2+|b|^2=1,
\\]
direct multiplication gives
\\[
\\varepsilon U\\varepsilon^{-1}
=\\begin{pmatrix}\\bar a&\\bar b\\\\-b&a\\end{pmatrix}
=\\bar U,
\\qquad
\\varepsilon U=\\bar U\\varepsilon.
\\]
Define the complex-linear map to the conjugate representation by the matrix
\\[
c:V\\longrightarrow\\bar V,\\qquad c(v)=\\varepsilon v
\\]
in the chosen coordinate conventions. Then
\\[
c\\rho(U)=\\bar\\rho(U)c.
\\]
This is the required SU(2)-equivariance. It is not the set-theoretic map (v\\mapsto\\bar v): that coordinate-conjugation map is antilinear and is not being used as (c).

If (c':V\\to\\bar V) is another complex-linear intertwiner, then (A=c^{-1}c') commutes with the irreducible defining representation of (SU(2)). Schur's lemma gives (A=\\lambda I), hence (c'=\\lambda c). If both intertwiners are unitary, \\(|\\lambda|=1\\), so the ambiguity is a phase.

For the doubled fibre \\(V\\oplus\\bar V\\), define
\\[
J_V(\\psi,\\bar\\chi)
=\\big(c^{-1}\\bar\\chi,c\\psi\\big).
\\]
With the direct-sum Hermitian inner product and unitary (c), this satisfies
\\[
J_V^*=J_V,\\qquad J_V^2=I,
\\]
and the parity projections \\(P_\\pm=(I\\pm J_V)/2\\) are orthogonal.

## 4. Covariant graph form and self-adjoint operator

Let \(U_e\in SU(2)\) be the first-sheet transport and \(\bar U_e\) its conjugate on the second sheet. For an oriented edge \(e=(g,h)\), define
\[
(D_U\Psi)_e=
\left(\psi_h-U_e\psi_g,\;
\bar\chi_h-\bar U_e\bar\chi_g\right).
\]
The sheet exchange and conjugate fibres induce an exchange \(J_E\) on edge fields, with
\[
D_UJ=J_ED_U.
\]

For positive edge coefficients \(\kappa_e>0\), set
\[
q_U[\Psi]=\sum_{e\in E}\kappa_e
\left(\|\psi_h-U_e\psi_g\|^2+
\|\bar\chi_h-\bar U_e\bar\chi_g\|^2\right).
\]
On the finite graph this form is non-negative and defined on the whole finite-dimensional Hilbert space. Its associated operator is
\[
L_U=D_U^*WD_U,\qquad W=\operatorname{diag}(\kappa_e I).
\]
Hence
\[
L_U=L_U^*\ge0,\qquad [L_U,J]=0.
\]
The even and odd sheet sectors reduce \(L_U\).

The baseline vertex measure is counting measure. The recovered embedding supplies edge lengths \(\ell_e=\|X_h-X_g\|\), but a canonical dual-cell vertex volume has not been fixed. No geometric volume measure is silently inserted into the Hilbert inner product. If positive vertex weights are later selected, the adjoint and operator must be recomputed in the weighted Hilbert space.

The convention \(\kappa_e=\ell_e^{-1}\) is the one-dimensional P1 edge-stiffness convention. It is a declared graph Dirichlet-energy model, not yet a uniquely derived four-dimensional PSU action.

## 5. Continuum self-adjointness template

For a compact smooth Riemannian manifold \((M,g)\) without boundary and a Hermitian bundle \(E\to M\) with unitary connection, the connection Laplacian on
\[
L^2(M,E\oplus\bar E;d\mu_g)
\]
is non-negative and self-adjoint on its standard \(H^2\) domain. If a unitary bundle involution preserves the metric, measure, connection and domain, it commutes with the operator and its parity sectors reduce it.

For a manifold with boundary, the boundary conditions and corresponding operator domain must be specified; self-adjointness cannot be claimed from the differential expression alone.

This is a standard continuum construction conditional on the existence of the smooth bundle, metric, connection and deck involution. Those data have not yet been derived as a limit of the Q6 seed.

## 6. Continuum-limit gate

The fixed 64-vertex graph has no continuum limit by itself. A theorem requires a refinement sequence \(G_h\), a fixed target Hilbert space, reconstruction maps, compatible measures, and closed forms \(q_h\). One must prove Mosco liminf and recovery, then deduce strong-resolvent convergence under the appropriate hypotheses. The discrete sheet involutions must also intertwine with the continuum involution.

The existing curved Ellis–RGS scalar refinement programme is a separate result and does not supply a spinorial two-sheeted Mosco theorem.

## 7. Separate two-field Q-ball model

The repository also contains a flat-space two-field gauged Q-ball action with a diagonal \(U(1)\) and a discrete field-exchange symmetry, documented in:
- docs/PSU_ANTITWIN_TWO_FIELD_QBALL_TECHNICAL_NOTE.md
- docs/PSU_ANTITWIN_TWO_FIELD_QBALL_BVP.md
- src/rgs_psu_antitwin_qball_bvp.py

That is a separate nonlinear field-theory/BVP model. It is not the finite Q6 fibre Laplacian above, and it does not by itself provide a self-adjoint Pellis resonance operator or a \(D_6\)-channel intertwiner.

## 8. Status of deliverable B

**Closed and ready to send**
- the finite scalar two-sheet Hilbert space and graph Laplacian;
- the fibre Hilbert space with the conjugate second-sheet bundle;
- the complex-linear sheet involution and orthogonal parity projections;
- the finite covariant graph form and its non-negative self-adjoint operator;
- the precise counting measure baseline and the status of geometric vertex measures;
- the conditional continuum self-adjointness template.

**Still open**
- a canonical map from all 64 Q6 states to any separate internal space;
- a derived geometric vertex measure and unique constitutive edge weights;
- a continuum refinement sequence and two-sheeted Mosco/strong-resolvent theorem;
- identification of this finite graph operator with the Pellis/Golden Resonance operator.

Thus B is closed as a finite operator-theoretic specification, not as a completed continuum physical theory.

## 9. Relation to PRISM

The comparison with PRISM must preserve the distinction between the 64-vertex Q6 seed and the periodic product graph
\[
C_{64}\square C_6^3\square C_5,
\]
which has \(69\,120\) vertices and degree \(10\). A map between their Hilbert spaces or operators must state which structures it preserves: vertices, edges, cells, group action, inner product, Laplacian, spectral projectors, or only a selected invariant subspace. No equality of operators follows from shared dimensions, labels or golden-ratio values.
