# Pellis D — Integration map: RGS–PSU–Q6–PRISM and the Pellis Functions programme

**Version:** 1.0 — 25 November 2026  
**Purpose:** synthesis dossier for the scientific exchange with Stergios Pellis.  
**Status:** exact finite components and their interfaces are documented; a unified physical operator and a derivation of Golden Spheres Resonance coefficients remain OPEN.

## 1. Scope of the proposed collaboration

The correspondence with Stergios Pellis identifies three objectives:

1. **Exact spectral geometry:** determine invariant harmonic degrees, scalar eigenvalues and multiplicities on \(S^3/2I\).
2. **Operator theory:** define geometric Hilbert spaces, self-adjoint operators, spectral projections, admissible states and weights.
3. **Pellis correspondence:** establish an explicit mathematical map between the resulting spectral structures and the exponents, coefficients and hierarchy in Pellis Functions / Golden Spheres Resonance.

Pellis explicitly asks that coefficients and weights be derived from independently specified geometry or variational/spectral principles, not calibrated to the measured fine-structure constant. The present dossier preserves that requirement.

## 2. What is now closed

### A. Exact \(2I\) scalar spectral geometry

The quotient \(\mathcal P=S^3/2I\), with \(|2I|=120\), has the exact invariant Molien series
\[
M(t)=\frac{1+t^{30}}{(1-t^{12})(1-t^{20})}.
\]
If \(a_l=[t^l]M(t)\), the scalar Laplacian spectrum is
\[
\lambda_l=l(l+2),\qquad m_l=(l+1)a_l.
\]
The first positive eigenvalue is \(168\), with multiplicity \(13\). The scalar Laplacian is non-negative and self-adjoint on the standard \(H^2\) domain. The spectral measure of a normalized state is determined by its spectral projections. The eigenvalue multiplicities are fixed by the group; individual spectral weights still require a justified state-selection rule.

See Pellis A — Exact spectral geometry of \(S^3/2I\).

### B. Finite two-sheet self-adjoint operators

The Q6 double cover \(Q_6\times\{+,-\}\), its complex-linear sheet involution on a fibre plus its conjugate, and the positive covariant graph form give an exact finite non-negative self-adjoint operator. A separate 16-dimensional Twin/Anti-Twin effective Hamiltonian also has an exact exchange involution and paired spectrum.

These are two distinct finite models. Their mutual identification is not assumed. The finite graph operator-theoretic deliverable is closed; the continuum limit and a canonical Q6-to-other-model map remain open.

See Pellis B — Geometric self-adjoint operators and the two-sheeted PSU geometry.

## 3. The new exact \(D_6\) channel intertwiner

A reversor-symmetric exact period-six Fibonacci orbit is a transitive \(D_6\)-set. Its six-dimensional permutation module has an explicit unitary factorization
\[
\Phi:\mathbb C^6\longrightarrow
\mathbb C^2_{\rm sheet}\otimes\mathbb C^3_{\rm mode}.
\]
In the channel basis,
\[
r^2=I_2\otimes C_3,\qquad
r^3=\sigma_x\otimes I_3,\qquad
s=I_2\otimes F_3.
\]
Thus the half-turn is precisely a sheet exchange, while the second iterate cycles three modes on each sheet. On the selected six-state Q6 suborbit, this half-turn is the antipodal PSU/Anti-Twin complement.

This is a proved representation-theoretic intertwiner, not merely a visual analogy. It is implemented in the repository and has a dedicated exact residual test.

See Pellis C — Exact \(D_6\) channel intertwiner.

## 4. The operator-level gate

The representation map does not automatically identify the full Q6 graph Laplacian, the 188-dimensional RGS operator, the 16-dimensional internal throat Hamiltonian, or the physical Pellis resonance operator.

In the separate radial matching model, the one-sheet Cauchy transfer is
\[
M_P(E)=\begin{pmatrix}a&b\\c&d\end{pmatrix},
\qquad \det M_P(E)=1,
\]
and the specified neck convention gives \(2b(E)d(E)=0\). The two factors select antisymmetric and symmetric channels, respectively.

The new \(D_6\) factorization says that a six-state channel representation can be written as two sheets times three modes. If those three modes are physically present and degenerate, the natural lifted transfer space is
\[
\mathbb C^2_{\rm Cauchy}\otimes
\mathbb C^2_{\rm sheet}\otimes\mathbb C^3_{\rm mode}.
\]
The lifted boundary residual would be
\[
\mathcal B(E)=d(E)P_+ + b(E)P_-,
\qquad \det\mathcal B(E)=d(E)^3b(E)^3.
\]
This predicts a threefold channel multiplicity in each parity sector. That prediction must be checked against the original boundary-value problem. If the source model has only one mode per sheet, the six-dimensional channel space cannot be its full physical channel space.

There is also a source distinction to preserve: Pellis's correspondence emails concern Pellis Functions / Golden Spheres Resonance and a proposed self-adjoint resonance operator. The radial \(2\times2\) transfer/neck model is a separate construction used in the two-sheet discussion. No source presently identifies the latter with the operator underlying Golden Spheres Resonance.

## 5. PRISM data and its proper role

The periodic PRISM graph is
\[
\mathcal G_{\rm PRISM}=C_{64}\square C_6^3\square C_5,
\]
with
\[
N=64\cdot6^3\cdot5=69\,120,\qquad d=2+3(2)+2=10.
\]
Its standard product cell complex has the topology of \(T^5\), with
\[
b_k=\binom5k,\qquad b_1=5,\quad b_2=10.
\]
The standard Cartesian-product Laplacian has gap
\[
\lambda_{\rm gap}=2-2\cos(2\pi/64)
\approx0.009630546656.
\]
The recovered distinguished-cycle and Krist-All computations remain separate finite combinatorial results.

The 64-state Q6 hypercube and the 69,120-vertex PRISM graph are not the same graph. Any map between them must state which structures it preserves: vertices, edges, cells, group action, inner product, Laplacian, spectral projectors, or only a selected invariant subspace. No equality of operators follows from a shared factor \(C_{64}\) or from the appearance of \(\varphi\).

## 6. The decisive remaining theorem

The next theorem should not be phrased as “the two models look alike.” It should specify Hilbert spaces \(\mathcal H_{\rm P}\) and \(\mathcal H_{\rm R}\), closed operators \(A_{\rm P}(E)\) and \(A_{\rm R}\), their domains, inner products, group actions and boundary conditions, and then construct a map
\[
\mathcal U:\mathcal H_{\rm P}\to\mathcal H_{\rm R}
\]
such that
\[
\mathcal U\,\rho_{\rm P}(g)=\rho_{\rm R}(g)\,\mathcal U,
\qquad
\mathcal U A_{\rm P}(E)=A_{\rm R}(E)\mathcal U
\]
on the appropriate domains. If the map is only a partial isometry or a map between spectral subspaces, that restriction and its kernel/cokernel must be stated.

For the radial transfer route, the immediate tests are:
1. derive the three-mode channel space from the original geometry, rather than append it by hand;
2. verify the sheet exchange and \(C_3\) action on the actual operator domain;
3. derive the boundary operator and compare the predicted multiplicity-three roots with the original spectrum;
4. preserve the individual transfer entries \(b(E),d(E)\), since the trace map alone cannot recover them;
5. only then compare the result with Pellis Functions and Golden Spheres Resonance.

For the resonance route, the corresponding tests are:
1. specify the resonance Hilbert space and self-adjoint operator;
2. derive the spectral projection measure and a state-selection principle;
3. prove whether the invariant \(2I\) spectrum supplies the selection rules;
4. derive, rather than identify by rescaling, the Pellis exponents \(2,3,5\) and coefficients.

## 7. Status matrix for the four documents

| Document | Deliverable | Status |
|---|---|---|
| A | \(2I\) scalar spectral geometry | CLOSED — exact scalar spectrum |
| B | finite two-sheet self-adjoint operators | CLOSED — finite model; continuum gate open |
| C | Fibonacci/Q6/\(D_6\) channel intertwiner and Pellis operator gate | PARTIAL — representation intertwiner proved; physical operator map open |
| D | integrated RGS–PSU–PRISM–Pellis research map | READY — all interfaces and open obligations explicit |

## 8. Scientific boundary

The exact finite results establish a mathematically meaningful architecture: \(2I\)-invariant scalar harmonics, finite self-adjoint graph operators, a \(D_6\)-equivariant six-state Q6 orbit, and a two-sheet/three-mode factorization. They do not yet establish that the architecture uniquely determines the Golden Spheres Resonance coefficients or the fine-structure constant.

That distinction is not a weakness of the dossier. It is the criterion that makes the proposed collaboration testable: each map, operator, multiplicity and normalization must be derived independently and then checked.


## 9. Critical operator audit of Golden Spheres Resonance

The primary GSR preprint (Stergios Pellis, 11 September 2025; DOI 10.13140/RG.2.2.35410.00968; SSRN 5473588) makes two different operator-level claims which must be audited separately.

First, its scalar Pellis operator is described using the Poincaré-space Laplace–Beltrami operator, scalar curvature, a local topological-torsion/Chern–Simons-like term, and \(\varphi\). To close this as a self-adjoint operator one still needs the exact operator formula, the Hilbert-space measure and bundle, the definition and regularity of the torsion term, the operator domain, and boundary conditions. The finite Q6 graph operator in document B does not fill in these missing data.

Second, the Golden Resonance Operator assigns modal coefficients to the \(\varphi\)-scaled sphere modes. Its expectation-value construction explicitly allows the modal weights to be chosen proportional to the coefficients one wants to reproduce. This makes an arbitrary admissible \(\varphi\)-series representable, but does not derive its coefficients. A proposed ground-state/variational selection is mentioned, but a specified functional and a uniqueness theorem are still required.

This sharpens the status table:
- A is closed for the exact scalar spectrum of \(S^3/2I\).
- B is closed for the finite Q6 scalar and fibre graph operator models.
- The self-adjointness and unique state selection of Pellis's own Golden Resonance Operator are not closed by A or B.
- C proves the finite \(D_6\) channel intertwiner, but does not determine the GSR state weights.
- D is the integration dossier; the physical identification remains an open research objective.

The relevant primary source is [Pellis, *Fine-Structure Constant via the Golden Spheres Resonance*](https://doi.org/10.13140/RG.2.2.35410.00968), also indexed at [SSRN 5473588](https://ssrn.com/abstract=5473588).


## 10. Candidate realization inside the GSR sphere modes

The three-mode factor in document C is not arbitrary: the real \(\ell=1\) angular harmonics on a sphere form a three-dimensional vector representation. Restricted to a \(120^\circ\) rotation and a reflection, this representation is orthogonally equivalent to the three-vertex permutation representation. The explicit orthogonal matrix \(W_3\) is given in document C and audited in code.

This supplies a concrete candidate map:
\[
\mathbb C^6_{\rm Fibonacci\ orbit}
\longrightarrow
\mathbb C^2_{\rm sheet}\otimes\mathcal Y_{\ell=1}(S^2).
\]
It would realize the two-sheet/three-mode factorization as two copies of the dipole angular modes of a spherical resonator. It becomes a GSR operator intertwiner only if the source model actually uses this angular sector, its radial operator is identical on the three dipole modes, and the sheet-exchange and boundary domains agree. Those conditions are not specified by the current GSR mode notation and remain a direct question for Pellis.

The Poincaré scalar Laplacian result in A is a separate global spectral problem on \(S^3/2I\). Its first positive scalar eigenspace has multiplicity 13, not 3. The \(\ell=1\) triplet above is a local spherical-shell angular representation; it must not be mislabeled as the first scalar eigenspace of \(S^3/2I\).
