# Pellis C — Exact \(D_6\) channel intertwiner and the operator-level correspondence gate

**Version:** 1.0 — 25 November 2026  
**Status:** the finite representation intertwiner is PROVED. The physical Pellis-operator intertwiner is CONDITIONAL/OPEN because the required three-mode channel space and its boundary action have not yet been identified in the source model.

## 1. The result: the six-orbit is a two-sheet, three-mode space

Let \(p\) be a reversor-symmetric point of exact period six for
\[
T(x,y,z)=(2xy-z,x,y),\qquad R(x,y,z)=(z,y,x).
\]
On the orbit \(\mathcal O(p)=\{T^kp:k=0,\ldots,5\}\), write \(e_k\) for the standard basis vector associated with \(T^kp\). The dihedral generators act by
\[
\rho(r)e_k=e_{k+1\pmod6},\qquad
\rho(s)e_k=e_{-k\pmod6}.
\]
They satisfy \(r^6=s^2=1\) and \(srs=r^{-1}\).

Define a six-dimensional channel space
\[
\mathcal K_{\rm ch}=\mathbb C^2_{\rm sheet}\otimes\mathbb C^3_{\rm mode}
\]
with orthonormal basis \(f_{\epsilon,j}=|\epsilon\rangle\otimes|j\rangle\), \(\epsilon=0,1\), \(j=0,1,2\). The explicit map is
\[
\boxed{
\Phi(e_{2j})=f_{0,j},\qquad
\Phi(e_{2j+3\ {\rm mod}\ 6})=f_{1,j}.
}
\]
Equivalently, the even orbit positions are sheet 0 and the odd positions are sheet 1, with the mode index chosen so that antipodal points have the same mode label. This is a permutation matrix, hence a unitary isomorphism.

Let \(C_3|j\rangle=|j+1\pmod3\rangle\) and \(F_3|j\rangle=|-j\pmod3\rangle\). In sheet-major ordering, the transported generators are
\[
\boxed{
\Phi\rho(r)\Phi^{-1}
=\begin{pmatrix}0&C_3^2\\ C_3^{-1}&0\end{pmatrix},
\qquad
\Phi\rho(s)\Phi^{-1}
=\begin{pmatrix}F_3&0\\0&F_3\end{pmatrix}.
}
\]
Consequently,
\[
\boxed{
\Phi\rho(r^2)\Phi^{-1}=I_2\otimes C_3,\qquad
\Phi\rho(r^3)\Phi^{-1}=\sigma_x\otimes I_3,\qquad
\Phi\rho(s)\Phi^{-1}=I_2\otimes F_3.
}
\]
These identities give the desired factorization without guessing a \(6\times6\) physical transfer matrix:

- \(r^2\) cycles the three modes within each sheet;
- \(r^3\) exchanges the two sheets, preserving the mode label;
- \(s\) reflects the three-mode triangle in each sheet.

The two eigenspaces of the sheet exchange are
\[
\mathcal K_\pm=\ker(\sigma_x\otimes I_3\mp I_6),
\qquad \dim\mathcal K_+=\dim\mathcal K_-=3,
\]
with orthogonal projections
\[
P_\pm=\frac12(I_6\pm\sigma_x\otimes I_3).
\]
These are symmetric/antisymmetric combinations across opposite vertices, not the coordinate spans of the two alternating triangles.

## 2. The three-mode characters

The cyclic permutation \(C_3\) is diagonalized by the unitary three-point Fourier transform. Its eigenvalues are
\[
1,\quad\omega,\quad\omega^2,\qquad \omega=e^{2\pi i/3}.
\]
The reflection \(F_3\) fixes the trivial character and exchanges the two nontrivial characters. Thus each sheet carries the regular three-point \(C_3\) representation, and the full six-state space carries the induced \(D_6\) action.

At the singular Fibonacci endpoint \(p_0=(0,0,0)\), the ambient matrix \(DT^2(p_0)\) has the same spectrum \(1,\omega,\omega^2\). This is a real \(C_3\)-representation comparison. The Fricke level is singular at \(p_0\), so it is not being described as a regular two-dimensional tangent Floquet plane.

## 3. Relation to the Q6 hexagram orbit

The six Q6 states
\[
000|111,\quad001|110,\quad011|100,\quad
111|000,\quad110|001,\quad100|011
\]
form a \(D_6\)-orbit under cyclic shifts and the chosen reflection of the six bit positions. The explicit orbit map in T150 identifies it with the Fibonacci six-orbit. On this suborbit the Q6 antipodal complement \(J\) is exactly \(r^3\). Therefore the channel involution \(\sigma_x\otimes I_3\) is the restriction of the Q6 PSU/Anti-Twin antipodal involution after the orbit identification.

This is an exact finite \(D_6\)-equivariant statement on six states. It does not make the six-state span invariant under the full Q6 hypercube Laplacian: the six-state set is an orbit under bit-position permutations, not a connected six-vertex subgraph of the Q6 1-skeleton.

## 4. The transfer-matrix question and the dimension obstruction

In the radial matching model considered separately, a one-sheet Cauchy transfer is
\[
M_P(E)=\begin{pmatrix}a(E)&b(E)\\c(E)&d(E)\end{pmatrix},
\qquad \det M_P(E)=1.
\]
The two-sheet neck convention used in that model gives
\[
2b(E)d(E)=0,
\]
with \(b(E)=0\) assigned to the antisymmetric channel and \(d(E)=0\) to the symmetric channel.

A raw \(2\times2\) Cauchy transfer cannot be unitarily equivalent to the six-dimensional \(D_6\) orbit representation: their dimensions and representation contents differ. The six-orbit supplies a channel factor \(\mathbb C^2_{\rm sheet}\otimes\mathbb C^3_{\rm mode}\), whereas \(M_P(E)\) acts on a two-dimensional Cauchy-data factor. The correct lifted space, if the three modes are genuine independent channels, is
\[
\mathcal H_{\rm lift}
=\mathbb C^2_{\rm Cauchy}\otimes
 \mathbb C^2_{\rm sheet}\otimes
 \mathbb C^3_{\rm mode},
\]
of dimension twelve.

If the same radial transfer acts identically on all six channels, its declared lift is
\[
\mathcal M_P(E)=M_P(E)\otimes I_6.
\]
Then
\[
\widehat\Phi=I_2\otimes\Phi
\]
is an exact unitary channel intertwiner:
\[
\widehat\Phi\,(M_P(E)\otimes I_6)
=(M_P(E)\otimes I_6)\,\widehat\Phi,
\]
where the right-hand side is written in the reordered Cauchy \(\otimes\) sheet \(\otimes\) mode basis. This is a mathematically exact tensor lift, conditional on the physical model actually possessing three degenerate modes per sheet. It is not evidence that such modes occur in Pellis's original resonance operator.

## 5. Boundary matching as a falsifiable multiplicity test

On the six-channel space, the natural parity-resolved boundary residual associated with the stated channel convention is
\[
\boxed{\mathcal B(E)=d(E)P_+ + b(E)P_-.}
\]
Its determinant is
\[
\det\mathcal B(E)=d(E)^3b(E)^3.
\]
Its zero set agrees with \(b(E)d(E)=0\), but each branch has a three-dimensional channel kernel when the corresponding scalar factor vanishes and the other does not. Thus the six-orbit lift predicts a triplication of the symmetric/antisymmetric channel multiplicity if the three modes are physically degenerate.

This is the decisive test, not a cosmetic resemblance:
1. derive the three-mode space from the physical Pellis geometry or operator;
2. derive its inner product and boundary domain;
3. verify that the boundary operator is \(d(E)P_+ + b(E)P_-\), or derive the correct replacement;
4. compare the predicted multiplicities with the original spectral problem.

If the physical model has only one mode per sheet, the six-dimensional representation cannot be its full channel space. One may retain the \(D_6\) orbit as a symmetry model, but must not call it a physical transfer intertwiner.

## 6. Do not conflate two Pellis research objects

The correspondence emails with Stergios Pellis concern his Pellis Functions / Golden Spheres Resonance programme, its coefficients, spectral weights, and a possible self-adjoint resonance operator. The radial \(2\times2\) transfer and neck condition above are a separate matching model used in the two-sheet discussion. No source examined here establishes that this transfer matrix is the operator underlying the Pellis Functions or Golden Spheres Resonance. The two lines must remain distinct until an explicit identification is supplied.

Pellis's 16 November correspondence asks for the two sheet spaces, their involution, the combined metric and measure, operator domains, a continuum-limit criterion, and the relation to PRISM. It does not itself specify the \(2\times2\) matrix \(M_P(E)\) or the three-mode channel space.

## 7. Evidence ledger

| Claim | Status |
|---|---|
| \(\Phi:\mathbb C^6\to\mathbb C^2\otimes\mathbb C^3\) is unitary | PROVED exactly |
| \(r^2\mapsto I_2\otimes C_3\) | PROVED exactly |
| \(r^3\mapsto\sigma_x\otimes I_3\) | PROVED exactly |
| \(s\mapsto I_2\otimes F_3\) | PROVED exactly |
| Q6 antipodal complement equals the half-turn on the selected six-state orbit | PROVED in T150 |
| The full Q6 Laplacian preserves that six-state span | FALSE in general; the selected orbit is not a Q6 induced connected component |
| A three-mode physical Pellis channel space exists | OPEN |
| The lifted transfer \(M_P(E)\otimes I_6\) is the actual Pellis operator | CONDITIONAL/OPEN |
| The threefold boundary multiplicity occurs physically | OPEN |
| The D6 construction derives Golden Spheres Resonance coefficients | OPEN |

## 8. Reproducibility

The exact permutation matrices and their intertwining identities are implemented in:
- src/fibonacci_d6_channel_intertwiner.py
- tests/test_fibonacci_d6_channel_intertwiner.py

Run:

    PYTHONPATH=. python -m pytest tests/test_fibonacci_d6_channel_intertwiner.py -q

The code audits unitarity, the \(D_6\) relations, and the \(r^2,r^3,s\) factorizations. It does not simulate or assert an unprovided physical Pellis transfer operator.


## 9. Critical audit of the published Golden Resonance Operator

The published preprint *Fine-Structure Constant via the Golden Spheres Resonance* (Stergios Pellis, 11 September 2025; DOI 10.13140/RG.2.2.35410.00968; SSRN 5473588) distinguishes a Pellis operator on scalar fields from a Golden Resonance Operator on a hierarchy of sphere modes.

The text describes the scalar Pellis operator using the Laplace–Beltrami operator on the Poincaré dodecahedral space, scalar curvature \(R(x)\), a term \(T(x)\) described as local topological torsion or Chern–Simons-like contribution, and \(\varphi\). However, the operator domain, precise bundle on which \(T\) acts, boundary conditions, and hypotheses guaranteeing self-adjointness must still be made explicit before this expression can be used as a closed self-adjoint spectral problem.

More decisively, the Golden Resonance Operator section introduces mode coefficients \(c_{n,m}\), and the subsequent expectation-value construction states that a convergent \(\varphi\)-series can be represented by choosing the modal weights \(|c_n|^2\) proportional to the desired coefficients and normalizing them. This proves representability, not unique spectral prediction: the target weights have been placed into the state. The paper proposes that a ground collective resonance might select the state by minimizing an energy functional with inter-shell coupling, but the functional, its domain, boundary conditions, and uniqueness/nondegeneracy conditions are not specified in the material examined here.

Therefore the operator-level Pellis correspondence requires an independent state-selection theorem. The D6 channel intertwiner derived above does not supply this missing variational principle and must not be presented as deriving the GSR coefficients.

Primary references:
- Stergios Pellis, *Fine-Structure Constant via the Golden Spheres Resonance*, 11 September 2025, DOI: https://doi.org/10.13140/RG.2.2.35410.00968
- SSRN record: https://ssrn.com/abstract=5473588


## 10. A concrete spherical-mode realization: the dipole sector \(\ell=1\)

There is a natural geometric realization of the three-mode representation if the standing-wave modes on each GSR sphere include the angular dipole sector. Let
\[
u_0=\frac{(1,1,1)}{\sqrt3},\qquad
u_1=\frac{(2,-1,-1)}{\sqrt6},\qquad
u_2=\frac{(0,1,-1)}{\sqrt2}.
\]
The orthogonal change of basis
\[
W_3=
\begin{pmatrix}
1/\sqrt3&1/\sqrt3&1/\sqrt3\\
2/\sqrt6&-1/\sqrt6&-1/\sqrt6\\
0&1/\sqrt2&-1/\sqrt2
\end{pmatrix}
\]
sends the three vertex modes to the invariant axis and the two-dimensional angular plane. It satisfies
\[
W_3C_3W_3^T=
\begin{pmatrix}
1&0&0\\
0&-\frac12&-\frac{\sqrt3}{2}\\
0&\frac{\sqrt3}{2}&-\frac12
\end{pmatrix},
\qquad
W_3F_3W_3^T=\operatorname{diag}(1,1,-1).
\]
This is precisely the real \(\ell=1\) spherical-harmonic representation restricted to a \(120^\circ\) rotation about an axis and a reflection through a plane containing that axis: the invariant coordinate corresponds to the zonal dipole, and the remaining two coordinates to the real \(m=\pm1\) pair.

Therefore
\[
\Phi_{\rm dipole}=(I_2\otimes W_3)\Phi
\]
is an explicit unitary map from the Fibonacci six-orbit to two copies of the \(\ell=1\) angular mode space of a spherical resonator. It is a concrete candidate bridge to the GSR nested-sphere modes.

**Physical condition still to verify:** the GSR source must identify its mode label \(m\) with angular spherical harmonics and must select the \(\ell=1\) sector (or another three-dimensional sector carrying this same \(D_3\) representation). If its mode label is radial, or if the relevant Poincaré-space operator has no such invariant three-dimensional subspace, this identification does not follow. The exact \(D_6\) representation map remains valid independently of that physical choice.

The orthogonal matrix \(W_3\) and the two generator identities are implemented and tested by the companion code.


## 11. Exact obstruction to restricting the full Q6 Laplacian to the six-state orbit

Write the selected six hexagrams as six-bit strings:
\[
000111,\ 001110,\ 011100,\ 111000,\ 110001,\ 100011.
\]
Every pair of distinct strings has Hamming distance \(2\), \(4\), or \(6\). Hence the induced subgraph on these six vertices has no edges, while each vertex still has all six hypercube neighbours outside the selected set.

Let \(P_{\mathcal O}\) be the coordinate projection onto their six-dimensional span and \(A_{Q_6}\) the full hypercube adjacency matrix. Then
\[
P_{\mathcal O}A_{Q_6}P_{\mathcal O}=0,\qquad
A_{Q_6}P_{\mathcal O}\ne0,
\]
and the exact Frobenius residuals are
\[
\|A_{Q_6}P_{\mathcal O}\|_F=6,\qquad
\|[A_{Q_6},P_{\mathcal O}]\|_F=6\sqrt2.
\]
Since \(L_{Q_6}=6I-A_{Q_6}\), the same commutator norm holds for the graph Laplacian. Thus the six-state orbit is invariant under the bit-position \(D_6\) action, but it is not a reducing subspace of the full Q6 adjacency or Laplacian. This is an exact obstruction to claiming that the orbit-level intertwiner already intertwines the full Q6 graph operator.


## 12. Finite operator-domain audit for the conditional radial lift

The representation intertwiner can now be carried one step further without assuming a physical GSR identification. Let \(M(E)\in\mathbb C^{2\times2}\) be any finite Cauchy transfer matrix. On the declared twelve-dimensional tensor product
\[
\mathcal H_{\rm lift}=\mathbb C^2_{\rm Cauchy}\otimes\mathbb C^2_{\rm sheet}\otimes\mathbb C^3_{\rm mode},
\]
define the *conditional* lift
\[
\widehat M(E)=M(E)\otimes I_6.
\]
The finite-dimensional domain is all of \(\mathcal H_{\rm lift}\), so it is invariant under the lifted D6 action. Since the D6 generators act on the channel factor only,
\[
[\widehat M(E),I_2\otimes\rho_{D_6}(r)]=
[\widehat M(E),I_2\otimes\rho_{D_6}(s)]=0
\]
identically, for arbitrary entries of \(M(E)\). This is an exact algebraic covariance statement for the tensor-product model; it is not a proof that the physical radial transfer is mode-independent.

On the six-channel factor, let \(S=\rho_{D_6}(r^3)=\sigma_x\otimes I_3\) and \(P_\pm=(I_6\pm S)/2\). The boundary residual
\[
B(E)=d(E)P_+ + b(E)P_-
\]
commutes with both D6 generators because the half-turn is central in D6. If \(b(E)=0\), \(d(E)\ne0\), its kernel is the three-dimensional antisymmetric sector; if \(d(E)=0\), \(b(E)\ne0\), its kernel is the three-dimensional symmetric sector. If both vanish, the entire six-dimensional channel space is in the kernel. For real \(b,d\), \(B(E)\) is self-adjoint.

The implementation is \`src/pellis_d6_operator_gate.py\`; tests are in \`tests/test_pellis_d6_operator_gate.py\`. Run:

    PYTHONPATH=. python -m pytest tests/test_pellis_d6_operator_gate.py tests/test_fibonacci_d6_channel_intertwiner.py -q

**Scope:** this closes the finite tensor-lift covariance and its algebraic kernel dimensions. The mode-degeneracy assumption, the physical boundary domain, and the identification of \(M(E)\) with Pellis/GSR remain OPEN.


## 13. The full D6 commutant: symmetry does not force a triplet degeneracy

The previous section's lift \(M(E)\otimes I_6\) is D6-covariant, but it is a *special isotropic case*, not the most general D6-covariant radial transfer. On the three-mode factor, the permutation representation decomposes as
\[
\mathbb C^3=\mathbb C u\oplus u^\perp,\qquad
u=(1,1,1)/\sqrt3,
\]
where \(\mathbb C u\) is the trivial one-dimensional representation and \(u^\perp\) is the irreducible two-dimensional representation of \(D_3\). The corresponding orthogonal projectors are
\[
P_0=uu^*,\qquad P_2=I_3-P_0.
\]
By Schur's lemma, the commutant on the mode factor is
\[
\operatorname{Comm}(D_3)=\{aP_0+bP_2:a,b\in\mathbb C\}.
\]
Consequently D6 covariance permits two independent radial Cauchy matrices:
\[
\mathbb M(E)=M_\parallel(E)\otimes I_{\rm sheet}\otimes P_0
+M_\perp(E)\otimes I_{\rm sheet}\otimes P_2.
\]
This operator commutes with the actual six-channel D6 generators. The equality \(M_\parallel(E)=M_\perp(E)\), which collapses it to \(M(E)\otimes I_6\), is an additional isotropy/degeneracy hypothesis; it does **not** follow from D6 symmetry alone.

Thus the representation-theoretic prediction is, in general, a singlet plus a doublet, not a triplet of identical modes. A physical threefold degeneracy requires an extra symmetry or a dynamical argument forcing equality of the two radial blocks. The new implementation is src/pellis_d6_commutant.py; its tests explicitly exhibit D6-covariant examples with unequal blocks. Run:

    PYTHONPATH=. python -m pytest tests/test_pellis_d6_commutant.py tests/test_pellis_d6_operator_gate.py tests/test_fibonacci_d6_channel_intertwiner.py -q

**Status correction:** the former M(E) tensor I_6 result remains a valid conditional special case. The general D6-covariant family is now the two-block expression above. The physical mode content and any additional degeneracy-producing symmetry remain OPEN.
