# Exact integer-lattice obstruction: Q6 Wilson fluxes versus the five D5 modes

## Theorem

Let (D_0\in\mathbb Z^{192\times64}) and
(D_1\in\mathbb Z^{240\times192}) be the oriented Q6 incidence and
square-boundary matrices. Put
\[
L_1=D_1^T D_1\in\mathbb Z^{192\times192}.
\]
Let (C\in\mathbb Z^{129\times192}) be the fundamental-cycle matrix:
its rows form a ℤ-basis of
\[
\mathcal Z_1=\ker(D_0^T)\cap\mathbb Z^{192}.
\]
Let (E_{12}\subset\mathbb R^{192}) be the top coexact eigenspace
\[
L_1v=12v,\qquad D_0^Tv=0,
\]
which has dimension five by the exact Q6 coexact 1-form spectrum
\(\lambda_k=2+2k\), with coexact multiplicities
\(6\binom5k-\binom6{k+1}\); at \(k=5\) this is five.

Then the five-mode map
\[
\Pi_5:\mathbb Z^{129}\longrightarrow E_{12},\qquad
c\longmapsto (cC)P_{12}
\]
has an integer kernel of rank exactly 124. Here (P_{12}) is the
orthogonal projector onto (E_{12}).

Moreover, for every \(\beta_t>0\), there exists a nonzero integer
cycle flux (n\in\mathcal Z_1\) such that
\[
nP_{12}=0,\qquad
\lambda_n(\beta_t)
=\prod_{e=1}^{192}\frac{I_{|n_e|}(\beta_t)}{I_0(\beta_t)}<1,
\]
whereas the vacuum has the same five-mode projection and transfer weight
\(\lambda_0=1\). Consequently, the full finite-\(\beta_t\) Wilson
transfer weight does **not** factor through the five D5 mode coordinates.

## Proof

1. **Rationality of the spectral projector.** The matrix (L_1) is
   integer symmetric. Its spectrum consists of (0,2,4,6,8,10,12)
   on this Q6 cochain sector. Therefore the projector onto eigenvalue
   12 is the rational polynomial
   \[
   P_{12}=p(L_1),\qquad
   p(x)=\frac{x(x-2)(x-4)(x-6)(x-8)(x-10)}{
   12\cdot10\cdot8\cdot6\cdot4\cdot2}.
   \]
   In particular (P_{12}\in\mathbb Q^{192\times192}). It has rank
   five. On the Coulomb/coexact subspace the (D_0D_0^T) term vanishes,
   so this is exactly the selected top coexact projector used by the
   five-mode reduction.

2. **The map has rational coefficients and rank five.** Since (C)
   is integral, (A=CP_{12}) is rational. The cycle rows span
   (\ker D_0^T\) over \(\mathbb R\), and (E_{12}\) is contained
   in that cycle space. Orthogonal projection from the cycle space onto
   (E_{12}) is therefore surjective: it fixes every vector in
   (E_{12}). Thus \(\operatorname{rank}_{\mathbb Q} A=5\).

3. **Apply the integer-kernel theorem.** Clear denominators in (A)
   to obtain an integer matrix (\widetilde A=dA) with the same
   kernel. For a rational matrix of rank five with 129 columns, its
   rational nullspace has dimension 124. A rational basis of that
   nullspace can be multiplied vector-by-vector by common denominators
   to give 124 linearly independent integer vectors. Hence
   \[
   \operatorname{rank}_{\mathbb Z}\ker(\widetilde A|_{\mathbb Z^{129}})
   =124.
   \]
   This is the standard lattice-kernel consequence of rational
   Gaussian elimination (equivalently, the structure theorem for
   subgroups of finitely generated free abelian groups).

4. **The invisible flux changes the Wilson weight.** Choose a nonzero
   integer vector (c\in\ker A\), and put (n=cC\). Since the
   fundamental-cycle rows are independent, (n\ne0), and by
   construction (nP_{12}=0). At least one component (n_e\ne0).
   For (x>0), (0<I_m(x)<I_0(x)) for every integer (m\ge1);
   this follows, for example, from the integral representation of
   (I_m) and strictness of \(|\cos(m\theta)|\le1\) on a set of
   positive measure. Every other normalized factor is at most one.
   Thus \(\lambda_n<1=\lambda_0\), although both fluxes project to
   the same five-mode point zero. The Wilson weight cannot be a
   function of those five coordinates alone. \(\square\)

## Consequence for the D5 Hua cutoff

This is stronger than a dimension-counting warning and stronger than a
numerical collision search: it proves the existence of infinitely many
nonzero integer Gauss-law fluxes invisible to the five-mode projection
but distinguishable by the exact finite-β Wilson transfer operator.

Therefore no scalar D5 generic norm (N(z)), and no cutoff
\(N(z)\ge\varepsilon\), can reproduce the full Wilson transfer
spectral selection using only the existing five-mode coordinates.
A relation \(\varepsilon=\varepsilon(\tau,\beta_t,a_t,\ldots)\)
requires extra data: a quotient retaining the discarded integer
directions, a conditional/averaged transfer weight on each fibre, or
a modified effective action derived by integrating those fibres out.
It cannot be obtained by simply identifying the two cutoff parameters.

This does **not** rule out a well-defined *effective* D5 measure after
a specified pushforward or fibre integration. It says that such a
measure is additional structure and must include the fibre prescription.

## Reproducibility and epistemic status

- **Exact theorem:** rational spectral projector; rank-five map; integer
  kernel rank 124; existence of a nonzero invisible integer flux with
  transfer weight different from vacuum.
- **Numerically audited in the repository:** orthonormal five-mode frame,
  projector/eigenspace residuals, and the 129-by-5 coordinate rank.
- **Still open:** explicit short integer kernel generators; canonical
  fibre measure/action; induced D5 density and its cutoff; physical
  continuum scaling and any fine-structure-constant derivation.

The proof uses the established Q6 coexact spectrum and exact incidence
matrices; it does not infer integer non-injectivity from real dimension
counting alone.
