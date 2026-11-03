# Phase 27 — Koopman bridge for RGS–PSU–QLG

This phase adds a reproducible Koopman/EDMD layer to the same RGS–PSU–PRISM–QLG
monograph. It does not introduce a separate physical theory.

## 1. Operator conventions

For a row-stochastic transition matrix \(P\), the finite Koopman operator acts
on observables by
\[
(Kf)(i)=\sum_jP_{ij}f(j),\qquad K=P.
\]
Probability columns evolve under the adjoint Markov action
\[
\mu_{n+1}=P^T\mu_n.
\]
For a strictly positive invariant law \(\mu\), the Hilbert-space adjoint in
\(L^2(\mu)\) is
\[
K^*=D_\mu^{-1}P^TD_\mu.
\]
The implementation audits row stochasticity, stationarity, preservation of
the constant observable, the weighted adjoint identity and the spectral radius.

## 2. Cantor–Parry–Q6 symbolic dynamics

The state space has 192 states \((i,g)\), with \(i\in\{0,1,2\}\) the
supergolden Parry automaton state and \(g\in\{0,\ldots,63\}\) the Q6 hexagram.
The canonical bit convention is \(g=\sum_{a=0}^5 2^a b_a\).
For each allowed Parry edge \(i\to j\), an emitted bit \(b\) is drawn using
the golden-weighted probabilities already defined in the unified operator;
it flips Q6 bit \((i+j)\bmod 6\). This gives a directed row-stochastic chain
before the stationary-flow symmetrization used by the static Laplacian.

The invariant distribution is \(\mu(i,g)=\pi_i/64\), where \(\pi\) is the
Parry stationary distribution. The uniform Q6 fibre is invariant because
each bit flip is a permutation of the 64 hexagrams. Tests check this identity
numerically to machine precision.

This finite Markov model is the six-symbol-cylinder truncation of the
symbolic construction; it is not itself the infinite Cantor shift.

## 3. Barycenter, covariance and Moiré observables

For particle/configuration positions \(X_i\) and normalized weights \(w_i\),
\[
X_B=\sum_iw_iX_i,\qquad
\Sigma_B=\sum_iw_i(X_i-X_B)(X_i-X_B)^T.
\]
The feature dictionary contains the barycenter, independent entries of the
covariance, cosine and sine of the barycentric relative phase, the real and
imaginary parts of the centered characteristic function, and the exact
averaged Moiré potential. The characteristic identity is retained exactly;
the covariance truncation is not substituted for it.

These features can be passed to EDMD using consecutive snapshots. For row-wise
feature matrices \(\Phi_0,\Phi_1\), EDMD solves
\[
\widehat K=\arg\min_K\|W^{1/2}(\Phi_0K-\Phi_1)\|_F^2
+\lambda\|K\|_F^2.
\]
The reported residual and rank are essential: EDMD is a finite Galerkin
approximation and an arbitrary feature dictionary is not guaranteed to be an
invariant Koopman subspace.

## 4. Doppler phase generator

For a constant relative angular frequency \(\Omega\), use the observable pair
\((\cos\theta,\sin\theta)\), with \(\dot\theta=\Omega\). In the row-feature
convention the generator is
\[
A_\Omega=\begin{pmatrix}0&\Omega\\-\Omega&0\end{pmatrix}.
\]
It is skew-symmetric and its semigroup is an exact rotation. This separates
the temporal Doppler generator from the spatial Moiré potential; dimensional
frequencies are not added to the dimensionless Q6 Laplacian.

## 5. RGS coupling and provenance boundary

The module intentionally does not invent a projection from the Q6 bits into
the 188D RGS tangent space. The existing affine estimator
derive_q6_projection_from_rgs_seed can be used when the original 64x3
\(O_h\) seed is present. The repository currently has the reconstruction
manifest and A14 report, but the binary \(X_{O_h}\) array has not been
verified in the tracked tree. Until that source is restored and the vierbein/
connection are supplied, no physical RGS Koopman generator or numerical
188D projection is claimed.

Likewise, the finite Parry x Q6 chain is a symbolic stochastic model, not a
claim that physical matter follows that Markov chain. The next analytic
obligation is to define the coupled state evolution from the RGS–PSU action,
then test convergence of finite Cantor cylinders and stability of Koopman
spectra as cylinder depth increases.

## Reproduction

~~~bash
PYTHONPATH=. python -m pytest tests/test_koopman_rgs_psu.py -q
~~~
