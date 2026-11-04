# Exact quadratic Wilson marginalization onto the Q6 five-mode sector

## Result

Let (D_0:mathbb R^{64}	omathbb R^{192}) and
(D_1:mathbb R^{192}	omathbb R^{240}) be the Q6 incidence maps, and
let (L_1=D_1^T D_1). On the Coulomb tangent
[
\mathcal C=\ker D_0^T,\qquad \dim\mathcal C=129,
]
let (B\in\mathbb R^{192\times5}) be an orthonormal basis of the
coexact eigenspace (L_1B=12B). Choose an orthonormal basis
(U\in\mathbb R^{192\times124}) for the orthogonal complement of
\(\operatorname{ran}B\) inside \(\mathcal C\).

The quadratic spatial Wilson action is
[
S_W^{(2)}(x)=\frac{\beta_s}{2}x^TL_1x,
\qquad x=Bq+U\eta.
]
Its Hessian in the split coordinates is
[
H=\beta_s
\begin{pmatrix}
B^TL_1B&B^TL_1U\\
U^TL_1B&U^TL_1U
\end{pmatrix}.
]
Since (L_1B=12B) and (B^TU=0),
[
B^TL_1B=12I_5,\qquad B^TL_1U=12B^TU=0.
]
Thus the Schur complement is exactly
[
H_{\rm eff}
=H_{qq}-H_{q\eta}H_{\eta\eta}^{-1}H_{\eta q}
=12\beta_s I_5.
]

The Gaussian fibre integral can therefore be done exactly:
[
\begin{aligned}
e^{-S_{\rm eff}^{(2)}(q)}
&=\int_{\mathbb R^{124}}
\exp\left[-\frac12\left(
12\beta_s|q|^2+
\beta_s\eta^T(U^TL_1U)\eta\right)\right]d^{124}\eta\\
&=\frac{(2\pi)^{62}}{
\sqrt{\det(\beta_s U^TL_1U)}}
\exp[-6\beta_s|q|^2].
\end{aligned}
]
Consequently,
[
\boxed{S_{\rm eff}^{(2)}(q)
=6\beta_s|q|^2+\frac12\log\det(\beta_s U^TL_1U)
-62\log(2\pi).}
]
The determinant term is independent of (q). The five-mode Maxwell
coefficient is unchanged by this Gaussian marginalization.

## Interpretation for Wyler–Bergman

This calculation establishes a useful negative result. At quadratic order,
integrating out the 124 real discarded Coulomb modes does **not** generate
a field-dependent Bergman factor, a Hua cutoff, or a factor
(V(D_5)^{-1/4}). It only produces a constant determinant. This is
consistent with the general fact that multiplying a measure by a
field-independent constant shifts the effective action by a constant.

The result does not rule out a field-dependent correction from the full
nonlinear compact Wilson action. For that problem one needs a specified
compact fibre prescription and the exact integer-flux sum. The existing
integer-lattice theorem proves that the full finite-\(\beta_t\) transfer
weight does not factor through the five retained coordinates: the discarded
integer fibre has rank 124 and contains nonzero fluxes invisible to the
five-mode projection but with transfer weight below the vacuum.

In particular, the Gaussian determinant computed here must not be
relabelled as the Bergman Jacobian. They are different objects.

## What remains open

1. Define the conditional Wilson fibre sum over each five-mode image,
   including its integer lattice and convergence.
2. Evaluate the resulting field-dependent pushforward weight, rather than
   assuming it is (N(z)^{-5}).
3. Establish a continuum limit and derive any cutoff from the transfer
   operator.
4. Only then compare the induced Maxwell normalization with Wyler's
   Euclidean volume expression.

No value of \(\alpha\), physical radius, or Wyler normalization is
predicted by this quadratic theorem.

## Reproducibility

Run:
```bash
pytest -q tests/test_rgs_q6_wilson_fiber_hessian.py
```

The executable audit is
`q6_wilson_fiber_hessian_report` in
`src/rgs_q6_wilson_fiber_hessian.py`.
