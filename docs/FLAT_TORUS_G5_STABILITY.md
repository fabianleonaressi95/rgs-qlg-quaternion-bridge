# Flat-torus G5: uniform stability, coercivity and compactness

**Status:** analytic G5 argument and finite-level diagnostics for the flat periodic P1 family.  
**Scope limit:** no curved RGS-Q6 refinement estimate is asserted.

## 1. Fixed family and exact norms

Let K_h be the periodic Freudenthal triangulation of T_L^d, d=2 or 3,
with h=L/N. Use the G4 embedding J_h into the common Hilbert space
L^2(T_L^d). For nodal data u_h, write U_h=J_h u_h. The consistent
mass and stiffness matrices satisfy the exact identities
[
|U_h|_{L^2}^2=u_h^*M_hu_h,qquad
|
abla U_h|_{L^2}^2=u_h^*K_hu_h.
]
Thus the discrete H1 norm is defined by the right-hand side sum; there
is no hidden mesh-dependent normalization.

## 2. Uniform shape regularity

Each simplex is a translate of h times one of the d! reference
Freudenthal simplices. There are finitely many reference matrices
A_sigma, independent of N. Consequently
[
sup_{h,sigma}kappa(A_{h,sigma})
=max_{sigmainmathcal R_d}kappa(A_sigma)<infty,
qquad
sup_{h,sigma}rac{operatorname{diam}(sigma)}h<infty.
]
The code reports the maximum element condition number and
max_edge_length/h on three or more levels. Their level-independence is
a numerical cross-check of this exact finite-reference-element argument.

The metric is the Euclidean metric used to build the element coordinates.
Element volumes, P1 gradients, M_h and K_h all derive from that same
metric. Hence metric consistency is exact for this flat benchmark.

## 3. Uniform coercivity

On the common-space finite-element subspace, define
[
|U_h|_{H^1}^2
=|U_h|_{L^2}^2+|
abla U_h|_{L^2}^2.
]
By the exact assembly identities,
[
u_h^*M_hu_h+u_h^*K_hu_h
=|J_hu_h|_{H^1}^2.
]
The coercivity estimate therefore holds with constant c=1 in this
definition of the discrete H1 norm. On the mean-zero subspace, the
continuum torus Poincare inequality gives the stronger estimate
[
u_h^*K_hu_hge
left(rac{2pi}{L}ight)^2u_h^*M_hu_h,
quad
mathbf 1^*M_hu_h=0.
]
This follows because J_hu_h is an ordinary mean-zero H1 function on
the flat torus. The smallest positive generalized eigenvalue of
K_h u=lambda M_h u must therefore be at least (2 pi/L)^2.
The numerical audit checks this lower bound; it is not the proof.

## 4. Compactness

Suppose a sequence of mesh functions U_h=J_hu_h has uniformly bounded
L2 norm and P1 energy. The identities above imply that U_h is bounded
in H1(T_L^d). Since T_L^d is compact, the Rellich-Kondrachov theorem
gives a subsequence converging strongly in L2(T_L^d).

This proves the compactness gate for the embedded flat P1 family.
It does not require a separate graph-oscillation argument: the
piecewise-linear reconstruction is an actual H1 function and its
gradient energy is exactly the assembled stiffness form.

## 5. Numerical audit

Run:

    PYTHONPATH=. python -m pytest tests/test_flat_torus_stability.py -q

The report checks at least three levels:
- shape-condition and scaled-diameter stability;
- positivity of the exact consistent mass;
- the first positive generalized stiffness eigenvalue against the
  continuum Poincare lower bound;
- the constant-mode stiffness residual.

The generalized eigenvalue computation is a finite-resolution
diagnostic. The analytic Poincare and Rellich arguments above supply
the corresponding mathematical statements for this flat family.

## 6. Boundary of the result

G5 is now established for the regular flat periodic P1 benchmark under
the stated construction. It is not transferred to RGS-Q6. The curved
problem still needs a refinement family whose coframe-derived element
maps have uniform bi-Lipschitz bounds, uniformly positive metric
eigenvalues, controlled dual volumes and compatible embeddings. Those
bounds cannot be inferred from a single 64-vertex graph or from the
flat reference family.

G6 remains open for curved RGS: liminf and recovery must be proved for
the actual coframe-derived forms. The present flat G5 result only
supplies the stability template and explicit constants to seek.
