# Curved Ellis–RGS G4: common-Hilbert-space P1 embedding

**Status:** explicit common L2 embedding and true curved-metric P1 mass/stiffness assembly implemented for the finite Ellis–RGS domain.  
**Scope:** M_L=[-L,L] x S2, with boundary at u=±L. This is not yet the complete noncompact two-ended wormhole or the full coframe-derived RGS-Q6 continuum.

## 1. One Hilbert space for every refinement

Fix
\[
\mathcal M_L=[-L,L]\times S^2,\qquad
g=du^2+R(u)^2d\Omega^2,\qquad
R(u)=\sqrt{u^2+r_0^2}.
\]
The common Hilbert space is
\[
\mathcal H=L^2(\mathcal M_L,d\mu_g),\qquad
d\mu_g=R(u)^2\,du\,d\Omega.
\]
All refinement levels approximate this same \(\mathcal M_L\), not different rescaled spaces.

Each angular seed face (n_a,n_b,n_c) is the radial projection of its planar chord triangle. In barycentric angular coordinates,
\[
w(\alpha,\beta)=(1-\alpha-\beta)n_a+\alpha n_b+\beta n_c,
\qquad
n(\alpha,\beta)=\frac{w(\alpha,\beta)}{\|w(\alpha,\beta)\|}.
\]
This gives a smooth chart on each spherical face. The charts agree on shared edges because the restriction to an edge depends only on its two endpoint directions.

The radial coordinate is affine on each interval. Every triangular prism is split into three reference tetrahedra with local coordinates (t,alpha,beta). The chart to the smooth manifold is
\[
\Phi_{h,\sigma}(t,\alpha,\beta)
=\big(u_j+\Delta u\,t,\ n(\alpha,\beta)\big).
\]
The exact R4 Ellis embedding F from G1 then maps this chart isometrically into the hypersurface.

## 2. Definition of J_h

Let V_h be the continuous nodal P1 space on the reference tetrahedral product complex. For nodal values u_h=(u_i), define
\[
(J_hu_h)\circ\Phi_{h,\sigma}(\xi)
=\sum_{a=0}^3u_{i_a}\lambda_a(\xi),
\]
where lambda_a are the four barycentric coordinates on the reference tetrahedron. The local maps agree on common faces, so J_hu_h is globally continuous and belongs to H1(M_L).

Thus
\[
J_h:\mathbb R^{V_h}\longrightarrow L^2(\mathcal M_L,d\mu_g)
\]
is an explicit common-space embedding. It is injective: a continuous P1 function that vanishes in L2 vanishes everywhere, hence at every mesh vertex.

The implementation exposes J_h(mesh, values) and the element evaluator evaluate_on_element(...) for pointwise evaluation in reference barycentric coordinates.

## 3. Exact curved metric in each reference element

Let y=(t,alpha,beta), and let A_sigma be the affine matrix from reference coordinates xi to y. The exact element map is
\[
\Psi_{h,\sigma}=F\circ\Phi_{h,\sigma}.
\]
Its Jacobian is evaluated analytically:
\[
D_\xi\Psi=D_y\Psi\,A_\sigma,\qquad
G_{h,\sigma}(\xi)=(D_\xi\Psi)^T D_\xi\Psi.
\]
The local common-space mass and Dirichlet matrices are
\[
(M_\sigma)_{ab}
=\int_{\widehat\sigma}\lambda_a\lambda_b
\sqrt{\det G_{h,\sigma}(\xi)}\,d\xi,
\]
\[
(K_\sigma)_{ab}
=\int_{\widehat\sigma}
(\nabla_\xi\lambda_a)^T G_{h,\sigma}(\xi)^{-1}
(\nabla_\xi\lambda_b)
\sqrt{\det G_{h,\sigma}(\xi)}\,d\xi.
\]
The global matrices are assembled from these curved integrals. They are not the chordal R4 simplex matrices from the G1 diagnostic; the latter remain available as a separate polyhedral approximation.

The code evaluates these integrals with tensor Gauss–Legendre quadrature after a Duffy map from the unit cube to the reference tetrahedron. The quadrature order is explicit and the total mass 1^T M_h 1 is checked against the exact manifold volume
\[
\operatorname{Vol}_g(\mathcal M_L)
=4\pi\left(\frac{2L^3}{3}+2r_0^2L\right).
\]
The mathematical forms above are defined by the exact integrals; quadrature is their reproducible numerical assembly.

## 4. G4 audit

Run:

    PYTHONPATH=. python -m pytest tests/test_rgs_ellis_refinement.py -v
    PYTHONPATH=. python examples/rgs_ellis_curved_refinement_audit.py

The G4 tests verify:
- the reference-tetrahedron quadrature weights sum to 1/6;
- nodal interpolation reproduces all four nodal values;
- the common-space mass matrix is symmetric positive definite;
- the curved stiffness is symmetric and annihilates constants;
- the integrated mass volume agrees with the analytic Ellis volume to the stated quadrature tolerance.

The G4 map and exact curved metric are established for this finite Ellis–RGS benchmark. The next proof gate is uniform G5 stability: uniform shape regularity of the chart maps and compactness in the common L2 space. Curved Mosco convergence is not inferred merely from a successful finite-level quadrature audit.
