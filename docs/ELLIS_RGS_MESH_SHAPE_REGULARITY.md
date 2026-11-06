# Conditional shape-regularity criterion for the Ellis-RGS refinement

**Status:** a sufficient-condition argument and numerical evidence, not yet a complete uniform-regularity proof for the generated family.

## 1. Smooth geometry

On the compact cylinder
\[
\mathcal M_L=[-L,L]\times S^2,\qquad
g=du^2+R(u)^2d\Omega^2,\quad R(u)=\sqrt{u^2+r_0^2},
\]
the embedding
\[
F(u,n)=(R(u)n,r_0\operatorname{arsinh}(u/r_0))
\]
is smooth. Its differential has singular values \(1,R(u),R(u)\) relative to the product reference metric \(du^2+d\Omega^2\). Therefore
\[
r_0\leq \sigma_{\min}(dF)\leq \sigma_{\max}(dF)
\leq \sqrt{L^2+r_0^2},
\]
and the metric distortion is uniformly bounded on the fixed compact domain. The same statement holds for the inverse chart differential.

## 2. Conditional mesh criterion

Let \(h_\ell\) denote the maximum spherical edge length of the level-\(\ell\) angular triangulation and let
\[
\Delta u_\ell=2L/(m_0 2^\ell)
\]
be the radial step. A sufficient route to uniform tetrahedral shape regularity is:

1. the angular triangulations are uniformly shape regular in spherical geodesic charts: every face has diameter comparable to \(h_\ell\), and its area is bounded below by \(c_\Omega h_\ell^2\);
2. radial/angular scales are comparable: \(0<c_-\leq \Delta u_\ell/h_\ell\leq c_+<\infty\);
3. the fixed prism staircase split has a nondegenerate reference tetrahedron family;
4. the smooth Ellis map and the angular normalization charts have uniformly bounded first derivatives and inverses on the compact domain.

Under these hypotheses, the parameter tetrahedra are uniformly shape regular. Composition with \(F\) preserves this property up to constants depending only on \(L,r_0,c_-,c_+\), the angular shape bound and the fixed staircase split. For sufficiently fine levels, the chordal tetrahedron Gram matrices inherit a uniform condition-number bound: the vertex-difference map is the derivative of the smooth parameterization plus an \(O(h_\ell^2)\) Taylor remainder, while its leading singular values are all comparable to \(h_\ell\).

This argument reduces the geometric part of the proof to verifying the angular subdivision estimate, the radial/angular comparability and the uniform Taylor remainder. It does **not** establish those estimates merely by observing three numerical levels.

## 3. What is now measured

The multi-level audit reports:

- minimum/maximum spherical angular edge lengths;
- minimum/maximum spherical face areas;
- the scale-free proxy \(\min_K |K|/\operatorname{diam}(K)^2\);
- minimum Gram eigenvalue and maximum tetrahedron Gram condition number;
- maximum chord diameter and relative P1 volume error.

Run:

\`\`\`bash
PYTHONPATH=. python examples/run_ellis_mesh_shape_regularity_audit.py --levels 0 1 2 --radial-base 8
\`\`\`

The current level-0/1/2 experiment has maximum tetrahedral Gram condition numbers 60.11, 57.22 and 59.61, respectively. The exact Ellis volume error decreases from 8.849% to 2.334% to 0.592%, with observed error-halving orders 1.923 and 1.980. These values are numerical evidence only.

## 4. Next proof obligations

- Prove the normalized-midpoint spherical subdivision has a uniform shape bound, starting from the recovered 56-direction/108-face seed.
- Bound \(\Delta u_\ell/h_\ell\) above and below using the initial angular edge extrema and the refinement rule.
- Convert the smooth-map Taylor estimate into a uniform bound for all four edge vectors of every tetrahedron.
- Prove equivalence of the discrete lumped P1 measure and the continuum Ellis \(L^2\) measure.
- Then prove a uniform discrete Poincare/Rellich estimate and construct recovery sequences for the curved energy.

Only after these steps is it appropriate to claim a uniform geometric hypothesis for the spinorial Mosco argument.


## Common Hilbert-space mass comparison

The helper \`audit_ellis_lumped_mass_equivalence\` assembles the curved-metric P1 mass matrix and compares it with the diagonal vertex measure used by the spin-Wilson operator. It computes the extreme eigenvalues of
\[
D_{\rm lump}^{-1/2}M_{\rm curved}D_{\rm lump}^{-1/2}
\]
using sparse extremal-eigenvalue routines, avoiding dense global matrices.

Run:
\`\`\`bash
PYTHONPATH=. python examples/run_ellis_lumped_mass_equivalence_audit.py --levels 0 1 --radial-base 8 --quadrature-order 4
\`\`\`

The finite-level eigenvalue bounds are diagnostics. A uniform bound over every refinement level still requires analytic control of the element Jacobian ratios and quadrature error.


### Conditional local proof of mass equivalence

On the reference tetrahedron \(\widehat K\), the affine P1 consistent mass matrix has eigenvalues in the ratio \([1/5,1]\) relative to the vertex-lumped mass matrix. Suppose the exact curved element Jacobian \(J_K(\xi)\) and the chordal element volume \(V_K^{\rm ch}\) satisfy the level-independent bounds
\[
0<c_J\leq \frac{J_K(\xi)}{6V_K^{\rm ch}}\leq C_J
\quad\text{for all }\xi\in\widehat K.
\]
Then, for every local P1 function,
\[
\frac{c_J}{5}\,m_{K,\rm lump}(v,v)
\leq m_{K,\rm curved}(v,v)
\leq C_J\,m_{K,\rm lump}(v,v).
\]
Summing over elements gives global spectral equivalence with the same constants. Thus the analytic issue is reduced to a uniform element-Jacobian/chordal-volume comparison; the finite eigenvalue computation is a diagnostic of that comparison, not its proof.

For levels 0 and 1 with radial_base=8, the computed generalized mass eigenvalue intervals are respectively [0.22510, 1.10338] and [0.22115, 1.02534]. The exact curved mass volume errors at quadrature order 4 are \(1.54\times10^{-6}\) and \(2.53\times10^{-8}\). The chordal lumped volume errors remain larger (8.849% and 2.334%), as expected for the piecewise-linear geometric measure. Uniformity over all levels is still OPEN.
