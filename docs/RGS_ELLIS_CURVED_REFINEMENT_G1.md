# Curved RGS refinement: Ellis–RGS geometric family (G1)

**Status:** curved refinement family constructed from the recovered RGS 64B seed and numerically audited on three levels.  
**Not claimed:** full RGS-EMS field-equation solution, Mosco convergence, or a continuum Dirac/Yang–Mills limit.

## 1. Recovered seed: no synthetic coordinates

The 64 coordinates are reconstructed exactly from the recovered
`X_Oh_reconstructed` array in
`A14_RGS_PSU_QLG_JACOBIAN_BUNDLE.zip` (SHA-256
`80442238349cd6ddc720dabae93d083075dea88decba8acf0a78c3ec03e0abcc`).
They consist of eight `O_h` orbits, with the six
absolute coordinate values recorded in
`data/rgs_188d/ellis_curved_seed_manifest.json`.
The orbit construction reproduces the archived 64x3 array entry by entry.

The 64 directions contain 56 distinct rays: eight pairs of source
vertices lie on repeated rays. The angular seed uses the 56 distinct
unit directions and retains the 64-to-56 inverse map. The convex hull
of these points on the unit sphere has 108 triangular faces.

The 192 Q6 edges and four-component edge coframe are reconstructed
from the A13 archive convention:
\[
e_{gh}^{\rm raw}=(X_h-X_g,\ \|X_h\|-\|X_g\|),\qquad
e_{gh}=e_{gh}^{\rm raw}/\|e_{gh}^{\rm raw}\|_{\mathbb R^4}.
\]
The archive is
`A13_RGS_PSU_QLG_KAGGLE_BUNDLE.zip` (SHA-256
`7ca0cdfeadbe1d693d4f776a5e22d957630b206b842f4b822f3afb578c488922`).
The Q6 edge coframe is preserved as its own graph-level data and audited;
it is not silently reinterpreted as the edge coframe of the new
tetrahedral mesh.

## 2. The curved target geometry

Use the RGS spherical branch
\[
R(u)=\sqrt{u^2+r_0^2},\qquad r_0=0.8292,
\]
with spatial metric
\[
g=du^2+R(u)^2d\Omega^2.
\]
The value \(r_0=0.8292\) is the throat parameter in the archived
Formulazione-4D diagnostic profile. That profile was not a solution
of the complete coupled field equations; here it fixes a geometric
benchmark, not a physical solution.

An exact isometric embedding into Euclidean \(\mathbb R^4\) is
\[
F(u,n)=\big(R(u)n,\ z(u)\big),\qquad
z(u)=r_0\,\operatorname{arsinh}(u/r_0),\quad n\in S^2.
\]
Indeed,
\[
R'(u)=u/R(u),\qquad z'(u)=r_0/R(u),\qquad
R'(u)^2+z'(u)^2=1,
\]
and hence
\[
F^*\delta_{\mathbb R^4}=du^2+R(u)^2d\Omega^2.
\]
This fixes the metric without free cube-sphere weights or fitted
anisotropy coefficients.

## 3. Refinement family

For a finite exhaustion interval \([-L,L]\), default \(L=3r_0\):

1. Start from the 56 distinct recovered RGS directions.
2. Triangulate their spherical convex hull (108 faces).
3. Refine every spherical edge by its normalized midpoint; each
   spherical triangle becomes four triangles.
4. Subdivide the radial interval into \(8\cdot2^\ell\) equal segments
   at level \(\ell\).
5. Form the product \(S^2\)-triangle x radial-interval prisms and split
   each prism into three tetrahedra by a globally consistent staircase
   rule.
6. Embed every vertex with the exact map \(F\) in \(\mathbb R^4\).

The first three levels have angular face counts 108, 432 and 1728.
The radial segment counts are 8, 16 and 32. The tetrahedral meshes
therefore form a nested geometric refinement family of the same
finite-domain manifold \([-L,L]\times S^2\).

Each tetrahedron is a flat 3-simplex in \(\mathbb R^4\). If its three
edge columns from one vertex are \(J_\sigma\), its induced Gram matrix
is
\[
G_\sigma=J_\sigma^T J_\sigma.
\]
The simplex volume is \(\sqrt{\det G_\sigma}/3!\). The P1 local
stiffness and consistent mass are assembled from this same Gram
matrix and volume; the mass is \(V_\sigma(I+\mathbf1\mathbf1^T)/20\).
Thus the element metric, volume and both operators share one geometry.

## 4. Structural and geometric audits

Run:

    python -m pip install -r requirements.txt
    PYTHONPATH=. python examples/rgs_ellis_curved_refinement_audit.py
    PYTHONPATH=. python -m pytest tests/test_rgs_ellis_refinement.py -v

The audit checks:
- exact reconstruction of the 64B seed and the 192-link Q6 coframe;
- the 56-direction spherical quotient and 108-face seed triangulation;
- nondegenerate 3D Gram matrices and element shape-condition bounds;
- exact hypersurface identity \(R'^2+z'^2=1\);
- exact simplicial chain identities \(B_1B_2=0\) and \(B_2B_3=0\);
- symmetry and constant-mode annihilation of the assembled P1 stiffness;
- consistency of the P1 mass volume with the sum of element volumes;
- convergence of the polyhedral volume toward
\[
\operatorname{Vol}_g([-L,L]\times S^2)
=4\pi\left(\frac{2L^3}{3}+2r_0^2L\right).
\]

The volume and shape checks are numerical evidence for the constructed
family. They are not, by themselves, a proof of uniform estimates for
all refinement levels.

## 5. Boundary and continuum status

The computational meshes approximate the compact manifold with
boundary \([-L,L]\times S^2\). The global Ellis spatial slice is
\(\mathbb R\times S^2\); it is noncompact. Passing from finite \(L\)
to the complete two-ended wormhole requires a separate exhaustion and
boundary-condition analysis. No periodic identification of the two
ends is made.

This phase establishes the **curved geometric family (G1)** and a
metric-consistent P1 discretization prototype. The curved common-space
embedding and the uniform G4/G5 estimates still need to be proved in
the chosen function spaces. Curved Mosco liminf/recovery (G6) remains
**OPEN**. Spin connection, SU(2) links, Dirac/overlap and Yang–Mills
limits remain separate later gates.

## 6. Evidence labels

- **RECOVERED:** original 64B coordinates and A13 Q6 edge/coframe
  conventions, with source archive hashes.
- **ANALYTICALLY EXACT:** the R4 embedding induces the stated Ellis
  warped-product metric.
- **NUMERICALLY AUDITED:** finite-level mesh nondegeneracy, chain
  identities, P1 matrix identities and volume refinement.
- **OPEN:** uniform all-level curved estimates, common-Hilbert-space
  Mosco convergence, and physical field-equation closure.
