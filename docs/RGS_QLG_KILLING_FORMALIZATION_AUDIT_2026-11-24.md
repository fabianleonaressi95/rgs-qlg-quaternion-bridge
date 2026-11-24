# RGS–QLG Killing and geometric formalization audit

**Date:** 24 November 2026  
**Scope:** independent mathematical review of the Killing/vierbein/Quaternion–QLG material proposed in PR #221  
**Classification:** exact algebraic identities separated from numerical evidence, incomplete discretizations, and open geometric claims.

## 1. Executive conclusion

The new repository material contains a useful reconstruction pipeline and a useful research direction, but its present Killing-vector conclusion is not established.

The following statements must not be promoted to theorems:

- “exactly three Killing vectors”;
- “the three modes are pure translations” as a classification of the isometry algebra;
- “homogeneous but not isotropic”;
- “the discrete Q6, 2I, or O_h symmetries extend partially to continuous isometries.”

The current src/vierbein/killing.py does not implement the Riemannian Killing equation. Its residual differentiates the components of the trial vector field and symmetrizes those partial derivatives, while it does not use the metric to lower indices, does not differentiate the metric, and does not include the Levi-Civita connection. Constant vector fields therefore have identically zero residual by construction, independently of whether they preserve the reconstructed metric.

The proper result is a formal correction and a reproducible verification programme, not a declaration that the three modes have been disproved numerically. The supplied JSON is evidence of what the current code computes; it is not evidence of the continuum Killing equation.

## 2. RGS metric and vierbein: exact algebraic layer

Let the embedded vertex positions be X_i in R^3, with
\[
d_{ij}=X_j-X_i,\qquad
w_{ij}=\exp(-\|d_{ij}\|^2/\ell^2).
\]
The covariance tensor used in the repository is
\[
C_i=\sum_{j\ne i}w_{ij}\,d_{ij}d_{ij}^{T}.
\]
For c_0>0, define
\[
g_i=c_0^{-1}C_i,\qquad e_i=g_i^{1/2},
\]
where the square root is the principal symmetric positive square root.

If C_i is positive definite, then g_i is positive definite and
\[
e_i^T e_i=e_i^2=g_i=c_0^{-1}C_i.
\]
This is an exact matrix identity. A numerical residual near machine precision verifies its implementation for the supplied finite data. It does not by itself establish a smooth Riemannian manifold, a continuum metric, or a Lorentzian spacetime.

For the reported finite archive, the JSON records 64/64 positive matrices, minimum eigenvalue approximately 0.591757, and maximum residual in e_i^T e_i=C_i/c_0 approximately 3.39e-15. These are numerical claims tied to the archive hash recorded in the JSON; they should remain labelled finite-data verification until the archive and script are rerun together from a clean checkout.

**Status:** g_i=C_i/c0 and e_i^T e_i=g_i are exact algebraic identities conditional on the definitions and positive definiteness. The reported finite values are numerical verification, not a continuum theorem.

## 3. Correct Killing equation

Let (M,g) be a smooth Riemannian 3-manifold and xi a smooth vector field. The Killing equation is
\[
\nabla_\mu\xi_\nu+\nabla_\nu\xi_\mu=0,
\qquad \xi_\nu=g_{\nu\lambda}\xi^\lambda.
\]
Equivalently,
\[
(\mathcal L_\xi g)_{\mu\nu}
=\xi^\lambda\partial_\lambda g_{\mu\nu}
+g_{\lambda\nu}\partial_\mu\xi^\lambda
+g_{\mu\lambda}\partial_\nu\xi^\lambda=0.
\]
The equivalence follows from metric compatibility and vanishing torsion of the Levi-Civita connection:
\[
(\mathcal L_\xi g)_{\mu\nu}
=\nabla_\mu\xi_\nu+\nabla_\nu\xi_\mu.
\]

The current code instead evaluates, in effect,
\[
\partial_\mu\xi^\nu+\partial_\nu\xi^\mu,
\]
with an ad hoc coordinate-wise difference quotient. This expression is not the Killing operator for a general nonconstant metric. In particular, it omits the term xi^lambda partial_lambda g_munu and the metric contractions. The g argument is stored by KillingOperator, but is not used in build_symmetric_deformation or killing_residual.

Consequently every constant trial field xi^a=constant has zero residual in this implementation. This is a tautology of the surrogate operator, not a test that translations are isometries of g.

The search is also not an exhaustive kernel computation: find_killing_vectors tests six prescribed ansatz fields (three constants and three rigid-rotation candidates). build_operator_matrix explicitly raises NotImplementedError. Thus “exactly three” is not established even for the implemented discrete surrogate, much less for the continuum metric.

## 4. A decisive consistency theorem for the claimed three translations

### Theorem 1 — Three coordinate translations force local flatness and six local Killing fields

Let g=g_munu(x) dx^mu dx^nu be a smooth positive-definite metric on a connected coordinate domain U in R^3. Suppose all three coordinate translation fields
\[
xi_{(1)}=\partial_{x^1},\quad
xi_{(2)}=\partial_{x^2},\quad
xi_{(3)}=\partial_{x^3}
\]
are Killing fields.

Then
\[
0=(\mathcal L_{\partial_{x^a}}g)_{\mu\nu}
=\partial_{x^a}g_{\mu\nu}
\]
for all a,mu,nu. Hence every component g_munu is constant on U. A constant positive-definite metric is flat. Moreover its local Killing algebra has dimension six: three translations and three independent linear fields
\[
xi_A(x)=A x,\qquad A^Tg+gA=0.
\]
The space of such A is isomorphic to so(3), hence has dimension three.

**Proof.** The coordinate expression for the Lie derivative gives the first implication. A constant positive-definite matrix is reduced to the Euclidean metric by a constant linear coordinate change, so its Riemann tensor vanishes. For every matrix A satisfying A^Tg+gA=0, the field Ax obeys L_(Ax)g=0. The skew-adjoint endomorphisms for a positive-definite inner product in dimension three form a three-dimensional Lie algebra. Together with the three translations this gives six independent local Killing fields. QED.

**Consequence.** The simultaneous claims “the three coordinate translations are exact Killing fields”, “there are exactly three Killing fields”, and “the metric is homogeneous but intrinsically anisotropic” cannot all hold for a smooth positive-definite 3D metric on a connected coordinate domain. Either the proposed translations are not true Killing fields, or the search for rotational Killing fields is incomplete. Component anisotropy in a coordinate matrix is not intrinsic anisotropy: a constant positive-definite metric is locally Euclidean.

This theorem does not rule out a three-dimensional isometry algebra of some other type, nor does it say that three arbitrary Killing fields on a 3-manifold always imply flatness. The conclusion above uses the specific claim that the fields are the three coordinate translations.

## 5. What a defensible discrete Killing audit requires

The 64 tensors g_i form a vertex-indexed tensor field. They do not alone specify derivatives, a smooth interpolation, or a Levi-Civita connection. Before claiming continuum Killing vectors, the programme must declare an interpolation or refinement family and prove its consistency.

A finite-level diagnostic can be defined explicitly. For each vertex i, let N(i) be its chosen neighbors and set
\[
A_i=\sum_{j\in N(i)}\omega_{ij}d_{ij}d_{ij}^T.
\]
If rank(A_i)=3, define the weighted least-squares gradient of a scalar vertex field f by
\[
\nabla_h f(i)=A_i^{-1}
\sum_{j\in N(i)}\omega_{ij}d_{ij}(f_j-f_i).
\]
Apply this componentwise to g_munu and to xi^mu. The discrete Lie derivative is then
\[
(\mathcal L^h_\xi g)_{i,\mu\nu}
=xi_i^\lambda(\nabla^h_\lambda g_{\mu\nu})_i
+g_{i,\lambda\nu}(\nabla^h_\mu xi^\lambda)_i
+g_{i,\mu\lambda}(\nabla^h_\nu xi^\lambda)_i.
\]
A normalized residual can be reported as
\[
rho_h(xi)=
\frac{\|\mathcal L^h_\xi g\|_{\omega}}
{\|g\|_{\omega}\,\|\nabla^h xi\|_{\omega}
+\|xi\|_{\omega}\,\|\nabla^h g\|_{\omega}+\varepsilon}.
\]
The weights, neighborhood, pseudoinverse policy, boundary treatment, tensor-component convention, and normalization must be fixed in code and in the report.

This is a discrete diagnostic, not automatically a continuum Killing theorem. Required gates:

1. Verify rank and conditioning of every A_i; report rank-deficient vertices.
2. Verify exact reproduction of constant and affine scalar/vector fields by the gradient reconstruction where the geometry permits it.
3. Use synthetic metrics with known isometry algebras: constant SPD metric (six fields), a generic anisotropic nonconstant metric (no assumed translations), and a known homogeneous curved metric (three or more known generators as appropriate).
4. Compute the full linearized discrete operator matrix and its numerical kernel; do not test only six hand-selected ansatz fields.
5. Compare the finite-difference Lie derivative with an independently assembled covariant Killing operator using reconstructed Christoffel symbols.
6. Repeat on a declared refinement family. Only a consistent, stable limit can support a continuum claim.

For a smooth metric, the Christoffel symbols are
\[
Gamma^rho_munu
=\tfrac12g^{rho sigma}
(partial_mu g_sigma nu+partial_nu g_sigma mu
-partial_sigma g_mu nu),
\]
and the covariant Killing operator must be assembled from these derivatives, not omitted.

## 6. Lie algebra and conservation laws

For genuine smooth Killing fields, their Lie brackets are Killing fields:
\[
[xi,eta]^\mu=xi^\nu\partial_\nu eta^\mu-eta^\nu\partial_\nu xi^\mu.
\]
Closure must be checked by solving for coefficients in a basis of the computed kernel, with a residual for the component orthogonal to that kernel. Testing only whether a bracket norm is small is not a proof of closure.

For a geodesic with tangent u^mu, a Killing field gives the conserved quantity
\[
Q_xi=g_{\mu\nu}xi^\mu u^\nu.
\]
The conservation law follows from the geodesic equation and nabla_(mu xi_nu)=0. A test that repeatedly evaluates the same fixed vertex is not a geodesic integration and cannot verify this law.

The current verify_killing_conservation explicitly uses a fixed vertex as a placeholder. Its zero variation must not be reported as geodesic conservation evidence.

## 7. Separate geometric structures in the RGS–QLG bridge

The following objects must remain distinct.

### 7.1 Metric / Levi-Civita geometry

- Spatial metric candidate: g_i=e_i^T e_i=C_i/c_0.
- A continuum metric requires interpolation/refinement and regularity.
- The Levi-Civita connection Gamma(g) is uniquely determined by a smooth metric; it is metric-compatible and torsion-free.
- The Killing equation is defined relative to this metric and connection.

### 7.2 Vierbein / local frame

A coframe e^a=e^a_mu dx^mu gives
\[
g_{\mu\nu}=\delta_{ab}e^a_\mu e^b_\nu
\]
for the spatial Riemannian sector. A Lorentzian extension requires a separately specified time coframe, lapse/shift or equivalent 4D construction, and internal signature eta_ab. A list of 64 spatial 3x3 matrices is not yet a 4D tetrad field.

### 7.3 Teleparallel connection and torsion

For a chosen spin connection omega^a_b,
\[
T^a=de^a+\omega^a_b\wedge e^b,\qquad
R^a_b(omega)=d omega^a_b+omega^a_c\wedge omega^c_b.
\]
Teleparallel geometry requires a specified flat spin connection R(omega)=0; torsion is then computed from the coframe and that connection. The Levi-Civita connection is not the teleparallel spin connection. Nonorientability, a varying vierbein, or a nontrivial quaternionic link field alone does not prove nonzero teleparallel torsion.

The current src/teleparallel/__init__.py returns an all-zero placeholder tensor and treats omitted curvature (None) as satisfying the flatness test. It is not a torsion computation or a TEGR audit.

### 7.4 QLG / Yang–Mills connection

A gauge connection A on a specified principal G-bundle has curvature
\[
F_A=dA+A\wedge A
\]
and holonomy given by path-ordered parallel transport. It is not automatically the Levi-Civita or teleparallel spin connection. An identification requires an explicit bundle map, representation, transformation law, and action.

The current src/yang_mills/__init__.py returns zero for Wilson loops and for every Yang–Mills action. It is a placeholder, not a numerical QLG result.

### 7.5 Quaternionic bridge

The proposed edge map
\[
M_{ij}=e_j e_i^{-1}
\]
and its polar decomposition M_ij=R_ij S_ij can be studied as a candidate discrete frame transport. Because e_i are invertible positive matrices, det(M_ij)>0, so the orthogonal polar factor has determinant +1. The lift R_ij to q_ij in SU(2) still requires a consistent sign/orientation convention.

Before calling this a gauge connection, verify:
- reverse-edge relation R_ji=R_ij^{-1} and compatible quaternion lift;
- local-frame transformation law and covariance of plaquette holonomy;
- path ordering and edge orientation on every face;
- gauge invariance of the Wilson action;
- independent comparison of matrix and quaternion holonomies.

The polar rotation discards the positive stretch S_ij. Therefore the full matrix M_ij, the polar rotation R_ij, and the lifted quaternion q_ij are distinct data; the latter is not a lossless encoding of the former.

## 8. Audit of PR #221 implementation scope

The source review of PR #221 identifies the following blockers to treating the PR as a completed RGS–QLG implementation:

| Component | Current state in proposed files | Correct classification |
|---|---|---|
| C_i, g_i=C_i/c0, principal e_i | Explicit finite construction | Algebraic definition; finite numerical checks reported |
| Killing residual | Omits g, partial g, and Christoffel terms | Not the Riemannian Killing operator |
| Killing search | Six fixed ansatz candidates; no assembled operator matrix | Candidate diagnostic only; not exhaustive |
| Lie bracket | Ad hoc finite differences and projection | Not a certified structure-constant computation |
| Geodesic conservation | Fixed-vertex placeholder | Not tested |
| Teleparallel torsion | Returns zero placeholder | Not implemented |
| Teleparallel flatness | Missing curvature treated as success | Not verified |
| Yang–Mills Wilson loop/action | Return zero placeholders | Not implemented |
| Quaternion bridge | Polar/lift pipeline proposed | Candidate discrete transport; covariance tests required |
| Tests | No dedicated test files appear in the PR #221 changed-file list | No regression suite for these claims is evidenced by this PR |
| Repository hygiene | PR includes tracked __pycache__/*.pyc files | Remove generated bytecode from source PR |

PR #221 also carries a large recovered .npz data artifact. Its provenance hash is recorded in the audit JSON; reviewers must ensure it is the canonical archive and not silently overwrite or duplicate an existing source of truth.

## 9. Corrected theorem/evidence ledger

| Statement | Status |
|---|---|
| g_i=C_i/c_0 and e_i^T e_i=g_i for the defined principal square root | **PROVED** (finite-dimensional algebra, under C_i>0,c_0>0) |
| 64 covariance matrices are positive definite with the reported minimum eigenvalue | **NUMERICALLY REPORTED**; rerun from canonical NPZ required |
| Three constant coordinate fields have zero residual in the current surrogate operator | **EXACT FOR THE CODE**, but not a Killing result |
| Three coordinate translations are Killing fields of a continuum RGS metric | **OPEN / NOT ESTABLISHED** |
| Exactly three Killing fields exist | **OPEN / NOT ESTABLISHED** |
| The metric is homogeneous but not isotropic | **NOT ESTABLISHED**; incompatible with the asserted three coordinate translations if all are exact Killing fields |
| Their Lie algebra is R^3 | **NOT ESTABLISHED** as an isometry algebra; the constant ansatz fields commute algebraically |
| Quaternionic polar link is a valid QLG gauge connection | **CANDIDATE** pending covariance, orientation, and holonomy tests |
| RGS data define a teleparallel torsion field / TEGR dynamics | **OPEN** |
| QLG Yang–Mills action and Wilson observables are implemented in the added placeholder module | **FALSE**; current functions return zero placeholders |

## 10. Required next implementation phase

Do not enlarge the current CI or add speculative physical claims. The minimal next phase is:

1. Keep the finite metric reconstruction and its exact identity as a separate, stable baseline.
2. Replace the current Killing surrogate with a metric-aware discrete Lie derivative, explicitly recording derivative reconstruction rank/conditioning.
3. Assemble the full linear operator and compute its kernel; compare against analytic synthetic controls.
4. Do not label the result “exactly three” until the full kernel and refinement behavior are independently checked.
5. Implement teleparallel torsion only after coframe, spin connection, orientation, and discrete exterior derivative have been specified.
6. Implement Wilson loops and Yang–Mills action from actual quaternion/matrix holonomies; remove zero placeholders from any completion claim.
7. Keep metric, Levi-Civita, teleparallel, and internal gauge connections separate in documentation and APIs.

**Bottom line:** the RGS covariance/vierbein identity is a legitimate finite geometric construction. The current Killing, torsion, geodesic-conservation, and Yang–Mills modules do not yet establish the advertised continuum or QLG conclusions. This audit preserves the useful construction while correcting the epistemic status of the claims.
