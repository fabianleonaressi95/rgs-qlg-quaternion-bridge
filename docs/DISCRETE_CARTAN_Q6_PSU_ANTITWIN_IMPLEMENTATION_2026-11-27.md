# Discrete Cartan implementation bridge — RGS / Q6 / PSU–Anti-Twin
## Companion implementation note — 27 November 2026

**Status:** implementation specification and falsifiable test plan. This note does not claim that a 4D Lorentzian complex or a physical continuum limit has already been constructed.

The primary mathematical specification is maintained in the theorem repository:
[Discrete Cartan construction proposal](https://github.com/fabianleonaressi95/theorem/blob/rgs-q6-antitwin-spacetime-bridge/research/DISCRETE_CARTAN_Q6_PSU_ANTITWIN_CONSTRUCTION_PROPOSAL_2026-11-27.md).

## 1. Data model to implement

The first executable prototype should use a small oriented simplicial 4-complex, initially one 4-simplex and its boundary, then a triangulated 4-cube. It must not pretend that the Q6 graph itself is a 4D spacetime mesh.

For each vertex \(v\), store a local Lorentz frame \(V_v\). For each oriented edge \(e=(v,w)\), store:

- \(L_{vw}:V_w\to V_v\), a proper orthochronous Lorentz link;
- \(E_{vw}\in V_v\), an edge coframe vector;
- inverse-edge constraints \(L_{wv}=L_{vw}^{-1}\), \(E_{wv}=-L_{vw}^{-1}E_{vw}\).

For a local Lorentz transformation \(h_v\), implement
\[
L_{vw}\mapsto h_vL_{vw}h_w^{-1},\qquad E_{vw}\mapsto h_vE_{vw}.
\]
Do not collapse \(L\) to a metric-derived Levi-Civita connection.

## 2. Face observables

Choose a base vertex and an ordered boundary path for every oriented triangular face \(f=(v_0,v_1,v_2)\). Define
\[
H_f=L_{v_0v_1}L_{v_1v_2}L_{v_2v_0},
\]
and
\[
\mathcal T_f=E_{v_0v_1}
+L_{v_0v_1}E_{v_1v_2}
+L_{v_0v_1}L_{v_1v_2}E_{v_2v_0}.
\]
The implementation must preserve the distinction:
- \(H_f\mapsto h_{v_0}H_fh_{v_0}^{-1}\): curvature holonomy;
- \(\mathcal T_f\mapsto h_{v_0}\mathcal T_f\): torsion/closure vector.

The logarithm \(\log H_f\) is allowed only on a declared local branch. Do not equate the raw holonomy matrix with a Lie-algebra curvature tensor at arbitrary lattice spacing.

## 3. Mandatory test matrix

| ID | Test | Expected result |
|---|---|---|
| DCart-01 | Boundary maps on the oriented complex | \(\partial^2=0\) exactly |
| DCart-02 | Random local Lorentz transformations | Edge and inverse-edge laws preserved |
| DCart-03 | Face gauge covariance | \(H_f'\!=h_0H_fh_0^{-1}\), \(\mathcal T_f'\!=h_0\mathcal T_f\) |
| DCart-04 | Reverse face orientation | \(H_{\bar f}=H_f^{-1}\), with basepoint transport handled |
| DCart-05 | Exact coframe \(E_{vw}=x_w-x_v\), \(L=I\) | \(\mathcal T_f=0\) |
| DCart-06 | Independent curvature/torsion configurations | Curvature and closure defect reported separately |
| DCart-07 | 3-cell Bianchi test | Ordered, basepoint-transported identity derived before numerical assertion |
| DCart-08 | Refinement sequence | Convergence to continuum curvature/torsion fluxes with stated norms |
| DCart-09 | PSU–Anti-Twin factor | \(J^2=I,\ JQJ^{-1}=-Q\); no identification with \(\mathcal T_f\) |
| DCart-10 | Spin lift | Lift consistency and \(\mathbb Z_2\) cycle signs verified |

DCart-07 and DCart-08 are mathematical gates. A numerical residual alone cannot promote either to a theorem.

## 4. Three geometry branches

Keep separate code paths or explicit constraint flags for:

1. Riemann–Cartan: independent coframe and Lorentz connection; curvature and torsion both allowed.
2. Teleparallel: curvature holonomy constrained to the flat branch; torsion may remain nonzero.
3. Torsion-free: closure/torsion constraint imposed; curvature may remain nonzero.

The teleparallel branch is not the general Cartan branch. The torsion-free branch is not a default assumption.

## 5. Anti-Twin and internal gauge factors

Keep Lorentz links \(L_{vw}\), internal gauge links \(U_{vw}\), and sheet exchange \(J_{\rm AT}\) in separate typed fields. For the two-sheet factor:
\[
J_{\rm AT}=I\otimes\tau_1,\quad Q=q_0\tau_3,\quad
J_{\rm AT}^2=I,\quad J_{\rm AT}QJ_{\rm AT}^{-1}=-Q.
\]
For conjugate U(1) links use \(U_e^{(-)}=(U_e^{(+)})^*\). These exact finite algebraic identities do not define Cartan torsion.

## 6. Implementation order

1. Build the oriented 4-simplex incidence data and verify \(\partial^2=0\).
2. Add Lorentz links and coframe edge data with inverse-edge constraints.
3. Implement local Lorentz gauge transformations and covariance tests.
4. Implement face holonomy and transported coframe closure.
5. Derive and implement the 3-cell Bianchi identities with explicit base paths.
6. Only then add spin lifts, matter operators, PSU–Anti-Twin couplings and a refinement study.

## 7. Evidence and literature guardrails

- Regge calculus supplies a simplicial framework for discrete geometry; it does not by itself establish the present torsion observable.
- DEC/FEEC supply cochain, Hodge and stability machinery; they do not choose a Lorentzian Cartan action.
- Braune et al. (2025), “Discrete Torsion of Connection Forms on Simplicial Meshes,” DOI 10.1145/3731197, is a direct comparison for discrete torsion constructions.
- S.-S. Xue (2009), “Quantum Regge Calculus of Einstein–Cartan theory,” *Physics Letters B* 682(3), 300–304, DOI 10.1016/j.physletb.2009.10.082, is Euclidean and uses a torsion-free Cartan structure equation in the stated construction. It must not be cited as proof of this independent torsionful proposal.
- Xue (2010), *Physical Review D* 82, 064039, DOI 10.1103/PhysRevD.82.064039, is a further comparison for 4D simplicial tetrad/connection variables.

**Promotion rule:** code passing DCart-01–06 establishes only implementation consistency for the selected finite complex. A theorem requires proofs of gauge covariance and Bianchi identities; a continuum claim additionally requires a refinement family and uniform consistency/stability estimates.
