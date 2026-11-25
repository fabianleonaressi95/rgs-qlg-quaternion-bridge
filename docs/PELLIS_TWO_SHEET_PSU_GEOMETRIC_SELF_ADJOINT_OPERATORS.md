# Geometric Self-Adjoint Operators and the Two-Sheeted PSU Geometry

**Technical note for the RGS–PSU–PRISM discussion with Stergios Pellis**  
**Version:** 0.1 — 25 November 2026  
**Status:** finite-model specification; continuum identification remains open

## Abstract

This note formalizes the two-sheeted PSU/Anti-Twin construction using only objects present in the recovered Q6/QLG archives. It distinguishes (i) the Q6 graph with an SU(2) fibre connection and its conjugate on the second sheet, (ii) the separate finite-dimensional Twin/Anti-Twin exchange algebra, and (iii) a proposed continuum geometry. The first two admit exact finite formulations and reproducible checks. The archives do not yet provide a canonical equivariant identification between the Q6 64-state space and the 8-dimensional internal Twin/Anti-Twin space, nor do they establish a continuum limit. Neither is assumed.

## Correspondence addressed by this note

In his reply of 16 November 2026, Pellis asked specifically for (i) the mathematical spaces corresponding to the two sheets, (ii) the involution or mapping relating them, and (iii) the metric and measure on the combined configuration. He then asked whether the two-sheeted structure has a natural relationship with the spectral or algebraic constructions in PRISM. His earlier reply of 5 November had separately emphasized the need for a precisely defined Hilbert space, operator domain, self-adjointness conditions, spectral measure, multiplicities and normalization. The operator-domain and continuum-limit sections below develop the mathematical consequences of those requests; they should not be misrepresented as verbatim bullet points from the 16 November email.

## 1. Base graph and two-sheeted lift

Let
\[
Q_6=(V,E),\qquad |V|=64,\quad |E|=192,
\]
be the recovered six-dimensional hypercube graph, with its canonical vertex and edge ordering. Its two-sheeted lift is
\[
\widetilde Q_6=Q_6\times\{+,-\}.
\]
The deck transformation is
\[
\tau(g,+)=(g,-),\qquad \tau(g,-)=(g,+),\qquad \tau^2=1.
\]
This is a two-fold covering of the finite graph. It is not, by itself, a wormhole, a non-orientable manifold, or a statement about spacetime topology.

For the scalar graph model, the Hilbert space is
\[
\mathcal H_0=\ell^2(V)\oplus\ell^2(V).
\]
For the archived rank-two fibre connection, use
\[
\mathcal H_{\rm fib}
=\ell^2(V;\mathbb C^2)\oplus
\ell^2(V;\overline{\mathbb C^2}).
\]
The bar denotes the conjugate vector space, not an unrelated second copy with an independently chosen complex structure. The direct-sum inner product is the sum of the standard fibre inner products. The rank-two model has 256 complex dimensions in total (128 per sheet); this is distinct from the scalar 128-dimensional vertex space and from the separate 16-dimensional throat model in Section 5.

## 2. Sheet involution

Let \(c:\mathbb C^2\to\overline{\mathbb C^2}\) be the canonical complex-linear identification with the conjugate vector space, and \(c^{-1}\) its inverse. Define
\[
J(\psi,\bar\chi)=
\big(c^{-1}(\bar\chi),c(\psi)\big).
\]
In coordinates this is often written \((\psi,\chi)\mapsto(\bar\chi,\bar\psi)\). The conjugate-space notation matters: \(J\) is complex-linear on \(\mathbb C^2\oplus\overline{\mathbb C^2}\), not an anti-linear map on one copy of \(\mathbb C^2\oplus\mathbb C^2\).

Directly,
\[
J^*=J,\qquad J^2=I.
\]
Thus
\[
P_{\rm even}=\tfrac12(I+J),\qquad
P_{\rm odd}=\tfrac12(I-J)
\]
are orthogonal complementary projections, and
\[
\mathcal H_{\rm fib}=\mathcal H_+\oplus\mathcal H_-,
\qquad \mathcal H_\pm=\ker(J\mp I).
\]
This even/odd decomposition is a theorem of the defined finite-dimensional involution.

## 3. Conjugate SU(2) connection and graph form

The recovered QLG archive contains edge transports \(U_e\in SU(2)\) on the first sheet and conjugate transports \(\bar U_e\) on the second. They are unitary:
\[
U_e^*U_e=I_2,\qquad \bar U_e^*\bar U_e=I_2.
\]
For an oriented edge \(e=(g,h)\), define the covariant edge difference
\[
(D_U\Psi)_e=
\left(\psi_h-U_e\psi_g,\;
\bar\chi_h-\bar U_e\bar\chi_g\right).
\]
The second component is interpreted in the conjugate fibre. Exchanging the sheets and conjugating the fibre identifies the two covariant-difference blocks; write this compatibility as \(D_UJ=J_ED_U\), where \(J_E\) is the corresponding exchange on edge fields.

Choose positive edge coefficients \(\kappa_e>0\). The quadratic form is
\[
q_U[\Psi]=\sum_{e\in E}\kappa_e
\left(\|\psi_h-U_e\psi_g\|^2+
\|\bar\chi_h-\bar U_e\bar\chi_g\|^2\right).
\]
On the finite graph it is non-negative and defined on the whole finite-dimensional Hilbert space. Its associated operator is
\[
L_U=D_U^*WD_U,\qquad W=\operatorname{diag}(\kappa_e I),
\]
so
\[
L_U=L_U^*\geq0,\qquad [L_U,J]=0.
\]
Therefore the even and odd sheet sectors are invariant under the covariant graph Laplacian.

**Metric and measure.** The current reproducible baseline uses counting measure on vertices. The recovered embedding \(X:V\to\mathbb R^3\) supplies edge lengths \(\ell_e=\|X_h-X_g\|\) when nonzero; their shortest-path metric is a metric on the finite connected graph. A canonical geometric vertex measure (for example, dual-cell volumes) has not yet been fixed in the cited archive. No such volume is silently inserted into the Hilbert inner product. If positive vertex weights \(\mu_g\) are later selected, the adjoint and operator must be recomputed using that measure.

The coefficient \(\kappa_e=\ell_e^{-1}\) is the one-dimensional P1 edge-stiffness convention used in the separate metric-dependent Q6 matter-sector benchmark. It is a declared graph Dirichlet-energy convention, not yet a derived four-dimensional continuum action or a unique PSU constitutive law.

## 4. What the QLG archive verifies

The A13 integrated audit records 64 embedded vertices, 192 Q6 links, a tangent frame \(Q\in\mathbb R^{192\times188}\), its orthonormality and four-dimensional normal complement, and an edge coframe of shape \(192\times4\) reconstructed from the embedding and oriented links. It checks the SO(3) link transports and SU(2) fibre transports, including orthogonality/determinant and unitarity. The second biverse sheet is represented by the complex-conjugate SU(2) connection.

The A13 audit explicitly does not identify the three Bell-correlation coordinates with the 64 Q6 states: no equivariant vertex-to-Q6 map is assumed. Its multifractal scan is a finite-graph diagnostic, not a proof of multifractality or continuum scaling.

The separate geometry-driven 64-to-188 feature-map script is a candidate construction-and-test pipeline. It warns that a projected-polar extension can commute with the archived kernel for a kernel-inclusion reason, and that a direct commutator is meaningful only if the tangent-coordinate bases coincide. It is not a proved canonical identification of the Q6 Anti-Twin action with the 188-dimensional operator.

## 5. Separate exact finite-dimensional Twin/Anti-Twin algebra

A different recovered archive contains an internal Hilbert space \(K\cong\mathbb C^8\), a Hermitian exchange involution \(E\), and Hermitian observables \(O_a\) satisfying
\[
E^*=E,\quad E^2=I_8,\qquad
O_a^*=O_a,\quad O_a^2=I_8,\quad EO_aE=-O_a.
\]
The exchange projections \(P_e=(I_8+E)/2\) and \(P_o=(I_8-E)/2\) are orthogonal and complementary. The archived checks verify these identities and the oddness of the three correlation observables.

On \(\mathbb C^2_{\rm sheet}\otimes K\), define the effective finite-dimensional model
\[
H(t,\delta)=-t\,\sigma_x\otimes I_8
             +\delta\,\sigma_z\otimes O_a,
\qquad t>0,\quad\delta\in\mathbb R.
\]
With \(\mathcal J=\sigma_x\otimes E\), one has
\[
\mathcal J^*=\mathcal J,\quad \mathcal J^2=I_{16},
\quad [H,\mathcal J]=0.
\]
Since \(\sigma_x\sigma_z+\sigma_z\sigma_x=0\) and \(O_a^2=I_8\),
\[
H^2=(t^2+\delta^2)I_{16}.
\]
Hence
\[
\operatorname{spec}(H)=
\{-\sqrt{t^2+\delta^2},+\sqrt{t^2+\delta^2}\},
\]
each eigenvalue having multiplicity eight. The lower branch is
\[
E_-=-t-\frac{\delta^2}{2t}
+O(\delta^4/t^3).
\]
The linear sheet-odd contribution cancels exactly in this declared model.

This is an exact finite-dimensional algebraic result for the archived 8-dimensional internal exchange and the effective hopping parameter \(t\). It is not a derivation of \(t\) from a geometric throat, an Einstein–Maxwell–scalar solution, or a continuum field equation. It does not establish Lorentz invariance or a physical PSU throat. The 16-dimensional model is not the rank-two Q6 connection Hilbert space of Sections 1–4.

## 6. Self-adjointness: finite graph and continuum candidate

For the finite graph, \(L_U=D_U^*WD_U\) is self-adjoint and non-negative on the entire finite-dimensional Hilbert space. No boundary-domain subtlety arises on a finite graph without external boundary vertices.

For a continuum candidate, let \((M,g)\) be a smooth Riemannian manifold and \(E\to M\) a Hermitian vector bundle with a unitary connection \(\nabla^E\). Set
\[
\widetilde E=E\oplus\overline E,\qquad
\mathcal H=L^2(M,\widetilde E;d\mu_g).
\]
On a compact manifold without boundary, the connection Laplacian \(\nabla^*\nabla\), with its standard Sobolev domain \(H^2(M,\widetilde E)\), is self-adjoint and non-negative. On a compact manifold with smooth boundary, boundary conditions (Dirichlet, or suitable Neumann/Robin conditions) and the corresponding operator domain must be specified; the domain cannot be omitted from a self-adjointness claim.

If \(J\) is a unitary bundle map preserving metric, measure, connection and boundary conditions, it commutes with the self-adjoint operator and the even/odd sectors reduce it. This is conditional continuum mathematics: the smooth bundle, connection and boundary data have not yet been derived from the finite Q6 seed.

## 7. Continuum-limit gate

The fixed 64-vertex Q6 graph has no continuum limit by itself. A continuum theorem requires a specified refinement sequence \(G_h\), not merely a larger numerical sample.

A defensible route is to: (1) specify a compact target \((M,g)\), measure \(d\mu_g\), bundle \(E\), connection and deck involution; (2) equip each refined graph with metric, vertex/edge measures and fibre transports; (3) define reconstruction maps \(J_h:\mathcal H_h\to L^2(M,E\oplus\bar E)\) and prove uniform norm/energy comparisons; (4) define closed non-negative forms \(q_h\) and prove Mosco liminf and recovery in the same Hilbert space; (5) deduce strong-resolvent convergence under the standard hypotheses; and (6) prove that the discrete involutions intertwine with \(J\) under reconstruction, so the even/odd sectors converge too.

The repository's curved Ellis–RGS P1 refinement programme supplies a separate scalar common-\(L^2\) framework. Its curved Mosco liminf/recovery theorem remains open. It is not a spinorial or two-sheeted Mosco theorem and cannot be cited as one.

## 8. Relationship with PRISM and the Pellis Functions

There are three distinct levels:

- **Exact finite graph algebra:** Q6 incidence, its graph Laplacian, and QLG SU(2)/conjugate-SU(2) transports are explicit finite objects.
- **Exact finite internal exchange algebra:** the 8-dimensional Twin/Anti-Twin involution and 16-dimensional effective Hamiltonian satisfy the identities in Section 5.
- **Cross-model and continuum correspondence:** a map identifying the Q6 sheet involution with the internal exchange \(E\), a PRISM-equivariant operator map, and a refinement limit have not been established.

The next mathematical test is to construct a specified linear map between the relevant Hilbert spaces and check its intertwining equations, for example
\[
T J_{Q_6}=\mathcal J T,\qquad
T L_{Q_6}=H_{\rm target}T,
\]
or the appropriate relation if the operators have different orders or spectral parameters. Dimensions, representations, domains and inner products must match. Until such a map is exhibited and verified, the relationship to the Pellis Functions is a research question, not a result.

The PRISM comparison must use its actual graph/cell complex and operator, preserving the distinction between the 69,120-vertex periodic PRISM graph and the 64-vertex Q6 seed. A shared appearance of the golden ratio or an icosahedral symmetry is not an intertwiner.

## 9. Evidence ledger

| Statement | Status | Scope |
|---|---|---|
| Q6 has 64 vertices and 192 edges | PROVED / exact finite combinatorics | Recovered canonical graph |
| Archived conjugate SU(2) transports are unitary | NUMERICALLY VERIFIED | Finite connection |
| The exchange \(J\) on \(E\oplus\bar E\) is a unitary involution | PROVED | Defined direct-sum construction |
| The covariant graph form gives a non-negative self-adjoint operator | PROVED | Once positive weights and transports are fixed |
| The archived 8D exchange and odd observables satisfy the stated identities | NUMERICALLY VERIFIED | Separate internal space |
| The 16D effective Hamiltonian has the exact paired spectrum | PROVED | Declared finite-dimensional model |
| A canonical map from Q6's 64 states to the internal 8D space exists | OPEN | No map in recovered archive |
| The graph operator is the unique physical PSU operator | OPEN | Constitutive/action derivation absent |
| The finite Q6 construction converges to a continuum two-sheeted manifold | OPEN | Refinement, reconstruction and Mosco proof absent |
| A PRISM/Pellis-Functions intertwiner exists | OPEN | Must be constructed and tested |

## 10. Source artifacts and reproducibility

Recovered artifacts used in preparing this note:

- A13_RGS_PSU_QLG_TWINE_ANTITWIN.py — integrated QLG, tangent/coframe, SU(2), Twin/Anti-Twin and finite-graph audit.
- build_psu_fiber_connection.py — geometric SU(2) fibre connection and conjugate biverse sheet.
- derive_RGS_PSU_AntiTwin_throat_symmetry_audit.py — exact finite-dimensional 16D throat-algebra audit.
- RGS_PSU_construct_geometric_intertwiner.py — candidate geometry-driven 64-to-188 feature/intertwiner pipeline, with explicit basis and kernel caveats.
- RGS_PSU_projector_constrained_operator.md — finite-dimensional constrained operator and distinction between finite-dimensional Mosco convergence and the still-open curved refinement limit.
- docs/recovered-rgs-psu-artifacts.md in the companion repository — recovered A13/A14/Q6 archive inventory and provenance.

This note gives Pellis an explicit mathematical object to inspect while identifying exactly which identifications and limits remain to be proved.
