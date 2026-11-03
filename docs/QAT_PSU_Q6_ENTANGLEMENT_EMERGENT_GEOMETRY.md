# QAT–PSU–Q6 entanglement bridge and emergent geometry

## Editorial status

This is Phase 20 of the same RGS–PSU–PRISM–QLG research programme. It is an internal technical insert in the existing 27-chapter monograph, not a new chapter or a separate theory. It carries forward the Phase 18 recovered \(O_h\) action and Phase 19 \(Q_6\to SU(2)^{192}\to U_p^{240}\to S_W\) construction.

## 1. Physical question and scope

QAT proposes a geometric picture involving a cube, an inscribed sphere, a null cone, a central event and two PSU–Anti-Twin sheets. The holographic-entanglement literature asks when quantum correlations between two sectors encode the connectedness of an emergent bulk geometry. The bridge to be tested is therefore:

\[
\text{QAT two-sheet kinematics}
\longrightarrow |\Psi\rangle\in\mathcal H_+\otimes\mathcal H_-
\longrightarrow \rho_+=\operatorname{Tr}_-|\Psi\rangle\langle\Psi|
\longrightarrow S_+=-\operatorname{Tr}(\rho_+\log\rho_+)
\longrightarrow \text{candidate geometric observable}.
\]

The last arrow is not automatic. It requires a reconstruction rule from quantum-information data to lengths, areas or a metric, plus a variational/dynamical principle. No equality between QAT ovalization and entanglement has yet been established.

## 2. Exact finite Q6 double-sheet model

Use the already fixed Q6 graph, with \(64\) vertices and combinatorial Laplacian \(L_6\). Its exact spectrum is

\[
\lambda_j=2j,\qquad m_j=\binom 6j,\qquad j=0,\ldots,6.
\]

Define a dimensionless heat-kernel amplitude

\[
K_\tau=e^{-\tau L_6/2},\qquad \tau\geq0,
\]

and the normalized pure state on two 64-dimensional sheets

\[
|\Psi_\tau\rangle=
\frac{1}{\|K_\tau\|_F}
\sum_{x,y=0}^{63}(K_\tau)_{xy}|x\rangle_+|y\rangle_-.
\]

The coefficient matrix is \(C_\tau=K_\tau/\|K_\tau\|_F\). The reduced state and entropy are

\[
\rho_+(\tau)=C_\tau C_\tau^\dagger
=\frac{e^{-\tau L_6}}{\operatorname{Tr}(e^{-\tau L_6})},
\qquad
S_+(\tau)=-\operatorname{Tr}(\rho_+\log\rho_+).
\]

Writing \(q=e^{-2\tau}\), the partition function is

\[
Z_6(\tau)=\sum_{j=0}^6\binom6j q^j=(1+q)^6.
\]

The eigenvalues of \(\rho_+\) are \(q^j/Z_6\), each with multiplicity \(\binom6j\). Therefore

\[
\boxed{S_+(\tau)=6h\!\left(\frac{1}{1+e^{2\tau}}\right)}
\]

where \(h(p)=-p\log p-(1-p)\log(1-p)\). This is an exact finite identity for this state ansatz.

Boundary checks:
- \(\tau=0\): \(K_0=I\), so \(C_0=I/8\). This is a maximally entangled state across the two sheets and \(S_+=\log64\).
- \(\tau\to\infty\): \(K_\tau\) tends to the projector onto the constant Q6 mode; the normalized coefficient matrix has rank one and \(S_+\to0\).

For a pure bipartite state, \(I(+:-)=2S_+\). The code computes the density matrix, entropy, purity, Schmidt spectrum and residuals, and compares them with the analytic formula.

## 3. Interpretation boundary: QAT and the parameter \(\tau\)

The parameter \(\tau\) is dimensionless and is a spectral correlation-scale parameter in the heat-kernel ansatz. It is **not** identified with proper time, cosmological time, QAT cone radius, the speed of reconnection, or an ovalization angle. A physical map must be derived from the QAT tetrad/coframe and a specified action, for example by defining a geometric scalar \(\chi[e]\) and proving/calibrating a relation \(\tau=f(\chi)\). Until then the relation is a proposed interface, not a result.

The QAT cube/sphere/null-cone construction should be represented by an explicit metric and causal structure. The PSU–Anti-Twin sheets become the two Hilbert factors only after their local degrees of freedom and inner products are specified. The current Q6 ansatz is a controlled finite testbed for that programme.

## 4. Connection to the holographic literature

- **Van Raamsdonk (2010), “Building up spacetime with quantum entanglement”, arXiv:1005.3035.** Motivates the hypothesis that connectedness of semiclassical bulk regions is tied to entanglement between corresponding quantum sectors. It does not supply the QAT-to-Q6 dictionary.
- **Ryu–Takayanagi (2006), arXiv:hep-th/0603001; extended treatment arXiv:hep-th/0605073.** In appropriate AdS/CFT settings, leading boundary entanglement entropy is related to the area of a bulk minimal surface. The finite Q6 model still needs a discrete minimal-surface definition and an independently fixed \(G_{\rm eff}\).
- **Maldacena–Susskind (2013), “Cool horizons for entangled black holes”, DOI 10.1002/prop.201300020.** ER=EPR is a conjectural correspondence in specified gravitational settings, not a general theorem that any entangled pair is a traversable wormhole. PSU–Anti-Twin must not be called a wormhole solely because its state is entangled.
- **Swingle (2012), “Entanglement renormalization and holography”, DOI 10.1103/PhysRevD.86.065007.** Motivates a scale direction built from multiscale entanglement renormalization. This is relevant to Q6/PRISM refinement, but requires explicit isometries/disentanglers and scale-dependent tensors.
- **Faulkner–Lewkowycz–Maldacena (2013), arXiv:1307.2892.** Shows that bulk quantum entanglement contributes beyond the leading area term in holographic entropy. This warns against identifying a finite graph cut with the full entropy.
- **Van Raamsdonk–Vilar López (2026), “Menagerie of Euclidean constructions for 3D holographic cosmologies”, arXiv:2601.10906, Phys. Rev. D 113, 106010.** Studies Euclidean wormhole saddles, heavy-particle closed cosmologies and AdS tubes. Its gluing and saddle-dominance conditions provide a later comparison for the PSU–Anti-Twin throat and the existing Euclidean wormhole calculations. They are not automatically the same construction.

## 5. Proposed falsifiable next tests

1. **QAT calibration:** construct \(\chi[e]\) from the existing tetrad, with invariance properties stated, and test whether one fixed map \(\tau=f(\chi)\) predicts the Q6 entropy curve across refinements.
2. **Discrete area:** define a cut/minimal surface \(\gamma_A\) in a specified Q6/PRISM bulk complex. Compute its area from edge/face weights derived from the coframe, not fitted to entropy.
3. **Area-law test:** compare \(S(A)\) with \(\operatorname{Area}(\gamma_A)/(4G_{\rm eff})\), fixing \(G_{\rm eff}\) independently and testing multiple partitions and refinements.
4. **Gauge coupling:** combine the state with the Phase 19 link variables only through a declared gauge-invariant coupling, such as a Wilson-loop-dressed kernel. Verify local gauge invariance before interpreting it.
5. **Continuum/refinement:** track entropy, cut area and their residual under Q6-to-PRISM refinement; a single 64-site fit is not a continuum result.
6. **Euclidean/Lorentzian bridge:** compare the Euclidean saddle action and boundary state with the QAT Lorentzian tetrad and PSU–Anti-Twin matching conditions. State the analytic continuation and boundary data.

## 6. Evidence classification

**Exact finite identity:** Q6 Laplacian spectrum and the entropy formula for the heat-kernel bipartite state above.

**Numerically verified:** matrix exponential, reduced density matrix, entropy and Schmidt spectrum agree with the analytic formula at the tested \(\tau\) values, subject to the CI report.

**Model choice:** the heat-kernel amplitude \(K_\tau\), the two-sheet Hilbert-space factorization and dimensionless \(\tau\).

**Open:** deriving \(K_\tau\) from QAT geometry; identifying the PSU sheets with holographic boundary/bulk sectors; a discrete RT surface; \(G_{\rm eff}\); an area law; a continuum holographic limit; physical wormhole equivalence.

## 7. Reproduction

\`\`\`bash
PYTHONPATH=. python -m pytest tests/test_q6_bisheet_entanglement.py -q
PYTHONPATH=. python examples/q6_bisheet_entanglement_audit.py
\`\`\`

The JSON report is written to \`results/json/q6_bisheet_entanglement.json\`.
