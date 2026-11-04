# Q6 hypercube law, Green–Kubo response and supergolden scaling

**Purpose.** Join the finite Q6 curvature identity, finite-dimensional Green–Kubo response and algebraic supergolden recurrence. This is a reproducible research bridge, not a derivation of physical Yang–Mills dynamics or a universal scaling law.

## 1. Hypercube / plaquette law

The canonical Q6 square 2-skeleton has 64 vertices, 192 edges and 240 square plaquettes. Under the original spectral/normalisation convention for the fourth-moment shift,
\[
\Delta\mu_4=8w^4\sum_{p=1}^{240}(\cos\Omega_p-1),\qquad
\mathcal F_{Q6}=\sum_{p=1}^{240}(1-\cos\Omega_p)
=-\frac{\Delta\mu_4}{8w^4}.
\]
At \(w=2\), \(\Delta\mu_4=-128\mathcal F_{Q6}\). For small fluxes,
\[
\mathcal F_{Q6}=\frac12\sum_p\Omega_p^2-\frac1{24}\sum_p\Omega_p^4
+O(\|\Omega\|^6).
\]
The cosine identity is exact once fluxes and spectral convention are fixed. It does not prove convergence to continuum Yang–Mills.

### Closed 2-cochain sector

For the full six-cube, \(B_3:C_3\to C_2\) has shape 240 by 160. The orthogonal projector onto \(\ker B_3^T\) is
\[
P_{\rm cl}=I-B_3(B_3^TB_3)^+B_3^T.
\]
The projected quadratic diagnostic is \(\frac12\Omega^TP_{\rm cl}\Omega\). It is a sector projection, not an algebraic rewrite of the unprojected cosine action. Use it only when the model explicitly restricts to closed 2-cochains. Its rank is 129.

The full six-cube has cell counts (64,192,240,160,60,12,1), boundary ranks (63,129,111,49,11,1), and Betti numbers (1,0,0,0,0,0,0). The Q6 2-skeleton alone has b2=111; higher cells fill those cycles.

## 2. Green–Kubo law

\[
\eta=\frac{1}{k_BT V}\int_0^\infty
\langle\delta\Pi_{xy}(t)\delta\Pi_{xy}(0)\rangle\,dt
=-\lim_{\omega\to0^+}\frac{\operatorname{Im}G^R_{\Pi\Pi}(\omega,0)}{\omega}.
\]
For \(C(t)=\sum_n A_ne^{-\gamma_nt}\), with positive rates,
\[
\eta=\frac{1}{k_BT V}\sum_n\frac{A_n}{\gamma_n}.
\]
For symmetric OU relaxation,
\[
C_\Pi(t)=b^Te^{-At}\Sigma b,\quad A=A^T\succeq0,\quad\Sigma=\Sigma^T\succeq0,
\]
and eigenpairs \((\gamma_n,v_n)\),
\[
A_n=(b^Tv_n)(v_n^T\Sigma b),\quad
\int_0^\infty C_\Pi(t)dt=b^TA^+\Sigma b,
\]
provided the zero-mode correlation component vanishes. Otherwise the integral diverges; code reports that rather than hiding it with a pseudoinverse. For \(C(t)=G_0e^{-2\kappa t}\), \(\eta=G_0/(2\kappa k_BT V)\). The shear observable is \(\Sigma_{xy}=-\partial_\gamma H|_0\).

A candidate Q6 relaxation operator is \(A=M^{1/2}P_{\rm cl}M^{1/2}\), but a physically justified mobility M, covariance Sigma, stress b, temperature and volume must be supplied. No physical parameters are silently selected.

## 3. Three supergolden results

Let \(\psi>1\) solve
\[
\psi^3-\psi^2-1=0,\qquad\psi=1.465571231876768\ldots
\]
with companion matrix
\[
A_\psi=\begin{pmatrix}1&0&1\\1&0&0\\0&1&0\end{pmatrix},
\quad\det(xI-A_\psi)=x^3-x^2-1.
\]

**Result I — shear response.** \(\Sigma_{\rm NS}=-\partial_\gamma H\); a geometric/gauge stress can be \(T_g=-\partial_\gamma S_W\) when a differentiable strained action is specified. This is a response definition; existence of thermodynamic/continuum limits requires regularity and convergence hypotheses.

**Result II — supergolden geometric series.** If \(a_{n+1}=a_n/\psi\), then
\[
\sum_{n=0}^{\infty}a_n=a_0\frac{\psi}{\psi-1}.
\]
This is an exact algebraic theorem. It does not prove that a physical spectrum, turbulent cascade, hypercube refinement or Yang–Mills coupling selects this ratio. The code offers \(\psi^{-n}\) as a testable refinement schedule and compares measured spectral ratios with \(\psi^2\) and \(\varphi^2\), without fitting or claiming universality.

**Result III — unitary equivalence of Green–Kubo data.** If
\[
H_f=UH_gU^\dagger,\quad\Pi_f=U\Pi_gU^\dagger,\quad
\rho_f=U\rho_gU^\dagger,
\]
then time-evolved stress correlators agree, hence so do Green–Kubo integrals when they exist and the same T,V normalisations are used. The implementation audits intertwining residuals. This conditional theorem does not establish that independently proposed RGS/PSU Hamiltonians are unitarily equivalent.

## 4. Evidence ledger

| Statement | Status |
|---|---|
| Q6 counts and full-cube boundary ranks | PROVED finite combinatorics |
| Cosine identity for supplied fluxes and fixed spectral convention | EXACT algebraic identity |
| Taylor expansion through quartic order | PROVED analytic expansion |
| Closed-sector projector formula | PROVED finite linear algebra |
| Finite OU Green–Kubo modal integral with zero-mode test | PROVED under stated matrix hypotheses |
| Supergolden polynomial and geometric-series sum | PROVED algebra |
| RGS refinement spectrum approaches \(\psi^2\) | OPEN / numerical test |
| Green–Kubo coefficient is physical PSU viscosity | OPEN / requires calibrated dynamics |
| Hypercube cosine action tends to continuum Yang–Mills | OPEN |
| Supergolden scaling is universal in turbulence/QCD/gravity | PHYSICAL HYPOTHESIS, not established |

## Reproduction

    PYTHONPATH=. python -m pytest tests/test_q6_hypercube_gk_supergolden.py -q
    PYTHONPATH=. python examples/q6_hypercube_gk_supergolden_audit.py
