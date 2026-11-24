# PSU–Anti-Twin two-field gauged Q-ball: symmetric stationary branch

**Status:** explicit variational model and symmetric stationary BVP implementation. The finite-radius branch is numerically testable; full two-sheet fluctuation stability, asymmetric branches and vierbein/TEGR backreaction remain open.

## 1. Recover the actual prior model

The potential is the one already specified in \`docs/PSU_ANTITWIN_TWO_FIELD_QBALL_TECHNICAL_NOTE.md\`, not a newly fitted interaction:
\[
I=|\Phi_+|^2+|\Phi_-|^2,\qquad
J=|\Phi_+|^2-|\Phi_-|^2,
\]
\[
\boxed{V=m^2I-\lambda I^2+gI^3+\alpha J^2
-\kappa(\Phi_+^*\Phi_-+\Phi_-^*\Phi_+).}
\]
It is invariant under the diagonal global/gauged \(U(1)\) and the exchange \(Z_2\). Here \(\alpha\geq0\) penalizes sheet-density imbalance and \(\kappa\geq0\) phase-locks the in-phase branch. They remain model parameters; no material calibration is claimed.

Use signature \((-+++)\), \(D_\mu=\partial_\mu-iqA_\mu\), \(\Phi_\pm=f_\pm(r)e^{i\omega t}\), \(A=A_0(r)dt\), and \(\Omega=\omega-qA_0\). The action is
\[
S=\int d^4x\left[-\frac14F_{\mu\nu}F^{\mu\nu}
-\sum_{s=\pm}(D_\mu\Phi_s)^*D^\mu\Phi_s-V\right].
\]

## 2. Coupled radial equations

For real stationary amplitudes, \(I=f_+^2+f_-^2\), \(J=f_+^2-f_-^2\). Variation gives
\[
f_+''+\frac2r f_+'=
(m^2-\Omega^2-2\lambda I+3gI^2+2\alpha J)f_+-\kappa f_-,
\]
\[
f_-''+\frac2r f_-'=
(m^2-\Omega^2-2\lambda I+3gI^2-2\alpha J)f_- -\kappa f_+,
\]
\[
A_0''+\frac2rA_0'=-2q\Omega I.
\]
The common Noether and electric charges are distinct:
\[
Q_N=8\pi\int_0^{R_{\max}}r^2\Omega I\,dr,\qquad
Q_e=qQ_N=-4\pi R_{\max}^2A_0'(R_{\max})
\]
(up to changing all linked signs if the convention for \(D_\mu\) is changed).

The energy is
\[
E=4\pi\int_0^{R_{\max}}r^2\left[
\Omega^2I+(f_+')^2+(f_-')^2+V(f_+,f_-)
+\frac12(A_0')^2\right]dr.
\]

## 3. Exact exchange-symmetric reduction

On the in-phase branch
\[
f_+=f_-=\frac{F}{\sqrt2},\qquad I=F^2,\qquad J=0,
\]
the \(\alpha J^2\) term vanishes on the background and the exchange term shifts the quadratic mass:
\[
V_{\rm sym}(F)=(m^2-\kappa)F^2-\lambda F^4+gF^6.
\]
Consequently,
\[
F''+\frac2rF'=[m_{\rm eff}^2-\Omega^2]F-2\lambda F^3+3gF^5,
\quad m_{\rm eff}^2=m^2-\kappa,
\]
\[
A_0''+\frac2rA_0'=-2q\Omega F^2.
\]
The necessary asymptotic frequency window is
\[
\boxed{\sqrt{m^2-\kappa-\frac{\lambda^2}{4g}}<\omega<
\sqrt{m^2-\kappa}}
\]
provided \(m^2-\kappa>\lambda^2/(4g)\). For the earlier \(m=1,\lambda=1.5,g=1\), the uncoupled \(\kappa=0\) interval is \(0.6614378<\omega<1\), as previously recorded. With \(\kappa=0.1\), the interval becomes approximately \(0.5809475<\omega<0.9486833\). These are necessary windows, not existence or stability theorems.

The implementation solves the two radial scalar equations and Maxwell equation together on the exchange-symmetric branch. It uses the one-field global BVP at \(m_{\rm eff}\) as a continuation seed, not as a substitute for the two-field equations.

## 4. Common and relative amplitude curvatures

At the symmetric background, the local amplitude Hessian diagonalizes in the exchange-even/odd basis. Its pointwise curvature coefficients are
\[
\boxed{h_+(r)=m^2-\kappa-\Omega^2-6\lambda F^2+15gF^4,}
\]
\[
\boxed{h_-(r)=m^2+\kappa-\Omega^2+
(-2\lambda+4\alpha)F^2+3gF^4.}
\]
The first is the common-amplitude direction; the second is the relative PSU/Anti-Twin amplitude direction. These are local potential-Hessian coefficients only. They are **not** eigenvalues of the full radial fluctuation operator: gradients, phase/gauge perturbations, Gauss constraints and charge-preserving boundary conditions still have to be included before any stability conclusion.

A useful exact consequence is that \(\alpha\) does not change the symmetric background, but it shifts the relative-amplitude curvature by \(4\alpha F^2\). Thus \(\alpha\) is invisible to the symmetric profile and directly controls the local cost of relative-sheet amplitude fluctuations.

## 5. Reproducible benchmark and classification

The test benchmark uses \(m=1,\lambda=1.5,g=1,\alpha=0.2,\kappa=0.1,q=0.05,F(0)=0.5,R_{\max}=20\). These \(\alpha,\kappa,q,F(0)\) values are explicit numerical test inputs, not values inferred from PSU physics.

The local development run gave approximately:
- \(\omega=0.89627868\);
- \(Q_N=67.60233\);
- \(E=65.53969\), \(E/Q_N=0.96948857\);
- Maxwell electric energy \(0.10847874\);
- Gauss-flux mismatch \(8.6\times10^{-9}\);
- exchange-profile mismatch at roundoff.

Independent residuals from the collocation interpolant were about \(4.4\times10^{-5}\) for each scalar equation and \(1.1\times10^{-5}\) for Maxwell. The repository tests set conservative tolerances around these values. These are finite-domain numerical checks, not a continuum existence theorem, and do not prove linear stability.

**Established by derivation:** common \(U(1)\), exchange symmetry, radial Euler–Lagrange equations, symmetric mass shift, necessary frequency window, charge/flux relation and local even/odd Hessian coefficients.

**Numerically testable in this implementation:** finite-radius symmetric branch, scalar and Maxwell residuals, Gauss law, \(q\to0\) recovery of the mass-shifted global Q-ball, and the \(\alpha\)-dependence of the relative curvature.

**Still open:** asymmetric stationary branches; full constrained fluctuation spectrum and stability; radial-domain and mesh convergence study; Maxwell magnetic sector; Q-ball/fluxon interaction; common vierbein/TEGR backreaction; radion dimensional reduction and stabilization; coupling to the Q6/T080/T088 discretizations through an action-preserving field/chain map.
