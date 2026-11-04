# PSU–Anti-Twin: two-field gauged Q-ball technical specification

**Status:** candidate mathematical model for independent verification, not a claim of new physical effects.  
**Conventions:** natural units c = hbar = 1; metric (-,+,+,+); D_mu = nabla_mu - i q A_mu; F_mu_nu = 2 partial_[mu A_nu].

## 1. Fields, symmetry and action

Use the minimal M0 model: one spacetime metric, one diagonal U(1) connection, and two complex scalar fields Phi_+, Phi_-. The discrete exchange X(Phi_+,Phi_-)=(Phi_-,Phi_+) obeys X^2=1 and is independent of the gauge symmetry. It does not mean antimatter, opposite electric charge, or a white-hole sector.

Define I=|Phi_+|^2+|Phi_-|^2 and J=|Phi_+|^2-|Phi_-|^2. A minimal exchange- and gauge-invariant potential is
[
V=m^2 I-\lambda I^2+gI^3+\alpha J^2
-\kappa(\Phi_+^*\Phi_-+\Phi_-^*\Phi_+).
]
The flat-space action is
[
S=\int d^4x\left[-\frac14F_{\mu\nu}F^{\mu\nu}
-\sum_{s=\pm}(D_\mu\Phi_s)^*D^\mu\Phi_s-V\right].
]
Gravity is omitted at this stage. For g>0 the sextic term controls large-field behaviour, but Q-ball existence still requires a separate potential/frequency analysis and a converged boundary-value solution.

## 2. Euler–Lagrange equations

The scalar equations are
[
D_\mu D^\mu\Phi_s-\frac{\partial V}{\partial\Phi_s^*}=0,
]
with
[
\frac{\partial V}{\partial\Phi_+^*}
=(m^2-2\lambda I+3gI^2+2\alpha J)\Phi_+-\kappa\Phi_-,
]
[
\frac{\partial V}{\partial\Phi_-^*}
=(m^2-2\lambda I+3gI^2-2\alpha J)\Phi_--\kappa\Phi_+.
]
The diagonal phase symmetry gives
[
j_N^\mu=i\sum_s[\Phi_s^*D^\mu\Phi_s-(D^\mu\Phi_s)^*\Phi_s],
\qquad \nabla_\mu j_N^\mu=0.
]
For kappa nonzero, the relative phase symmetry is broken; only the diagonal U(1) charge is conserved. For kappa=0 the symmetry enlarges to U(1)_+ x U(1)_-.

Varying the stated action with D_mu=nabla_mu-iqA_mu gives
[
\nabla_\mu F^{\mu\nu}=+qj_N^\nu.
]
The sign is convention-dependent: changing the sign in D_mu changes the assigned electric-charge sign. The Maxwell equation, current, and Gauss-flux diagnostic must be changed together.

For Phi=f(r) exp(i omega t), A=A_0(r)dt, Omega=omega-qA_0, one has j_N^0=2 Omega f^2 and
[
A_0''+\frac2rA_0'=-2q\Omega f^2.
]
The outward electric flux is Q_flux=-4 pi R^2 A_0'(R)=+q Q_N, where Q_N=4 pi integral_0^R 2 Omega f^2 r^2 dr.

**Repository audit:** an earlier patch changed the radial source to the positive sign while leaving the independent residual formula at the negative sign. That mismatch caused the reported residual of about 0.1826. The solver is now being aligned with the action-derived negative radial source, positive flux charge Q_e=+q Q_N, and matching initial iterate. The CI must confirm the correction.

## 3. Exchange-symmetric stationary branch

Choose Phi_+=Phi_-=f(r)exp(i omega t)/sqrt(2). Then I=f^2, J=0, and the relative phase is zero. The coupling shifts the quadratic coefficient to m_eff^2=m^2-kappa. The radial equations are
[
f''+\frac2r f'=(m_{\rm eff}^2-\Omega^2)f-2\lambda f^3+3gf^5,
]
[
A_0''+\frac2rA_0'=-2q\Omega f^2,\qquad \Omega=\omega-qA_0.
]
Regularity requires f'(0)=A_0'(0)=0. At a finite outer radius use f(Rmax) approximately zero and A_0(Rmax)=0; the latter fixes residual constant-gauge freedom. Convergence in Rmax is mandatory.

A necessary asymptotic frequency window for this symmetric branch is
[
\sqrt{m_{\rm eff}^2-\frac{\lambda^2}{4g}}<\omega<m_{\rm eff},
]
when m_eff^2>lambda^2/(4g). This is necessary, not sufficient for existence or stability. The relative-phase-pi branch has a different effective mass shift and must be derived separately.

## 4. Stability and superfluid limit

A converged stationary BVP is not a stability proof. Linearize the full time-dependent two-field Maxwell system, retain amplitude, phase and gauge perturbations, remove pure-gauge modes, and enforce Gauss' constraint. Report physical eigenvalues, negative/zero modes, mesh and domain convergence, and the limits q->0, kappa->0, alpha->0. Compare with global and gauged one-field Q-balls.

Where Phi_s is nonzero, write Phi_s=sqrt(n_s) exp(i theta_s). In the nonrelativistic regime, after removing the rest-mass phase, v_s=(grad theta_s-q A)/m_eff. Continuity and phase equations contain the quantum-pressure term Q_s=-(1/(2m_eff)) Delta sqrt(n_s)/sqrt(n_s). This is a hydrodynamic rewriting of the same condensate, not an additional fluid. A normal component requires independent degrees of freedom and a finite-temperature closure.

## 5. Gravity and horizon criteria

Only after the flat-space system passes its field-equation, Gauss-law, convergence and linear-stability tests should its stress tensor be coupled to a dynamical vierbein/metric. In spherical symmetry use future-directed radial null expansions
[
\theta_\pm=\frac2R k_\pm^\mu\nabla_\mu R.
]
With a declared orientation, trapped regions have theta_+<0 and theta_-<0; anti-trapped regions have theta_+>0 and theta_->0. A PSU–Anti-Twin bounce requires a regular dynamical solution connecting the relevant regions, while preserving gravitational constraints, Gauss law, energy/charge balance and perturbative stability. None of this follows merely from the existence of a Q-ball.

## 6. Reproducibility gates

1. Freeze conventions and derive all equations from the same action.
2. Verify scalar residuals, Maxwell residuals and integrated Gauss law independently.
3. Recover the global Q-ball at q=0.
4. Demonstrate radial-domain and mesh convergence.
5. Construct two-component stationary branches and their linearized spectra.
6. Only then add vierbein/TEGR backreaction and horizon evolution.

**Classification:** the action and Euler–Lagrange equations are a specified candidate model; the symmetric radial reduction follows from it. Existence and stability of two-component gauged branches, a controlled superfluid regime, and any black-to-white-hole transition remain open questions.
