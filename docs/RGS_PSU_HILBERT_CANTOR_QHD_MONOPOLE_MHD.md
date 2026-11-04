# RGS–PSU–Q6 integration: Hilbert, Cantor, QHD, monopoles, MHD and gauge

This extension joins the rank-three discrete bundle/multifractal diagnostic
to finite symbolic geometry, quantum-fluid equations and gauge/fluid
benchmarks. It deliberately distinguishes exact finite identities from
continuum and physical claims.

## 1. Hilbert ordering and Cantor addresses

Each six-bit Q6 word $g\in\{0,1\}^6$ is used as a level-3 Hilbert distance
$d=g$, producing a bijection to the $8\times8$ grid. The same bits define
the finite middle-third Cantor coordinate

\[
x_C(g)=\sum_{k=0}^{5}\frac{2b_k}{3^{k+1}}.
\]

For branch probabilities $p_0,p_1>0$, $p_0+p_1=1$, cylinder masses are
$\mu([b_1\ldots b_n])=\prod_k p_{b_k}$. The 64 words are only finite
cylinders; the infinite Cantor measure requires all levels and a consistent
branch process. Hilbert ordering supplies locality of the finite address
curve, not a proof that the Cantor embedding is geometrically local.

## 2. Quantum hydrodynamics

For $\Psi=\sqrt{\rho}\exp(iS/\hbar)$ and $v=\nabla S/m$, the Madelung
equations (with external potential $V$) are

\[
\partial_t\rho+\nabla\cdot(\rho v)=0,
\qquad
m(\partial_t v+v\cdot\nabla v)=-\nabla(V+Q),
\]
\[
Q=-\frac{\hbar^2}{2m}\frac{\Delta\sqrt{\rho}}{\sqrt{\rho}}.
\]

The code supplies a periodic 1-D centered-difference benchmark. A gauge
coupled version replaces derivatives by covariant derivatives
$D_i=\partial_i-iqA_i/\hbar$; the exact discrete U(1) covariance test is
implemented separately on link variables.

## 3. Monopoles and topology

For a closed oriented 2-cycle,

\[
c_1=\frac{q}{2\pi\hbar}\int_\Sigma F\in\mathbb Z
\]

in the convention where the flux has been normalized to include the
coupling. The audit accepts oriented face fluxes and checks integer
quantization. It does **not** construct a global vector potential: a
nonzero Chern class obstructs a single global potential and requires
overlapping gauge patches / transition functions. This is the discrete
Dirac-monopole gate to be completed on the actual Q6/RGS 2-complex.

## 4. Magnetohydrodynamics

The finite benchmark uses the reduced 1-D transverse system

\[
\partial_t u+u\partial_xu
=-\frac{B\partial_xB}{\mu_0\rho}+\nu\partial_x^2u,
\qquad
\partial_t B+\partial_x(uB)=\eta\partial_x^2B.
\]

It includes the current $J=\mu_0^{-1}\partial_xB$ and Lorentz-force
density $-B\partial_xB/\mu_0$. This is a reduced periodic model, not the
full three-dimensional incompressible MHD system, nor helium-II HVBK.
For a superfluid/normal mixture the appropriate extension is two-fluid
hydrodynamics with quantized vortex filaments and mutual friction.

## 5. Gauge–fluid coupling and unified action target

The link field $U_{vw}\in U(1)$ transports the charged order parameter.
A covariant edge difference is tested under
$\psi_v\mapsto g_v\psi_v$ and
$U_{v\leftarrow w}\mapsto g_vU_{v\leftarrow w}g_w^{-1}$.
The next action-level model should vary

\[
S= S_{\rm RGS}[e,\omega]+S_{\rm YM}[U]
  +S_{\rm QHD}[\rho,S,A]+S_{\rm MHD}[u,B]
  +S_{\rm int}[e,\omega,U,\rho,S,u,B].
\]

The interaction terms cannot be selected by analogy: they must be fixed by
symmetry, dimensions, and a variational principle. The U(1) monopole sector
is topological; SU(2)/Spin and SU(3) colour require their own associated
bundles and representations. The canonical RGS frame connection found in
the previous audit is flat; nonzero gauge curvature must come from
independent dynamical link variables.

## 6. Evidence ledger

- **Exact finite construction:** Q6-to-Hilbert bijection; Cantor finite
  addresses; Bernoulli cylinder normalization; gauge covariance of U(1)
  link differences.
- **Numerical benchmark:** periodic 1-D Madelung and reduced MHD right-hand
  sides; discrete flux integer check.
- **Not established:** an RGS-derived invariant multifractal measure; a
  nonzero monopole Chern class on the actual RGS/Q6 complex with patch data;
  full 3-D MHD/HVBK; coupled QHD–Yang–Mills equations; TEGR/YM/QHD/MHD
  continuum limit; Mosco convergence of the coupled system.


## 7. Fluxons: compact phase winding and sine-Gordon solitons

Here *fluxons* means quantized magnetic-flux kinks of a compact U(1)
phase field, as in the long-junction sine-Gordon limit. On a periodic
link lattice, with $U_{i\leftarrow i+1}=e^{ia_i}$ and
$\psi_i=e^{i\theta_i}$, define the gauge-invariant compact phase drop

\[
\delta_i=\operatorname{Arg}
\exp\{i(\theta_{i+1}-\theta_i+a_i)\}.
\]

It is invariant under $\theta_i\mapsto\theta_i+\chi_i$ and
$a_i\mapsto a_i+\chi_i-\chi_{i+1}$. The fluxon winding audit sums
$\sum_i\delta_i/(2\pi)$ and declares it quantized only when the
background link holonomy is trivial. A nontrivial background holonomy
changes the global sector; the code deliberately refuses to label a
noninteger result as an integer fluxon number.

In the continuum long-junction limit the static phase obeys

\[
\lambda_J^2\,\partial_x^2\theta=\sin\theta,
\qquad
\theta_K(x)=4\arctan e^{(x-X)/\lambda_J}.
\]

The normalized magnetic profile is
$B/\Phi_0=(2\pi)^{-1}\partial_x\theta_K$ and its integral is one
flux quantum $\Phi_0$. The implementation supplies the analytic kink,
its magnetic profile, a finite-difference sine-Gordon residual, and
gauge-invariant winding tests.

This is a controlled U(1)/Josephson sine-Gordon sector. It is not yet a
derivation of flux tubes in SU(3) QCD, electroweak strings, RGS torsion
defects, or a three-dimensional Abrikosov vortex lattice. Those require
their own gauge group, energy functional, boundary conditions, and
topological classification.


## 8. Q-balls: localized charged scalar configurations

Introduce a complex scalar field with global U(1) symmetry,
$\\Phi(t,r)=f(r)e^{i\\omega t}$, and the sextic potential

\\[
V(f)=m^2f^2-\\lambda f^4+g f^6,
\\qquad m,\\lambda,g>0.
\\]

With the convention $\\mathcal L=\\partial_\\mu\\Phi^*
\\partial^\\mu\\Phi-V(|\\Phi|)$, the three-dimensional radial field
equation is

\\[
f''+\\frac{2}{r}f'
=\\frac12 V'(f)-\\omega^2f
=(m^2-\\omega^2)f-2\\lambda f^3+3g f^5.
\\]

The Noether charge and energy are

\\[
Q=2\\omega\\,4\\pi\\int_0^\\infty r^2f^2\\,dr,
\\]
\\[
E=4\\pi\\int_0^\\infty r^2
\\left[\\omega^2f^2+(f')^2+V(f)\\right]dr.
\\]

A necessary existence window follows from
$\\min_{f>0}V(f)/f^2=m^2-\\lambda^2/(4g)$:
\\[
\\sqrt{m^2-\\frac{\\lambda^2}{4g}}<\\omega<m,
\\]
with a positive lower endpoint requiring $\\lambda^2<4gm^2$.
This is an existence criterion, not a proof of linear or nonlinear stability.

The implementation provides the potential, its derivative, this frequency
window, the radial equation residual, and radial energy/Noether-charge
integrals. A Gaussian profile is included strictly as a variational trial
for exercising the diagnostics: it is **not** labelled a stationary
Q-ball solution. A true solution requires solving the nonlinear boundary
value problem with $f'(0)=0$, $f(\\infty)=0$, then checking convergence,
$E(Q)$ and the fluctuation spectrum.

The current Q-ball sector is the ungauged global-U(1) model. Gauging this
U(1) introduces Maxwell dynamics and Coulomb energy and can alter or
destroy the stable branch; it must not be obtained by merely replacing
$\\omega$ with a covariant frequency in the scalar equation. The coupling
to the existing fluxon sector is therefore an action-level task still open.
No identification is made between Q-balls, sine-Gordon fluxons, QCD
chromoelectric tubes, or RGS torsion defects.
