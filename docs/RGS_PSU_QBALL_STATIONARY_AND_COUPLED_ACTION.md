# RGS–PSU–Q6: stationary Q-ball boundary-value problem and shared gauge/geometry action

## 1. What changed

`src/rgs_psu_qball_bvp.py` solves the nonlinear, three-dimensional, flat-space,
ungauged global-(U(1)) Q-ball radial boundary-value problem. Unlike the
previous Gaussian trial, the returned profile satisfies the nonlinear
Euler–Lagrange equation to the collocation tolerance on a finite radial
interval.

For
[
\mathcal L=\partial_\mu\Phi^*\partial^\mu\Phi-V(|\Phi|),\qquad
\Phi(t,r)=f(r)e^{i\omega t},
]
and
[
V(f)=m^2f^2-\lambda f^4+g f^6,
]
the radial equation is
[
f''+\frac{2}{r}f'
=(m^2-\omega^2)f-2\lambda f^3+3g f^5.
]
The solver fixes (f(0)=f_0>0), imposes (f'(0)=0) and (f(R_{\max})=0),
and treats (omega) as an unknown eigenparameter. The regular-singular
radial term is handled by the collocation solver's singular matrix, rather
than by dividing by (r=0).

The energy and charge use
[
Q=8\pi\omega\int_0^{R_{\max}}r^2f^2\,dr,
\qquad
E=4\pi\int_0^{R_{\max}}r^2
\left[\omega^2f^2+(f')^2+V(f)\right]dr.
]

## 2. Numerical validation and interpretation

The regression tests check:
- nontrivial regular profile and both boundary conditions;
- frequency inside the necessary open interval
  (\sqrt{m^2-\lambda^2/(4g)}<\omega<m);
- collocation and independently evaluated radial-equation residuals;
- positive finite energy and Noether charge;
- mesh-refinement convergence of (omega), (E/Q), and (Q).

For the reference test (m=\lambda=1), (g=0.5), (f(0)=1), the
solver converges to (omega\simeq0.905607329) and
(E/Q\simeq0.988291616) on (0\le r\le20). These are numerical
finite-domain results, not universal constants. (E/Q<m) is an energetic
binding diagnostic; it does not establish linear or nonlinear stability.
A fluctuation eigenvalue analysis and radial-domain truncation study remain
necessary before making a stability claim.

## 3. Shared (U(1)) gauge and RGS geometry: action-level structure

The common gauge/geometry backbone can be written, with signature
((-+++)), as the candidate action
[
S=\int d^4x\,e\left[
\frac{M_{\rm Pl}^2}{2}\,T
-\frac14 F_{\mu\nu}F^{\mu\nu}
-g^{\mu\nu}(D_\mu\Phi)^*D_\nu\Phi
-V(|\Phi|)
\right]+S_{\rm matter},
]
where
[
e=\det(e^a{}_{\mu}),\quad
g_{\mu\nu}=\eta_{ab}e^a{}_{\mu}e^b{}_{\nu},\quad
D_\mu=\partial_\mu-iqA_\mu,\quad
F=dA.
]
Here (T) denotes the TEGR torsion scalar in a fixed convention. The
vierbein defines the metric and volume form; the same (A_\mu) supplies
the (U(1)) connection for the charged scalar and the electromagnetic
field strength used by the magnetofluid sector. This avoids introducing
two unrelated electromagnetic fields for Q-ball/fluxon and MHD.

For a separate Madelung fluid ((\rho,S)), a gauge-covariant
nonrelativistic matter action may be taken as
[
S_{\rm QHD}=\int dt\,d^3x\left[
-\rho(\partial_tS+qA_0)
-\frac{\rho}{2m}|\nabla S-q\mathbf A|^2
-\frac{\hbar^2}{8m\rho}|\nabla\rho|^2
-\mathcal U(\rho)
\right].
]
Its phase variation gives charge/particle continuity, and its density
variation gives the quantum Hamilton–Jacobi equation. This sector should
only be added if it represents matter independent of (Phi); using both
for the same condensate would double-count degrees of freedom.

The fluxon sector is not an extra arbitrary field: it is a topological
configuration of the same compact (U(1)) connection/phase, with winding
defined by gauge-invariant link holonomies. The MHD variables are a
fluid velocity and the electric/magnetic decomposition of (F_{\mu\nu}).
Ideal MHD imposes the conducting-fluid constraint (F_{\mu\nu}u^\nu=0);
resistive MHD is dissipative and cannot be inferred from the conservative
action above without an explicit open-system or constitutive prescription.

## 4. Curved, gauged stationary Q-ball equation

For the static spherical metric
[
ds^2=-N(r)^2dt^2+A(r)^2dr^2+r^2d\Omega^2
]
and ansatz
[
\Phi=f(r)e^{i\omega t},\qquad A=A_0(r)dt,
\qquad \Omega(r)=\omega-qA_0(r),
]
the scalar equation becomes
[
\frac{1}{NAr^2}\frac{d}{dr}
\left(\frac{Nr^2}{A}f'\right)
+\frac{\Omega(r)^2}{N^2}f
-\frac12V'(f)=0.
]
Equivalently,
[
f''+\left(\frac{N'}N+\frac2r-\frac{A'}A\right)f'
=A^2\left[\frac12V'(f)-\frac{\Omega^2}{N^2}f\right].
]
This equation exhibits the required coupling: the vierbein/metric changes
the radial measure and kinetic operator, while the gauge potential enters
through the invariant local frequency (Omega). The Maxwell equation
must be solved simultaneously with the scalar equation and sourced by its
Noether current. The current and sign conventions must be fixed from the
same action before coding the coupled solver.

## 5. What is and is not established

**Implemented and numerically tested:** the flat, ungauged stationary
global-(U(1)) Q-ball BVP; regularity; boundary conditions; residuals;
mesh convergence; energy and charge diagnostics.

**Derived action-level target:** common (U(1)) connection for charged
scalar/fluxon and electromagnetic MHD field; vierbein metric and volume
coupling; independent QHD Madelung action; curved/gauged radial scalar
equation.

**Still open:** simultaneous nonlinear Maxwell–Q-ball BVP; curved
RGS/teleparallel backreaction; fluctuation spectrum and stability;
fluxon–Q-ball interaction solutions; full 3-D MHD or HVBK coupling;
a Mosco/continuum limit for the coupled operators. No torsion–scalar
coupling constant is inserted by analogy: any such term must be derived
from the chosen symmetry and variational principles.
